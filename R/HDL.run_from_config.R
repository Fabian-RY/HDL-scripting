##########################################################################
#                                                                        #
#   Método HDL para obtención de correlaciones genicas en GWAS summaries #
#   Autor: Fabián Robledo a partir del paquete HDL                       #
#   Automatización del script para N muestras                            #
#                                                                        #
##########################################################################

library(HDL)
library(data.table)

# Constantes
version <- "0.9" # Versión actual del paquete
## Data.table funciona en el script original con 1.12.1, si no busca esa versión instalada
## No he visto motivos suficientes para ello, en este script lo he puesto como warning
data.table_minimum_version <- "1.12.1"
package <- "data.table"

# Solo necesitamos un argumento: el config que contiene todos los parámetros:
# Dicho config es un csv, disponible un ejemplo template_test.csv en el paquete
process_args <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  # Solo permitimos un argumento, el config
  # Ninguno o más de uno finaliza el script
  if (len(args) == 0){ # Si no tiene un csv como argumento no hay nada que procesar
    stop("Argument for the csv config file requirerd")
  }
  else if (len(args) > 1){ # De momento solo aceptará un argumento.
    stop("Only one config allowed")
  }
  else{
    return(read.csv(args[1]))
  }
}

# Función derivada de smart.read del archivo HDL.run.R
# Ahora mismo es equivalente, pero se puede adaptar a más tipos de archivos
read_GWAS_summary <- function(path){
  path.split <- unlist(strsplit(path, split = "\\."))
  file.type <- path.split[length(path.split)]
  if(file.type == "rds"){
    return(readRDS(path))
  } else if(file.type == "gz" | file.type == "bgz"){
    options(datatable.fread.input.cmd.message=FALSE)
    return(fread(input = paste("zcat < ",path)))
  } else{
    try_error <- try(return(fread(path)))
    if(!is.null(try_error)){
      error.message <- "This file type is not supported by fread function in data.table package. Please reformat it to .txt, .csv or .tsv."
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      stop(error.message)
    }
  }
}
# A veces es posible que no queramos reanalizar un registro concreto
# La variable ignore permite decidir si se reanaliza o si se salta este registro
# Si NO existe resultado, lo ejecuta siempre
# Si existe resultado e ignore es FALSO, lo salta
# Si existe resultado e ignore es TRUE, lo reanaliza y sobreescribe
# La función solo determina si se debe saltar o no, el hecho en sí de saltarse un paso
# se hace en la función principal
check_to_skip_calc <- function(x, ignore=F){
  if(x == ""){
     # Si está vacio, no guardamos el archivo y no hay nada que comprobar
  }
  else if(file.exists(x) && !(ignore)){ # No sobreescribe si no lo ignoras explicitamente
   message(paste0("Output file ", x, " already exists. Not overwriting. Skipping..."))
   return(TRUE) # True s
  }
  else if(file.exists(x) && ignore){
    message(paste0("Output file ", x, "already exists. Overwriting existing value"))
  }
  # Si no existe el archivo. o decidimos ignorarlo, avanzamos y lo sobreescribimos o creamos
  return(FALSE)
}

# Necesita datatable 1.12.1 o superior. Comprobamos que exista y si es una versión anterios avisamos
check_package_version <- function(x, package){
  package.version <- packageVersion(package)
  if(package.version < x){
    message(paste0(package, " detected with an old version: (", package.version ," . Recommended is at least ",x))
  }
}

# HDL utiliza el mínimo de SNP de ambos gwas
get_min_SNP <- function(gwas1, gwas2){
  return(min(length(gwas1$N), length(gwas2$N)))
}

# Función principar para realizar el HDL. Utiliza los argumentos recibidos en el config
# Después, escribe el resultado
perform_HDL <- function(gwas1, gwas2, LD.path, Nref, N0, output.file, eigen.cut, jackknife.df, fill.missing, intercept){
  res.HDL <- HDL.rg(    gwas1, gwas2, LD.path, Nref, N0, output.file, eigen.cut, jackknife.df, fill.missing, intercept.output = intercept)
  # Escribimos el resultado (he hecho un apaño porque arg.lines es un valor por defecto no explicado)
  if(output.file != ""){
    fConn <- file(output.file)
    Lines <- readLines(fConn)
    writeLines(c(output.file,Lines), con = fConn)
    close(fConn)
    if(jackknife.df == TRUE){
      write.table(res.HDL$jackknife.df, file = paste0(output.file,"_jackknife.df",".txt"))
    }
  }
}

# Realiza todo el proceso para la linea del conf en cuestión
#x <- info
calc_correlation <- function(x){
  print(x)
  output_file <- x[[4]] # Archivo de salida
  output.file <- ifelse(length(output_file)==0, "", output_file)
  ignore <- as.logical(x[[11]])
  check <- check_to_skip_calc(output.file, ignore) # Comprobamos si el archivo existe
  # Y si es el caso, decidimos si lo reanalizamos. Si no lo reanalizamos, termina el procesamiento
  # Lo miramos al principio para que sea rapido en decidir si lo reanaliza o no
  if(check){
    return("Existe, no se reanalizará")
  }
  if(x[[1]] == x[[2]]){return("")}
  gwas1 <- read_GWAS_summary(x[[1]]) # Path del archivo de gwas 1. Lo cargamos
  gwas1$N <- as.numeric(gwas1$N)
  gwas2 <- read_GWAS_summary(x[[2]])# Path del archivo de gwas 2. Lo cargamos
  gwas2$N <- as.numeric(gwas2$N)
  LD_path <- x[[3]] # Path del archivo que contiene la información de LD para HDL
  log_file <- x[[5]] # Archivo de log
  N0 <- as.numeric(x[[7]]) # Mínimo número común de SNPs entre ambos datasets. Está como argumento, pero se infiere de los datos si es 0.
  if(as.numeric(N0) == 0){
    # Si el número es 0, lo imputamos a partir de los datos. Funciona como comodín "No lo sabemos de antemano".
    N0 <- get_min_SNP(gwas1, gwas2)
  }
  Nref <- x[[6]] # Número de SNP en la referencia. 0 asumen que es la referencia por defecto
               # de HDL que tiene 335265
  Nref <- ifelse(length(Nref)==0 || is.na(Nref), 335265, as.numeric(Nref))
  eigen_cut <- x[[8]] # Parámetro que tengo que dilucidar todavía
  jackknife <- x[[9]] # Idem
  eigen.cut <- ifelse(length(eigen_cut)==0 || is.na(eigen_cut), "automatic", as.numeric(eigen_cut))
  jackknife.df <- ifelse(length(jackknife)==0 || is.na(jackknife), FALSE, as.logical(jackknife))
  intercept.output <-FALSE
  fill_missing <- x[[10]] # Decide cómo imputar datos faltantes:
                                            # con la mediana del dato, el minimo, el máximo
                                            # Cualquier otra cosa hace que no impute y descarte ese registro
  if(!any(fill_missing==c("median", "min", "max"))) fill_missing <- NULL # Si no es una de las opciones disponibles, que sea NULL por defecto
  perform_HDL(gwas1, gwas2, LD_path, Nref, N0, output_file, eigen.cut, jackknife.df, fill_missing, intercept.output)
}

args <- commandArgs(trailingOnly = TRUE)
# Todos los argumentos para automatizar deben estár en el archivo csv acompañante.

check_package_version(data.table_minimum_version, package) # Necesita data.table 1.12.1 o superior

#info <- process_args() # -- TO DO
info <- read.csv("/home/fabian/Escritorio/PhD_data/WP1/datasets_WP1/config.txt")
apply(info, 1, calc_correlation)
