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
version <- "0.9"
data.table_minimum_version <- "1.12.1"
package <- "data.table"

# Lee los datos guardados. Es complatible con RDSs y archivos gz y bgz
#
read_GWAS_summary <- function(path){
  path.split <- unlist(strsplit(path, split = "\\."))
  file.type <- path.split[length(path.split)]
  print(file.type)
  print("===============")
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

get_min_SNP <- function(gwas1, gwas2){
  return(min(gwas1$N, gwas2$N))
}

perform_HDL <- function(gwas1, gwas2, LD.path, Nref, N0, output.file, eigen.cut, jackknife.df, fill.missing, intercept){
  res.HDL <- HDL.rg(    gwas1, gwas2, LD.path, Nref, N0, output.file, eigen.cut, jackknife.df, fill.missing, intercept.output = intercept)
  print("done")
  if(output.file != ""){
    fConn <- file(output.file)
    Lines <- readLines(fConn)
    #writeLines(c(args.print,Lines), con = fConn)
    close(fConn)
    if(jackknife.df == TRUE){
      write.table(res.HDL$jackknife.df, file = paste0(output.file,"_jackknife.df",".txt"))
    }
  }
}

# Realiza todo el proceso para la linea del conf en cuestión
calc_correlation <- function(x){
  gwas1 <- read_GWAS_summary(x[1]) # Path del archivo de gwas 1
  gwas2 <- read_GWAS_summary(x[2])# Path del archivo de gwas 2
  LD_path <- x[3] # Path del archivo que contiene la información de LD para HDL
  output_file <- x[4] # Archivo de salida
  log_file <- x[5] # Archivo de log
  N0 <- x[7] # Mínimo número común de SNPs entre ambos datasets
  if(as.numeric(N0) == 0){
    # Si el número es 0, lo imputamos a partir de los datos. 0 funciona como comodín.
    N0 <- get_min_SNP(gwas1, gwas2)
  }
  Nref <- x[6]
  Nref <- ifelse(length(Nref)==0 || is.na(Nref), 335265, as.numeric(Nref))
  eigen_cut <- x[8]
  jackknife <- x[9]
  eigen.cut <- ifelse(length(eigen_cut)==0 || is.na(eigen_cut), "automatic", as.numeric(eigen_cut))
  jackknife.df <- ifelse(length(jackknife)==0 || is.na(jackknife), FALSE, as.logical(jackknife))
  intercept.output <-FALSE
  fill_missing <- x[10][["fill.missing.N"]]
  if(!any(fill_missing==c("median", "min", "max"))) fill_missing <- NULL
  output.file <- ifelse(length(output_file)==0, "", output_file)
  ignore <- as.logical(x[11])
  check <- check_to_skip_calc(output_file, ignore)
  if(check){
    return("")
  }
  print(fill_missing)
  perform_HDL(gwas1, gwas2, LD_path, Nref, N0, output_file, eigen.cut, jackknife.df, fill_missing, intercept.output)
}

args <- commandArgs(trailingOnly = TRUE)
# Todos los argumentos para automatizar deben estár en el archivo csv acompañante.

check_package_version(data.table_minimum_version, package) # Necesita data.table 1.12.1 o superior

#info <- process_args() # -- TO DO
info <- read.csv("R/template_test.csv")
apply(info, 1, calc_correlation)
