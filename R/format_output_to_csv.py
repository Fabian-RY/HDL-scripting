#! /usr/bin/env python3

## WIP: script para parsear el output de HDL y generar csvs con la información de heredabilidad
## Y correlaciones

import argparse
import gzip
from collections import defaultdict

# Variables que vamos a utilizar como header en el archivo final (en este orden)
filename="file"
heritability_1_header="Heritability_1"
heritability_2_header="Heritability_2"
covariance_header="Covariance"
correlation_header="Correlation"
pvalue_header="pvalue"

# Strings que hay en el archivo y que queremos reconocer
heritability_1_startline="Heritability of phenotype 1:"
heritability_2_startline="Heritability of phenotype 2:"
covariance_startline="Genetic Covariance:"
correlation_startline="Genetic Correlation:"
pvalue_startline="P:"

# Modos de abrir de los archivos para leer o escribir
read_mode="rt"
write_mode="wt"
append_mode="at"

def parseargs():
  """
    parseargs -> Parsea los argumentos usados en linea de comandos
      outfile -> Archivo de salida, opcional, si no se especifica, lo escribe en stdout.
      gz -> Si el archivo está comprimido en gz, entonces se debe indicar True. Por defecto, es texto plano (False)
      separator -> El separador a usar en el archivo de salida. Opcional. Por defecto es coma ",".
      infiles -> Lista de archivos a procesar: Necesario uno o más archivos, se pueden dar usando * en la terminal o uno después de otro
  """
  parser = argparse.ArgumentParser()
  parser.add_argument("--outfile","-o", type=str, action="store",dest="outfile", default="/dev/stdout")
  parser.add_argument("--separator","-sep", type=str, action="store",dest="sep", default=",")
  parser.add_argument("--gz","-gz", type=bool, action="store",dest="gz", default=False)
  parser.add_argument("--files", nargs="+", action="store", type=str, dest="infiles")
  return parser.parse_args()

def parse_HDL_results(file:str, sep:str, gzipped:bool=False) -> str:
  """
    Abre, parsea el archivo y devuelve una string con la información ya procesada.

    Para ello, lee linea a linea y se queda con el elemento que forma el valor deseado: heredabilidad 1 o 2, correlacion,
    covarianza y pvalor

    Utiliza las variables declaradas al inicio del archivo para detectar la linea correcta
  """
  #if gzipped: open=gzip.open() # Sobreescribimos la funcion open solo si está comprimido, para poder ofrecer compatibilidad con gzip
  heritability1 = heritability2 = covariance = correlation  = pvalue = "" # Las declaramos vacías para que si no se encuentra algún valor continúe y no falle
  with open(file, read_mode) as fhandler:
    for line in fhandler:
      if line.startswith(heritability_1_startline): 
        heritability1:str = line.strip("\n").split()[4]        
      elif line.startswith(heritability_2_startline): 
        heritability2:str = line.strip("\n").split()[4]
      elif line.startswith(covariance_startline): 
        covariance:str = line.strip("\n").split()[2]
      elif line.startswith(correlation_startline): 
        correlation:str = line.strip("\n").split()[2]
      elif line.startswith(pvalue_startline): 
        pvalue:str=line.strip("\n").split()[1]
  return sep.join([file, heritability1, heritability2, covariance, correlation, pvalue])
      
def main():
  """
    Función principal:
      - Parsea argumentos
      - Ejecuta los parseos
      - Genera el documento final y lo escribe donde toca

      Map permite ejecutar todas las lineas a la vez
  """
  args = parseargs()
  results_parsed = list(map(lambda x: parse_HDL_results(x, args.sep, args.gz), args.infiles)) # Parseamos los N archivos
  with open(args.outfile, write_mode) as filewriter: # Escribimos el header
      print(",".join([filename, heritability_1_header, heritability_2_header, covariance_header, correlation_header, pvalue_header]), file=filewriter)
  with open(args.outfile, append_mode) as fileappender: # Escribimos el contenido
      print("\n".join(results_parsed), file=fileappender)
  pass

if __name__ == "__main__":
  main()

