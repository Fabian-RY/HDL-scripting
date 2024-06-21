#! /usr/bin/env python3

import os.path
import os
import itertools
import math
from urllib.request import urlretrieve
from tqdm import tqdm

CONFIG_SEP:str = ","
READ_MODE:str = "rt"
#MODES:str = ["download", "munge","ldsc","HDL"]
#MODES:str = ["munge","ldsc","HDL"]
#MODES:str = ["ldsc",]
#MODES:str = [ "unify-results"]
MODES:str = ["HR"]
CONFIG_FILE:str = "/home/fabian/Escritorio/PhD/data/WP1/analisis/Datasets-info.csv"
RAWDATA_FOLDER:str = "/home/fabian/Escritorio/PhD/data/WP1/analisis/rawdata/"
MUNGE_FOLDER:str = "/home/fabian/Escritorio/PhD/data/WP1/analisis/munged-data/"
LDSC_FOLDER:str="/home/fabian/Escritorio/PhD/data/WP1/analisis/tmp_ldsc/"
RESULTS:str="/home/fabian/Escritorio/PhD/data/WP1/analisis/results/"
PROCESS_ALL:bool = False
HDL_FOLDER:str = "/home/fabian/Escritorio/PhD/data/WP1/analisis/results/HDL"
GWAS_TYPE="UKB.Neale"
MUNGE_PY_PATH:str="singularity exec /home/fabian/Escritorio/PhD/Software/WP1/singularity/lsdc.sif python2 /home/fabian/Escritorio/PhD/Software/WP1/ldsc/munge_sumstats.py"
LDSC_PY_PATH:str="singularity exec /home/fabian/Escritorio/PhD/Software/WP1/singularity/lsdc.sif python2 /home/fabian/Escritorio/PhD/Software/WP1/ldsc/ldsc.py"
HDL_R_PATH:str="Rscript /home/fabian/Escritorio/PhD/Software/WP1/HDL/HDL.run.R"
R_SCRIPT_PATH:str="Rscript"
HEATMAP_SCRIPT_PATH:str="/home/fabian/Escritorio/PhD/Software/WP1/HDL-scripting/R/results.R"
SNP_LIST_PATH:str="/home/fabian/Descargas/w_hm3.snplist"
REF_LD_PATH:str="/home/fabian/Escritorio/PhD/data/WP1/datasets/ldsc_input/eur_w_ld_chr/"
LD_PATH:str="/home/fabian/Escritorio/PhD/data/WP1/HDL/UKB_imputed_SVD_eigen99_extraction/"
WILD_LD_PATH:str="/home/fabian/Escritorio/PhD/data/WP1/datasets/ldsc_input/eur_w_ld_chr/"
CHUNKSIZE:int=50000
LDSC_RESULTS:str="/home/fabian/Escritorio/PhD/data/WP1/analisis/results/"
LDSC_RESULTS_FILE:str="LDSC-results.csv"
HEATMAP_WIDTH:int = 1000 # Píxeles
HEATMAP_HEIGHT:int = 1000 # Píxeles
IMAGE_PATH:str = "/home/fabian/Escritorio/PhD/data/WP1/analisis/results/"
MR_SCRIPT_PATH:str = "/home/fabian/Escritorio/PhD/Software/WP1/HDL-scripting/R/MR.R"
MR_results:str="/home/fabian/Escritorio/PhD/data/WP1/analisis/MR/"

def parse_config(csv_path:str, sep) -> dict:
    """
        Parses the config file.

        Returns a dict in which every key is the ID and every value
        is a list with all the collumns of the file

        Asumes first row is header that has the colnames
    """
    parsed_data = dict()
    with open(csv_path, READ_MODE) as csv:
        for id, line in enumerate(csv):
            if (id == 0):
                header:list = line.strip("\n").split(sep)
                continue
            data:list = line.strip("\n").split(sep)
            id:str = data[0]
            parsed_data[id] = dict(zip(header[1:], data[1:]))
    return parsed_data

def download_datasets(config, sep, outfolder, verbose=False, process_all=False) -> None:
    """
        Downloads datasets from the link indicated in the Sumstats-link column of the
        config file. Datasets with "discarded" flag in config won't be descarded unless
        the process_all flag is used. Datasets will be saved in the outfolder folder

        Has some verbose options, including the reason why it was ignored.
    """
    config_data = parse_config(config, sep)
    ids = list(config_data.keys())
    for idn in tqdm(range(len(ids))):
        id = ids[idn]
        dataset_info = config_data[id]
        url = dataset_info["Sumstats-link"]
        file_name = dataset_info["rawdata"]
        dest_path = os.path.join(outfolder, file_name)
        if(os.path.exists(dest_path) and not process_all): 
            print("Ignoring {}".format(id))
            continue # Si ya existe y no le hemos dicho que reprocese todo, simplemente pasa de ello
            # Ejecutando con process all = true reinicias todo desde el principio
        if(verbose): print("Downloading dataset with ID: ", id)
        if(process_all or int(dataset_info["discarded"])== 0 ):  
            # Si no está descartado y o si forzamos, entonces 
            # Se descarga el archivo
            urlretrieve(url, dest_path)
            if (verbose): print("Downloaded dataset with ID: ", id)
        elif (verbose): print("Dataset with ID discarded for reason: {}".format(dataset_info["Reason"]) )

def munge_datasets(config, sep, rawdata_folder, munge_folder, verbose=True, process_all=False) -> None:
    """
        Procesess original datasets to munge it into an LDSC-compatible input

        Uses info from the config and executes it according to the global config values
    """
    config_data = parse_config(config, sep)
    ids = list(config_data.keys())
    for idn in tqdm(range(len(ids))):
        id = ids[idn]
        dataset_info = config_data[id]
        if(process_all or int(dataset_info["discarded"])== 0 ):
            if (verbose): print("Munging dataset with ID: {}".format(id))
            rawdata = os.path.join(rawdata_folder, dataset_info["rawdata"])
            outfile = os.path.join(munge_folder, dataset_info["munged-dataset"])
            snpid = dataset_info["variant-id"]
            Ncas = dataset_info["Ncas"] 
            Ncon = dataset_info["Ncol"]
            N = dataset_info["N"]
            if(os.path.exists(outfile+".sumstats.gz") and not process_all):
                print("Ignoring dataset")
                continue
            elif (Ncas != "-" and Ncon != "-" ):
                # In case N cas and N con are available, uses them
                Ncommand = "--N-cas-col {} --N-con-col {}".format(Ncas, Ncon)
            elif (not N == "-" ):
                # In case N is available, uses N
                Ncommand = "--N-col {}".format(dataset_info["N"])
            elif (dataset_info["sample"] != "-"): 
                # Else uses the number of individuals as an aproximation
                print("Unable to determine N. Aproximating using number of patients")
                N = int(dataset_info["sample"])
                Ncommand = "--N {}".format(N)
            else:
                # If nothing else works, the dataset is ignored
                # Not enough info to process it
                print("N unavailable and cannot be estimated.Ignored")
                continue
            command = "{} --sumstats {} --out {} --snp {} --merge-alleles {}  {} --chunksize {} ".format(MUNGE_PY_PATH, 
                                                                                                     rawdata,
                                                                                                     outfile,
                                                                                                     snpid,
                                                                                                     SNP_LIST_PATH,
                                                                                                     Ncommand,
                                                                                                     CHUNKSIZE)
            print("Executing: ", command)
            os.system(command)
    pass

def ldsc(config, sep, munge_folder, ldsc_folder, verbose=True, process_all=False) -> None:
    if(not os.path.exists(LDSC_FOLDER)): os.mkdir(ldsc_folder)
    config_data = parse_config(config, sep)
    inputs = [config_data[key]["munged-dataset"] for key in config_data.keys() 
              if config_data[key]["discarded"]]
    permutations = list(itertools.product(inputs, inputs))
    for permutation in tqdm(permutations):
        input_a:str = os.path.join(munge_folder, permutation[0]+".sumstats.gz")
        input_b:str = os.path.join(munge_folder, permutation[1]+".sumstats.gz")
        outfile:str = os.path.join(ldsc_folder, "{}@{}.txt".format(permutation[0], permutation[1]))
        outfile_alt:str = os.path.join(ldsc_folder, "{}@{}.txt".format(permutation[1], permutation[0]))
        # or os.path.exists(outfile_alt+".log")
        if (not process_all and (os.path.exists(outfile+".log") )): 
            continue
        command = "{} --rg {},{} --ref-ld-chr {} --w-ld-chr {}  --out {} > /dev/null".format(
                                    LDSC_PY_PATH,
                                    input_a,
                                    input_b,
                                    REF_LD_PATH,
                                    WILD_LD_PATH,
                                    outfile
                                )
        os.system(command)
    pass

def filter_by_results(pvalue:float, rg:float, se:float,
                      limpval:float, limrg:float, limse: float) -> bool:
    """
        If all three conditions are met, the correlation is good, and returns true

        if not, returns false
    """
    if (pvalue < limpval and math.abs(rg) < limrg and se < limse):
        return True
    return False

def HDL(config, sep, munge_folder, hdl_folder, verbose=True, process_all=True) -> None:
    if (not os.path.exists(hdl_folder)): os.mkdir(hdl_folder)
    config_data = parse_config(config, sep)
    inputs = [config_data[key]["munged-dataset"] for key in config_data.keys() 
              if config_data[key]["discarded"]]
    permutations = list(itertools.product(inputs, inputs))
    for permutation in tqdm(permutations):
        input_a:str = os.path.join(munge_folder, permutation[0]+".sumstats.gz")
        input_b:str = os.path.join(munge_folder, permutation[1]+".sumstats.gz")
        outfile:str = os.path.join(hdl_folder, "{}@{}.txt".format(permutation[0], permutation[1]))
        logfile:str = os.path.join(hdl_folder, "{}@{}.log".format(permutation[0], permutation[1]))
        outfile_alt:str = os.path.join(hdl_folder, "{}@{}.txt".format(permutation[1], permutation[0]))
        # or os.path.exists(outfile_alt+".log")
        if (not process_all and (os.path.exists(outfile+".log"))): 
            continue
        command = "{} gwas1.df={} gwas2.df={} LD.path={} GWAS.type={} output.file={} fill.missing.N={} > /dev/null".format(
                                          HDL_R_PATH,
                                          input_a,
                                          input_b,
                                          LD_PATH,
                                          GWAS_TYPE,
                                          outfile,
                                          0)
        os.system(command)
    pass
        
def parse_ldsc_results(config, sep,  ldsc_folder, result_folder, results_file_name) -> None:
    print(result_folder)
    config_data = parse_config(config, sep)
    inputs = [config_data[key]["munged-dataset"] for key in config_data.keys() 
              if config_data[key]["discarded"]]
    permutations = list(itertools.product(inputs, inputs))
    header = []
    results = []
    for permutation in tqdm(permutations):
        outfile:str = os.path.join(ldsc_folder, "{}@{}.txt.log".format(permutation[0], permutation[1]))
        with open(outfile, "rt") as fhand:
            for line in fhand:
                if not line.startswith("p1"): continue
                if not header: 
                    header = line.strip("\n").split()
                    print(header)
                result = next(fhand)
                results.append(result.strip("\n").split())
        with open(os.path.join(result_folder, results_file_name), "wt") as whand:
            whand.write(sep.join(header)+"\n")
            for result in results:
                whand.write(sep.join(result)+"\n")

def heatmap_ldsc(sep, results_path, outfile_prefix) -> None:
    #if not os.path.exists(results_path): os.mkdir(results_path)
    print(outfile_prefix)
    command = "{} {} {} {} {} {} {}".format(
        R_SCRIPT_PATH,
        HEATMAP_SCRIPT_PATH,
        results_path,
        sep,
        outfile_prefix+"heatmap",
        HEATMAP_WIDTH,
        HEATMAP_HEIGHT,
    )
    os.system(command)

def unify_results(config, sep, ldsc_folder, results_folder, ldsc_results, hdl_folder, verbose=True):
    parse_ldsc_results(config, sep, ldsc_folder, results_folder, ldsc_results)
    heatmap_ldsc(sep, os.path.join(results_folder, ldsc_results), IMAGE_PATH)

def HR(config, sep, rawdata_path, mr_results):
    if (not os.path.exists(mr_results)): os.mkdir(mr_results)
    config_data = parse_config(config, sep)
    inputs = [config_data[key]["rawdata"] for key in config_data.keys() 
              if config_data[key]["discarded"]]
    variants_ids = [config_data[key]["variant-id"] for key in config_data.keys() 
              if config_data[key]["discarded"]]
    print(len(inputs))
    sumstats_combinations = list(itertools.product(zip(inputs, variants_ids), zip(inputs,variants_ids)))

    for sumstats1, sumstats2 in sumstats_combinations:
        outfile = "{}/{}@{}.mr".format(mr_results, sumstats1[0], sumstats2[0])
        if os.path.exists(outfile): continue
        command = "{} {} {} {} {} {} > {}".format(R_SCRIPT_PATH, MR_SCRIPT_PATH, rawdata_path+sumstats1[0], 
                                                   rawdata_path+sumstats2[0], sumstats1[1], sumstats2[1], outfile)
        print(command)
        os.system(command)
    pass

if __name__ == "__main__":
    if "download" in MODES:
        download_datasets(CONFIG_FILE, CONFIG_SEP, RAWDATA_FOLDER, verbose=True, process_all=PROCESS_ALL)
    if "munge" in MODES:
        munge_datasets(CONFIG_FILE, CONFIG_SEP, RAWDATA_FOLDER, MUNGE_FOLDER, verbose=True, process_all=PROCESS_ALL)
    if "ldsc" in MODES:
        ldsc(CONFIG_FILE, CONFIG_SEP, MUNGE_FOLDER, LDSC_FOLDER, verbose=True, process_all=PROCESS_ALL)
    if "HDL" in MODES:
        HDL(CONFIG_FILE, CONFIG_SEP, MUNGE_FOLDER, HDL_FOLDER, verbose=True, process_all=True)
    if "unify-results" in MODES:
        unify_results(CONFIG_FILE, CONFIG_SEP, LDSC_FOLDER, LDSC_RESULTS, LDSC_RESULTS_FILE, HDL_FOLDER)
    if "HR" in MODES:
        HR(CONFIG_FILE, CONFIG_SEP, RAWDATA_FOLDER, MR_results)