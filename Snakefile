from pathlib import Path
import snakemake.utils  # Load snakemake API
import os
import datetime
import re
import glob
import zipfile
import json
import numpy

__author__ = "Marine AGLAVE"


#using: snakemake -j5 --configfile /home/m_aglave/Bureau/scRNAseq_10X/Param_snakfile_alignment.yaml --use-conda --printshellcmds -s ./Snakefile



### parameters ###################################################################################################################################

#### Pipeline ####
STEPS = config['Steps']
PIPELINE_FOLDER = workflow.snakefile
PIPELINE_FOLDER = PIPELINE_FOLDER.replace("/Snakefile", "")
#CONDA_ENV_SING =  PIPELINE_FOLDER + "/envs/conda/Singularity.yml"

if "Alignment_countTable_GE" in STEPS:
    ### Sample/Project
    if 'Alignment_countTable_GE' in config and 'sample.name.ge' not in config['Alignment_countTable_GE']: sys.exit("Error: No sample.name.ge in configfile (Alignment_countTable_GE)!")
    if 'Alignment_countTable_GE' in config and 'input.dir.ge' not in config['Alignment_countTable_GE']: sys.exit("Error: No input.dir.ge in configfile (Alignment_countTable_GE)!")
    if 'Alignment_countTable_GE' in config and 'output.dir.ge' not in config['Alignment_countTable_GE']: sys.exit("Error: No output.dir.ge in configfile (Alignment_countTable_GE)!")
    ALIGN_SAMPLE_NAME_GE_RAW = config['Alignment_countTable_GE']['sample.name.ge']
    ALIGN_INPUT_DIR_GE_RAW = os.path.normpath(config['Alignment_countTable_GE']['input.dir.ge'])
    ALIGN_OUTPUT_DIR_GE = os.path.normpath(config['Alignment_countTable_GE']['output.dir.ge'])
    ALIGN_INPUT_DIR_GE = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/fastq/")
    ### Index
    KINDEX_GE = config['Alignment_countTable_GE']['kindex.ge'] if 'Alignment_countTable_GE' in config and 'kindex.ge' in config['Alignment_countTable_GE'] else sys.exit("Error: No kindex.ge in configfile (Alignment_countTable_GE)!")
    TR2GFILE_GE = config['Alignment_countTable_GE']['tr2g.file.ge'] if 'Alignment_countTable_GE' in config and 'tr2g.file.ge' in config['Alignment_countTable_GE'] else sys.exit("Error: No tr2g.file.ge in configfile (Alignment_countTable_GE)!")
    REF_TXT_GE = config['Alignment_countTable_GE']['reference.txt'] if 'reference.txt' in config['Alignment_countTable_GE'] else "<insert_you_reference_here>"
    ### File names
    ALIGN_SAMPLE_NAME_GE = []
    ALIGN_SYMLINK_FILES_GE = []
    ALIGN_SYMLINK_FILES_NAME_GE = []
    for i in range(0,len(ALIGN_SAMPLE_NAME_GE_RAW),1):
        #check samples names and add "_GE" if needed
        ALIGN_SAMPLE_NAME_GE.append(ALIGN_SAMPLE_NAME_GE_RAW[i] + "_GE") if (ALIGN_SAMPLE_NAME_GE_RAW[i][len(ALIGN_SAMPLE_NAME_GE_RAW[i])-3:] != "_GE") else ALIGN_SAMPLE_NAME_GE.append(ALIGN_SAMPLE_NAME_GE_RAW[i])
        ORIG_FILES = glob.glob(os.path.join(ALIGN_INPUT_DIR_GE_RAW, str(ALIGN_SAMPLE_NAME_GE_RAW[i]) + "*_R1_*.f*q*")) + glob.glob(os.path.join(ALIGN_INPUT_DIR_GE_RAW, str(ALIGN_SAMPLE_NAME_GE_RAW[i]) + "*_R2_*.f*q*"))
        #files with path and extention
        ALIGN_SYMLINK_FILES_GE = ALIGN_SYMLINK_FILES_GE + [ os.path.normpath(ALIGN_INPUT_DIR_GE + "/" + os.path.basename(file).replace(ALIGN_SAMPLE_NAME_GE_RAW[i], ALIGN_SAMPLE_NAME_GE[i])) for file in ORIG_FILES]
    #files without path and extention
    ALIGN_SYMLINK_FILES_NAME_GE = [os.path.splitext(os.path.splitext(os.path.basename(x))[0])[0] for x in ALIGN_SYMLINK_FILES_GE]

if "Alignment_countTable_ADT" in STEPS:
    ### Sample/Project
    if 'Alignment_countTable_ADT' in config and 'sample.name.adt' not in config['Alignment_countTable_ADT']: sys.exit("Error: No sample.name.adt in configfile (Alignment_countTable_ADT)!")
    if 'Alignment_countTable_ADT' in config and 'input.dir.adt' not in config['Alignment_countTable_ADT']: sys.exit("Error: No input.dir.adt in configfile (Alignment_countTable_ADT)!")
    if 'Alignment_countTable_ADT' in config and 'output.dir.adt' not in config['Alignment_countTable_ADT']: sys.exit("Error: No output.dir.adt in configfile (Alignment_countTable_ADT)!")
    ALIGN_SAMPLE_NAME_ADT_RAW = config['Alignment_countTable_ADT']['sample.name.adt']
    ALIGN_INPUT_DIR_ADT_RAW = os.path.normpath(config['Alignment_countTable_ADT']['input.dir.adt'])
    ALIGN_OUTPUT_DIR_ADT = os.path.normpath(config['Alignment_countTable_ADT']['output.dir.adt'])
    ALIGN_INPUT_DIR_ADT = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/fastq/")
    ### Index
    KINDEX_ADT = config['Alignment_countTable_ADT']['kindex.adt'] if 'Alignment_countTable_ADT' in config and 'kindex.adt' in config['Alignment_countTable_ADT'] else sys.exit("Error: No kindex.adt in configfile (Alignment_countTable_ADT)!")
    TR2GFILE_ADT = config['Alignment_countTable_ADT']['tr2g.file.adt'] if 'Alignment_countTable_ADT' in config and 'tr2g.file.adt' in config['Alignment_countTable_ADT'] else sys.exit("Error: No tr2g.file.adt in configfile (Alignment_countTable_ADT)!")
    ### File names
    ALIGN_SAMPLE_NAME_ADT = []
    ALIGN_SYMLINK_FILES_ADT = []
    ALIGN_SYMLINK_FILES_NAME_ADT = []
    for i in range(0,len(ALIGN_SAMPLE_NAME_ADT_RAW),1):
        #check samples names and add "_ADT" if needed
        ALIGN_SAMPLE_NAME_ADT.append(ALIGN_SAMPLE_NAME_ADT_RAW[i] + "_ADT") if (ALIGN_SAMPLE_NAME_ADT_RAW[i][len(ALIGN_SAMPLE_NAME_ADT_RAW[i])-4:] != "_ADT") else ALIGN_SAMPLE_NAME_ADT.append(ALIGN_SAMPLE_NAME_ADT_RAW[i])
        ORIG_FILES = glob.glob(os.path.join(ALIGN_INPUT_DIR_ADT_RAW, str(ALIGN_SAMPLE_NAME_ADT_RAW[i]) + "*_R1_*.f*q*")) + glob.glob(os.path.join(ALIGN_INPUT_DIR_ADT_RAW, str(ALIGN_SAMPLE_NAME_ADT_RAW[i]) + "*_R2_*.f*q*"))
        #files with path and extention
        ALIGN_SYMLINK_FILES_ADT = ALIGN_SYMLINK_FILES_ADT + [ os.path.normpath(ALIGN_INPUT_DIR_ADT + "/" + os.path.basename(file).replace(ALIGN_SAMPLE_NAME_ADT_RAW[i], ALIGN_SAMPLE_NAME_ADT[i])) for file in ORIG_FILES]
    #files without path and extention
    ALIGN_SYMLINK_FILES_NAME_ADT = [os.path.splitext(os.path.splitext(os.path.basename(x))[0])[0] for x in ALIGN_SYMLINK_FILES_ADT]

if "Alignment_annotations_TCR_BCR" in STEPS:
    ### Sample/Project
    if 'Alignment_annotations_TCR_BCR' in config and 'sample.name.tcr' not in config['Alignment_annotations_TCR_BCR'] and 'sample.name.bcr' not in config['Alignment_annotations_TCR_BCR']: sys.exit("Error: No sample.name.tcr or sample.name.bcr in configfile (Alignment_annotations_TCR_BCR)!")
    if 'Alignment_annotations_TCR_BCR' in config and 'input.dir.tcr' not in config['Alignment_annotations_TCR_BCR'] and 'input.dir.bcr' not in config['Alignment_annotations_TCR_BCR']: sys.exit("Error: No input.dir.tcr or input.dir.bcr in configfile (Alignment_annotations_TCR_BCR)!")
    if 'Alignment_annotations_TCR_BCR' in config and 'output.dir.tcr_bcr' not in config['Alignment_annotations_TCR_BCR']: sys.exit("Error: No output.dir.tcr_bcr in configfile (Alignment_annotations_TCR_BCR)!")
    ALIGN_SAMPLE_NAME_TCR_RAW = config['Alignment_annotations_TCR_BCR']['sample.name.tcr'] if 'sample.name.tcr' in config['Alignment_annotations_TCR_BCR'] else None
    ALIGN_SAMPLE_NAME_BCR_RAW = config['Alignment_annotations_TCR_BCR']['sample.name.bcr'] if 'sample.name.bcr' in config['Alignment_annotations_TCR_BCR'] else None
    ALIGN_INPUT_DIR_TCR_RAW = os.path.normpath(config['Alignment_annotations_TCR_BCR']['input.dir.tcr'] + "/") if 'input.dir.tcr' in config['Alignment_annotations_TCR_BCR'] else None
    ALIGN_INPUT_DIR_BCR_RAW = os.path.normpath(config['Alignment_annotations_TCR_BCR']['input.dir.bcr'] + "/") if 'input.dir.bcr' in config['Alignment_annotations_TCR_BCR'] else None
    ALIGN_OUTPUT_DIR_TCR_BCR = os.path.normpath(config['Alignment_annotations_TCR_BCR']['output.dir.tcr_bcr'])
    ALIGN_INPUT_DIR_TCR_BCR = os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/fastq/")
    ### Index
    CRINDEX_TCR_BCR=config['Alignment_annotations_TCR_BCR']['crindex.tcr_bcr'] if ('Alignment_annotations_TCR_BCR' in config and 'crindex.tcr_bcr' in config['Alignment_annotations_TCR_BCR']) else "/mnt/beegfs/database/bioinfo/single-cell/TCR_REFERENCES/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0"
    ### File names
    #TCR
    ALIGN_SAMPLE_NAME_TCR = []
    ALIGN_ORIG_FILES_TCR = []
    ALIGN_SYMLINK_FILES_TCR = []
    ALIGN_SYMLINK_FILES_NAME_TCR = []
    if ALIGN_SAMPLE_NAME_TCR_RAW is not None:
        for i in range(0,len(ALIGN_SAMPLE_NAME_TCR_RAW),1):
            #check samples names and add "_TCR" if needed
            ALIGN_SAMPLE_NAME_TCR.append(ALIGN_SAMPLE_NAME_TCR_RAW[i] + "_TCR") if (ALIGN_SAMPLE_NAME_TCR_RAW[i][len(ALIGN_SAMPLE_NAME_TCR_RAW[i])-4:] != "_TCR") else ALIGN_SAMPLE_NAME_TCR.append(ALIGN_SAMPLE_NAME_TCR_RAW[i])
            ORIG_FILES = glob.glob(os.path.join(ALIGN_INPUT_DIR_TCR_RAW, str(ALIGN_SAMPLE_NAME_TCR_RAW[i]) + "*_R1_*.f*q*")) + glob.glob(os.path.join(ALIGN_INPUT_DIR_TCR_RAW, str(ALIGN_SAMPLE_NAME_TCR_RAW[i]) + "*_R2_*.f*q*"))
             #files with path and extention
            ALIGN_ORIG_FILES_TCR = ALIGN_ORIG_FILES_TCR + ORIG_FILES
            for file in ORIG_FILES:
                if re.match(str(ALIGN_SAMPLE_NAME_TCR_RAW[i] + "_S[0-9]+_L00[0-9]{1}_R[1-2]{1}_.*"), os.path.basename(file)) is not None: #good name format
                    ALIGN_SYMLINK_FILES_TCR = ALIGN_SYMLINK_FILES_TCR + [ os.path.normpath(ALIGN_INPUT_DIR_TCR_BCR + "/" + os.path.basename(file).replace(ALIGN_SAMPLE_NAME_TCR_RAW[i], ALIGN_SAMPLE_NAME_TCR[i])) ]
                elif re.match(str(ALIGN_SAMPLE_NAME_TCR_RAW[i] + "_[0-9]{1}_S[0-9]+_R[1-2]{1}_.*"), os.path.basename(file)) is not None: # => reformat
                    res_match = re.match(str(ALIGN_SAMPLE_NAME_TCR_RAW[i] + "_(?P<nb>[0-9]{1})_(?P<S>S[0-9]+)_(?P<R_compl>R[1-2]{1}_.*)"), os.path.basename(file))
                    ALIGN_SYMLINK_FILES_TCR = ALIGN_SYMLINK_FILES_TCR + [ os.path.normpath(str(ALIGN_INPUT_DIR_TCR_BCR + "/" + ALIGN_SAMPLE_NAME_TCR[i] + "_" + res_match.group('S') + "_L00" + res_match.group('nb') + "_" + res_match.group('R_compl'))) ]
                else:
                    sys.exit("File names for TCR not recognized. It must be like mysample_2_S1_R1_001.fastq.gz or mysample_S1_L002_R1_001.fastq.gz")
        #files without path and extention
        ALIGN_SYMLINK_FILES_NAME_TCR = [os.path.splitext(os.path.splitext(os.path.basename(x))[0])[0] for x in ALIGN_SYMLINK_FILES_TCR]
    #BCR
    ALIGN_SAMPLE_NAME_BCR = []
    ALIGN_ORIG_FILES_BCR = []
    ALIGN_SYMLINK_FILES_BCR = []
    ALIGN_SYMLINK_FILES_NAME_BCR = []
    if ALIGN_SAMPLE_NAME_BCR_RAW is not None:
        for i in range(0,len(ALIGN_SAMPLE_NAME_BCR_RAW),1):
            #check samples names and add "_BCR" if needed
            ALIGN_SAMPLE_NAME_BCR.append(ALIGN_SAMPLE_NAME_BCR_RAW[i] + "_BCR") if (ALIGN_SAMPLE_NAME_BCR_RAW[i][len(ALIGN_SAMPLE_NAME_BCR_RAW[i])-4:] != "_BCR") else ALIGN_SAMPLE_NAME_BCR.append(ALIGN_SAMPLE_NAME_BCR_RAW[i])
            ORIG_FILES = glob.glob(os.path.join(ALIGN_INPUT_DIR_BCR_RAW, str(ALIGN_SAMPLE_NAME_BCR_RAW[i]) + "*_R1_*.f*q*")) + glob.glob(os.path.join(ALIGN_INPUT_DIR_BCR_RAW, str(ALIGN_SAMPLE_NAME_BCR_RAW[i]) + "*_R2_*.f*q*"))
            #files with path and extention
            ALIGN_ORIG_FILES_BCR = ALIGN_ORIG_FILES_BCR + ORIG_FILES
            for file in ORIG_FILES:
                if re.match(str(ALIGN_SAMPLE_NAME_BCR_RAW[i] + "_S[0-9]+_L00[0-9]{1}_R[1-2]{1}_.*"), os.path.basename(file)) is not None: #good name format
                    ALIGN_SYMLINK_FILES_BCR = ALIGN_SYMLINK_FILES_BCR + [ os.path.normpath(ALIGN_INPUT_DIR_TCR_BCR + "/" + os.path.basename(file).replace(ALIGN_SAMPLE_NAME_BCR_RAW[i], ALIGN_SAMPLE_NAME_BCR[i])) ]
                elif re.match(str(ALIGN_SAMPLE_NAME_BCR_RAW[i] + "_[0-9]{1}_S[0-9]+_R[1-2]{1}_.*"), os.path.basename(file)) is not None: # => reformat
                    res_match = re.match(str(ALIGN_SAMPLE_NAME_BCR_RAW[i] + "_(?P<nb>[0-9]{1})_(?P<S>S[0-9]+)_(?P<R_compl>R[1-2]{1}_.*)"), os.path.basename(file))
                    ALIGN_SYMLINK_FILES_BCR = ALIGN_SYMLINK_FILES_BCR + [ os.path.normpath(str(ALIGN_INPUT_DIR_TCR_BCR + "/" + ALIGN_SAMPLE_NAME_BCR[i] + "_" + res_match.group('S') + "_L00" + res_match.group('nb') + "_" + res_match.group('R_compl'))) ]
                else:
                    sys.exit("File names for BCR not recognized. It must be like mysample_2_S1_R1_001.fastq.gz or mysample_S1_L002_R1_001.fastq.gz")
        #files without path and extention
        ALIGN_SYMLINK_FILES_NAME_BCR = [os.path.splitext(os.path.splitext(os.path.basename(x))[0])[0] for x in ALIGN_SYMLINK_FILES_BCR]
    #Fusion TCR/BCR
    ALIGN_SAMPLE_NAME_TCR_BCR_RAW = ALIGN_SAMPLE_NAME_TCR_RAW  + ALIGN_SAMPLE_NAME_BCR_RAW
    ALIGN_SAMPLE_NAME_TCR_BCR = ALIGN_SAMPLE_NAME_TCR + ALIGN_SAMPLE_NAME_BCR
    ALIGN_ORIG_FILES_TCR_BCR = ALIGN_ORIG_FILES_TCR + ALIGN_ORIG_FILES_BCR
    ALIGN_SYMLINK_FILES_TCR_BCR = ALIGN_SYMLINK_FILES_TCR + ALIGN_SYMLINK_FILES_BCR
    ALIGN_SYMLINK_FILES_NAME_TCR_BCR = ALIGN_SYMLINK_FILES_NAME_TCR + ALIGN_SYMLINK_FILES_NAME_BCR

if "Alignment_countTable_GE" in STEPS or "Alignment_countTable_ADT" in STEPS:
    # 10X Technology
    if 'Alignment_countTable_GE' in config and 'sctech' in config['Alignment_countTable_GE']:
        SCTECH = config['Alignment_countTable_GE']['sctech']
    elif 'Alignment_countTable_ADT' in config and 'sctech' in config['Alignment_countTable_ADT']:
        SCTECH = config['Alignment_countTable_ADT']['sctech']
    else:
        SCTECH = '10xv3' # '10xv2' '10xv3'
    if SCTECH == '10xv3' :
        WHITELISTNAME = PIPELINE_FOLDER + '/resources/WHITELISTS/3M-february-2018.txt' # '737K-august-2016.txt' '3M-february-2018.txt'
    elif SCTECH == '10xv2' :
        WHITELISTNAME = PIPELINE_FOLDER + '/resources/WHITELISTS/737K-august-2016.txt'
    else :
        sys.exit("Error: sctech doesn't exist! Only '10xv2' and '10xv3' are available.\n")

if "Alignment_countTable_GE" in STEPS or "Alignment_annotations_TCR_BCR" in STEPS:
    # Fastq-screen Index
    if 'Alignment_countTable_GE' in config and 'fastqscreen_index' in config['Alignment_countTable_GE']:
        FASTQSCREEN_INDEX = config['Alignment_countTable_GE']['fastqscreen_index']
    elif 'Alignment_annotations_TCR_BCR' in config and 'fastqscreen_index' in config['Alignment_annotations_TCR_BCR']:
        FASTQSCREEN_INDEX = config['Alignment_annotations_TCR_BCR']['fastqscreen_index']
    else :
        FASTQSCREEN_INDEX = "/mnt/beegfs/database/bioinfo/single-cell/INDEX/FASTQ_SCREEN/0.14.0/fastq_screen.conf"

if "Alignment_countTable_GE" in STEPS or "Alignment_countTable_ADT" in STEPS or "Alignment_annotations_TCR_BCR" in STEPS:
    # Cutadapt parameters
    ADAPTERSEQ='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
    MINBASEQ=28
    # Name of conda environment
    CONDA_ENV_QC_ALIGN_GE_ADT = PIPELINE_FOLDER + "/envs/conda/QC_Alignment.yml"

if "Droplets_QC_GE" in STEPS:
    ### Sample/Project
    if 'Droplets_QC_GE' in config and 'sample.name.ge' in config['Droplets_QC_GE'] and 'input.dir.ge' in config['Droplets_QC_GE']:
        QC_SAMPLE_NAME_GE_RAW = config['Droplets_QC_GE']['sample.name.ge']
        QC_INPUT_DIR_GE = config['Droplets_QC_GE']['input.dir.ge']
        #check samples names and add "_GE" if needed
        QC_SAMPLE_NAME_GE = []
        for i in range(0,len(QC_SAMPLE_NAME_GE_RAW),1):
            QC_SAMPLE_NAME_GE.append(QC_SAMPLE_NAME_GE_RAW[i] + "_GE") if (QC_SAMPLE_NAME_GE_RAW[i][len(QC_SAMPLE_NAME_GE_RAW[i])-3:] != "_GE") else QC_SAMPLE_NAME_GE.append(QC_SAMPLE_NAME_GE_RAW[i])
    elif 'sample.name.ge' in config['Alignment_countTable_GE'] and 'input.dir.ge' in config['Alignment_countTable_GE']  and "Alignment_countTable_GE" in STEPS:
        sys.stderr.write("Note: No sample.name.ge or input.dir.ge find in Droplets_QC_GE section of configfile; sample.name.ge and input.dir.ge will be determine from Alignment_countTable_GE step for Droplets_QC_GE step!\n")
        QC_SAMPLE_NAME_GE = copy.deepcopy(ALIGN_SAMPLE_NAME_GE)
        QC_INPUT_DIR_GE = [os.path.join(ALIGN_OUTPUT_DIR_GE, str(x), "KALLISTOBUS") for x in ALIGN_SAMPLE_NAME_GE]
    else:
        sys.exit("Error: No sample.name.ge or/and input.dir.ge in configfile!\n")
    if 'Droplets_QC_GE' in config and 'output.dir.ge' in config['Droplets_QC_GE'] :
        QC_OUTPUT_DIR_GE = config['Droplets_QC_GE']['output.dir.ge']
    elif 'output.dir.ge' in config['Alignment_countTable_GE'] :
        QC_OUTPUT_DIR_GE = [os.path.join(ALIGN_OUTPUT_DIR_GE, str(x)) for x in ALIGN_SAMPLE_NAME_GE]
        sys.stderr.write("Note: No output.dir.ge find in Droplets_QC_GE section of configfile; output.dir.ge will be determine from Alignment_countTable_GE step for Droplets_QC_GE step!\n")
    else :
        sys.exit("Error: No output.dir.ge find in configfile!\n")
    QC_SPECIES = config['Droplets_QC_GE']['species'] if ('Droplets_QC_GE' in config and 'species' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['species'] != None) else "NULL"
    QC_AUTHOR_NAME = config['Droplets_QC_GE']['author.name'].replace(", ", ",").replace(" ", "_") if ('Droplets_QC_GE' in config and 'author.name' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['author.name'] != None) else "NULL"
    QC_AUTHOR_MAIL = config['Droplets_QC_GE']['author.mail'].replace(", ", ",") if ('Droplets_QC_GE' in config and 'author.mail' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['author.mail'] != None) else "NULL"
    ### Analysis Parameters
    # Emptydrops
    QC_EMPTYDROPS_FDR = config['Droplets_QC_GE']['emptydrops.fdr'] if ('Droplets_QC_GE' in config and 'emptydrops.fdr' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['emptydrops.fdr'] != None) else "NULL"
    QC_DROPLETS_LIMIT = config['Droplets_QC_GE']['droplets.limit'] if ('Droplets_QC_GE' in config and 'droplets.limit' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['droplets.limit'] != None) else "NULL"
    QC_EMPTYDROPS_RETAIN = config['Droplets_QC_GE']['emptydrops.retain'] if ('Droplets_QC_GE' in config and 'emptydrops.retain' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['emptydrops.retain'] != None) else "NULL"
    # Translate ENSG into Gene Symbol
    QC_TRANSLATION_BOOL = config['Droplets_QC_GE']['translation'] if ('Droplets_QC_GE' in config and 'translation' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['translation'] != None) else "NULL"
    # QC cell
    QC_PCMITO_MIN = config['Droplets_QC_GE']['pcmito.min'] if ('Droplets_QC_GE' in config and 'pcmito.min' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['pcmito.min'] != None) else "NULL"
    QC_PCMITO_MAX = config['Droplets_QC_GE']['pcmito.max'] if ('Droplets_QC_GE' in config and 'pcmito.max' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['pcmito.max'] != None) else "NULL"
    QC_PCRIBO_MIN = config['Droplets_QC_GE']['pcribo.min'] if ('Droplets_QC_GE' in config and 'pcribo.min' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['pcribo.min'] != None) else "NULL"
    QC_PC_RIBO_MAX = config['Droplets_QC_GE']['pcribo.max'] if ('Droplets_QC_GE' in config and 'pcribo.max' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['pcribo.max'] != None) else "NULL"
    QC_MIN_FEATURES = config['Droplets_QC_GE']['min.features'] if ('Droplets_QC_GE' in config and 'min.features' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['min.features'] != None) else "NULL"
    QC_MIN_COUNTS = config['Droplets_QC_GE']['min.counts'] if ('Droplets_QC_GE' in config and 'min.counts' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['min.counts'] != None) else "NULL"
    # QC gene
    QC_MIN_CELLS = config['Droplets_QC_GE']['min.cells'] if ('Droplets_QC_GE' in config and 'min.cells' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['min.cells'] != None) else "NULL"
    ### Databases
    # QC
    QC_MT_FILE = config['Droplets_QC_GE']['mt.genes.file'] if ('Droplets_QC_GE' in config and 'mt.genes.file' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['mt.genes.file'] != None) else "NULL"
    QC_RB_FILE = config['Droplets_QC_GE']['crb.genes.file'] if ('Droplets_QC_GE' in config and 'crb.genes.file' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['crb.genes.file'] != None) else "NULL"
    QC_ST_FILE = config['Droplets_QC_GE']['str.genes.file'] if ('Droplets_QC_GE' in config and 'str.genes.file' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['str.genes.file'] != None) else "NULL"
    # Translation into gene Symbols
    QC_TRANSLATION_FILE = config['Droplets_QC_GE']['translation.file'] if ('Droplets_QC_GE' in config and 'translation.file' in config['Droplets_QC_GE'] and config['Droplets_QC_GE']['translation.file'] != None) else "NULL"
    ### Snakefile parameters
    #check end paths (del "/" if necessary)
    for i in range(0,len(QC_INPUT_DIR_GE),1):
        QC_INPUT_DIR_GE[i] = os.path.normpath(QC_INPUT_DIR_GE[i])
        QC_OUTPUT_DIR_GE[i] = os.path.normpath(QC_OUTPUT_DIR_GE[i])
    #Correspondance sample/input/output
    dic_SAMPLE_NAME_GE_INFO = {}
    for i in range(0,len(QC_SAMPLE_NAME_GE),1):
        dic_SAMPLE_NAME_GE_INFO[QC_SAMPLE_NAME_GE[i]] = {}
        dic_SAMPLE_NAME_GE_INFO[QC_SAMPLE_NAME_GE[i]]['QC_INPUT_DIR'] = QC_INPUT_DIR_GE[i]
        dic_SAMPLE_NAME_GE_INFO[QC_SAMPLE_NAME_GE[i]]['QC_OUTPUT_DIR'] = QC_OUTPUT_DIR_GE[i]

if "Filtering_GE" in STEPS:
    ### Sample/Project
    if 'Filtering_GE' in config and 'sample.name.ge' in config['Filtering_GE'] and 'input.rda.ge' in config['Filtering_GE'] :
        FILERING_SAMPLE_NAME_GE_RAW = config['Filtering_GE']['sample.name.ge']
        FILERING_INPUT_RDA_GE = config['Filtering_GE']['input.rda.ge']
        #check samples names and add "_GE" if needed
        FILERING_SAMPLE_NAME_GE = []
        for i in range(0,len(FILERING_SAMPLE_NAME_GE_RAW),1):
            FILERING_SAMPLE_NAME_GE.append(FILERING_SAMPLE_NAME_GE_RAW[i] + "_GE") if (FILERING_SAMPLE_NAME_GE_RAW[i][len(FILERING_SAMPLE_NAME_GE_RAW[i])-3:] != "_GE") else FILERING_SAMPLE_NAME_GE.append(FILERING_SAMPLE_NAME_GE_RAW[i])
    elif "Droplets_QC_GE" in STEPS:
        sys.stderr.write("Note: No input.rda.ge find in Filtering_GE section of configfile; input.rda.ge will be determine from Droplets_QC_GE step for Filtering_GE step!\n")
        FILERING_SAMPLE_NAME_GE = copy.deepcopy(QC_SAMPLE_NAME_GE)
        FILERING_INPUT_RDA_GE = [os.path.normpath(dic_SAMPLE_NAME_GE_INFO[x]['QC_OUTPUT_DIR'] + ("/QC_droplets/" if str(QC_EMPTYDROPS_RETAIN) == "NULL" else "/QC_droplets_retain" + str(QC_EMPTYDROPS_RETAIN) + "/") + str(x) + "_QC_NON-NORMALIZED.rda") for x in QC_SAMPLE_NAME_GE]
    else:
        sys.exit("Error: No sample.name.ge or/and input.rda.ge in configfile!\n")
    if 'Filtering_GE' in config and 'output.dir.ge' in config['Filtering_GE'] :
        FILERING_OUTPUT_DIR_GE = config['Filtering_GE']['output.dir.ge']
    elif "Droplets_QC_GE" in STEPS:
        FILERING_OUTPUT_DIR_GE = copy.deepcopy(QC_OUTPUT_DIR_GE)
        sys.stderr.write("Note: No output.dir.ge find in Filtering_GE section of configfile; output.dir.ge will be determine from Droplets_QC_GE step for Filtering_GE step!\n")
    else :
        sys.exit("Error: No output.dir.ge find in configfile!\n")
    FILERING_AUTHOR_NAME = config['Filtering_GE']['author.name'].replace(", ", ",").replace(" ", "_") if ('Filtering_GE' in config and 'author.name' in config['Filtering_GE'] and config['Filtering_GE']['author.name'] != None) else "NULL"
    FILERING_AUTHOR_MAIL = config['Filtering_GE']['author.mail'].replace(", ", ",") if ('Filtering_GE' in config and 'author.mail' in config['Filtering_GE'] and config['Filtering_GE']['author.mail'] != None) else "NULL"
    ### Analysis Parameters
    # QC cell
    FILERING_PCMITO_MIN = config['Filtering_GE']['pcmito.min'] if ('Filtering_GE' in config and 'pcmito.min' in config['Filtering_GE'] and config['Filtering_GE']['pcmito.min'] != None) else "0"
    FILERING_PCMITO_MAX = config['Filtering_GE']['pcmito.max'] if ('Filtering_GE' in config and 'pcmito.max' in config['Filtering_GE'] and config['Filtering_GE']['pcmito.max'] != None) else "0.2"
    FILERING_PCRIBO_MIN = config['Filtering_GE']['pcribo.min'] if ('Filtering_GE' in config and 'pcribo.min' in config['Filtering_GE'] and config['Filtering_GE']['pcribo.min'] != None) else "0"
    FILERING_PC_RIBO_MAX = config['Filtering_GE']['pcribo.max'] if ('Filtering_GE' in config and 'pcribo.max' in config['Filtering_GE'] and config['Filtering_GE']['pcribo.max'] != None) else "1"
    FILERING_MIN_FEATURES = config['Filtering_GE']['min.features'] if ('Filtering_GE' in config and 'min.features' in config['Filtering_GE'] and config['Filtering_GE']['min.features'] != None) else "200"
    FILERING_MIN_COUNTS = config['Filtering_GE']['min.counts'] if ('Filtering_GE' in config and 'min.counts' in config['Filtering_GE'] and config['Filtering_GE']['min.counts'] != None) else "1000"
    # QC gene
    FILERING_MIN_CELLS = config['Filtering_GE']['min.cells'] if ('Filtering_GE' in config and 'min.cells' in config['Filtering_GE'] and config['Filtering_GE']['min.cells'] != None) else "5"
    # Doublets
    FILERING_DOUBLET_FILTER_METHOD= config['Filtering_GE']['doublets.filter.method'] if ('Filtering_GE' in config and 'doublets.filter.method' in config['Filtering_GE'] and config['Filtering_GE']['doublets.filter.method'] != None) else "NULL"
    ### Databases
    # QC
    FILERING_CC_SEURAT_FILE = config['Filtering_GE']['cc.seurat.file'] if ('Filtering_GE' in config and 'cc.seurat.file' in config['Filtering_GE'] and config['Filtering_GE']['cc.seurat.file'] != None) else "NULL"
    FILERING_CC_CYCLONE_FILE = config['Filtering_GE']['cc.cyclone.file'] if ('Filtering_GE' in config and 'cc.cyclone.file' in config['Filtering_GE'] and config['Filtering_GE']['cc.cyclone.file'] != None) else "NULL"
    ### Snakefile parameters
    #Correspondance sample/input/output
    dic_FILTER_INFO = {}
    for i in range(0,len(FILERING_SAMPLE_NAME_GE),1):
        dic_FILTER_INFO[FILERING_SAMPLE_NAME_GE[i]] = {}
        dic_FILTER_INFO[FILERING_SAMPLE_NAME_GE[i]]['FILTER_INPUT_RDA'] = FILERING_INPUT_RDA_GE[i]
        dic_FILTER_INFO[FILERING_SAMPLE_NAME_GE[i]]['FILTER_OUTPUT_DIR'] = FILERING_OUTPUT_DIR_GE[i]
    FILTERS_FOLDER = "F" + str(FILERING_MIN_FEATURES) + "_C" + str(FILERING_MIN_COUNTS) + "_M" + str(FILERING_PCMITO_MIN) + "-" + str(FILERING_PCMITO_MAX) + "_R" + str(FILERING_PCRIBO_MIN) + "-" + str(FILERING_PC_RIBO_MAX) + "_G" + str(FILERING_MIN_CELLS)
    #name of the doublets identification method
    FILERING_DOUBLET_FILTER_METHOD_NAME = "all" if FILERING_DOUBLET_FILTER_METHOD == "NULL" else FILERING_DOUBLET_FILTER_METHOD

if "Norm_DimRed_Eval_GE" in STEPS: #alias NDRE_
    ### Sample/Project
    if ('Norm_DimRed_Eval_GE' in config) and ('sample.name.ge' in config['Norm_DimRed_Eval_GE']) and ('input.rda.ge' in config['Norm_DimRed_Eval_GE']) :
        NDRE_SAMPLE_NAME_GE_RAW = config['Norm_DimRed_Eval_GE']['sample.name.ge']
        NDRE_INPUT_RDA_GE = config['Norm_DimRed_Eval_GE']['input.rda.ge']
        #check samples names and add "_GE" if needed
        NDRE_SAMPLE_NAME_GE = []
        for i in range(0,len(NDRE_SAMPLE_NAME_GE_RAW),1):
            NDRE_SAMPLE_NAME_GE.append(NDRE_SAMPLE_NAME_GE_RAW[i] + "_GE") if (NDRE_SAMPLE_NAME_GE_RAW[i][len(NDRE_SAMPLE_NAME_GE_RAW[i])-3:] != "_GE") else NDRE_SAMPLE_NAME_GE.append(NDRE_SAMPLE_NAME_GE_RAW[i])
    elif "Filtering_GE" in STEPS:
        sys.stderr.write("Note: No input.rda.ge find in Norm_DimRed_Eval_GE section of configfile; input.rda.ge will be determine from Filtering_GE step for Norm_DimRed_Eval_GE step!\n")
        NDRE_SAMPLE_NAME_GE = copy.deepcopy(FILERING_SAMPLE_NAME_GE)
        if FILERING_DOUBLET_FILTER_METHOD_NAME == "none":
            NDRE_INPUT_RDA_GE = [os.path.normpath(dic_FILTER_INFO[x]['FILTER_OUTPUT_DIR'] + "/" + FILTERS_FOLDER + "/DOUBLETSKEPT/" + x + "_DOUBLETSKEPT_NON-NORMALIZED.rda") for x in FILERING_SAMPLE_NAME_GE]
        else:
            NDRE_INPUT_RDA_GE = [os.path.normpath(dic_FILTER_INFO[x]['FILTER_OUTPUT_DIR'] + "/" + FILTERS_FOLDER + "/DOUBLETSFILTER_" + FILERING_DOUBLET_FILTER_METHOD_NAME + "/" + x + "_FILTERED_NON-NORMALIZED.rda") for x in FILERING_SAMPLE_NAME_GE]
    else:
        sys.exit("Error: No sample.name.ge or/and input.rda.ge in configfile!\n")
    if ('Norm_DimRed_Eval_GE' in config) and ('output.dir.ge' in config['Norm_DimRed_Eval_GE']) :
        NDRE_OUTPUT_DIR_GE = config['Norm_DimRed_Eval_GE']['output.dir.ge']
    elif "Filtering_GE" in STEPS:
        NDRE_OUTPUT_DIR_GE = [dic_FILTER_INFO[x]['FILTER_OUTPUT_DIR'] + "/" + FILTERS_FOLDER + ("/DOUBLETSKEPT" if FILERING_DOUBLET_FILTER_METHOD_NAME == "none" else ("/DOUBLETSFILTER_" + FILERING_DOUBLET_FILTER_METHOD_NAME)) for x in FILERING_SAMPLE_NAME_GE]
        sys.stderr.write("Note: No output.dir.ge find in Norm_DimRed_Eval_GE section of configfile; output.dir.ge will be determine from Filtering_GE step for Norm_DimRed_Eval_GE step!\n")
    else :
        sys.exit("Error: No output.dir.ge find in configfile!\n")
    ### Analysis Parameters
    NDRE_AUTHOR_NAME = config['Norm_DimRed_Eval_GE']['author.name'].replace(", ", ",").replace(" ", "_") if ('Norm_DimRed_Eval_GE' in config and 'author.name' in config['Norm_DimRed_Eval_GE'] and config['Norm_DimRed_Eval_GE']['author.name'] != None) else "NULL"
    NDRE_AUTHOR_MAIL = config['Norm_DimRed_Eval_GE']['author.mail'].replace(", ", ",") if ('Norm_DimRed_Eval_GE' in config and 'author.mail' in config['Norm_DimRed_Eval_GE'] and config['Norm_DimRed_Eval_GE']['author.mail'] != None) else "NULL"
    NDRE_EVAL_MARKERS = config['Norm_DimRed_Eval_GE']['eval.markers'].replace(", ", ",") if ('Norm_DimRed_Eval_GE' in config and 'eval.markers' in config['Norm_DimRed_Eval_GE'] and config['Norm_DimRed_Eval_GE']['eval.markers'] != None) else "NULL" # formated "GAPDH,Actin,other"
    # Normalization and dimension reduction
    NDRE_FEATURES_N = config['Norm_DimRed_Eval_GE']['features.n'] if ('Norm_DimRed_Eval_GE' in config and 'features.n' in config['Norm_DimRed_Eval_GE'] and config['Norm_DimRed_Eval_GE']['features.n'] != None) else "NULL"
    NDRE_NORM_METHOD = config['Norm_DimRed_Eval_GE']['norm.method'] if ('Norm_DimRed_Eval_GE' in config and 'norm.method' in config['Norm_DimRed_Eval_GE'] and config['Norm_DimRed_Eval_GE']['norm.method'] != None) else 'SCTransform'
    NDRE_DIMRED_METHOD = config['Norm_DimRed_Eval_GE']['dimred.method'] if ('Norm_DimRed_Eval_GE' in config and 'dimred.method' in config['Norm_DimRed_Eval_GE'] and config['Norm_DimRed_Eval_GE']['dimred.method'] != None) else "pca"
    NDRE_VTR = config['Norm_DimRed_Eval_GE']['vtr'].replace(", ", ",") if ('Norm_DimRed_Eval_GE' in config and 'vtr' in config['Norm_DimRed_Eval_GE'] and config['Norm_DimRed_Eval_GE']['vtr'] != None) else "NULL"
    NDRE_VTR_SCALE = config['Norm_DimRed_Eval_GE']['vtr.scale'] if ('Norm_DimRed_Eval_GE' in config and 'vtr.scale' in config['Norm_DimRed_Eval_GE'] and config['Norm_DimRed_Eval_GE']['vtr.scale'] != None) else "NULL"
    NDRE_DIM_MAX = config['Norm_DimRed_Eval_GE']['dims.max'] if ('Norm_DimRed_Eval_GE' in config and 'dims.max' in config['Norm_DimRed_Eval_GE'] and config['Norm_DimRed_Eval_GE']['dims.max'] != None) else 49
    NDRE_DIM_MIN = config['Norm_DimRed_Eval_GE']['dims.min'] if ('Norm_DimRed_Eval_GE' in config and 'dims.min' in config['Norm_DimRed_Eval_GE'] and config['Norm_DimRed_Eval_GE']['dims.min'] != None) else 3
    NDRE_DIM_STEPS = config['Norm_DimRed_Eval_GE']['dims.steps'] if ('Norm_DimRed_Eval_GE' in config and 'dims.steps' in config['Norm_DimRed_Eval_GE'] and config['Norm_DimRed_Eval_GE']['dims.steps'] != None) else 2
    NDRE_RES_MAX = config['Norm_DimRed_Eval_GE']['res.max'] if ('Norm_DimRed_Eval_GE' in config and 'res.max' in config['Norm_DimRed_Eval_GE'] and config['Norm_DimRed_Eval_GE']['res.max'] != None) else 1.2
    NDRE_RES_MIN = config['Norm_DimRed_Eval_GE']['res.min'] if ('Norm_DimRed_Eval_GE' in config and 'res.min' in config['Norm_DimRed_Eval_GE'] and config['Norm_DimRed_Eval_GE']['res.min'] != None) else 0.1
    NDRE_RES_STEPS = config['Norm_DimRed_Eval_GE']['res.steps'] if ('Norm_DimRed_Eval_GE' in config and 'res.steps' in config['Norm_DimRed_Eval_GE'] and config['Norm_DimRed_Eval_GE']['res.steps'] != None) else 0.1
    ### Snakefile parameters
    #Correspondance sample/input/output
    dic_NDRE_INFO = {}
    for i in range(0,len(NDRE_SAMPLE_NAME_GE),1):
        dic_NDRE_INFO[NDRE_SAMPLE_NAME_GE[i]] = {}
        dic_NDRE_INFO[NDRE_SAMPLE_NAME_GE[i]]['NDRE_INPUT_RDA'] = NDRE_INPUT_RDA_GE[i]
        dic_NDRE_INFO[NDRE_SAMPLE_NAME_GE[i]]['NDRE_OUTPUT_DIR'] = NDRE_OUTPUT_DIR_GE[i]
    #Names
    NDRE_NORM_VTR = NDRE_NORM_METHOD if (NDRE_NORM_METHOD == "LogNormalize" or NDRE_VTR == "NULL") else (NDRE_NORM_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(NDRE_VTR.split(","))))))
    NDRE_DIMRED_VTR = NDRE_DIMRED_METHOD if (NDRE_DIMRED_METHOD == "pca" or NDRE_VTR == "NULL") else (NDRE_DIMRED_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(NDRE_VTR.split(","))))))
    POSSIBLE_DIM = ["%.0f" % number for number in numpy.arange(NDRE_DIM_MIN,NDRE_DIM_MAX+1,NDRE_DIM_STEPS)]
    POSSIBLE_RES = ["%.1f" % number for number in numpy.arange(NDRE_RES_MIN,NDRE_RES_MAX+0.1,NDRE_RES_STEPS)]
    ASSAY = "RNA" if NDRE_NORM_METHOD == "LogNormalize" else "SCT"

if "Clust_Markers_Annot_GE" in STEPS:
    ### Sample/Project
    if ('Clust_Markers_Annot_GE' in config) and ('sample.name.ge' in config['Clust_Markers_Annot_GE']) and ('input.rda.ge' in config['Clust_Markers_Annot_GE']) :
        CMA_SAMPLE_NAME_GE_RAW = config['Clust_Markers_Annot_GE']['sample.name.ge']
        CMA_INPUT_RDA_GE = config['Clust_Markers_Annot_GE']['input.rda.ge']
        #check samples names and add "_GE" if needed
        CMA_SAMPLE_NAME_GE = []
        for i in range(0,len(CMA_SAMPLE_NAME_GE_RAW),1):
            CMA_SAMPLE_NAME_GE.append(CMA_SAMPLE_NAME_GE_RAW[i] + "_GE") if (CMA_SAMPLE_NAME_GE_RAW[i][len(CMA_SAMPLE_NAME_GE_RAW[i])-3:] != "_GE") else CMA_SAMPLE_NAME_GE.append(CMA_SAMPLE_NAME_GE_RAW[i])
    elif "Norm_DimRed_Eval_GE" in STEPS:
        sys.stderr.write("Note: No input.rda.ge and sample.name.ge find in Clust_Markers_Annot_GE section of configfile; input.rda.ge and sample.name.ge will be determine from Norm_DimRed_Eval_GE step for Clust_Markers_Annot_GE step!\n")
        CMA_SAMPLE_NAME_GE = copy.deepcopy(NDRE_SAMPLE_NAME_GE)
        CMA_INPUT_RDA_GE = [os.path.normpath(dic_NDRE_INFO[x]['NDRE_OUTPUT_DIR'] + "/" + NDRE_NORM_VTR + "/" + NDRE_DIMRED_VTR + "/" + x + "_" + NDRE_NORM_VTR + "_" + NDRE_DIMRED_VTR + ".rda") for x in NDRE_SAMPLE_NAME_GE]
    else:
        sys.exit("Error: No sample.name.ge or/and input.rda.ge in configfile!\n")
    if ('Clust_Markers_Annot_GE' in config) and ('output.dir.ge' in config['Clust_Markers_Annot_GE']) :
        CMA_OUTPUT_DIR_GE = config['Clust_Markers_Annot_GE']['output.dir.ge']
    elif "Norm_DimRed_Eval_GE" in STEPS:
        CMA_OUTPUT_DIR_GE = [os.path.normpath(dic_NDRE_INFO[x]['NDRE_OUTPUT_DIR'] + "/" + NDRE_NORM_VTR + "/" + NDRE_DIMRED_VTR) for x in NDRE_SAMPLE_NAME_GE]
        sys.stderr.write("Note: No output.dir.ge find in Clust_Markers_Annot_GE section of configfile; output.dir.ge will be determine from Norm_DimRed_Eval_GE step for Clust_Markers_Annot_GE step!\n")
    else :
        sys.exit("Error: No output.dir.ge find in configfile!\n")
    ### Analysis Parameters
    CMA_AUTHOR_NAME = config['Clust_Markers_Annot_GE']['author.name'].replace(", ", ",").replace(" ", "_") if ('Clust_Markers_Annot_GE' in config and 'author.name' in config['Clust_Markers_Annot_GE'] and config['Clust_Markers_Annot_GE']['author.name'] != None) else "NULL"
    CMA_AUTHOR_MAIL = config['Clust_Markers_Annot_GE']['author.mail'].replace(", ", ",") if ('Clust_Markers_Annot_GE' in config and 'author.mail' in config['Clust_Markers_Annot_GE'] and config['Clust_Markers_Annot_GE']['author.mail'] != None) else "NULL"
    CMA_MARKFILE = config['Clust_Markers_Annot_GE']['markfile'].replace(", ", ",") if ('Clust_Markers_Annot_GE' in config and 'markfile' in config['Clust_Markers_Annot_GE'] and config['Clust_Markers_Annot_GE']['markfile'] != None) else "NULL" # formated "file1,file2,file3"
    # Normalization and dimension reduction
    CMA_KEEP_DIM = config['Clust_Markers_Annot_GE']['keep.dims'] if ('Clust_Markers_Annot_GE' in config and 'keep.dims' in config['Clust_Markers_Annot_GE'] and config['Clust_Markers_Annot_GE']['keep.dims'] != None) else "NULL"
    CMA_KEEP_RES = config['Clust_Markers_Annot_GE']['keep.res'] if ('Clust_Markers_Annot_GE' in config and 'keep.res' in config['Clust_Markers_Annot_GE'] and config['Clust_Markers_Annot_GE']['keep.res'] != None) else "NULL"
    CMA_CFR_MINSCORE = config['Clust_Markers_Annot_GE']['cfr.minscore'] if ('Clust_Markers_Annot_GE' in config and 'cfr.minscore' in config['Clust_Markers_Annot_GE'] and config['Clust_Markers_Annot_GE']['cfr.minscore'] != None) else "NULL"
    CMA_SR_MINSCORE = config['Clust_Markers_Annot_GE']['sr.minscore'] if ('Clust_Markers_Annot_GE' in config and 'sr.minscore' in config['Clust_Markers_Annot_GE'] and config['Clust_Markers_Annot_GE']['sr.minscore'] != None) else "NULL"
    ### Snakefile parameters
    #check end paths (add "/" if necessary)
    for i in range(0,len(CMA_OUTPUT_DIR_GE),1):
        CMA_OUTPUT_DIR_GE[i] = os.path.normpath(CMA_OUTPUT_DIR_GE[i])
    #Correspondance sample/input/output
    dic_CMA_INFO = {}
    CMA_COMPLEMENT = []
    for i in range(0,len(CMA_SAMPLE_NAME_GE),1):
        dic_CMA_INFO[CMA_SAMPLE_NAME_GE[i]] = {}
        dic_CMA_INFO[CMA_SAMPLE_NAME_GE[i]]['CMA_INPUT_RDA'] = CMA_INPUT_RDA_GE[i]
        dic_CMA_INFO[CMA_SAMPLE_NAME_GE[i]]['CMA_OUTPUT_DIR'] = CMA_OUTPUT_DIR_GE[i]
        compl = os.path.splitext(os.path.basename(CMA_INPUT_RDA_GE[i]))[0]
        if compl.startswith(CMA_SAMPLE_NAME_GE[i]):
            compl = compl[len(CMA_SAMPLE_NAME_GE[i]):]
        CMA_COMPLEMENT.append(compl)
    #Names
    CMA_CLUST_FOLDER = "dims" + str(CMA_KEEP_DIM) + "_res" + str(CMA_KEEP_RES)

if "Droplets_QC_GE" in STEPS or "Filtering_GE" in STEPS or "Norm_DimRed_Eval_GE" in STEPS or "Clust_Markers_Annot_GE" in STEPS or "Adding_ADT" in STEPS:
    SINGULARITY_ENV = PIPELINE_FOLDER + "/envs/singularity/single_cell.simg"

if "Adding_ADT" in STEPS:
    ### Sample/Project
    if 'Adding_ADT' in config and 'input.rda.ge' in config['Adding_ADT'] :
        ADD_ADT_INPUT_RDA_GE = config['Adding_ADT']['input.rda.ge']
    elif "Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("Note: No input.rda.ge find in Adding_ADT section of configfile; input.rda.ge will be determine from Clust_Markers_Annot_GE step for Adding_ADT step!\n")
        ADD_ADT_INPUT_RDA_GE = [os.path.normpath(os.path.dirname(dic_CMA_INFO[CMA_SAMPLE_NAME_GE[x]]['CMA_INPUT_RDA']) + "/" + CMA_CLUST_FOLDER + "/" + CMA_SAMPLE_NAME_GE[x] + CMA_COMPLEMENT[x] + "_" + str(CMA_KEEP_DIM) + "_" + str(CMA_KEEP_RES) + ".rda") for x in range(len(CMA_SAMPLE_NAME_GE))]
    else:
        sys.exit("Error: No input.rda.ge in configfile!\n")
    if 'Adding_ADT' in config and 'sample.name.adt' in config['Adding_ADT'] and 'input.dir.adt' in config['Adding_ADT'] :
        ADD_ADT_INPUT_DIR_ADT = config['Adding_ADT']['input.dir.adt']
        ADD_ADT_SAMPLE_NAME_ADT_RAW = config['Adding_ADT']['sample.name.adt']
        #check samples names and add "_ADT" if needed
        ADD_ADT_SAMPLE_NAME_ADT = []
        for i in range(0,len(ADD_ADT_SAMPLE_NAME_GE_RAW),1):
            ADD_ADT_SAMPLE_NAME_ADT.append(ADD_ADT_SAMPLE_NAME_ADT_RAW[i] + "_ADT") if (ADD_ADT_SAMPLE_NAME_ADT_RAW[i][len(ADD_ADT_SAMPLE_NAME_ADT_RAW[i])-4:] != "_ADT") else ADD_ADT_SAMPLE_NAME_ADT.append(ADD_ADT_SAMPLE_NAME_ADT_RAW[i])
    elif "Alignment_countTable_ADT" in STEPS:
        sys.stderr.write("Note: No sample.name.adt or input.dir.adt find in Adding_ADT section of configfile; sample.name.adt and input.dir.adt will be determine from Alignment_countTable_ADT step for Adding_ADT step!\n")
        ADD_ADT_INPUT_DIR_ADT = [ os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/" + str(x) + "/KALLISTOBUS") for x in ALIGN_SAMPLE_NAME_ADT]
        ADD_ADT_SAMPLE_NAME_ADT = copy.deepcopy(ALIGN_SAMPLE_NAME_ADT)
    else:
        sys.exit("Error: No sample.name.adt or input.dir.adt in configfile!\n")
    ### Analysis Parameters
    ADD_ADT_AUTHOR_NAME = config['Adding_ADT']['author.name'].replace(", ", ",").replace(" ", "_") if ('Adding_ADT' in config and 'author.name' in config['Adding_ADT'] and config['Adding_ADT']['author.name'] != None) else "NULL"
    ADD_ADT_AUTHOR_MAIL = config['Adding_ADT']['author.mail'].replace(", ", ",") if ('Adding_ADT' in config and 'author.mail' in config['Adding_ADT'] and config['Adding_ADT']['author.mail'] != None) else "NULL"
    ADD_ADT_GENE_NAMES = config['Adding_ADT']['gene.names'].replace(", ", ",") if ('Adding_ADT' in config and 'gene.names' in config['Adding_ADT'] and config['Adding_ADT']['gene.names'] != None) else "NULL"
    sys.stderr.write(ADD_ADT_GENE_NAMES)
    ADD_ADT_MAX_CUTOFF = config['Adding_ADT']['ADT.max.cutoff'].replace(", ", ",") if ('Adding_ADT' in config and 'ADT.max.cutoff' in config['Adding_ADT'] and config['Adding_ADT']['ADT.max.cutoff'] != None) else "NULL"
    ADD_ADT_MIN_CUTOFF = config['Adding_ADT']['ADT.min.cutoff'].replace(", ", ",") if ('Adding_ADT' in config and 'ADT.min.cutoff' in config['Adding_ADT'] and config['Adding_ADT']['ADT.min.cutoff'] != None) else "NULL"
    ### Snakefile parameters
    #Correspondance input/output
    ADD_ADT_OUTPUT = [os.path.splitext(x)[0] for x in ADD_ADT_INPUT_RDA_GE]
    dic_ADD_ADT_INFO = {}
    for i in range(0,len(ADD_ADT_OUTPUT),1):
        dic_ADD_ADT_INFO[ADD_ADT_OUTPUT[i]] = {}
        dic_ADD_ADT_INFO[ADD_ADT_OUTPUT[i]]['ADD_ADT_INPUT_RDA_GE'] = ADD_ADT_INPUT_RDA_GE[i]
        dic_ADD_ADT_INFO[ADD_ADT_OUTPUT[i]]['ADD_ADT_INPUT_DIR_ADT'] = ADD_ADT_INPUT_DIR_ADT[i]
        dic_ADD_ADT_INFO[ADD_ADT_OUTPUT[i]]['ADD_ADT_SAMPLE_NAME_ADT'] = ADD_ADT_SAMPLE_NAME_ADT[i]

if "Adding_TCR" in STEPS:
    ### Sample/Project
    if 'Adding_TCR' in config and 'input.rda' in config['Adding_TCR'] :
        ADD_TCR_INPUT_RDA_GE = config['Adding_TCR']['input.rda']
    elif "Adding_ADT" in STEPS:
        sys.stderr.write("Note: No input.rda find in Adding_TCR section of configfile; input.rda will be determine from Adding_ADT step for Adding_TCR step!\n")
        ADD_TCR_INPUT_RDA_GE = [ x + "_ADT.rda" for x in ADD_ADT_OUTPUT]
    elif "Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("Note: No input.rda find in Adding_TCR section of configfile; input.rda will be determine from Clust_Markers_Annot_GE step for Adding_TCR step!\n")
        ADD_TCR_INPUT_RDA_GE = [os.path.normpath(os.path.dirname(dic_CMA_INFO[CMA_SAMPLE_NAME_GE[x]]['CMA_INPUT_RDA']) + "/" + CMA_CLUST_FOLDER + "/" + CMA_SAMPLE_NAME_GE[x] + CMA_COMPLEMENT[x] + "_" + str(CMA_KEEP_DIM) + "_" + str(CMA_KEEP_RES) + ".rda") for x in range(len(CMA_SAMPLE_NAME_GE))]
    else:
        sys.exit("Error: No input.rda in configfile!\n")
    if 'Adding_TCR' in config and 'vdj.input.file.tcr' in config['Adding_TCR'] :
        ADD_TCR_INPUT_CSV_TCR = config['Adding_TCR']['vdj.input.file.tcr']
    elif "Alignment_annotations_TCR_BCR" in STEPS:
        sys.stderr.write("Note: No vdj.input.file.tcr find in Adding_TCR section of configfile; vdj.input.file.tcr will be determine from Alignment_annotations_TCR_BCR step for Adding_TCR step!\n")
        ALIGN_SAMPLE_NAME_TCR = [sample for sample in ALIGN_SAMPLE_NAME_TCR_BCR if bool(re.match(".+_TCR", sample))]
        ADD_TCR_INPUT_CSV_TCR = [ os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/" + x + "/" + x + "_CellRanger/outs/filtered_contig_annotations.csv") for x in ALIGN_SAMPLE_NAME_TCR]
    else:
        sys.exit("Error: No vdj.input.file.tcr in configfile!\n")
    ### Analysis Parameters
    ADD_TCR_AUTHOR_NAME = config['Adding_TCR']['author.name'].replace(", ", ",").replace(" ", "_") if ('Adding_TCR' in config and 'author.name' in config['Adding_TCR'] and config['Adding_TCR']['author.name'] != None) else "NULL"
    ADD_TCR_AUTHOR_MAIL = config['Adding_TCR']['author.mail'].replace(", ", ",") if ('Adding_TCR' in config and 'author.mail' in config['Adding_TCR'] and config['Adding_TCR']['author.mail'] != None) else "NULL"
    ### Snakefile parameters
    #Correspondance input/output
    ADD_TCR_OUTPUT = [os.path.splitext(x)[0] for x in ADD_TCR_INPUT_RDA_GE]
    dic_ADD_TCR_INFO = {}
    for i in range(0,len(ADD_TCR_OUTPUT),1):
        dic_ADD_TCR_INFO[ADD_TCR_OUTPUT[i]] = {}
        dic_ADD_TCR_INFO[ADD_TCR_OUTPUT[i]]['ADD_TCR_INPUT_RDA_GE'] = ADD_TCR_INPUT_RDA_GE[i]
        dic_ADD_TCR_INFO[ADD_TCR_OUTPUT[i]]['ADD_TCR_INPUT_CSV_TCR'] = ADD_TCR_INPUT_CSV_TCR[i]

if "Adding_BCR" in STEPS:
    ### Sample/Project
    if 'Adding_BCR' in config and 'input.rda' in config['Adding_BCR'] :
        ADD_BCR_INPUT_RDA_GE = config['Adding_BCR']['input.rda']
    elif "Adding_TCR" in STEPS:
        sys.stderr.write("Note: No input.rda.ge find in Cerebro section of configfile; input.rda.ge will be determine from Adding_TCR step for Adding_BCR step!\n")
        ADD_BCR_INPUT_RDA_GE = [ x + "_TCR.rda" for x in ADD_TCR_OUTPUT]
    elif "Adding_ADT" in STEPS:
        sys.stderr.write("Note: No input.rda find in Adding_BCR section of configfile; input.rda will be determine from Adding_ADT step for Adding_BCR step!\n")
        ADD_BCR_INPUT_RDA_GE = [ x + "_ADT.rda" for x in ADD_ADT_OUTPUT]
    elif "Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("Note: No input.rda find in Adding_BCR section of configfile; input.rda will be determine from Clust_Markers_Annot_GE step for Adding_BCR step!\n")
        ADD_BCR_INPUT_RDA_GE = [os.path.normpath(os.path.dirname(dic_CMA_INFO[CMA_SAMPLE_NAME_GE[x]]['CMA_INPUT_RDA']) + "/" + CMA_CLUST_FOLDER + "/" + CMA_SAMPLE_NAME_GE[x] + CMA_COMPLEMENT[x] + "_" + str(CMA_KEEP_DIM) + "_" + str(CMA_KEEP_RES) + ".rda") for x in range(len(CMA_SAMPLE_NAME_GE))]
    else:
        sys.exit("Error: No input.rda in configfile!\n")
    if 'Adding_BCR' in config and 'vdj.input.file.bcr' in config['Adding_BCR'] :
        ADD_BCR_INPUT_CSV_BCR = config['Adding_BCR']['vdj.input.file.bcr']
    elif "Alignment_annotations_TCR_BCR" in STEPS:
        sys.stderr.write("Note: No vdj.input.file.bcr find in Adding_BCR section of configfile; vdj.input.file.bcr will be determine from Alignment_annotations_BCR_BCR step for Adding_BCR step!\n")
        ALIGN_SAMPLE_NAME_BCR = [sample for sample in ALIGN_SAMPLE_NAME_TCR_BCR if bool(re.match(".+_BCR", sample))]
        ADD_BCR_INPUT_CSV_BCR = [ os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/" + x + "/" + x + "_CellRanger/outs/filtered_contig_annotations.csv") for x in ALIGN_SAMPLE_NAME_BCR]
    else:
        sys.exit("Error: No vdj.input.file.bcr in configfile!\n")
    ### Analysis Parameters
    ADD_BCR_AUTHOR_NAME = config['Adding_BCR']['author.name'].replace(", ", ",").replace(" ", "_") if ('Adding_BCR' in config and 'author.name' in config['Adding_BCR'] and config['Adding_BCR']['author.name'] != None) else "NULL"
    ADD_BCR_AUTHOR_MAIL = config['Adding_BCR']['author.mail'].replace(", ", ",") if ('Adding_BCR' in config and 'author.mail' in config['Adding_BCR'] and config['Adding_BCR']['author.mail'] != None) else "NULL"
    ### Snakefile parameters
    #Correspondance input/output
    ADD_BCR_OUTPUT = [os.path.splitext(x)[0] for x in ADD_BCR_INPUT_RDA_GE]
    dic_ADD_BCR_INFO = {}
    for i in range(0,len(ADD_BCR_OUTPUT),1):
        dic_ADD_BCR_INFO[ADD_BCR_OUTPUT[i]] = {}
        dic_ADD_BCR_INFO[ADD_BCR_OUTPUT[i]]['ADD_BCR_INPUT_RDA_GE'] = ADD_BCR_INPUT_RDA_GE[i]
        dic_ADD_BCR_INFO[ADD_BCR_OUTPUT[i]]['ADD_BCR_INPUT_CSV_BCR'] = ADD_BCR_INPUT_CSV_BCR[i]

if "Alignment_annotations_TCR_BCR" in STEPS or "Adding_TCR" in STEPS or "Adding_BCR" in STEPS:
    SINGULARITY_ENV_TCR_BCR = PIPELINE_FOLDER + "/envs/singularity/single_cell_TCR_BCR.simg"

if "Cerebro" in STEPS:
    ### Sample/Project
    if 'Cerebro' in config and 'input.rda.ge' in config['Cerebro'] :
        CEREBRO_INPUT_RDA = config['Cerebro']['input.rda']
    elif "Adding_BCR" in STEPS:
        sys.stderr.write("Note: No input.rda.ge find in Cerebro section of configfile; input.rda.ge will be determine from Adding_BCR_GE step for Cerebro step!\n")
        CEREBRO_INPUT_RDA = [ x + "_BCR.rda" for x in ADD_BCR_OUTPUT]
    elif "Adding_TCR" in STEPS:
        sys.stderr.write("Note: No input.rda.ge find in Cerebro section of configfile; input.rda.ge will be determine from Adding_TCR_GE step for Cerebro step!\n")
        CEREBRO_INPUT_RDA = [ x + "_TCR.rda" for x in ADD_TCR_OUTPUT]
    elif "Adding_ADT" in STEPS:
        sys.stderr.write("Note: No input.rda find in Cerebro section of configfile; input.rda will be determine from Adding_ADT step for Cerebro step!\n")
        CEREBRO_INPUT_RDA = [ x + "_ADT.rda" for x in ADD_ADT_OUTPUT]
    elif "Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("Note: No input.rda find in Cerebro section of configfile; input.rda will be determine from Clust_Markers_Annot_GE step for Cerebro step!\n")
        CEREBRO_INPUT_RDA = [os.path.normpath(os.path.dirname(dic_CMA_INFO[CMA_SAMPLE_NAME_GE[x]]['CMA_INPUT_RDA']) + "/" + CMA_CLUST_FOLDER + "/" + CMA_SAMPLE_NAME_GE[x] + CMA_COMPLEMENT[x] + "_" + str(CMA_KEEP_DIM) + "_" + str(CMA_KEEP_RES) + ".rda") for x in range(len(CMA_SAMPLE_NAME_GE))]
    else:
        sys.exit("Error: No input.rda in configfile!\n")
    ### Analysis Parameters
    CEREBRO_AUTHOR_NAME = config['Cerebro']['author.name'].replace(", ", ",").replace(" ", "_") if ('Cerebro' in config and 'author.name' in config['Cerebro'] and config['Cerebro']['author.name'] != None) else "NULL"
    CEREBRO_AUTHOR_MAIL = config['Cerebro']['author.mail'].replace(", ", ",") if ('Cerebro' in config and 'author.mail' in config['Cerebro'] and config['Cerebro']['author.mail'] != None) else "NULL"
    # Normalization and dimension reduction
    CEREBRO_VERSION = config['Cerebro']['version'] if ('Cerebro' in config and 'version' in config['Cerebro'] and config['Cerebro']['version'] != None) else "v1.3"
    CEREBRO_GROUPS = config['Cerebro']['groups'].replace(", ", ",") if ('Cerebro' in config and 'groups' in config['Cerebro'] and config['Cerebro']['groups'] != None) else "NULL"
    CEREBRO_REMOVE_OTHER_RED = config['Cerebro']['remove.other.reductions'] if ('Cerebro' in config and 'remove.other.reductions' in config['Cerebro'] and config['Cerebro']['remove.other.reductions'] != None) else "NULL"
    CEREBRO_REMOVE_OTHER_IDENT = config['Cerebro']['remove.other.idents'] if ('Cerebro' in config and 'remove.other.idents' in config['Cerebro'] and config['Cerebro']['remove.other.idents'] != None) else "NULL"
    CEREBRO_REMOVE_MT = config['Cerebro']['remove.mt.genes'] if ('Cerebro' in config and 'remove.mt.genes' in config['Cerebro'] and config['Cerebro']['remove.mt.genes'] != None) else "NULL"
    CEREBRO_REMOVE_CRB = config['Cerebro']['remove.crb.genes'] if ('Cerebro' in config and 'remove.crb.genes' in config['Cerebro'] and config['Cerebro']['remove.crb.genes'] != None) else "NULL"
    CEREBRO_REMOVE_STR = config['Cerebro']['remove.str.genes'] if ('Cerebro' in config and 'remove.str.genes' in config['Cerebro'] and config['Cerebro']['remove.str.genes'] != None) else "NULL"
    CEREBRO_ONLY_POS_DE = config['Cerebro']['only.pos.DE'] if ('Cerebro' in config and 'only.pos.DE' in config['Cerebro'] and config['Cerebro']['only.pos.DE'] != None) else "NULL"
    CEREBRO_REMOVE_CUSTOM_DE = config['Cerebro']['remove.custom.DE'] if ('Cerebro' in config and 'remove.custom.DE' in config['Cerebro'] and config['Cerebro']['remove.custom.DE'] != None) else "NULL"
    CEREBRO_GMT_FILE = config['Cerebro']['gmt.file'] if ('Cerebro' in config and 'gmt.file' in config['Cerebro'] and config['Cerebro']['gmt.file'] != None) else "NULL"
    ### Snakefile parameters
    #Creation output complement + extention:
    CEREBRO_COMPLEMENT = ""
    if CEREBRO_REMOVE_MT == "TRUE" or CEREBRO_REMOVE_MT == "True": CEREBRO_COMPLEMENT = CEREBRO_COMPLEMENT + "_noMT"
    if CEREBRO_REMOVE_CRB == "TRUE" or CEREBRO_REMOVE_CRB == "True": CEREBRO_COMPLEMENT = CEREBRO_COMPLEMENT + "_noRB"
    if CEREBRO_REMOVE_STR == "TRUE" or CEREBRO_REMOVE_STR == "True": CEREBRO_COMPLEMENT = CEREBRO_COMPLEMENT + "_noSTR"
    CEREBRO_COMPLEMENT_CRB = [CEREBRO_COMPLEMENT + ".crb"]
    if (CEREBRO_GROUPS != None or CEREBRO_GROUPS != "NULL") and CEREBRO_VERSION == "v1.2":
        for group in CEREBRO_GROUPS.split(','):
            CEREBRO_COMPLEMENT_CRB.append(CEREBRO_COMPLEMENT +  '_clusterIs_' + group + ".crb")
    #Correspondance sample/input/output
    CEREBRO_INPUT_RDA_NO_EXTENTION = [os.path.splitext(x)[0] for x in CEREBRO_INPUT_RDA]
    #Singularity environnement
    if CEREBRO_VERSION == "v1.2":
        SINGULARITY_ENV_CEREBRO = PIPELINE_FOLDER + "/envs/singularity/single_cell_oldcerebro.simg"
    elif CEREBRO_VERSION == "v1.3":
        SINGULARITY_ENV_CEREBRO = PIPELINE_FOLDER + "/envs/singularity/single_cell.simg"
    else:
        sys.exit("Error: Unknown version of cerebro in configfile!\n")


onstart:
    sys.stderr.write("\n############################################################# \n")
    sys.stderr.write("\n\n\t Single-cell RNA-seq pipeline \n\n")
    sys.stderr.write("\n############################################################# \n\n")
    sys.stderr.write("***************** PARAMETERS ******************\n")
    if "Alignment_countTable_GE" in STEPS:
        sys.stderr.write("\n" + "Alignment_countTable_GE:" + "\n")
        sys.stderr.write("SAMPLE(S):\n")
        if isinstance(ALIGN_SAMPLE_NAME_GE_RAW, list):
            for sample in ALIGN_SAMPLE_NAME_GE_RAW:
                sys.stderr.write("\t" + str(sample) + "\n")
        else:
            sys.stderr.write("\t" + str(ALIGN_SAMPLE_NAME_GE_RAW) + "\n")

    if "Alignment_countTable_ADT" in STEPS:
        sys.stderr.write("\n" + "Alignment_countTable_ADT:" + "\n")
        sys.stderr.write("SAMPLE(S):\n")
        if isinstance(ALIGN_SAMPLE_NAME_ADT_RAW, list):
            for sample in ALIGN_SAMPLE_NAME_ADT_RAW:
                sys.stderr.write("\t" + str(sample) + "\n")
        else:
            sys.stderr.write("\t" + str(ALIGN_SAMPLE_NAME_ADT_RAW) + "\n")

    if "Alignment_annotations_TCR_BCR" in STEPS:
        sys.stderr.write("\n" + "Alignment_annotations_TCR_BCR:" + "\n")
        sys.stderr.write("SAMPLE(S):\n")
        if isinstance(ALIGN_SAMPLE_NAME_TCR_BCR_RAW, list):
            for sample in ALIGN_SAMPLE_NAME_TCR_BCR_RAW:
                sys.stderr.write("\t" + str(sample) + "\n")
        else:
            sys.stderr.write("\t" + str(ALIGN_SAMPLE_NAME_TCR_BCR_RAW) + "\n")

    if "Droplets_QC_GE" in STEPS:
        sys.stderr.write("\n" + "Droplets_QC_GE:" + "\n")
        sys.stderr.write("SAMPLE(S):\n")
        if isinstance(QC_SAMPLE_NAME_GE, list):
            for sample in QC_SAMPLE_NAME_GE:
                sys.stderr.write("\t" + str(sample) + "\n")
        else:
            sys.stderr.write("\t" + str(QC_SAMPLE_NAME_GE) + "\n")

    if "Filtering_GE" in STEPS:
        sys.stderr.write("\n" + "Filtering_GE:" + "\n")
        sys.stderr.write("SAMPLE(S):\n")
        if isinstance(FILERING_SAMPLE_NAME_GE, list):
            for sample in FILERING_SAMPLE_NAME_GE:
                sys.stderr.write("\t" + str(sample) + "\n")
        else:
            sys.stderr.write("\t" + str(FILERING_SAMPLE_NAME_GE) + "\n")

    if "Norm_DimRed_Eval_GE" in STEPS:
        sys.stderr.write("\n" + "Norm_DimRed_Eval_GE:" + "\n")
        sys.stderr.write("SAMPLE(S):\n")
        if isinstance(NDRE_SAMPLE_NAME_GE, list):
            for sample in NDRE_SAMPLE_NAME_GE:
                sys.stderr.write("\t" + str(sample) + "\n")
        else:
            sys.stderr.write("\t" + str(NDRE_SAMPLE_NAME_GE) + "\n")
    if "Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("\n" + "Clust_Markers_Annot_GE:" + "\n")
        sys.stderr.write("SAMPLE(S):\n")
        if isinstance(CMA_SAMPLE_NAME_GE, list):
            for sample in CMA_SAMPLE_NAME_GE:
                sys.stderr.write("\t" + str(sample) + "\n")
        else:
            sys.stderr.write("\t" + str(CMA_SAMPLE_NAME_GE) + "\n")
    if "Adding_ADT" in STEPS:
        sys.stderr.write("\n" + "Adding_ADT:" + "\n")
        sys.stderr.write("RDA FILE(S):\n")
        if isinstance(ADD_ADT_INPUT_RDA_GE, list):
            for sample in ADD_ADT_INPUT_RDA_GE:
                sys.stderr.write("\t" + str(sample) + "\n")
        else:
            sys.stderr.write("\t" + str(ADD_ADT_INPUT_RDA_GE) + "\n")
    if "Adding_TCR" in STEPS:
        sys.stderr.write("\n" + "Adding_TCR:" + "\n")
        sys.stderr.write("RDA FILE(S):\n")
        if isinstance(ADD_TCR_INPUT_RDA_GE, list):
            for sample in ADD_TCR_INPUT_RDA_GE:
                sys.stderr.write("\t" + str(sample) + "\n")
        else:
            sys.stderr.write("\t" + str(ADD_TCR_INPUT_RDA_GE) + "\n")
    if "Adding_BCR" in STEPS:
        sys.stderr.write("\n" + "Adding_BCR:" + "\n")
        sys.stderr.write("RDA FILE(S):\n")
        if isinstance(ADD_BCR_INPUT_RDA_GE, list):
            for sample in ADD_BCR_INPUT_RDA_GE:
                sys.stderr.write("\t" + str(sample) + "\n")
        else:
            sys.stderr.write("\t" + str(ADD_BCR_INPUT_RDA_GE) + "\n")
    if "Cerebro" in STEPS:
        sys.stderr.write("\n" + "Cerebro:" + "\n")
        sys.stderr.write("RDA FILE(S):\n")
        if isinstance(CEREBRO_INPUT_RDA, list):
            for sample in CEREBRO_INPUT_RDA:
                sys.stderr.write("\t" + str(sample) + "\n")
        else:
            sys.stderr.write("\t" + str(CEREBRO_INPUT_RDA) + "\n")

    sys.stderr.write("\n***************** RUN ******************\n")
    return []


### rule all ###################################################################################################################################
include: "rules/Rule_all.smk"
rule all:
    input:
        **get_targets(STEPS)
    message:
        "Single-cell RNA-seq pipeline done!"


### real rules ###################################################################################################################################
if "Alignment_countTable_GE" in STEPS:
    include: "rules/Alignment_countTable_GE.smk"

if "Alignment_countTable_ADT" in STEPS:
    include: "rules/Alignment_countTable_ADT.smk"

if "Alignment_annotations_TCR_BCR" in STEPS:
    include: "rules/Alignment_annotations_TCR_BCR.smk"

if "Droplets_QC_GE" in STEPS:
    include: "rules/Droplets_QC_GE.smk"

if "Filtering_GE" in STEPS:
    include: "rules/Filtering_GE.smk"

if "Norm_DimRed_Eval_GE" in STEPS:
    include: "rules/Norm_DimRed_Eval_GE.smk"

if "Clust_Markers_Annot_GE" in STEPS:
    include: "rules/Clust_Markers_Annot_GE.smk"

if "Adding_ADT" in STEPS:
    include: "rules/Adding_ADT.smk"

if "Adding_TCR" in STEPS:
    include: "rules/Adding_TCR.smk"

if "Adding_BCR" in STEPS:
    include: "rules/Adding_BCR.smk"

if "Cerebro" in STEPS:
    include: "rules/Cerebro.smk"
