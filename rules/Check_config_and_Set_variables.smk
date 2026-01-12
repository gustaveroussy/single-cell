### parameters ###################################################################################################################################
sys.stderr.write("\n#################### Setting Parameters ####################\n\n")

#function to get value from config file and format if None or "" or "NULL"
def getelseDefault(step, parameter, default_value):
    if step in config and parameter in config[step] and config[step][parameter] != None and config[step][parameter] != "" and config[step][parameter] != "NULL":
        return config[step][parameter]
    else:
        return default_value

#Manage steps
STEPS = config['Steps']

#Manage analysis author
AUTHOR_NAME = config['Author.name'].replace(", ", ",") if 'Author.name' in config else "NULL"
AUTHOR_MAIL = config['Author.mail'].replace(", ", ",") if 'Author.mail' in config else "NULL"

#Manage species
if "Alignment_countTable_GE" in STEPS or "Alignment_annotations_TCR_BCR" in STEPS or "Droplets_QC_GE" in STEPS:
    SPECIES = config['Species'] if 'Species' in config else sys.exit("Error: No 'Species' found in configfile! Please set 'Species:' under 'Steps:'. Possibles choices are 'homo_sapiens' or 'mus_musculus'.\n")


if "Alignment_countTable_GE" in STEPS or "Alignment_countTable_ADT" in STEPS:
    if 'Sctech' in config :
        SCTECH = config['Sctech']
    else :
        SCTECH = "10xv4_3p"
        sys.stderr.write("Single-cell technology is set to '10xv4_3p'. If you would set another technology, please set 'Sctech:' under 'Steps:' \n")
    if 'WHITELISTNAME' in config :
        WHITELISTNAME = config["WHITELISTNAME"]
    elif SCTECH == '10xv4_3p' :
        WHITELISTNAME = PIPELINE_FOLDER + '/resources/WHITELISTS/3M-3pgex-may-2023_TRU.txt'
    elif SCTECH == '10xv3_5p' :
        WHITELISTNAME = PIPELINE_FOLDER + '/resources/WHITELISTS/3M-5pgex-jan-2023.txt'
    elif SCTECH == '10xv3_3p' :
        WHITELISTNAME = PIPELINE_FOLDER + '/resources/WHITELISTS/3M-february-2018.txt'
    elif SCTECH == '10xv2_5p' :
        WHITELISTNAME = PIPELINE_FOLDER + '/resources/WHITELISTS/737K-august-2016.txt'
    else :
        sys.exit("Error: Sctech doesn't exist! Only '10xv2_5p', '10xv3_3p', '10xv3_5p' and '10xv4_3p' are available.\n")

if "Alignment_countTable_GE" in STEPS or "Alignment_annotations_TCR_BCR" in STEPS:
    # Fastq-screen Index
    if 'Alignment_countTable_GE' in config and 'fastqscreen_index' in config['Alignment_countTable_GE']:
        FASTQSCREEN_INDEX = config['Alignment_countTable_GE']['fastqscreen_index']
    elif 'Alignment_annotations_TCR_BCR' in config and 'fastqscreen_index' in config['Alignment_annotations_TCR_BCR']:
        FASTQSCREEN_INDEX = config['Alignment_annotations_TCR_BCR']['fastqscreen_index']
    else :
        FASTQSCREEN_INDEX = "/mnt/beegfs02/database/bioinfo/single-cell/INDEX/FASTQ_SCREEN/0.14.0/fastq_screen.conf"

if "Alignment_countTable_GE" in STEPS:
    ### Sample/Project
    if 'Alignment_countTable_GE' in config and 'sample.name.ge' not in config['Alignment_countTable_GE']: sys.exit("Error: No sample.name.ge in configfile (Alignment_countTable_GE)!")
    if 'Alignment_countTable_GE' in config and 'input.dir.ge' not in config['Alignment_countTable_GE']: sys.exit("Error: No input.dir.ge in configfile (Alignment_countTable_GE)!")
    if 'Alignment_countTable_GE' in config and 'output.dir.ge' not in config['Alignment_countTable_GE']: sys.exit("Error: No output.dir.ge in configfile (Alignment_countTable_GE)!")
    ALIGN_SAMPLE_NAME_GE_RAW = config['Alignment_countTable_GE']['sample.name.ge']
    ALIGN_INPUT_DIR_GE_RAW = os.path.normpath(config['Alignment_countTable_GE']['input.dir.ge'])
    if not os.path.exists(ALIGN_INPUT_DIR_GE_RAW) :
        sys.exit("Error: input.dir.ge (" + ALIGN_INPUT_DIR_GE_RAW +") not found! Wrong path? \n")
    ALIGN_OUTPUT_DIR_GE = os.path.normpath(config['Alignment_countTable_GE']['output.dir.ge'])
    ALIGN_INPUT_DIR_GE = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/fastq/")
    ### Index
    if 'Alignment_countTable_GE' in config and 'kindex.ge' in config['Alignment_countTable_GE'] and 'tr2g.file.ge' in config['Alignment_countTable_GE']:
        KINDEX_GE = config['Alignment_countTable_GE']['kindex.ge']
        TR2GFILE_GE = config['Alignment_countTable_GE']['tr2g.file.ge']
        REF_TXT_GE = config['Alignment_countTable_GE']['reference.txt'] if 'reference.txt' in config['Alignment_countTable_GE'] else "<insert_your_reference_here>"
    else:
        if SPECIES == "homo_sapiens":
            KINDEX_GE = "/mnt/beegfs02/database/bioinfo/single-cell/INDEX/KALLISTO/0.51.1/homo_sapiens/GRCh38/Ensembl/r99/cDNA_LINCs_MIRs/GRCH38_r99_cDNA_linc_mir.idx"
            TR2GFILE_GE = "/mnt/beegfs02/database/bioinfo/single-cell/INDEX/KALLISTO/0.51.1/homo_sapiens/GRCh38/Ensembl/r99/cDNA_LINCs_MIRs/GRCH38_r99_cDNA_linc_mir_tr2gs.txt"
            REF_TXT_GE = "Ensembl reference transcriptome v99 corresponding to the homo sapiens GRCH38 build,"
        elif SPECIES == "mus_musculus":
            KINDEX_GE = "/mnt/beegfs02/database/bioinfo/single-cell/INDEX/KALLISTO/0.51.1/mus_musculus/GRCm38/Ensembl/r99/cDNA_LINCs_MIRs/GRCm38_r99_cDNA_linc_mir.idx"
            TR2GFILE_GE = "/mnt/beegfs02/database/bioinfo/single-cell/INDEX/KALLISTO/0.51.1/mus_musculus/GRCm38/Ensembl/r99/cDNA_LINCs_MIRs/GRCm38_r99_cDNA_linc_mir_tr2gs.txt"
            REF_TXT_GE = "Ensembl reference transcriptome v99 corresponding to the mus musculus GRCm38 build,"
        else:
            sys.exit("Error: Only 'homo_sapiens' and 'mus_musculus' are possible for Alignment_countTable_GE! If you have another species, please set 'kindex.ge' and 'tr2g.file.ge' parameters.\n")
    ### File names
    ALIGN_SAMPLE_NAME_GE = []
    ALIGN_SYMLINK_FILES_GE = []
    ALIGN_SYMLINK_FILES_NAME_GE = []
    for i in range(0,len(ALIGN_SAMPLE_NAME_GE_RAW),1):
        #check samples names and add "_GE" if needed
        ALIGN_SAMPLE_NAME_GE.append(ALIGN_SAMPLE_NAME_GE_RAW[i] + "_GE") if (ALIGN_SAMPLE_NAME_GE_RAW[i][len(ALIGN_SAMPLE_NAME_GE_RAW[i])-3:] != "_GE") else ALIGN_SAMPLE_NAME_GE.append(ALIGN_SAMPLE_NAME_GE_RAW[i])
        ORIG_FILES = glob.glob(os.path.join(ALIGN_INPUT_DIR_GE_RAW, str(ALIGN_SAMPLE_NAME_GE_RAW[i]) + "_[1-4]_S*_R1_*.f*q*")) + glob.glob(os.path.join(ALIGN_INPUT_DIR_GE_RAW, str(ALIGN_SAMPLE_NAME_GE_RAW[i]) + "_S[0-9]*_R1_*.f*q*")) + glob.glob(os.path.join(ALIGN_INPUT_DIR_GE_RAW, str(ALIGN_SAMPLE_NAME_GE_RAW[i]) + "_[1-4]_S*_R2_*.f*q*")) + glob.glob(os.path.join(ALIGN_INPUT_DIR_GE_RAW, str(ALIGN_SAMPLE_NAME_GE_RAW[i]) + "_S[0-9]*_R2_*.f*q*"))
        if (ORIG_FILES is None):
            sys.exit("Error: No fastq files found for " + ALIGN_SAMPLE_NAME_GE_RAW[i] + " sample. Wrong name?\n")
        #files with path and extention
        ALIGN_SYMLINK_FILES_GE = ALIGN_SYMLINK_FILES_GE + [ os.path.normpath(ALIGN_INPUT_DIR_GE + "/" + os.path.basename(file).replace(ALIGN_SAMPLE_NAME_GE_RAW[i], ALIGN_SAMPLE_NAME_GE[i],1)) for file in ORIG_FILES]
    #files without path and extention
    ALIGN_SYMLINK_FILES_NAME_GE = [os.path.splitext(os.path.splitext(os.path.basename(x))[0])[0] for x in ALIGN_SYMLINK_FILES_GE]

if "Alignment_countTable_ADT" in STEPS:
    ### Sample/Project
    if 'Alignment_countTable_ADT' in config and 'sample.name.adt' not in config['Alignment_countTable_ADT']: sys.exit("Error: No sample.name.adt in configfile (Alignment_countTable_ADT)!")
    if 'Alignment_countTable_ADT' in config and 'input.dir.adt' not in config['Alignment_countTable_ADT']: sys.exit("Error: No input.dir.adt in configfile (Alignment_countTable_ADT)!")
    if 'Alignment_countTable_ADT' in config and 'output.dir.adt' not in config['Alignment_countTable_ADT']: sys.exit("Error: No output.dir.adt in configfile (Alignment_countTable_ADT)!")
    ALIGN_SAMPLE_NAME_ADT_RAW = config['Alignment_countTable_ADT']['sample.name.adt']
    ALIGN_INPUT_DIR_ADT_RAW = os.path.normpath(config['Alignment_countTable_ADT']['input.dir.adt'])
    if not os.path.exists(ALIGN_INPUT_DIR_ADT_RAW) :
        sys.exit("Error: input.dir.adt not found! Wrong path? \n")
    ALIGN_OUTPUT_DIR_ADT = os.path.normpath(config['Alignment_countTable_ADT']['output.dir.adt'])
    ALIGN_INPUT_DIR_ADT = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/fastq/")
    ### Index
    FEATURES_FA_ADT = config['Alignment_countTable_ADT']['features.fa.adt'] if 'Alignment_countTable_ADT' in config and 'features.fa.adt' in config['Alignment_countTable_ADT'] else sys.exit("Error: No features.fa.adt in configfile (Alignment_countTable_ADT)!")
    ### File names
    ALIGN_SAMPLE_NAME_ADT = []
    ALIGN_SYMLINK_FILES_ADT = []
    ALIGN_SYMLINK_FILES_NAME_ADT = []
    for i in range(0,len(ALIGN_SAMPLE_NAME_ADT_RAW),1):
        #check samples names and add "_ADT" if needed
        ALIGN_SAMPLE_NAME_ADT.append(ALIGN_SAMPLE_NAME_ADT_RAW[i] + "_ADT") if (ALIGN_SAMPLE_NAME_ADT_RAW[i][len(ALIGN_SAMPLE_NAME_ADT_RAW[i])-4:] != "_ADT") else ALIGN_SAMPLE_NAME_ADT.append(ALIGN_SAMPLE_NAME_ADT_RAW[i])
        ORIG_FILES = glob.glob(os.path.join(ALIGN_INPUT_DIR_ADT_RAW, str(ALIGN_SAMPLE_NAME_ADT_RAW[i]) + "*_R1_*.f*q*")) + glob.glob(os.path.join(ALIGN_INPUT_DIR_ADT_RAW, str(ALIGN_SAMPLE_NAME_ADT_RAW[i]) + "*_R2_*.f*q*"))
        if (ORIG_FILES is None):
            sys.exit("Error: No fastq files found for " + ALIGN_SAMPLE_NAME_ADT_RAW[i] + " sample. Wrong name?\n")
        #files with path and extention
        ALIGN_SYMLINK_FILES_ADT = ALIGN_SYMLINK_FILES_ADT + [ os.path.normpath(ALIGN_INPUT_DIR_ADT + "/" + os.path.basename(file).replace(ALIGN_SAMPLE_NAME_ADT_RAW[i], ALIGN_SAMPLE_NAME_ADT[i],1)) for file in ORIG_FILES]
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
    if ALIGN_INPUT_DIR_TCR_RAW is not None and not os.path.exists(ALIGN_INPUT_DIR_TCR_RAW) :
        sys.exit("Error: input.dir.tcr not found! Wrong path? \n")
    ALIGN_INPUT_DIR_BCR_RAW = os.path.normpath(config['Alignment_annotations_TCR_BCR']['input.dir.bcr'] + "/") if 'input.dir.bcr' in config['Alignment_annotations_TCR_BCR'] else None
    if ALIGN_INPUT_DIR_BCR_RAW is not None and not os.path.exists(ALIGN_INPUT_DIR_BCR_RAW) :
        sys.exit("Error: input.dir.bcr not found! Wrong path? \n")
    ALIGN_OUTPUT_DIR_TCR_BCR = os.path.normpath(config['Alignment_annotations_TCR_BCR']['output.dir.tcr_bcr'])
    ALIGN_INPUT_DIR_TCR_BCR = os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/fastq/")
    ### Index
    if SPECIES == "homo_sapiens":
        CRINDEX_TCR_BCR="/references/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0"
    elif SPECIES == "mus_musculus":
        CRINDEX_TCR_BCR="/references/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0"
    else:
        sys.exit("Error: Only 'homo_sapiens' and 'mus_musculus' are possible for Alignment_annotations_TCR_BCR!\n")
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
            if (ORIG_FILES is None):
                sys.exit("Error: No fastq files found for " + ALIGN_SAMPLE_NAME_TCR_RAW[i] + " sample. Wrong name?\n")
            #files with path and extention
            ALIGN_ORIG_FILES_TCR = ALIGN_ORIG_FILES_TCR + ORIG_FILES
            for file in ORIG_FILES:
                if re.match(str(ALIGN_SAMPLE_NAME_TCR_RAW[i] + "_S[0-9]+_L00[0-9]{1}_R[1-2]{1}_.*"), os.path.basename(file)) is not None: #good name format
                    ALIGN_SYMLINK_FILES_TCR = ALIGN_SYMLINK_FILES_TCR + [ os.path.normpath(ALIGN_INPUT_DIR_TCR_BCR + "/" + os.path.basename(file).replace(ALIGN_SAMPLE_NAME_TCR_RAW[i], ALIGN_SAMPLE_NAME_TCR[i],1)) ]
                elif re.match(str(ALIGN_SAMPLE_NAME_TCR_RAW[i] + "_[0-9]{1}_S[0-9]+_R[1-2]{1}_.*"), os.path.basename(file)) is not None: # => reformat
                    res_match = re.match(str(ALIGN_SAMPLE_NAME_TCR_RAW[i] + "_(?P<nb>[0-9]{1})_(?P<S>S[0-9]+)_(?P<R_compl>R[1-2]{1}_.*)"), os.path.basename(file))
                    ALIGN_SYMLINK_FILES_TCR = ALIGN_SYMLINK_FILES_TCR + [ os.path.normpath(str(ALIGN_INPUT_DIR_TCR_BCR + "/" + ALIGN_SAMPLE_NAME_TCR[i] + "_" + res_match.group('S') + "_L00" + res_match.group('nb') + "_" + res_match.group('R_compl'))) ]
                else:
                    sys.exit("File names for TCR not recognized. It must be like mysample_2_S1_R1_001.fastq.gz or mysample_S1_L002_R1_001.fastq.gz")
        #files without path and extention
        ALIGN_SYMLINK_FILES_NAME_TCR = [os.path.splitext(os.path.splitext(os.path.basename(x))[0])[0] for x in ALIGN_SYMLINK_FILES_TCR]
    else:
        ALIGN_SAMPLE_NAME_TCR_RAW = []
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
            if (ORIG_FILES is None):
                sys.exit("Error: No fastq files found for " + ALIGN_SAMPLE_NAME_BCR_RAW[i] + " sample. Wrong name?\n")
            #files with path and extention
            ALIGN_ORIG_FILES_BCR = ALIGN_ORIG_FILES_BCR + ORIG_FILES
            for file in ORIG_FILES:
                if re.match(str(ALIGN_SAMPLE_NAME_BCR_RAW[i] + "_S[0-9]+_L00[0-9]{1}_R[1-2]{1}_.*"), os.path.basename(file)) is not None: #good name format
                    ALIGN_SYMLINK_FILES_BCR = ALIGN_SYMLINK_FILES_BCR + [ os.path.normpath(ALIGN_INPUT_DIR_TCR_BCR + "/" + os.path.basename(file).replace(ALIGN_SAMPLE_NAME_BCR_RAW[i], ALIGN_SAMPLE_NAME_BCR[i],1)) ]
                elif re.match(str(ALIGN_SAMPLE_NAME_BCR_RAW[i] + "_[0-9]{1}_S[0-9]+_R[1-2]{1}_.*"), os.path.basename(file)) is not None: # => reformat
                    res_match = re.match(str(ALIGN_SAMPLE_NAME_BCR_RAW[i] + "_(?P<nb>[0-9]{1})_(?P<S>S[0-9]+)_(?P<R_compl>R[1-2]{1}_.*)"), os.path.basename(file))
                    ALIGN_SYMLINK_FILES_BCR = ALIGN_SYMLINK_FILES_BCR + [ os.path.normpath(str(ALIGN_INPUT_DIR_TCR_BCR + "/" + ALIGN_SAMPLE_NAME_BCR[i] + "_" + res_match.group('S') + "_L00" + res_match.group('nb') + "_" + res_match.group('R_compl'))) ]
                else:
                    sys.exit("File names for BCR not recognized. It must be like mysample_2_S1_R1_001.fastq.gz or mysample_S1_L002_R1_001.fastq.gz")
        #files without path and extention
        ALIGN_SYMLINK_FILES_NAME_BCR = [os.path.splitext(os.path.splitext(os.path.basename(x))[0])[0] for x in ALIGN_SYMLINK_FILES_BCR]
    else:
        ALIGN_SAMPLE_NAME_BCR_RAW = []
    #Fusion TCR/BCR
    ALIGN_SAMPLE_NAME_TCR_BCR_RAW = ALIGN_SAMPLE_NAME_TCR_RAW  + ALIGN_SAMPLE_NAME_BCR_RAW
    ALIGN_SAMPLE_NAME_TCR_BCR = ALIGN_SAMPLE_NAME_TCR + ALIGN_SAMPLE_NAME_BCR
    ALIGN_ORIG_FILES_TCR_BCR = ALIGN_ORIG_FILES_TCR + ALIGN_ORIG_FILES_BCR
    ALIGN_SYMLINK_FILES_TCR_BCR = ALIGN_SYMLINK_FILES_TCR + ALIGN_SYMLINK_FILES_BCR
    ALIGN_SYMLINK_FILES_NAME_TCR_BCR = ALIGN_SYMLINK_FILES_NAME_TCR + ALIGN_SYMLINK_FILES_NAME_BCR

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
    ### Analysis Parameters
    # Emptydrops
    QC_EMPTYDROPS_FDR = getelseDefault('Droplets_QC_GE', 'emptydrops.fdr', 'NULL')
    QC_DROPLETS_LIMIT = getelseDefault('Droplets_QC_GE', 'droplets.limit', 'NULL')
    QC_EMPTYDROPS_RETAIN = getelseDefault('Droplets_QC_GE', 'emptydrops.retain', 'NULL')
    # QC cell
    QC_PCMITO_MIN = getelseDefault('Droplets_QC_GE', 'pcmito.min', 'NULL')
    QC_PCMITO_MAX = getelseDefault('Droplets_QC_GE', 'pcmito.max', 'NULL')
    QC_PCRIBO_MIN = getelseDefault('Droplets_QC_GE', 'pcribo.min', 'NULL')
    QC_PC_RIBO_MAX = getelseDefault('Droplets_QC_GE', 'pcribo.max', 'NULL')
    QC_MIN_FEATURES = getelseDefault('Droplets_QC_GE', 'min.features', 'NULL')
    QC_MIN_COUNTS = getelseDefault('Droplets_QC_GE', 'min.counts', 'NULL')
    # QC gene
    QC_MIN_CELLS = getelseDefault('Droplets_QC_GE', 'min.cells', 'NULL')
    ### Databases
    # Metadata file
    QC_METADATA_FILE = getelseDefault('Droplets_QC_GE', 'metadata.file', 'NULL').replace(", ", ",")
    # QC
    QC_MT_FILE = getelseDefault('Droplets_QC_GE', 'mt.genes.file', 'NULL')
    QC_RB_FILE = getelseDefault('Droplets_QC_GE', 'crb.genes.file', 'NULL')
    QC_ST_FILE = getelseDefault('Droplets_QC_GE', 'str.genes.file', 'NULL')
    ### Snakefile parameters
    #check end paths (del "/" if necessary)
    for i in range(0,len(QC_INPUT_DIR_GE),1):
        QC_INPUT_DIR_GE[i] = os.path.normpath(QC_INPUT_DIR_GE[i])
        QC_OUTPUT_DIR_GE[i] = os.path.normpath(QC_OUTPUT_DIR_GE[i])
    #Correspondance sample/input/output
    dic_SAMPLE_NAME_GE_INFO = {}
    for i in range(0,len(QC_INPUT_DIR_GE),1):
            dic_SAMPLE_NAME_GE_INFO[QC_INPUT_DIR_GE[i]] = {}
            dic_SAMPLE_NAME_GE_INFO[QC_INPUT_DIR_GE[i]]['QC_OUTPUT_DIR'] = QC_OUTPUT_DIR_GE[i]
            dic_SAMPLE_NAME_GE_INFO[QC_INPUT_DIR_GE[i]]['SAMPLE_NAME_GE'] = QC_SAMPLE_NAME_GE[i]

if "Filtering_GE" in STEPS:
    ### Sample/Project
    if 'Filtering_GE' in config and 'sample.name.ge' in config['Filtering_GE'] and 'input.rda.ge' in config['Filtering_GE'] :
        FILTERING_SAMPLE_NAME_GE_RAW = config['Filtering_GE']['sample.name.ge']
        FILTERING_INPUT_RDA_GE = config['Filtering_GE']['input.rda.ge']
        #check samples names and add "_GE" if needed
        FILTERING_SAMPLE_NAME_GE = []
        for i in range(0,len(FILTERING_SAMPLE_NAME_GE_RAW),1):
            FILTERING_SAMPLE_NAME_GE.append(FILTERING_SAMPLE_NAME_GE_RAW[i] + "_GE") if (FILTERING_SAMPLE_NAME_GE_RAW[i][len(FILTERING_SAMPLE_NAME_GE_RAW[i])-3:] != "_GE") else FILTERING_SAMPLE_NAME_GE.append(FILTERING_SAMPLE_NAME_GE_RAW[i])
    elif "Droplets_QC_GE" in STEPS:
        sys.stderr.write("Note: No input.rda.ge find in Filtering_GE section of configfile; input.rda.ge will be determine from Droplets_QC_GE step for Filtering_GE step!\n")
        FILTERING_SAMPLE_NAME_GE = copy.deepcopy(QC_SAMPLE_NAME_GE)
        FILTERING_INPUT_RDA_GE = [os.path.normpath(dic_SAMPLE_NAME_GE_INFO[x]['QC_OUTPUT_DIR'] + ("/QC_droplets/" if str(QC_EMPTYDROPS_RETAIN) == "NULL" else "/QC_droplets_retain" + str(QC_EMPTYDROPS_RETAIN) + "/") + str(dic_SAMPLE_NAME_GE_INFO[x]['SAMPLE_NAME_GE']) + "_QC_NON-NORMALIZED.rda") for x in QC_INPUT_DIR_GE]
    else:
        sys.exit("Error: No sample.name.ge or/and input.rda.ge in configfile!\n")
    if 'Filtering_GE' in config and 'output.dir.ge' in config['Filtering_GE'] :
        FILTERING_OUTPUT_DIR_GE = [os.path.normpath(x) for x in config['Filtering_GE']['output.dir.ge']]
    elif "Droplets_QC_GE" in STEPS:
        FILTERING_OUTPUT_DIR_GE = copy.deepcopy(QC_OUTPUT_DIR_GE)
        sys.stderr.write("Note: No output.dir.ge find in Filtering_GE section of configfile; output.dir.ge will be determine from Droplets_QC_GE step for Filtering_GE step!\n")
    else :
        sys.exit("Error: No output.dir.ge find in configfile!\n")
    ### Analysis Parameters
    # QC cell
    FILTERING_PCMITO_MIN = getelseDefault('Filtering_GE', 'pcmito.min', '0')
    FILTERING_PCMITO_MAX = getelseDefault('Filtering_GE', 'pcmito.max', '20')
    FILTERING_PCRIBO_MIN = getelseDefault('Filtering_GE', 'pcribo.min', '0')
    FILTERING_PC_RIBO_MAX = getelseDefault('Filtering_GE', 'pcribo.max', '100')
    FILTERING_MIN_FEATURES = getelseDefault('Filtering_GE', 'min.features', '200')
    FILTERING_MIN_COUNTS = getelseDefault('Filtering_GE', 'min.counts', '1000')
    # QC gene
    FILTERING_MIN_CELLS = getelseDefault('Filtering_GE', 'min.cells', '5')
    # Doublets
    FILTERING_DOUBLET_FILTER_METHOD = getelseDefault('Filtering_GE', 'doublets.filter.method', 'NULL')
    ### Databases
    # Metadata file
    FILTERING_METADATA_FILE = getelseDefault('Filtering_GE', 'metadata.file', 'NULL')
    # QC
    FILTERING_CC_SEURAT_FILE = getelseDefault('Filtering_GE', 'cc.seurat.file', 'NULL')
    FILTERING_CC_CYCLONE_FILE = getelseDefault('Filtering_GE', 'cc.cyclone.file', 'NULL')
    ### Snakefile parameters
    #Correspondance sample/input/output
    dic_FILTER_INFO = {}
    for i in range(0,len(FILTERING_INPUT_RDA_GE),1):
        dic_FILTER_INFO[FILTERING_INPUT_RDA_GE[i]] = {}
        dic_FILTER_INFO[FILTERING_INPUT_RDA_GE[i]]['FILTERING_SAMPLE_NAME_GE'] = FILTERING_SAMPLE_NAME_GE[i]
        dic_FILTER_INFO[FILTERING_INPUT_RDA_GE[i]]['FILTER_OUTPUT_DIR'] = FILTERING_OUTPUT_DIR_GE[i]
    FILTERS_FOLDER = "F" + str(FILTERING_MIN_FEATURES) + "_C" + str(FILTERING_MIN_COUNTS) + "_M" + str(FILTERING_PCMITO_MIN) + "-" + str(FILTERING_PCMITO_MAX) + "_R" + str(FILTERING_PCRIBO_MIN) + "-" + str(FILTERING_PC_RIBO_MAX) + "_G" + str(FILTERING_MIN_CELLS)
    #name of the doublets identification method
    FILTERING_DOUBLET_FILTER_METHOD_NAME = "all" if FILTERING_DOUBLET_FILTER_METHOD == "NULL" else FILTERING_DOUBLET_FILTER_METHOD

if "Norm_DimRed_Eval_GE" in STEPS:
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
        NDRE_SAMPLE_NAME_GE = copy.deepcopy(FILTERING_SAMPLE_NAME_GE)
        if FILTERING_DOUBLET_FILTER_METHOD_NAME == "none":
            NDRE_INPUT_RDA_GE = [os.path.normpath(dic_FILTER_INFO[x]['FILTER_OUTPUT_DIR'] + "/" + FILTERS_FOLDER + "/DOUBLETSKEPT/" + str(dic_FILTER_INFO[x]['FILTERING_SAMPLE_NAME_GE']) + "_DOUBLETSKEPT_NON-NORMALIZED.rda") for x in FILTERING_INPUT_RDA_GE]
        else:
            NDRE_INPUT_RDA_GE = [os.path.normpath(dic_FILTER_INFO[x]['FILTER_OUTPUT_DIR'] + "/" + FILTERS_FOLDER + "/DOUBLETSFILTER_" + FILTERING_DOUBLET_FILTER_METHOD_NAME + "/" + str(dic_FILTER_INFO[x]['FILTERING_SAMPLE_NAME_GE']) + "_FILTERED_NON-NORMALIZED.rda") for x in FILTERING_INPUT_RDA_GE]
    else:
        sys.exit("Error: No sample.name.ge or/and input.rda.ge in configfile!\n")
    if ('Norm_DimRed_Eval_GE' in config) and ('output.dir.ge' in config['Norm_DimRed_Eval_GE']) :
        NDRE_OUTPUT_DIR_GE = [os.path.normpath(x) for x in config['Norm_DimRed_Eval_GE']['output.dir.ge']]
    elif "Filtering_GE" in STEPS:
        NDRE_OUTPUT_DIR_GE = [dic_FILTER_INFO[x]['FILTER_OUTPUT_DIR'] + "/" + FILTERS_FOLDER + ("/DOUBLETSKEPT" if FILTERING_DOUBLET_FILTER_METHOD_NAME == "none" else ("/DOUBLETSFILTER_" + FILTERING_DOUBLET_FILTER_METHOD_NAME)) for x in FILTERING_INPUT_RDA_GE]
        sys.stderr.write("Note: No output.dir.ge find in Norm_DimRed_Eval_GE section of configfile; output.dir.ge will be determine from Filtering_GE step for Norm_DimRed_Eval_GE step!\n")
    else :
        sys.exit("Error: No output.dir.ge find in configfile!\n")
    ### Analysis Parameters
    NDRE_EVAL_MARKERS = getelseDefault('Norm_DimRed_Eval_GE', 'eval.markers', 'NULL').replace(", ", ",")
    # Normalization and dimension reduction
    NDRE_FEATURES_N = getelseDefault('Norm_DimRed_Eval_GE', 'features.n', 'NULL')
    NDRE_NORM_METHOD = getelseDefault('Norm_DimRed_Eval_GE', 'norm.method', 'SCTransform')
    NDRE_HVG_FINDVARIABLEFEATURESMIX = getelseDefault('Norm_DimRed_Eval_GE', 'HVG.FindVariableFeaturesMix', 'FALSE')
    NDRE_DIMRED_METHOD = getelseDefault('Norm_DimRed_Eval_GE', 'dimred.method', 'pca')
    #NDRE_VTR_BIASES_NORM = [ x.replace(", ", ",") for x in config['Norm_DimRed_Eval_GE']['vtr.biases.norm']] if ('Norm_DimRed_Eval_GE' in config and 'vtr.biases.norm' in config['Norm_DimRed_Eval_GE'] and config['Norm_DimRed_Eval_GE']['vtr.biases.norm'] != None) else ["NULL"]
    NDRE_VTR_BIASES_NORM = getelseDefault('Norm_DimRed_Eval_GE', 'vtr.biases.norm', ['NULL'])
    if not isinstance(NDRE_VTR_BIASES_NORM, list): NDRE_VTR_BIASES_NORM = list(NDRE_VTR_BIASES_NORM)
    NDRE_VTR_BIASES_NORM = ["NULL" if x == "" else x.replace(", ", ",") for x in NDRE_VTR_BIASES_NORM]
    #NDRE_VTR_BIASES_DIMRED = [ x.replace(", ", ",") for x in config['Norm_DimRed_Eval_GE']['vtr.biases.dimred']] if ('Norm_DimRed_Eval_GE' in config and 'vtr.biases.dimred' in config['Norm_DimRed_Eval_GE'] and config['Norm_DimRed_Eval_GE']['vtr.biases.dimred'] != None) else ["NULL"]
    NDRE_VTR_BIASES_DIMRED = getelseDefault('Norm_DimRed_Eval_GE', 'vtr.biases.dimred', ['NULL'])
    if not isinstance(NDRE_VTR_BIASES_DIMRED, list): NDRE_VTR_BIASES_DIMRED = list(NDRE_VTR_BIASES_DIMRED)
    NDRE_VTR_BIASES_DIMRED = ["NULL" if x == "" else x.replace(", ", ",") for x in NDRE_VTR_BIASES_DIMRED]
    NDRE_VTR_SCALE = getelseDefault('Norm_DimRed_Eval_GE', 'vtr.scale', 'NULL')
    NDRE_REGEX_REMOVE_HVG = getelseDefault('Norm_DimRed_Eval_GE', 'regex.genes.to.remove.from.HVG', 'NULL')
    NDRE_DIM_MAX = getelseDefault('Norm_DimRed_Eval_GE', 'dims.max', 49)
    NDRE_SKIP_EVALDIMRES = getelseDefault('Norm_DimRed_Eval_GE', 'skip.eval_dims_res', 'NULL')
    NDRE_EVAL_DIM_MAX = getelseDefault('Norm_DimRed_Eval_GE', 'eval.dims.max', 49)
    NDRE_EVAL_DIM_MIN = getelseDefault('Norm_DimRed_Eval_GE', 'eval.dims.min', 9)
    NDRE_EVAL_DIM_STEPS = getelseDefault('Norm_DimRed_Eval_GE', 'eval.dims.steps', 2)
    NDRE_EVAL_RES_MAX = getelseDefault('Norm_DimRed_Eval_GE', 'eval.res.max', 1.2)
    NDRE_EVAL_RES_MIN = getelseDefault('Norm_DimRed_Eval_GE', 'eval.res.min', 0.1)
    NDRE_EVAL_RES_STEPS = getelseDefault('Norm_DimRed_Eval_GE', 'eval.res.steps', 0.1)
    NDRE_EVAL_PTSIZE = getelseDefault('Norm_DimRed_Eval_GE', 'eval.pt.size', 'NULL').replace(", ", ",")
    # Metadata file
    NDRE_METADATA_FILE = getelseDefault('Norm_DimRed_Eval_GE', 'metadata.file', 'NULL')
    ### Snakefile parameters
    #Correspondance sample/input/output
    dic_NDRE_INFO = {}
    for i in range(0,len(NDRE_INPUT_RDA_GE),1):
        dic_NDRE_INFO[NDRE_INPUT_RDA_GE[i]] = {}
        dic_NDRE_INFO[NDRE_INPUT_RDA_GE[i]]['NDRE_SAMPLE_NAME_GE'] = NDRE_SAMPLE_NAME_GE[i]
        dic_NDRE_INFO[NDRE_INPUT_RDA_GE[i]]['NDRE_OUTPUT_DIR'] = NDRE_OUTPUT_DIR_GE[i]
    #Names
    #NDRE_NORM_VTR = [NDRE_NORM_METHOD] if (NDRE_VTR_BIASES_NORM == "NULL") else ([NDRE_NORM_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(x.split(","))))) for x in NDRE_VTR_BIASES_NORM])
    NDRE_NORM_VTR = [NDRE_NORM_METHOD if (x == "NULL") else (NDRE_NORM_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(x.split(",")))))) for x in NDRE_VTR_BIASES_NORM]
    #NDRE_DIMRED_VTR = [NDRE_DIMRED_METHOD] if (NDRE_DIMRED_METHOD == "pca" or NDRE_DIMRED_METHOD == "ica" or NDRE_DIMRED_METHOD == "mds" or NDRE_VTR_BIASES_DIMRED == "NULL") else ([NDRE_DIMRED_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(x.split(","))))) for x in NDRE_VTR_BIASES_DIMRED])
    NDRE_DIMRED_VTR = [NDRE_DIMRED_METHOD if (x == "NULL") else (NDRE_DIMRED_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(x.split(",")))))) for x in NDRE_VTR_BIASES_DIMRED]

if "Clust_Markers_Annot_GE" in STEPS: #alias CMA
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
        CMA_SAMPLE_NAME_GE = []
        CMA_INPUT_RDA_GE = []
        for i in NDRE_INPUT_RDA_GE:
            for j in NDRE_NORM_VTR:
                for k in NDRE_DIMRED_VTR:
                    CMA_SAMPLE_NAME_GE.append(dic_NDRE_INFO[i]['NDRE_SAMPLE_NAME_GE'])
                    CMA_INPUT_RDA_GE.append(os.path.normpath(dic_NDRE_INFO[i]['NDRE_OUTPUT_DIR'] + "/" + j + "/" + k + "/" + dic_NDRE_INFO[i]['NDRE_SAMPLE_NAME_GE'] + "_" + j + "_" + k + ".rda"))
    else:
        sys.exit("Error: No sample.name.ge or/and input.rda.ge in configfile!\n")
    if ('Clust_Markers_Annot_GE' in config) and ('output.dir.ge' in config['Clust_Markers_Annot_GE']) :
        CMA_OUTPUT_DIR_GE = config['Clust_Markers_Annot_GE']['output.dir.ge']
    elif "Norm_DimRed_Eval_GE" in STEPS:
        CMA_OUTPUT_DIR_GE = [os.path.dirname(x) for x in CMA_INPUT_RDA_GE]
        sys.stderr.write("Note: No output.dir.ge find in Clust_Markers_Annot_GE section of configfile; output.dir.ge will be determine from Norm_DimRed_Eval_GE step for Clust_Markers_Annot_GE step!\n")
    else :
        sys.exit("Error: No output.dir.ge find in configfile!\n")
    ### Analysis Parameters
    CMA_MARKFILE = getelseDefault('Clust_Markers_Annot_GE', 'markfile', 'NULL').replace(", ", ",")
    CMA_MARKERS_PTSIZE = getelseDefault('Clust_Markers_Annot_GE', 'markers.pt.size', 'NULL').replace(", ", ",")
    CMA_MARKERS_ORDER = getelseDefault('Clust_Markers_Annot_GE', 'markers.order', 'NULL').replace(", ", ",")
    # Normalization and dimension reduction
    CMA_KEEP_DIM_RES = getelseDefault('Clust_Markers_Annot_GE', 'keep.dim.res', 'NULL')
    CMA_KEEP_DIM_RES = [str(dim_res).replace(".0", "").replace(",0","") for dim_res in CMA_KEEP_DIM_RES]
    # Annotation
    CMA_CFR_MINSCORE = getelseDefault('Clust_Markers_Annot_GE', 'cfr.minscore', 'NULL')
    CMA_SR_MINSCORE = getelseDefault('Clust_Markers_Annot_GE', 'sr.minscore', 'NULL')
    CMA_CUSTOM_SCE_REF = getelseDefault('Clust_Markers_Annot_GE', 'custom.sce.ref', 'NULL').replace(", ", ",")
    CMA_CUSTOM_MARKERS_REF = getelseDefault('Clust_Markers_Annot_GE', 'custom.markers.ref', 'NULL').replace(", ", ",")
    # SKIP
    CMA_SKIP_TECHNICAL = getelseDefault('Clust_Markers_Annot_GE', 'skip.technical_plots', 'NULL')
    CMA_SKIP_ANNOT = getelseDefault('Clust_Markers_Annot_GE', 'skip.annotation', 'NULL')
    CMA_SKIP_MARKERS_IDENT = getelseDefault('Clust_Markers_Annot_GE', 'skip.markers_identification', 'NULL')
    # Metadata file
    CMA_METADATA_FILE = getelseDefault('Clust_Markers_Annot_GE', 'metadata.file', 'NULL').replace(", ", ",")
    ### Snakefile parameters
    #check end paths (add "/" if necessary)
    for i in range(0,len(CMA_OUTPUT_DIR_GE),1):
        CMA_OUTPUT_DIR_GE[i] = os.path.normpath(CMA_OUTPUT_DIR_GE[i])
    #Correspondance sample/input/output
    dic_CMA_INFO = {}
    CMA_COMPLEMENT = []
    for i in range(0,len(CMA_INPUT_RDA_GE),1):
        dic_CMA_INFO[CMA_INPUT_RDA_GE[i]] = {}
        dic_CMA_INFO[CMA_INPUT_RDA_GE[i]]['CMA_SAMPLE_NAME_GE'] = CMA_SAMPLE_NAME_GE[i]
        dic_CMA_INFO[CMA_INPUT_RDA_GE[i]]['CMA_OUTPUT_DIR'] = CMA_OUTPUT_DIR_GE[i]
        compl = os.path.splitext(os.path.basename(CMA_INPUT_RDA_GE[i]))[0]
        if compl.startswith(CMA_SAMPLE_NAME_GE[i]):
            compl = compl[len(CMA_SAMPLE_NAME_GE[i]):]
        CMA_COMPLEMENT.append(compl)
    #Names
    CMA_CLUST_FOLDERS = [("dims" + str(dim_res).replace("_", "_res")) for dim_res in CMA_KEEP_DIM_RES]

if "Adding_ADT" in STEPS:
    ### Sample/Project
    if 'Adding_ADT' in config and 'input.rda.ge' in config['Adding_ADT'] :
        ADD_ADT_INPUT_RDA_GE = config['Adding_ADT']['input.rda.ge']
    elif "Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("Note: No input.rda.ge find in Adding_ADT section of configfile; input.rda.ge will be determine from Clust_Markers_Annot_GE step for Adding_ADT step!\n")
        ADD_ADT_INPUT_RDA_GE = []
        ADD_ADT_SAMPLE_NAME_GE = []
        for i in range(len(CMA_INPUT_RDA_GE)):
            for j in range(len(CMA_CLUST_FOLDERS)):
                ADD_ADT_SAMPLE_NAME_GE.append(dic_CMA_INFO[CMA_INPUT_RDA_GE[i]]['CMA_SAMPLE_NAME_GE'])
                ADD_ADT_INPUT_RDA_GE.append(os.path.normpath(dic_CMA_INFO[CMA_INPUT_RDA_GE[i]]['CMA_OUTPUT_DIR'] + "/" + CMA_CLUST_FOLDERS[j] + "/" + dic_CMA_INFO[CMA_INPUT_RDA_GE[i]]['CMA_SAMPLE_NAME_GE'] + CMA_COMPLEMENT[i] + "_" + CMA_KEEP_DIM_RES[j] + ".rda"))
    else:
        sys.exit("Error: No input.rda.ge in configfile!\n")
    if 'Adding_ADT' in config and 'sample.name.adt' in config['Adding_ADT'] and 'input.dir.adt' in config['Adding_ADT'] :
        ADD_ADT_INPUT_DIR_ADT = config['Adding_ADT']['input.dir.adt']
        ADD_ADT_SAMPLE_NAME_ADT_RAW = config['Adding_ADT']['sample.name.adt']
        #check samples names and add "_ADT" if needed
        ADD_ADT_SAMPLE_NAME_ADT = []
        for i in range(0,len(ADD_ADT_SAMPLE_NAME_ADT_RAW),1):
            ADD_ADT_SAMPLE_NAME_ADT.append(ADD_ADT_SAMPLE_NAME_ADT_RAW[i] + "_ADT") if (ADD_ADT_SAMPLE_NAME_ADT_RAW[i][len(ADD_ADT_SAMPLE_NAME_ADT_RAW[i])-4:] != "_ADT") else ADD_ADT_SAMPLE_NAME_ADT.append(ADD_ADT_SAMPLE_NAME_ADT_RAW[i])
    elif "Alignment_countTable_ADT" in STEPS:
        sys.stderr.write("Note: No sample.name.adt or input.dir.adt find in Adding_ADT section of configfile; sample.name.adt and input.dir.adt will be determine from Alignment_countTable_ADT step for Adding_ADT step!\n")
        ADD_ADT_INPUT_DIR_ADT = [ os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/" + str(x) + "/KALLISTOBUS") for x in ALIGN_SAMPLE_NAME_ADT]
        ADD_ADT_SAMPLE_NAME_ADT = copy.deepcopy(ALIGN_SAMPLE_NAME_ADT)
    else:
        sys.exit("Error: No sample.name.adt or input.dir.adt in configfile!\n")
    ### Analysis Parameters
    ADD_ADT_GENE_NAMES = getelseDefault('Adding_ADT', 'gene.names', 'NULL').replace(", ", ",")
    ### Snakefile parameters
    #Correspondance between ADT sample names, input folder of ADT and GE sample names in yaml
    if "Alignment_countTable_GE" in STEPS:
        TMP_SAMPLE_NAME_GE = copy.deepcopy(ALIGN_SAMPLE_NAME_GE)
    elif "Droplets_QC_GE" in STEPS:
        TMP_SAMPLE_NAME_GE = copy.deepcopy(QC_SAMPLE_NAME_GE)
    elif "Filtering_GE" in STEPS:
        TMP_SAMPLE_NAME_GE = copy.deepcopy(FILTERING_SAMPLE_NAME_GE)
    elif "Norm_DimRed_Eval_GE" in STEPS:
        TMP_SAMPLE_NAME_GE = copy.deepcopy(NDRE_SAMPLE_NAME_GE)
    elif "Clust_Markers_Annot_GE" in STEPS:
        TMP_SAMPLE_NAME_GE = copy.deepcopy(CMA_SAMPLE_NAME_GE)
    if "Alignment_countTable_GE" in STEPS or "Droplets_QC_GE" in STEPS or "Filtering_GE" in STEPS or "Norm_DimRed_Eval_GE" in STEPS or "Clust_Markers_Annot_GE" in STEPS:
        dic_ADD_ADT_CORRES_SAMPLE_NAMES = {}
        for i in range(len(TMP_SAMPLE_NAME_GE)):
            dic_ADD_ADT_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_GE[i]] = {}
            dic_ADD_ADT_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_GE[i]]['ADD_ADT_INPUT_DIR_ADT'] = ADD_ADT_INPUT_DIR_ADT[i]
            dic_ADD_ADT_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_GE[i]]['ADD_ADT_SAMPLE_NAME_ADT'] = ADD_ADT_SAMPLE_NAME_ADT[i]
        #Correspondance input/output (output = folder and file name without the extention "_ADT.rda")
        ADD_ADT_OUTPUT = [os.path.splitext(x)[0] for x in ADD_ADT_INPUT_RDA_GE]
        dic_ADD_ADT_INFO = {}
        for i in range(len(ADD_ADT_OUTPUT)):
            dic_ADD_ADT_INFO[ADD_ADT_OUTPUT[i]] = {}
            dic_ADD_ADT_INFO[ADD_ADT_OUTPUT[i]]['ADD_ADT_INPUT_RDA_GE'] = ADD_ADT_INPUT_RDA_GE[i]
            dic_ADD_ADT_INFO[ADD_ADT_OUTPUT[i]]['ADD_ADT_SAMPLE_NAME_GE'] = ADD_ADT_SAMPLE_NAME_GE[i] #
            dic_ADD_ADT_INFO[ADD_ADT_OUTPUT[i]]['ADD_ADT_INPUT_DIR_ADT'] = dic_ADD_ADT_CORRES_SAMPLE_NAMES[ADD_ADT_SAMPLE_NAME_GE[i]]['ADD_ADT_INPUT_DIR_ADT']
            dic_ADD_ADT_INFO[ADD_ADT_OUTPUT[i]]['ADD_ADT_SAMPLE_NAME_ADT'] = dic_ADD_ADT_CORRES_SAMPLE_NAMES[ADD_ADT_SAMPLE_NAME_GE[i]]['ADD_ADT_SAMPLE_NAME_ADT']
    else:
        ADD_ADT_OUTPUT = [os.path.splitext(x)[0] for x in ADD_ADT_INPUT_RDA_GE]
        dic_ADD_ADT_INFO = {}
        for i in range(len(ADD_ADT_OUTPUT)):
            dic_ADD_ADT_INFO[ADD_ADT_OUTPUT[i]] = {}
            dic_ADD_ADT_INFO[ADD_ADT_OUTPUT[i]]['ADD_ADT_INPUT_RDA_GE'] = ADD_ADT_INPUT_RDA_GE[i]
            dic_ADD_ADT_INFO[ADD_ADT_OUTPUT[i]]['ADD_ADT_INPUT_DIR_ADT'] = ADD_ADT_INPUT_DIR_ADT[i]
            dic_ADD_ADT_INFO[ADD_ADT_OUTPUT[i]]['ADD_ADT_SAMPLE_NAME_ADT'] = ADD_ADT_SAMPLE_NAME_ADT[i]

if "Adding_TCR" in STEPS:
    ### Sample/Project
    if 'Adding_TCR' in config and 'input.rda' in config['Adding_TCR'] :
        ADD_TCR_INPUT_RDA_GE = config['Adding_TCR']['input.rda']
    elif "Adding_ADT" in STEPS and "Clust_Markers_Annot_GE" not in STEPS:
        sys.stderr.write("Note: No input.rda find in Adding_TCR section of configfile; input.rda will be determine from Adding_ADT step for Adding_TCR step!\n")
        ADD_TCR_INPUT_RDA_GE = [ x + "_ADT.rda" for x in ADD_ADT_OUTPUT]
    elif "Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("Note: No input.rda find in Adding_TCR section of configfile; input.rda will be determine from Clust_Markers_Annot_GE step for Adding_TCR step!\n")
        #OLD: ADD_TCR_INPUT_RDA_GE = [os.path.normpath(os.path.dirname(dic_CMA_INFO[CMA_SAMPLE_NAME_GE[x]]['CMA_INPUT_RDA']) + "/" + CMA_CLUST_FOLDER + "/" + CMA_SAMPLE_NAME_GE[x] + CMA_COMPLEMENT[x] + "_" + str(CMA_KEEP_DIM) + "_" + str(CMA_KEEP_RES) + ".rda") for x in range(len(CMA_SAMPLE_NAME_GE))]
        ADD_TCR_INPUT_RDA_GE = []
        ADD_TCR_SAMPLE_NAME_GE = []
        for i in range(len(CMA_INPUT_RDA_GE)):
            for j in range(len(CMA_CLUST_FOLDERS)):
                ADD_TCR_SAMPLE_NAME_GE.append(dic_CMA_INFO[CMA_INPUT_RDA_GE[i]]['CMA_SAMPLE_NAME_GE'])
                ADD_TCR_INPUT_RDA_GE.append(os.path.normpath(dic_CMA_INFO[CMA_INPUT_RDA_GE[i]]['CMA_OUTPUT_DIR'] + "/" + CMA_CLUST_FOLDERS[j] + "/" + dic_CMA_INFO[CMA_INPUT_RDA_GE[i]]['CMA_SAMPLE_NAME_GE'] + CMA_COMPLEMENT[i] + "_" + CMA_KEEP_DIM_RES[j] + ".rda"))
        if "Adding_ADT" in STEPS:
            ADD_TCR_INPUT_RDA_GE = [ x + "_ADT.rda" for x in ADD_ADT_OUTPUT]
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
    ### Snakefile parameters
    #Correspondance between TCR sample names, input folder of TCR and GE sample names in yaml
    if "Alignment_countTable_GE" in STEPS:
        TMP_SAMPLE_NAME_GE = copy.deepcopy(ALIGN_SAMPLE_NAME_GE)
    elif "Droplets_QC_GE" in STEPS:
        TMP_SAMPLE_NAME_GE = copy.deepcopy(QC_SAMPLE_NAME_GE)
    elif "Filtering_GE" in STEPS:
        TMP_SAMPLE_NAME_GE = copy.deepcopy(FILTERING_SAMPLE_NAME_GE)
    elif "Norm_DimRed_Eval_GE" in STEPS:
        TMP_SAMPLE_NAME_GE = copy.deepcopy(NDRE_SAMPLE_NAME_GE)
    elif "Clust_Markers_Annot_GE" in STEPS:
        TMP_SAMPLE_NAME_GE = copy.deepcopy(CMA_SAMPLE_NAME_GE)
    elif "Adding_ADT" in STEPS:
        TMP_SAMPLE_NAME_GE = copy.deepcopy(ADD_ADT_SAMPLE_NAME_GE)
    if "Clust_Markers_Annot_GE" in STEPS:
        dic_ADD_TCR_CORRES_SAMPLE_NAMES = {}
        for i in range(len(TMP_SAMPLE_NAME_GE)):
            dic_ADD_TCR_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_GE[i]] = {}
            dic_ADD_TCR_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_GE[i]]['ADD_TCR_INPUT_CSV_TCR'] = ADD_TCR_INPUT_CSV_TCR[i]
            #dic_ADD_TCR_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_GE[i]]['ALIGN_SAMPLE_NAME_TCR'] = ALIGN_SAMPLE_NAME_TCR[i] #
        #Correspondance input/output (output = folder and file name without the extention "_ADT.rda")
        ADD_TCR_OUTPUT = [os.path.splitext(x)[0] for x in ADD_TCR_INPUT_RDA_GE]
        dic_ADD_TCR_INFO = {}
        for i in range(len(ADD_TCR_OUTPUT)):
            dic_ADD_TCR_INFO[ADD_TCR_OUTPUT[i]] = {}
            dic_ADD_TCR_INFO[ADD_TCR_OUTPUT[i]]['ADD_TCR_INPUT_RDA_GE'] = ADD_TCR_INPUT_RDA_GE[i]
            #dic_ADD_TCR_INFO[ADD_TCR_OUTPUT[i]]['ADD_TCR_SAMPLE_NAME_GE'] = ADD_TCR_SAMPLE_NAME_GE[i] #
            dic_ADD_TCR_INFO[ADD_TCR_OUTPUT[i]]['ADD_TCR_INPUT_CSV_TCR'] = dic_ADD_TCR_CORRES_SAMPLE_NAMES[ADD_TCR_SAMPLE_NAME_GE[i]]['ADD_TCR_INPUT_CSV_TCR']
            #dic_ADD_TCR_INFO[ADD_TCR_OUTPUT[i]]['ALIGN_SAMPLE_NAME_TCR'] = dic_ADD_TCR_CORRES_SAMPLE_NAMES[ADD_TCR_SAMPLE_NAME_GE[i]]['ALIGN_SAMPLE_NAME_TCR'] #
    else: # from "Adding_TCR" step or "Adding_ADT"(wo CMA)
        ADD_TCR_OUTPUT = [os.path.splitext(x)[0] for x in ADD_TCR_INPUT_RDA_GE]
        dic_ADD_TCR_INFO = {}
        for i in range(len(ADD_TCR_OUTPUT)):
            dic_ADD_TCR_INFO[ADD_TCR_OUTPUT[i]] = {}
            dic_ADD_TCR_INFO[ADD_TCR_OUTPUT[i]]['ADD_TCR_INPUT_RDA_GE'] = ADD_TCR_INPUT_RDA_GE[i]
            dic_ADD_TCR_INFO[ADD_TCR_OUTPUT[i]]['ADD_TCR_INPUT_CSV_TCR'] = ADD_TCR_INPUT_CSV_TCR[i]

    #sys.exit("Stop")

if "Adding_BCR" in STEPS:
    ### Sample/Project
    if 'Adding_BCR' in config and 'input.rda' in config['Adding_BCR'] :
        ADD_BCR_INPUT_RDA_GE = config['Adding_BCR']['input.rda']
    elif "Adding_TCR" in STEPS and "Clust_Markers_Annot_GE" not in STEPS:
        sys.stderr.write("Note: No input.rda.ge find in Adding_BCR section of configfile; input.rda.ge will be determine from Adding_TCR step for Adding_BCR step!\n")
        ADD_BCR_INPUT_RDA_GE = [ x + "_TCR.rda" for x in ADD_TCR_OUTPUT]
    elif "Adding_ADT" in STEPS and "Clust_Markers_Annot_GE" not in STEPS:
        sys.stderr.write("Note: No input.rda find in Adding_BCR section of configfile; input.rda will be determine from Adding_ADT step for Adding_BCR step!\n")
        ADD_BCR_INPUT_RDA_GE = [ x + "_ADT.rda" for x in ADD_ADT_OUTPUT]
    elif "Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("Note: No input.rda find in Adding_BCR section of configfile; input.rda will be determine from Clust_Markers_Annot_GE step for Adding_BCR step!\n")
        #OLD: ADD_BCR_INPUT_RDA_GE = [os.path.normpath(os.path.dirname(dic_CMA_INFO[CMA_SAMPLE_NAME_GE[x]]['CMA_INPUT_RDA']) + "/" + CMA_CLUST_FOLDER + "/" + CMA_SAMPLE_NAME_GE[x] + CMA_COMPLEMENT[x] + "_" + str(CMA_KEEP_DIM) + "_" + str(CMA_KEEP_RES) + ".rda") for x in range(len(CMA_SAMPLE_NAME_GE))]
        ADD_BCR_INPUT_RDA_GE = []
        ADD_BCR_SAMPLE_NAME_GE = []
        for i in range(len(CMA_INPUT_RDA_GE)):
            for j in range(len(CMA_CLUST_FOLDERS)):
                ADD_BCR_SAMPLE_NAME_GE.append(dic_CMA_INFO[CMA_INPUT_RDA_GE[i]]['CMA_SAMPLE_NAME_GE'])
                ADD_BCR_INPUT_RDA_GE.append(os.path.normpath(dic_CMA_INFO[CMA_INPUT_RDA_GE[i]]['CMA_OUTPUT_DIR'] + "/" + CMA_CLUST_FOLDERS[j] + "/" + dic_CMA_INFO[CMA_INPUT_RDA_GE[i]]['CMA_SAMPLE_NAME_GE'] + CMA_COMPLEMENT[i] + "_" + CMA_KEEP_DIM_RES[j] + ".rda"))
        if "Adding_TCR" in STEPS:
            ADD_BCR_INPUT_RDA_GE = [ x + "_TCR.rda" for x in ADD_TCR_OUTPUT]
        elif "Adding_ADT" in STEPS:
            ADD_TCR_INPUT_RDA_GE = [ x + "_ADT.rda" for x in ADD_ADT_OUTPUT]
    else:
        sys.exit("Error: No input.rda in configfile!\n")
    if 'Adding_BCR' in config and 'vdj.input.file.bcr' in config['Adding_BCR'] :
        ADD_BCR_INPUT_CSV_BCR = config['Adding_BCR']['vdj.input.file.bcr']
    elif "Alignment_annotations_TCR_BCR" in STEPS:
        sys.stderr.write("Note: No vdj.input.file.bcr find in Adding_BCR section of configfile; vdj.input.file.bcr will be determine from Alignment_annotations_TCR_BCR step for Adding_BCR step!\n")
        ALIGN_SAMPLE_NAME_BCR = [sample for sample in ALIGN_SAMPLE_NAME_TCR_BCR if bool(re.match(".+_BCR", sample))]
        ADD_BCR_INPUT_CSV_BCR = [ os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/" + x + "/" + x + "_CellRanger/outs/filtered_contig_annotations.csv") for x in ALIGN_SAMPLE_NAME_BCR]
    else:
        sys.exit("Error: No vdj.input.file.bcr in configfile!\n")
    ### Analysis Parameters
    ### Snakefile parameters
    #Correspondance between BCR sample names, input folder of BCR and GE sample names in yaml
    if "Alignment_countTable_GE" in STEPS:
        TMP_SAMPLE_NAME_GE = copy.deepcopy(ALIGN_SAMPLE_NAME_GE)
    elif "Droplets_QC_GE" in STEPS:
        TMP_SAMPLE_NAME_GE = copy.deepcopy(QC_SAMPLE_NAME_GE)
    elif "Filtering_GE" in STEPS:
        TMP_SAMPLE_NAME_GE = copy.deepcopy(FILTERING_SAMPLE_NAME_GE)
    elif "Norm_DimRed_Eval_GE" in STEPS:
        TMP_SAMPLE_NAME_GE = copy.deepcopy(NDRE_SAMPLE_NAME_GE)
    elif "Clust_Markers_Annot_GE" in STEPS:
        TMP_SAMPLE_NAME_GE = copy.deepcopy(CMA_SAMPLE_NAME_GE)
    if "Clust_Markers_Annot_GE" in STEPS:
        dic_ADD_BCR_CORRES_SAMPLE_NAMES = {}
        for i in range(len(TMP_SAMPLE_NAME_GE)):
            dic_ADD_BCR_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_GE[i]] = {}
            dic_ADD_BCR_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_GE[i]]['ADD_BCR_INPUT_CSV_BCR'] = ADD_BCR_INPUT_CSV_BCR[i]
            #dic_ADD_BCR_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_GE[i]]['ALIGN_SAMPLE_NAME_BCR'] = ALIGN_SAMPLE_NAME_BCR[i]
        #Correspondance input/output (output = folder and file name without the extention "_ADT.rda")
        ADD_BCR_OUTPUT = [os.path.splitext(x)[0] for x in ADD_BCR_INPUT_RDA_GE]
        dic_ADD_BCR_INFO = {}
        for i in range(len(ADD_BCR_OUTPUT)):
            dic_ADD_BCR_INFO[ADD_BCR_OUTPUT[i]] = {}
            dic_ADD_BCR_INFO[ADD_BCR_OUTPUT[i]]['ADD_BCR_INPUT_RDA_GE'] = ADD_BCR_INPUT_RDA_GE[i]
            #dic_ADD_BCR_INFO[ADD_BCR_OUTPUT[i]]['ADD_BCR_SAMPLE_NAME_GE'] = ADD_BCR_SAMPLE_NAME_GE[i]
            dic_ADD_BCR_INFO[ADD_BCR_OUTPUT[i]]['ADD_BCR_INPUT_CSV_BCR'] = dic_ADD_BCR_CORRES_SAMPLE_NAMES[ADD_BCR_SAMPLE_NAME_GE[i]]['ADD_BCR_INPUT_CSV_BCR']
            #dic_ADD_BCR_INFO[ADD_BCR_OUTPUT[i]]['ALIGN_SAMPLE_NAME_BCR'] = dic_ADD_BCR_CORRES_SAMPLE_NAMES[ADD_BCR_SAMPLE_NAME_GE[i]]['ALIGN_SAMPLE_NAME_BCR']
    else:
        ADD_BCR_OUTPUT = [os.path.splitext(x)[0] for x in ADD_BCR_INPUT_RDA_GE]
        dic_ADD_BCR_INFO = {}
        for i in range(len(ADD_BCR_OUTPUT)):
            dic_ADD_BCR_INFO[ADD_BCR_OUTPUT[i]] = {}
            dic_ADD_BCR_INFO[ADD_BCR_OUTPUT[i]]['ADD_BCR_INPUT_RDA_GE'] = ADD_BCR_INPUT_RDA_GE[i]
            dic_ADD_BCR_INFO[ADD_BCR_OUTPUT[i]]['ADD_BCR_INPUT_CSV_BCR'] = ADD_BCR_INPUT_CSV_BCR[i]

if "Int_Norm_DimRed_Eval_GE" in STEPS:
    ### Sample/Project
    if ('Int_Norm_DimRed_Eval_GE' in config) and ('name.int' in config['Int_Norm_DimRed_Eval_GE']) and ('input.list.rda' in config['Int_Norm_DimRed_Eval_GE']) :
        INT_NDRE_NAME_INT = config['Int_Norm_DimRed_Eval_GE']['name.int']
        INT_NDRE_INPUT_LIST_RDA_GE = config['Int_Norm_DimRed_Eval_GE']['input.list.rda']
    else:
        sys.exit("Error: No name.int or/and input.list.rda in configfile!\n")
    if ('Int_Norm_DimRed_Eval_GE' in config) and ('output.dir.int' in config['Int_Norm_DimRed_Eval_GE']) :
        INT_NDRE_OUTPUT_DIR_GE = [os.path.normpath(x) for x in config['Int_Norm_DimRed_Eval_GE']['output.dir.int']]
    else :
        sys.exit("Error: No output.dir.int find in configfile!\n")
    ### Analysis Parameters
    INT_NDRE_EVAL_MARKERS = getelseDefault('Int_Norm_DimRed_Eval_GE', 'eval.markers', 'NULL').replace(", ", ",")
    # Load data
    INT_NDRE_MIN_CELLS = getelseDefault('Int_Norm_DimRed_Eval_GE', 'min.cells', 'NULL')
    # Integration
    if ('Int_Norm_DimRed_Eval_GE' in config) and ('integration.method' in config['Int_Norm_DimRed_Eval_GE']):
        INT_NDRE_INT_METHOD = getelseDefault('Int_Norm_DimRed_Eval_GE', 'integration.method', 'NULL')
    else :
        sys.exit("Error: No integration.method found in configfile!\n")
    INT_NDRE_VTR_BATCH = getelseDefault('Int_Norm_DimRed_Eval_GE', 'vtr.batch', 'orig.ident')
    # Normalization and dimension reduction
    INT_NDRE_FEATURES_N = getelseDefault('Int_Norm_DimRed_Eval_GE', 'features.n', 'NULL')
    INT_NDRE_NORM_METHOD = getelseDefault('Int_Norm_DimRed_Eval_GE', 'norm.method', 'SCTransform')
    INT_NDRE_HVG_FINDVARIABLEFEATURESMIX = getelseDefault('Int_Norm_DimRed_Eval_GE', 'HVG.FindVariableFeaturesMix', 'FALSE')
    INT_NDRE_DIMRED_METHOD = getelseDefault('Int_Norm_DimRed_Eval_GE', 'dimred.method', 'pca')
    #INT_NDRE_VTR_BIASES_NORM = [ x.replace(", ", ",") for x in config['Int_Norm_DimRed_Eval_GE']['vtr.biases.norm']] if ('Int_Norm_DimRed_Eval_GE' in config and 'vtr.biases.norm' in config['Int_Norm_DimRed_Eval_GE'] and config['Int_Norm_DimRed_Eval_GE']['vtr.biases.norm'] != None) else "NULL"
    INT_NDRE_VTR_BIASES_NORM = getelseDefault('Int_Norm_DimRed_Eval_GE', 'vtr.biases.norm', ['NULL'])
    if not isinstance(INT_NDRE_VTR_BIASES_NORM, list): INT_NDRE_VTR_BIASES_NORM = list(INT_NDRE_VTR_BIASES_NORM)
    INT_NDRE_VTR_BIASES_NORM = ["NULL" if x == "" else x.replace(", ", ",") for x in INT_NDRE_VTR_BIASES_NORM]
    #INT_NDRE_VTR_BIASES_DIMRED = [ x.replace(", ", ",") for x in config['Int_Norm_DimRed_Eval_GE']['vtr.biases.dimred']] if ('Int_Norm_DimRed_Eval_GE' in config and 'vtr.biases.dimred' in config['Int_Norm_DimRed_Eval_GE'] and config['Int_Norm_DimRed_Eval_GE']['vtr.biases.dimred'] != None) else "NULL"
    INT_NDRE_VTR_BIASES_DIMRED = getelseDefault('Int_Norm_DimRed_Eval_GE', 'vtr.biases.dimred', ['NULL'])
    if not isinstance(INT_NDRE_VTR_BIASES_DIMRED, list): INT_NDRE_VTR_BIASES_DIMRED = list(INT_NDRE_VTR_BIASES_DIMRED)
    INT_NDRE_VTR_BIASES_DIMRED = ["NULL" if x == "" else x.replace(", ", ",") for x in INT_NDRE_VTR_BIASES_DIMRED]
    INT_NDRE_VTR_SCALE = getelseDefault('Int_Norm_DimRed_Eval_GE', 'vtr.scale', 'NULL')
    INT_NDRE_REGEX_REMOVE_HVG = getelseDefault('Int_Norm_DimRed_Eval_GE', 'regex.genes.to.remove.from.HVG', 'NULL')
    INT_NDRE_DIM_MAX = getelseDefault('Int_Norm_DimRed_Eval_GE', 'dims.max', 49)
    INT_NDRE_SKIP_EVALDIMRES = getelseDefault('Int_Norm_DimRed_Eval_GE', 'skip.eval_dims_res', 'NULL')
    INT_NDRE_EVAL_DIM_MAX = getelseDefault('Int_Norm_DimRed_Eval_GE', 'eval.dims.max', 49)
    INT_NDRE_EVAL_DIM_MIN = getelseDefault('Int_Norm_DimRed_Eval_GE', 'eval.dims.min', 9)
    INT_NDRE_EVAL_DIM_STEPS = getelseDefault('Int_Norm_DimRed_Eval_GE', 'eval.dims.steps', 2)
    INT_NDRE_EVAL_RES_MAX = getelseDefault('Int_Norm_DimRed_Eval_GE', 'eval.res.max', 1.2)
    INT_NDRE_EVAL_RES_MIN = getelseDefault('Int_Norm_DimRed_Eval_GE', 'eval.res.min', 0.1)
    INT_NDRE_EVAL_RES_STEPS = getelseDefault('Int_Norm_DimRed_Eval_GE', 'eval.res.steps', 0.1)
    INT_NDRE_EVAL_PTSIZE = getelseDefault('Int_Norm_DimRed_Eval_GE', 'eval.pt.size', 'NULL').replace(", ", ",")
    # Metadata file
    INT_NDRE_METADATA_FILE = getelseDefault('Int_Norm_DimRed_Eval_GE', 'metadata.file', 'NULL').replace(", ", ",")
    ### Snakefile parameters
    #Correspondance sample/input/output
    dic_INT_NDRE_INFO = {}
    for i in range(0,len(INT_NDRE_INPUT_LIST_RDA_GE),1):
        dic_INT_NDRE_INFO[INT_NDRE_INPUT_LIST_RDA_GE[i]] = {}
        dic_INT_NDRE_INFO[INT_NDRE_INPUT_LIST_RDA_GE[i]]['INT_NDRE_NAME_INT'] = INT_NDRE_NAME_INT[i]
        dic_INT_NDRE_INFO[INT_NDRE_INPUT_LIST_RDA_GE[i]]['INT_NDRE_OUTPUT_DIR'] = INT_NDRE_OUTPUT_DIR_GE[i] + "/GROUPED_ANALYSIS/INTEGRATED/" + INT_NDRE_NAME_INT[i]
    #Names
    #if (INT_NDRE_VTR_BIASES_NORM == ["NULL"]):
    #    INT_NDRE_NORM_VTR = [INT_NDRE_NORM_METHOD]
    #else :
    #    INT_NDRE_NORM_VTR = [INT_NDRE_NORM_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(x.split(","))))) for x in INT_NDRE_VTR_BIASES_NORM]
    INT_NDRE_NORM_VTR = [INT_NDRE_NORM_METHOD if (x == "NULL") else (INT_NDRE_NORM_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(x.split(",")))))) for x in INT_NDRE_VTR_BIASES_NORM]

    if (INT_NDRE_INT_METHOD == "scbfa" or INT_NDRE_INT_METHOD == "bpca" or INT_NDRE_INT_METHOD == "Harmony" or INT_NDRE_INT_METHOD == "Seurat_CCAIntegration" or INT_NDRE_INT_METHOD == "Seurat_RPCAIntegration" or INT_NDRE_INT_METHOD == "Seurat_HarmonyIntegration" or INT_NDRE_INT_METHOD == "Seurat_scVIIntegration"):
        if (INT_NDRE_VTR_BATCH == "NULL") : sys.exit("Error: No vtr.batch can't be empty with scbfa, bpca, Harmony or Seurat integration!\n")
    #if (INT_NDRE_DIMRED_METHOD == "pca" or INT_NDRE_DIMRED_METHOD == "mds" or INT_NDRE_DIMRED_METHOD == "ica" or INT_NDRE_VTR_BIASES_DIMRED == ["NULL"]):
    #    INT_NDRE_DIMRED_VTR = [INT_NDRE_DIMRED_METHOD]
    #else:
    #    INT_NDRE_DIMRED_VTR = [INT_NDRE_DIMRED_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(x.split(","))))) for x in INT_NDRE_VTR_BIASES_DIMRED]
    if (INT_NDRE_DIMRED_METHOD == "pca" or INT_NDRE_DIMRED_METHOD == "mds" or INT_NDRE_DIMRED_METHOD == "ica"):
        INT_NDRE_DIMRED_VTR = [INT_NDRE_DIMRED_METHOD if (x == "NULL") else (INT_NDRE_DIMRED_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(x.split(",")))))) for x in INT_NDRE_VTR_BIASES_DIMRED]
    #if (INT_NDRE_INT_METHOD == "scbfa" or INT_NDRE_INT_METHOD == "bpca"):
        #if (INT_NDRE_VTR_BIASES_DIMRED == ["NULL"]) :
        #    INT_NDRE_DIMRED_VTR = [INT_NDRE_INT_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(INT_NDRE_VTR_BATCH.split(",")))))]
        #else :
        #    INT_NDRE_DIMRED_VTR = [INT_NDRE_INT_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(x.split(","))) + list(dict.fromkeys(INT_NDRE_VTR_BATCH.split(","))))) for x in INT_NDRE_VTR_BIASES_DIMRED]
    if (INT_NDRE_INT_METHOD == "scbfa" or INT_NDRE_INT_METHOD == "bpca"):
        INT_NDRE_DIMRED_VTR = [INT_NDRE_INT_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(INT_NDRE_VTR_BATCH.split(","))))) if (x == "NULL") else (INT_NDRE_INT_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(x.split(","))) + list(dict.fromkeys(INT_NDRE_VTR_BATCH.split(",")))))) for x in INT_NDRE_VTR_BIASES_DIMRED]

if "Int_Clust_Markers_Annot_GE" in STEPS:
    ### Sample/Project
    if ('Int_Clust_Markers_Annot_GE' in config) and ('name.int' in config['Int_Clust_Markers_Annot_GE']) and ('input.rda.int' in config['Int_Clust_Markers_Annot_GE']) :
        INT_CMA_NAME_INT = config['Int_Clust_Markers_Annot_GE']['name.int']
        INT_CMA_INPUT_RDA_GE = config['Int_Clust_Markers_Annot_GE']['input.rda.int']
    elif "Int_Norm_DimRed_Eval_GE" in STEPS:
        sys.stderr.write("Note: No input.rda.int and name.int find in Int_Clust_Markers_Annot_GE section of configfile; input.rda.int and name.int will be determine from Int_Norm_DimRed_Eval_GE step for Int_Clust_Markers_Annot_GE step!\n")
        INT_CMA_NAME_INT = []
        INT_CMA_INPUT_RDA_GE = []
        for i in INT_NDRE_INPUT_LIST_RDA_GE:
            for j in INT_NDRE_NORM_VTR:
                for k in INT_NDRE_DIMRED_VTR:
                    INT_CMA_NAME_INT.append(dic_INT_NDRE_INFO[i]['INT_NDRE_NAME_INT'])
                    INT_CMA_INPUT_RDA_GE.append(os.path.normpath(dic_INT_NDRE_INFO[i]['INT_NDRE_OUTPUT_DIR'] + "/" + j + "/" + k + "/" + dic_INT_NDRE_INFO[i]['INT_NDRE_NAME_INT'] + "_" + j + "_" + k + ".rda"))
    else:
        sys.exit("Error: No name.int or/and input.rda.int in configfile!\n")
    if ('Int_Clust_Markers_Annot_GE' in config) and ('output.dir.int' in config['Int_Clust_Markers_Annot_GE']) :
        INT_CMA_OUTPUT_DIR_GE = config['Int_Clust_Markers_Annot_GE']['output.dir.int']
    elif "Int_Norm_DimRed_Eval_GE" in STEPS:
        #INT_CMA_OUTPUT_DIR_GE = [os.path.normpath(dic_INT_NDRE_INFO[x]['INT_NDRE_OUTPUT_DIR'] + "/" + INT_NDRE_NORM_VTR + "/" + INT_NDRE_DIMRED_VTR) for x in INT_NDRE_NAME_INT]
        INT_CMA_OUTPUT_DIR_GE = [os.path.dirname(x) for x in INT_CMA_INPUT_RDA_GE]
        sys.stderr.write("Note: No output.dir.int find in Int_Clust_Markers_Annot_GE section of configfile; output.dir.int will be determine from Int_Norm_DimRed_Eval_GE step for Int_Clust_Markers_Annot_GE step!\n")
    else :
        sys.exit("Error: No output.dir.int find in configfile!\n")
    ### Analysis Parameters
    # Normalization and dimension reduction
    INT_CMA_KEEP_DIM_RES = getelseDefault('Int_Clust_Markers_Annot_GE', 'keep.dim.res', 'NULL')
    INT_CMA_KEEP_DIM_RES = [str(dim_res).replace(".0", "").replace(",0","") for dim_res in INT_CMA_KEEP_DIM_RES]
    # Annotation
    INT_CMA_CFR_MINSCORE = getelseDefault('Int_Clust_Markers_Annot_GE', 'cfr.minscore', 'NULL')
    INT_CMA_SR_MINSCORE = getelseDefault('Int_Clust_Markers_Annot_GE', 'sr.minscore', 'NULL')
    INT_CMA_CUSTOM_SCE_REF = getelseDefault('Int_Clust_Markers_Annot_GE', 'custom.sce.ref', 'NULL').replace(", ", ",")
    INT_CMA_CUSTOM_MARKERS_REF = getelseDefault('Int_Clust_Markers_Annot_GE', 'custom.markers.ref', 'NULL').replace(", ", ",")    
    # Markers
    INT_CMA_MARKFILE = getelseDefault('Int_Clust_Markers_Annot_GE', 'markfile', 'NULL').replace(", ", ",")
    INT_CMA_MARKERS_PTSIZE = getelseDefault('Int_Clust_Markers_Annot_GE', 'markers.pt.size', 'NULL').replace(", ", ",")
    INT_CMA_MARKERS_ORDER = getelseDefault('Int_Clust_Markers_Annot_GE', 'markers.order', 'NULL').replace(", ", ",")
    # SKIP
    INT_CMA_SKIP_TECHNICAL = getelseDefault('Int_Clust_Markers_Annot_GE', 'skip.technical_plots', 'NULL')
    INT_CMA_SKIP_ANNOT = getelseDefault('Int_Clust_Markers_Annot_GE', 'skip.annotation', 'NULL')
    INT_CMA_SKIP_MARKERS_IDENT = getelseDefault('Int_Clust_Markers_Annot_GE', 'skip.markers_identification', 'NULL')
    # Metadata file
    INT_CMA_METADATA_FILE = getelseDefault('Int_Clust_Markers_Annot_GE', 'metadata.file', 'NULL').replace(", ", ",")
    ### Snakefile parameters
    #check end paths (del "/" if necessary)
    for i in range(0,len(INT_CMA_OUTPUT_DIR_GE),1):
        INT_CMA_OUTPUT_DIR_GE[i] = os.path.normpath(INT_CMA_OUTPUT_DIR_GE[i])
    #Correspondance sample/input/output
    dic_INT_CMA_INFO = {}
    INT_CMA_COMPLEMENT = []
    for i in range(0,len(INT_CMA_INPUT_RDA_GE),1):
        dic_INT_CMA_INFO[INT_CMA_INPUT_RDA_GE[i]] = {}
        dic_INT_CMA_INFO[INT_CMA_INPUT_RDA_GE[i]]['INT_CMA_NAME_INT'] = INT_CMA_NAME_INT[i]
        dic_INT_CMA_INFO[INT_CMA_INPUT_RDA_GE[i]]['INT_CMA_OUTPUT_DIR'] = INT_CMA_OUTPUT_DIR_GE[i]
        compl = os.path.splitext(os.path.basename(INT_CMA_INPUT_RDA_GE[i]))[0]
        if compl.startswith(INT_CMA_NAME_INT[i]):
            compl = compl[len(INT_CMA_NAME_INT[i]):]
        INT_CMA_COMPLEMENT.append(compl)
    #Names
    INT_CMA_CLUST_FOLDERS = [("dims" + str(dim_res).replace("_", "_res")) for dim_res in INT_CMA_KEEP_DIM_RES]

if "Int_Adding_ADT" in STEPS:
    ### Sample/Project
    if 'Int_Adding_ADT' in config and 'input.rda' in config['Int_Adding_ADT']:
        INT_ADD_ADT_INPUT_RDA = config['Int_Adding_ADT']['input.rda']
    elif "Int_Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("Note: No input.rda find in Int_Adding_ADT section of configfile; input.rda will be determine from Int_Clust_Markers_Annot_GE step for Int_Adding_ADT step!\n")
        INT_ADD_ADT_INPUT_RDA = []
        INT_ADD_ADT_NAME_INT = []
        for i in range(len(INT_CMA_INPUT_RDA_GE)):
            for j in range(len(INT_CMA_CLUST_FOLDERS)):
                INT_ADD_ADT_NAME_INT.append(dic_INT_CMA_INFO[INT_CMA_INPUT_RDA_GE[i]]['INT_CMA_NAME_INT'])
                INT_ADD_ADT_INPUT_RDA.append(os.path.normpath(dic_INT_CMA_INFO[INT_CMA_INPUT_RDA_GE[i]]['INT_CMA_OUTPUT_DIR'] + "/" + INT_CMA_CLUST_FOLDERS[j] + "/" + dic_INT_CMA_INFO[INT_CMA_INPUT_RDA_GE[i]]['INT_CMA_NAME_INT'] + INT_CMA_COMPLEMENT[i] + "_" + INT_CMA_KEEP_DIM_RES[j] + ".rda"))
    else:
        sys.exit("Error: No input.rda in configfile!\n")
    if 'Int_Adding_ADT' in config and 'samples.name.adt' in config['Int_Adding_ADT'] and 'input.dirs.adt' in config['Int_Adding_ADT'] :
        INT_ADD_ADT_INPUT_DIR_ADT = [ x.replace(", ", ",") for x in config['Int_Adding_ADT']['input.dirs.adt']]
        INT_ADD_ADT_SAMPLE_NAME_ADT_RAW = [ x.replace(", ", ",") for x in config['Int_Adding_ADT']['samples.name.adt']]
        #check samples names and add "_ADT" if needed
        INT_ADD_ADT_SAMPLE_NAME_ADT = []
        for i in range(0,len(INT_ADD_ADT_SAMPLE_NAME_ADT_RAW),1):
            list_sample_tmp = []
            for sample in INT_ADD_ADT_SAMPLE_NAME_ADT_RAW[i].split(","):
                if sample[len(sample)-4:] != "_ADT" :
                    list_sample_tmp.append(sample + "_ADT")
                else :
                    list_sample_tmp.append(sample)
            INT_ADD_ADT_SAMPLE_NAME_ADT.append(",".join(list_sample_tmp))
    else:
        sys.exit("Error: No samples.name.adt or input.dirs.adt in configfile!\n")
    ### Analysis Parameters
    INT_ADD_ADT_GENE_NAMES = getelseDefault('Int_Adding_ADT', 'gene.names', 'NULL').replace(", ", ",")
    ### Snakefile parameters
    #Correspondance input/output
    if "Int_Norm_DimRed_Eval_GE" in STEPS:
        TMP_SAMPLE_NAME_INT = copy.deepcopy(INT_NDRE_NAME_INT)
    elif "Clust_Markers_Annot_GE" in STEPS:
        TMP_SAMPLE_NAME_INT = copy.deepcopy(INT_CMA_NAME_INT)
    if "Int_Norm_DimRed_Eval_GE" in STEPS or "Int_Clust_Markers_Annot_GE" in STEPS:
        dic_INT_ADD_ADT_CORRES_SAMPLE_NAMES = {}
        for i in range(len(TMP_SAMPLE_NAME_INT)):
            dic_INT_ADD_ADT_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_INT[i]] = {}
            dic_INT_ADD_ADT_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_INT[i]]['INT_ADD_ADT_INPUT_DIR_ADT'] = INT_ADD_ADT_INPUT_DIR_ADT[i]
            dic_INT_ADD_ADT_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_INT[i]]['INT_ADD_ADT_SAMPLE_NAME_ADT'] = INT_ADD_ADT_SAMPLE_NAME_ADT[i]
        #Correspondance input/output (output = folder+file name without the extention "_ADT.rda")
        INT_ADD_ADT_OUTPUT = [os.path.splitext(x)[0] for x in INT_ADD_ADT_INPUT_RDA]
        dic_INT_ADD_ADT_INFO = {}
        for i in range(len(INT_ADD_ADT_OUTPUT)):
            dic_INT_ADD_ADT_INFO[INT_ADD_ADT_OUTPUT[i]] = {}
            dic_INT_ADD_ADT_INFO[INT_ADD_ADT_OUTPUT[i]]['INT_ADD_ADT_INPUT_RDA'] = INT_ADD_ADT_INPUT_RDA[i]
            dic_INT_ADD_ADT_INFO[INT_ADD_ADT_OUTPUT[i]]['INT_ADD_ADT_INPUT_DIR_ADT'] = dic_INT_ADD_ADT_CORRES_SAMPLE_NAMES[INT_ADD_ADT_NAME_INT[i]]['INT_ADD_ADT_INPUT_DIR_ADT']
            dic_INT_ADD_ADT_INFO[INT_ADD_ADT_OUTPUT[i]]['INT_ADD_ADT_SAMPLE_NAME_ADT'] = dic_INT_ADD_ADT_CORRES_SAMPLE_NAMES[INT_ADD_ADT_NAME_INT[i]]['INT_ADD_ADT_SAMPLE_NAME_ADT']
    else: #Int_Adding_ADT
        INT_ADD_ADT_OUTPUT = [os.path.splitext(x)[0] for x in INT_ADD_ADT_INPUT_RDA]
        dic_INT_ADD_ADT_INFO = {}
        for i in range(len(INT_ADD_ADT_OUTPUT)):
            dic_INT_ADD_ADT_INFO[INT_ADD_ADT_OUTPUT[i]] = {}
            dic_INT_ADD_ADT_INFO[INT_ADD_ADT_OUTPUT[i]]['INT_ADD_ADT_INPUT_RDA'] = INT_ADD_ADT_INPUT_RDA[i]
            dic_INT_ADD_ADT_INFO[INT_ADD_ADT_OUTPUT[i]]['INT_ADD_ADT_INPUT_DIR_ADT'] = INT_ADD_ADT_INPUT_DIR_ADT[i]
            dic_INT_ADD_ADT_INFO[INT_ADD_ADT_OUTPUT[i]]['INT_ADD_ADT_SAMPLE_NAME_ADT'] = INT_ADD_ADT_SAMPLE_NAME_ADT[i]

if "Int_Adding_TCR" in STEPS:
    ### Sample/Project
    if 'Int_Adding_TCR' in config and 'input.rda' in config['Int_Adding_TCR'] :
        INT_ADD_TCR_INPUT_RDA = config['Int_Adding_TCR']['input.rda']
    elif "Int_Adding_ADT" in STEPS and not "Int_Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("Note: No input.rda find in Int_Adding_TCR section of configfile; input.rda will be determine from Int_Adding_ADT step for Int_Adding_TCR step!\n")
        INT_ADD_TCR_INPUT_RDA = [ x + "_ADT.rda" for x in INT_ADD_ADT_OUTPUT]
    elif "Int_Clust_Markers_Annot_GE" in STEPS:
        #identification of the end of rda files
        rda_ending=".rda"
        message = True
        if "Int_Adding_ADT" in STEPS:
            rda_ending = "_ADT" + rda_ending
            if message: sys.stderr.write("Note: No input.rda find in Int_Adding_TCR section of configfile; input.rda will be determine from Int_Adding_ADT step for Int_Adding_TCR step!\n")
            message = False
        if message:
            sys.stderr.write("Note: No input.rda find in Int_Adding_TCR section of configfile; input.rda will be determine from Clust_Markers_Annot_GE step for Int_Adding_TCR step!\n")
        #identification of name and rda files.
        INT_ADD_TCR_INPUT_RDA = []
        INT_ADD_TCR_SAMPLE_NAME_INT = []
        for i in range(len(INT_CMA_INPUT_RDA_GE)):
            for j in range(len(INT_CMA_CLUST_FOLDERS)):
                INT_ADD_TCR_INPUT_RDA.append(os.path.normpath(dic_INT_CMA_INFO[INT_CMA_INPUT_RDA_GE[i]]['INT_CMA_OUTPUT_DIR'] + "/" + INT_CMA_CLUST_FOLDERS[j] + "/" + dic_INT_CMA_INFO[INT_CMA_INPUT_RDA_GE[i]]['INT_CMA_NAME_INT'] + INT_CMA_COMPLEMENT[i] + "_" + INT_CMA_KEEP_DIM_RES[j] + rda_ending))
                INT_ADD_TCR_SAMPLE_NAME_INT.append(dic_INT_CMA_INFO[INT_CMA_INPUT_RDA_GE[i]]['INT_CMA_NAME_INT'])
    else:
        sys.exit("Error: No input.rda in configfile!\n")

    if 'Int_Adding_TCR' in config and 'vdj.input.files.tcr' in config['Int_Adding_TCR'] :
        INT_ADD_TCR_INPUT_CSV_TCR = [ x.replace(", ", ",") for x in config['Int_Adding_TCR']['vdj.input.files.tcr']]
    else:
        sys.exit("Error: No vdj.input.files.tcr in configfile!\n")
    ### Snakefile parameters:
    #Correspondance between input folders of TCR (of this step) and Int sample names (in previous step) in yaml (because some rda input can be duplicated with previous steps, but not folder of TCR from this step)
    if "Int_Norm_DimRed_Eval_GE" in STEPS:
        TMP_SAMPLE_NAME_INT = copy.deepcopy(INT_NDRE_NAME_INT)
    elif "Int_Clust_Markers_Annot_GE" in STEPS:
        TMP_SAMPLE_NAME_INT = copy.deepcopy(INT_CMA_NAME_INT)
    if "Int_Norm_DimRed_Eval_GE" in STEPS or "Int_Clust_Markers_Annot_GE" in STEPS:
        dic_INT_ADD_TCR_CORRES_SAMPLE_NAMES = {}
        for i in range(len(TMP_SAMPLE_NAME_INT)):
            dic_INT_ADD_TCR_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_INT[i]] = {}
            dic_INT_ADD_TCR_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_INT[i]]['INT_ADD_TCR_INPUT_CSV_TCR'] = INT_ADD_TCR_INPUT_CSV_TCR[i]
        #Correspondance input/output (output = folder and file name without the extention "_ADT.rda")
        INT_ADD_TCR_OUTPUT = [os.path.splitext(x)[0] for x in INT_ADD_TCR_INPUT_RDA]
        dic_INT_ADD_TCR_INFO = {}
        for i in range(len(INT_ADD_TCR_OUTPUT)):
            dic_INT_ADD_TCR_INFO[INT_ADD_TCR_OUTPUT[i]] = {}
            dic_INT_ADD_TCR_INFO[INT_ADD_TCR_OUTPUT[i]]['INT_ADD_TCR_INPUT_RDA'] = INT_ADD_TCR_INPUT_RDA[i]
            dic_INT_ADD_TCR_INFO[INT_ADD_TCR_OUTPUT[i]]['INT_ADD_TCR_INPUT_CSV_TCR'] = dic_INT_ADD_TCR_CORRES_SAMPLE_NAMES[INT_ADD_TCR_SAMPLE_NAME_INT[i]]['INT_ADD_TCR_INPUT_CSV_TCR']
    else: #the pipeline begun directly by ADT or TCR, so not need to make the correspondance with previous steps
        INT_ADD_TCR_OUTPUT = [os.path.splitext(x)[0] for x in INT_ADD_TCR_INPUT_RDA]
        dic_INT_ADD_TCR_INFO = {}
        for i in range(len(INT_ADD_TCR_OUTPUT)):
            dic_INT_ADD_TCR_INFO[INT_ADD_TCR_OUTPUT[i]] = {}
            dic_INT_ADD_TCR_INFO[INT_ADD_TCR_OUTPUT[i]]['INT_ADD_TCR_INPUT_RDA'] = INT_ADD_TCR_INPUT_RDA[i]
            dic_INT_ADD_TCR_INFO[INT_ADD_TCR_OUTPUT[i]]['INT_ADD_TCR_INPUT_CSV_TCR'] = INT_ADD_TCR_INPUT_CSV_TCR[i]

if "Int_Adding_BCR" in STEPS:
    ### Sample/Project
    if 'Int_Adding_BCR' in config and 'input.rda' in config['Int_Adding_BCR'] :
        INT_ADD_BCR_INPUT_RDA = config['Int_Adding_BCR']['input.rda']
    elif "Int_Adding_TCR" in STEPS and not "Int_Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("Note: No input.rda find in Int_Adding_BCR section of configfile; input.rda will be determine from Int_Adding_TCR step for Int_Adding_BCR step!\n")
        INT_ADD_BCR_INPUT_RDA = [ x + "_TCR.rda" for x in INT_ADD_TCR_OUTPUT]
    elif "Int_Adding_ADT" in STEPS and not "Int_Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("Note: No input.rda find in Int_Adding_BCR section of configfile; input.rda will be determine from Int_Adding_ADT step for Int_Adding_BCR step!\n")
        INT_ADD_BCR_INPUT_RDA = [ x + "_ADT.rda" for x in INT_ADD_ADT_OUTPUT]
    elif "Int_Clust_Markers_Annot_GE" in STEPS:
        #identification of the end of rda files
        rda_ending=".rda"
        message = True
        if "Int_Adding_TCR" in STEPS:
            rda_ending = "_TCR" + rda_ending
            sys.stderr.write("Note: No input.rda find in Int_Adding_BCR section of configfile; input.rda will be determine from Int_Adding_TCR step for Int_Adding_BCR step!\n")
            message = False
        if "Int_Adding_ADT" in STEPS:
            rda_ending = "_ADT" + rda_ending
            if message: sys.stderr.write("Note: No input.rda find in Int_Adding_BCR section of configfile; input.rda will be determine from Int_Adding_ADT step for Int_Adding_BCR step!\n")
            message = False
        if message:
            sys.stderr.write("Note: No input.rda find in Int_Adding_BCR section of configfile; input.rda will be determine from Clust_Markers_Annot_GE step for Int_Adding_BCR step!\n")
        #identification of name and rda files.
        INT_ADD_BCR_INPUT_RDA = []
        INT_ADD_BCR_SAMPLE_NAME_INT = []
        for i in range(len(INT_CMA_INPUT_RDA_GE)):
            for j in range(len(INT_CMA_CLUST_FOLDERS)):
                INT_ADD_BCR_INPUT_RDA.append(os.path.normpath(dic_INT_CMA_INFO[INT_CMA_INPUT_RDA_GE[i]]['INT_CMA_OUTPUT_DIR'] + "/" + INT_CMA_CLUST_FOLDERS[j] + "/" + dic_INT_CMA_INFO[INT_CMA_INPUT_RDA_GE[i]]['INT_CMA_NAME_INT'] + INT_CMA_COMPLEMENT[i] + "_" + INT_CMA_KEEP_DIM_RES[j] + rda_ending))
                INT_ADD_BCR_SAMPLE_NAME_INT.append(dic_INT_CMA_INFO[INT_CMA_INPUT_RDA_GE[i]]['INT_CMA_NAME_INT'])
    else:
        sys.exit("Error: No input.rda in configfile!\n")
    if 'Int_Adding_BCR' in config and 'vdj.input.files.bcr' in config['Int_Adding_BCR'] :
        INT_ADD_BCR_INPUT_CSV_BCR = [ x.replace(", ", ",") for x in config['Int_Adding_BCR']['vdj.input.files.bcr']]
    else:
        sys.exit("Error: No vdj.input.files.bcr in configfile!\n")
    ### Snakefile parameters
    #Correspondance between input folders of BCR (of this step) and Int sample names (in previous step) in yaml (because some rda input can be duplicated with previous steps, but not folder of BCR from this step)
    if "Int_Norm_DimRed_Eval_GE" in STEPS:
        TMP_SAMPLE_NAME_INT = copy.deepcopy(INT_NDRE_NAME_INT)
    elif "Int_Clust_Markers_Annot_GE" in STEPS:
        TMP_SAMPLE_NAME_INT = copy.deepcopy(INT_CMA_NAME_INT)
    if "Int_Norm_DimRed_Eval_GE" in STEPS or "Int_Clust_Markers_Annot_GE" in STEPS:
        dic_INT_ADD_BCR_CORRES_SAMPLE_NAMES = {}
        for i in range(len(TMP_SAMPLE_NAME_INT)):
            dic_INT_ADD_BCR_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_INT[i]] = {}
            dic_INT_ADD_BCR_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_INT[i]]['INT_ADD_BCR_INPUT_CSV_BCR'] = INT_ADD_BCR_INPUT_CSV_BCR[i]
        #Correspondance input/output (output = folder and file name without the extention "_ADT.rda")
        INT_ADD_BCR_OUTPUT = [os.path.splitext(x)[0] for x in INT_ADD_BCR_INPUT_RDA]
        dic_INT_ADD_BCR_INFO = {}
        for i in range(len(INT_ADD_BCR_OUTPUT)):
            dic_INT_ADD_BCR_INFO[INT_ADD_BCR_OUTPUT[i]] = {}
            dic_INT_ADD_BCR_INFO[INT_ADD_BCR_OUTPUT[i]]['INT_ADD_BCR_INPUT_RDA'] = INT_ADD_BCR_INPUT_RDA[i]
            dic_INT_ADD_BCR_INFO[INT_ADD_BCR_OUTPUT[i]]['INT_ADD_BCR_INPUT_CSV_BCR'] = dic_INT_ADD_BCR_CORRES_SAMPLE_NAMES[INT_ADD_BCR_SAMPLE_NAME_INT[i]]['INT_ADD_BCR_INPUT_CSV_BCR']
    else: #the pipeline begun directly by ADT or TCR or BCR, so not need to make the correspondance with previous steps
        INT_ADD_BCR_OUTPUT = [os.path.splitext(x)[0] for x in INT_ADD_BCR_INPUT_RDA]
        dic_INT_ADD_BCR_INFO = {}
        for i in range(len(INT_ADD_BCR_OUTPUT)):
            dic_INT_ADD_BCR_INFO[INT_ADD_BCR_OUTPUT[i]] = {}
            dic_INT_ADD_BCR_INFO[INT_ADD_BCR_OUTPUT[i]]['INT_ADD_BCR_INPUT_RDA'] = INT_ADD_BCR_INPUT_RDA[i]
            dic_INT_ADD_BCR_INFO[INT_ADD_BCR_OUTPUT[i]]['INT_ADD_BCR_INPUT_CSV_BCR'] = INT_ADD_BCR_INPUT_CSV_BCR[i]

if "Grp_Norm_DimRed_Eval_GE" in STEPS:
    ### Sample/Project
    if ('Grp_Norm_DimRed_Eval_GE' in config) and ('name.grp' in config['Grp_Norm_DimRed_Eval_GE']) and ('input.list.rda' in config['Grp_Norm_DimRed_Eval_GE']) :
        GRP_NDRE_NAME_GRP = config['Grp_Norm_DimRed_Eval_GE']['name.grp']
        GRP_NDRE_INPUT_LIST_RDA_GE = config['Grp_Norm_DimRed_Eval_GE']['input.list.rda']
    else:
        sys.exit("Error: No name.grp or/and input.list.rda in configfile!\n")
    if ('Grp_Norm_DimRed_Eval_GE' in config) and ('output.dir.grp' in config['Grp_Norm_DimRed_Eval_GE']) :
        GRP_NDRE_OUTPUT_DIR_GE = [os.path.normpath(x) for x in config['Grp_Norm_DimRed_Eval_GE']['output.dir.grp']]
    else :
        sys.exit("Error: No output.dir.grp find in configfile!\n")
    ### Analysis Parameters
    GRP_NDRE_EVAL_MARKERS = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'eval.markers', 'NULL').replace(", ", ",")
    # Load data
    GRP_NDRE_MIN_CELLS = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'min.cells', 'NULL')
    # Normalization and dimension reduction
    GRP_NDRE_INDIV_NORM = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'individual.norm', 'FALSE')
    GRP_NDRE_FEATURES_N = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'features.n', 'NULL')
    GRP_NDRE_NORM_METHOD = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'norm.method', 'SCTransform')
    GRP_NDRE_HVG_FINDVARIABLEFEATURESMIX = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'HVG.FindVariableFeaturesMix', 'FALSE')
    GRP_NDRE_DIMRED_METHOD = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'dimred.method', 'pca')
    #GRP_NDRE_VTR_BIASES_NORM = [ x.replace(", ", ",") for x in config['Grp_Norm_DimRed_Eval_GE']['vtr.biases.norm']] if ('Grp_Norm_DimRed_Eval_GE' in config and 'vtr.biases.norm' in config['Grp_Norm_DimRed_Eval_GE'] and config['Grp_Norm_DimRed_Eval_GE']['vtr.biases.norm'] != None) else "NULL"
    GRP_NDRE_VTR_BIASES_NORM = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'vtr.biases.norm', ['NULL'])
    if not isinstance(GRP_NDRE_VTR_BIASES_NORM, list): GRP_NDRE_VTR_BIASES_NORM = list(GRP_NDRE_VTR_BIASES_NORM)
    GRP_NDRE_VTR_BIASES_NORM = ["NULL" if x == "" else x.replace(", ", ",") for x in GRP_NDRE_VTR_BIASES_NORM]
    #GRP_NDRE_VTR_BIASES_DIMRED = [ x.replace(", ", ",") for x in config['Grp_Norm_DimRed_Eval_GE']['vtr.biases.dimred']] if ('Grp_Norm_DimRed_Eval_GE' in config and 'vtr.biases.dimred' in config['Grp_Norm_DimRed_Eval_GE'] and config['Grp_Norm_DimRed_Eval_GE']['vtr.biases.dimred'] != None) else "NULL"
    GRP_NDRE_VTR_BIASES_DIMRED = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'vtr.biases.dimred', ['NULL'])
    if not isinstance(GRP_NDRE_VTR_BIASES_DIMRED, list): GRP_NDRE_VTR_BIASES_DIMRED = list(GRP_NDRE_VTR_BIASES_DIMRED)
    GRP_NDRE_VTR_BIASES_DIMRED = ["NULL" if x == "" else x.replace(", ", ",") for x in GRP_NDRE_VTR_BIASES_DIMRED]
    GRP_NDRE_VTR_SCALE = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'vtr.scale', 'NULL')
    GRP_NDRE_REGEX_REMOVE_HVG = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'regex.genes.to.remove.from.HVG', 'NULL')
    GRP_NDRE_DIM_MAX = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'dims.max', 49)
    GRP_NDRE_SKIP_EVALDIMRES = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'skip.eval_dims_res', 'NULL')
    GRP_NDRE_EVAL_DIM_MAX = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'eval.dims.max', 49)
    GRP_NDRE_EVAL_DIM_MIN = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'eval.dims.min', 9)
    GRP_NDRE_EVAL_DIM_STEPS = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'eval.dims.step', 2)
    GRP_NDRE_EVAL_RES_MAX = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'eval.res.max', 1.2)
    GRP_NDRE_EVAL_RES_MIN = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'eval.res.min', 0.1)
    GRP_NDRE_EVAL_RES_STEPS = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'eval.res.steps', 0.1)
    GRP_NDRE_EVAL_PTSIZE = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'eval.pt.size', 'NULL')
    # Metadata file
    GRP_NDRE_METADATA_FILE = getelseDefault('Grp_Norm_DimRed_Eval_GE', 'metadata.file', 'NULL').replace(", ", ",")
    ### Snakefile parameters
    #Correspondance sample/input/output
    dic_GRP_NDRE_INFO = {}
    for i in range(0,len(GRP_NDRE_INPUT_LIST_RDA_GE),1):
        dic_GRP_NDRE_INFO[GRP_NDRE_INPUT_LIST_RDA_GE[i]] = {}
        dic_GRP_NDRE_INFO[GRP_NDRE_INPUT_LIST_RDA_GE[i]]['GRP_NDRE_NAME_GRP'] = GRP_NDRE_NAME_GRP[i]
        dic_GRP_NDRE_INFO[GRP_NDRE_INPUT_LIST_RDA_GE[i]]['GRP_NDRE_OUTPUT_DIR'] = GRP_NDRE_OUTPUT_DIR_GE[i] + "/GROUPED_ANALYSIS/NO_INTEGRATED/" + GRP_NDRE_NAME_GRP[i]
    #Names
    GRP_NDRE_NORM_VTR = [GRP_NDRE_NORM_METHOD if (x == "NULL") else (GRP_NDRE_NORM_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(x.split(",")))))) for x in GRP_NDRE_VTR_BIASES_NORM]
    #elif (GRP_NDRE_VTR_BIASES_NORM == "NULL"):
    #    GRP_NDRE_NORM_VTR = [GRP_NDRE_NORM_METHOD]
    #else :
    #    #GRP_NDRE_NORM_VTR =  GRP_NDRE_NORM_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(GRP_NDRE_VTR_BIASES_NORM.split(",")))))
    #    GRP_NDRE_NORM_VTR = [GRP_NDRE_NORM_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(x.split(","))))) for x in GRP_NDRE_VTR_BIASES_NORM]
    #if (GRP_NDRE_DIMRED_METHOD == "pca" or GRP_NDRE_DIMRED_METHOD == "mds" or GRP_NDRE_DIMRED_METHOD == "ica" or GRP_NDRE_VTR_BIASES_DIMRED == "NULL"):
    #    GRP_NDRE_DIMRED_VTR = [GRP_NDRE_DIMRED_METHOD]
    #else:
    #    GRP_NDRE_DIMRED_VTR = [GRP_NDRE_DIMRED_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(x.split(","))))) for x in GRP_NDRE_VTR_BIASES_DIMRED]
    GRP_NDRE_DIMRED_VTR = [GRP_NDRE_DIMRED_METHOD if (x == "NULL") else (GRP_NDRE_DIMRED_METHOD + "_" + "_".join(sorted(list(dict.fromkeys(x.split(",")))))) for x in GRP_NDRE_VTR_BIASES_DIMRED]

    
if "Grp_Clust_Markers_Annot_GE" in STEPS:
    ### Sample/Project
    if ('Grp_Clust_Markers_Annot_GE' in config) and ('name.grp' in config['Grp_Clust_Markers_Annot_GE']) and ('input.rda.grp' in config['Grp_Clust_Markers_Annot_GE']) :
        GRP_CMA_NAME_GRP = config['Grp_Clust_Markers_Annot_GE']['name.grp']
        GRP_CMA_INPUT_RDA_GE = config['Grp_Clust_Markers_Annot_GE']['input.rda.grp']
    elif "Grp_Norm_DimRed_Eval_GE" in STEPS:
        sys.stderr.write("Note: No input.rda.grp and name.grp find in Grp_Clust_Markers_Annot_GE section of configfile; input.rda.grp and name.grp will be determine from Grp_Norm_DimRed_Eval_GE step for Grp_Clust_Markers_Annot_GE step!\n")
        #GRP_CMA_NAME_GRP = copy.deepcopy(GRP_NDRE_NAME_GRP)
        #GRP_CMA_INPUT_RDA_GE = [os.path.normpath(dic_GRP_NDRE_INFO[x]['GRP_NDRE_OUTPUT_DIR'] + "/" + GRP_NDRE_NORM_VTR + "/" + GRP_NDRE_DIMRED_VTR + "/" + x + "_" + GRP_NDRE_NORM_VTR + "_" + GRP_NDRE_DIMRED_VTR + ".rda") for x in GRP_NDRE_NAME_GRP]
        GRP_CMA_NAME_GRP = []
        GRP_CMA_INPUT_RDA_GE = []
        for i in GRP_NDRE_INPUT_LIST_RDA_GE:
            for j in GRP_NDRE_NORM_VTR:
                for k in GRP_NDRE_DIMRED_VTR:
                    GRP_CMA_NAME_GRP.append(dic_GRP_NDRE_INFO[i]['GRP_NDRE_NAME_GRP'])
                    GRP_CMA_INPUT_RDA_GE.append(os.path.normpath(dic_GRP_NDRE_INFO[i]['GRP_NDRE_OUTPUT_DIR'] + "/" + j + "/" + k + "/" + dic_GRP_NDRE_INFO[i]['GRP_NDRE_NAME_GRP'] + "_" + j + "_" + k + ".rda"))
    else:
        sys.exit("Error: No name.grp or/and input.rda.grp in configfile!\n")

    if ('Grp_Clust_Markers_Annot_GE' in config) and ('output.dir.grp' in config['Grp_Clust_Markers_Annot_GE']) :
        GRP_CMA_OUTPUT_DIR_GE = config['Grp_Clust_Markers_Annot_GE']['output.dir.grp']
    elif "Grp_Norm_DimRed_Eval_GE" in STEPS:
        GRP_CMA_OUTPUT_DIR_GE = [os.path.dirname(x) for x in GRP_CMA_INPUT_RDA_GE]
        sys.stderr.write("Note: No output.dir.grp find in Grp_Clust_Markers_Annot_GE section of configfile; output.dir.grp will be determine from Grp_Norm_DimRed_Eval_GE step for Grp_Clust_Markers_Annot_GE step!\n")
    else :
        sys.exit("Error: No output.dir.grp find in configfile!\n")
    ### Analysis Parameters
    GRP_CMA_MARKFILE = getelseDefault('Grp_Clust_Markers_Annot_GE', 'markfile', 'NULL').replace(", ", ",")
    GRP_CMA_MARKERS_PTSIZE = getelseDefault('Grp_Clust_Markers_Annot_GE', 'markers.pt.size', 'NULL').replace(", ", ",")
    GRP_CMA_MARKERS_ORDER = getelseDefault('Grp_Clust_Markers_Annot_GE', 'markers.order', 'NULL').replace(", ", ",")
    GRP_CMA_CUSTOM_SCE_REF = getelseDefault('Grp_Clust_Markers_Annot_GE', 'custom.sce.ref', 'NULL').replace(", ", ",")
    GRP_CMA_CUSTOM_MARKERS_REF = getelseDefault('Grp_Clust_Markers_Annot_GE', 'custom.markers.ref', 'NULL').replace(", ", ",")
    # Normalization and dimension reduction
    GRP_CMA_KEEP_DIM_RES = getelseDefault('Grp_Clust_Markers_Annot_GE', 'keep.dim.res', 'NULL')
    GRP_CMA_KEEP_DIM_RES = [str(dim_res).replace(".0", "").replace(",0","") for dim_res in GRP_CMA_KEEP_DIM_RES]
    GRP_CMA_CFR_MINSCORE =  getelseDefault('Grp_Clust_Markers_Annot_GE', 'cfr.minscores', 'NULL')
    GRP_CMA_SR_MINSCORE =  getelseDefault('Grp_Clust_Markers_Annot_GE', 'sr.minscore', 'NULL')
    # SKIP
    GRP_CMA_SKIP_TECHNICAL = getelseDefault('Grp_Clust_Markers_Annot_GE', 'skip.technical_plots', 'NULL')
    GRP_CMA_SKIP_ANNOT = getelseDefault('Grp_Clust_Markers_Annot_GE', 'skip.annotation', 'NULL')
    GRP_CMA_SKIP_MARKERS_IDENT = getelseDefault('Grp_Clust_Markers_Annot_GE', 'skip.markers_identification', 'NULL')
    # Metadata file
    GRP_CMA_METADATA_FILE = getelseDefault('Grp_Clust_Markers_Annot_GE', 'metadata.file', 'NULL').replace(", ", ",")
    ### Snakefile parameters
    #check end paths (del "/" if necessary)
    for i in range(0,len(GRP_CMA_OUTPUT_DIR_GE),1):
        GRP_CMA_OUTPUT_DIR_GE[i] = os.path.normpath(GRP_CMA_OUTPUT_DIR_GE[i])
    #Correspondance sample/input/output
    dic_GRP_CMA_INFO = {}
    GRP_CMA_COMPLEMENT = []
    for i in range(0,len(GRP_CMA_INPUT_RDA_GE),1):
        dic_GRP_CMA_INFO[GRP_CMA_INPUT_RDA_GE[i]] = {}
        dic_GRP_CMA_INFO[GRP_CMA_INPUT_RDA_GE[i]]['GRP_CMA_NAME_GRP'] = GRP_CMA_NAME_GRP[i]
        dic_GRP_CMA_INFO[GRP_CMA_INPUT_RDA_GE[i]]['GRP_CMA_OUTPUT_DIR'] = GRP_CMA_OUTPUT_DIR_GE[i]
        compl = os.path.splitext(os.path.basename(GRP_CMA_INPUT_RDA_GE[i]))[0]
        if compl.startswith(GRP_CMA_NAME_GRP[i]):
            compl = compl[len(GRP_CMA_NAME_GRP[i]):]
        GRP_CMA_COMPLEMENT.append(compl)
    #Names
    GRP_CMA_CLUST_FOLDERS = [("dims" + str(dim_res).replace("_", "_res")) for dim_res in GRP_CMA_KEEP_DIM_RES]

if "Grp_Adding_ADT" in STEPS:
    ### Sample/Project
    if 'Grp_Adding_ADT' in config and 'input.rda' in config['Grp_Adding_ADT']:
        GRP_ADD_ADT_INPUT_RDA = config['Grp_Adding_ADT']['input.rda']
    elif "Grp_Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("Note: No input.rda find in Grp_Adding_ADT section of configfile; input.rda will be determine from Grp_Clust_Markers_Annot_GE step for Grp_Adding_ADT step!\n")
        #GRP_ADD_ADT_INPUT_RDA = [os.path.normpath(os.path.dirname(dic_GRP_CMA_INFO[GRP_CMA_NAME_GRP[x]]['GRP_CMA_INPUT_RDA']) + "/" + GRP_CMA_CLUST_FOLDER + "/" + GRP_CMA_NAME_GRP[x] + GRP_CMA_COMPLEMENT[x] + "_" + str(GRP_CMA_KEEP_DIM) + "_" + str(GRP_CMA_KEEP_RES) + ".rda") for x in range(len(GRP_CMA_NAME_GRP))]
        GRP_ADD_ADT_INPUT_RDA = []
        GRP_ADD_ADT_NAME_GRP = []
        for i in range(len(GRP_CMA_INPUT_RDA_GE)):
            for j in range(len(GRP_CMA_CLUST_FOLDERS)):
                GRP_ADD_ADT_NAME_GRP.append(dic_GRP_CMA_INFO[GRP_CMA_INPUT_RDA_GE[i]]['GRP_CMA_NAME_GRP'])
                GRP_ADD_ADT_INPUT_RDA.append(os.path.normpath(dic_GRP_CMA_INFO[GRP_CMA_INPUT_RDA_GE[i]]['GRP_CMA_OUTPUT_DIR'] + "/" + GRP_CMA_CLUST_FOLDERS[j] + "/" + dic_GRP_CMA_INFO[GRP_CMA_INPUT_RDA_GE[i]]['GRP_CMA_NAME_GRP'] + GRP_CMA_COMPLEMENT[i] + "_" + GRP_CMA_KEEP_DIM_RES[j] + ".rda"))
    else:
        sys.exit("Error: No input.rda in configfile!\n")
    if 'Grp_Adding_ADT' in config and 'samples.name.adt' in config['Grp_Adding_ADT'] and 'input.dirs.adt' in config['Grp_Adding_ADT'] :
        GRP_ADD_ADT_INPUT_DIR_ADT = [ x.replace(", ", ",") for x in config['Grp_Adding_ADT']['input.dirs.adt']]
        GRP_ADD_ADT_SAMPLE_NAME_ADT_RAW = [ x.replace(", ", ",") for x in config['Grp_Adding_ADT']['samples.name.adt']]
        #check samples names and add "_ADT" if needed
        GRP_ADD_ADT_SAMPLE_NAME_ADT = []
        for i in range(0,len(GRP_ADD_ADT_SAMPLE_NAME_ADT_RAW),1):
            list_sample_tmp = []
            for sample in GRP_ADD_ADT_SAMPLE_NAME_ADT_RAW[i].split(","):
                if sample[len(sample)-4:] != "_ADT" :
                    list_sample_tmp.append(sample + "_ADT")
                else :
                    list_sample_tmp.append(sample)
            GRP_ADD_ADT_SAMPLE_NAME_ADT.append(",".join(list_sample_tmp))
    else:
        sys.exit("Error: No samples.name.adt or input.dirs.adt in configfile!\n")
    ### Analysis Parameters
    GRP_ADD_ADT_GENE_NAMES = getelseDefault('Grp_Adding_ADT', 'gene.names', 'NULL').replace(", ", ",")
    ### Snakefile parameters
    #Correspondance input/output
    if "Grp_Norm_DimRed_Eval_GE" in STEPS:
        TMP_SAMPLE_NAME_GRP = copy.deepcopy(GRP_NDRE_NAME_GRP)
    elif "Clust_Markers_Annot_GE" in STEPS:
        TMP_SAMPLE_NAME_GRP = copy.deepcopy(GRP_ADD_ADT_NAME_GRP)
    if "Grp_Norm_DimRed_Eval_GE" in STEPS or "Grp_Clust_Markers_Annot_GE" in STEPS:
        dic_GRP_ADD_ADT_CORRES_SAMPLE_NAMES = {}
        for i in range(len(TMP_SAMPLE_NAME_GRP)):
            dic_GRP_ADD_ADT_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_GRP[i]] = {}
            dic_GRP_ADD_ADT_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_GRP[i]]['GRP_ADD_ADT_INPUT_DIR_ADT'] = GRP_ADD_ADT_INPUT_DIR_ADT[i]
            dic_GRP_ADD_ADT_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_GRP[i]]['GRP_ADD_ADT_SAMPLE_NAME_ADT'] = GRP_ADD_ADT_SAMPLE_NAME_ADT[i]
        #Correspondance input/output (output = folder+file name without the extention "_ADT.rda")
        GRP_ADD_ADT_OUTPUT = [os.path.splitext(x)[0] for x in GRP_ADD_ADT_INPUT_RDA]
        dic_GRP_ADD_ADT_INFO = {}
        for i in range(len(GRP_ADD_ADT_OUTPUT)):
            dic_GRP_ADD_ADT_INFO[GRP_ADD_ADT_OUTPUT[i]] = {}
            dic_GRP_ADD_ADT_INFO[GRP_ADD_ADT_OUTPUT[i]]['GRP_ADD_ADT_INPUT_RDA'] = GRP_ADD_ADT_INPUT_RDA[i]
            dic_GRP_ADD_ADT_INFO[GRP_ADD_ADT_OUTPUT[i]]['GRP_ADD_ADT_INPUT_DIR_ADT'] = dic_GRP_ADD_ADT_CORRES_SAMPLE_NAMES[GRP_ADD_ADT_NAME_GRP[i]]['GRP_ADD_ADT_INPUT_DIR_ADT']
            dic_GRP_ADD_ADT_INFO[GRP_ADD_ADT_OUTPUT[i]]['GRP_ADD_ADT_SAMPLE_NAME_ADT'] = dic_GRP_ADD_ADT_CORRES_SAMPLE_NAMES[GRP_ADD_ADT_NAME_GRP[i]]['GRP_ADD_ADT_SAMPLE_NAME_ADT']
    else:
        GRP_ADD_ADT_OUTPUT = [os.path.splitext(x)[0] for x in GRP_ADD_ADT_INPUT_RDA]
        dic_GRP_ADD_ADT_INFO = {}
        for i in range(len(GRP_ADD_ADT_OUTPUT)):
            dic_GRP_ADD_ADT_INFO[GRP_ADD_ADT_OUTPUT[i]] = {}
            dic_GRP_ADD_ADT_INFO[GRP_ADD_ADT_OUTPUT[i]]['GRP_ADD_ADT_INPUT_RDA'] = GRP_ADD_ADT_INPUT_RDA[i]
            dic_GRP_ADD_ADT_INFO[GRP_ADD_ADT_OUTPUT[i]]['GRP_ADD_ADT_INPUT_DIR_ADT'] = GRP_ADD_ADT_INPUT_DIR_ADT[i]
            dic_GRP_ADD_ADT_INFO[GRP_ADD_ADT_OUTPUT[i]]['GRP_ADD_ADT_SAMPLE_NAME_ADT'] = GRP_ADD_ADT_SAMPLE_NAME_ADT[i]

if "Grp_Adding_TCR" in STEPS:
    ### Sample/Project
    if 'Grp_Adding_TCR' in config and 'input.rda' in config['Grp_Adding_TCR'] :
        GRP_ADD_TCR_INPUT_RDA = config['Grp_Adding_TCR']['input.rda']
    elif "Grp_Adding_ADT" in STEPS and not "Grp_Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("Note: No input.rda find in Grp_Adding_TCR section of configfile; input.rda will be determine from Grp_Adding_ADT step for Grp_Adding_TCR step!\n")
        GRP_ADD_TCR_INPUT_RDA = [ x + "_ADT.rda" for x in GRP_ADD_ADT_OUTPUT]
    elif "Grp_Clust_Markers_Annot_GE" in STEPS:
        #identification of the end of rda files
        rda_ending=".rda"
        message = True
        if "Grp_Adding_ADT" in STEPS:
            rda_ending = "_ADT" + rda_ending
            if message: sys.stderr.write("Note: No input.rda find in Grp_Adding_TCR section of configfile; input.rda will be determine from Grp_Adding_ADT step for Grp_Adding_TCR step!\n")
            message = False
        if message:
            sys.stderr.write("Note: No input.rda find in Grp_Adding_TCR section of configfile; input.rda will be determine from Clust_Markers_Annot_GE step for Grp_Adding_TCR step!\n")
        #identification of name and rda files.
        GRP_ADD_TCR_INPUT_RDA = []
        GRP_ADD_TCR_SAMPLE_NAME_GRP = []
        for i in range(len(GRP_CMA_INPUT_RDA_GE)):
            for j in range(len(GRP_CMA_CLUST_FOLDERS)):
                GRP_ADD_TCR_INPUT_RDA.append(os.path.normpath(dic_GRP_CMA_INFO[GRP_CMA_INPUT_RDA_GE[i]]['GRP_CMA_OUTPUT_DIR'] + "/" + GRP_CMA_CLUST_FOLDERS[j] + "/" + dic_GRP_CMA_INFO[GRP_CMA_INPUT_RDA_GE[i]]['GRP_CMA_NAME_GRP'] + GRP_CMA_COMPLEMENT[i] + "_" + GRP_CMA_KEEP_DIM_RES[j] + rda_ending))
                GRP_ADD_TCR_SAMPLE_NAME_GRP.append(dic_GRP_CMA_INFO[GRP_CMA_INPUT_RDA_GE[i]]['GRP_CMA_NAME_GRP'])
    else:
        sys.exit("Error: No input.rda in configfile!\n")
    if 'Grp_Adding_TCR' in config and 'vdj.input.files.tcr' in config['Grp_Adding_TCR'] :
        GRP_ADD_TCR_INPUT_CSV_TCR = [ x.replace(", ", ",") for x in config['Grp_Adding_TCR']['vdj.input.files.tcr']]
    else:
        sys.exit("Error: No vdj.input.files.tcr in configfile!\n")
    ### Snakefile parameters
    #Correspondance between ADT sample names, input folder of ADT and GE sample names in yaml
    if "Grp_Norm_DimRed_Eval_GE" in STEPS:
        TMP_SAMPLE_NAME_GRP = copy.deepcopy(GRP_NDRE_NAME_GRP)
    elif "Grp_Clust_Markers_Annot_GE" in STEPS:
        TMP_SAMPLE_NAME_GRP = copy.deepcopy(GRP_CMA_NAME_GRP)
    elif "Grp_Adding_ADT" in STEPS:
        TMP_SAMPLE_NAME_GRP = copy.deepcopy(GRP_ADD_ADT_SAMPLE_NAME_ADT)
    if "Grp_Norm_DimRed_Eval_GE" in STEPS or "Grp_Clust_Markers_Annot_GE" in STEPS:
        dic_GRP_ADD_TCR_CORRES_SAMPLE_NAMES = {}
        for i in range(len(TMP_SAMPLE_NAME_GRP)):
            dic_GRP_ADD_TCR_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_GRP[i]] = {}
            dic_GRP_ADD_TCR_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_GRP[i]]['GRP_ADD_TCR_INPUT_CSV_TCR'] = GRP_ADD_TCR_INPUT_CSV_TCR[i]
        #Correspondance input/output (output = folder and file name without the extention "_ADT.rda")
        GRP_ADD_TCR_OUTPUT = [os.path.splitext(x)[0] for x in GRP_ADD_TCR_INPUT_RDA]
        dic_GRP_ADD_TCR_INFO = {}
        for i in range(len(GRP_ADD_TCR_OUTPUT)):
            dic_GRP_ADD_TCR_INFO[GRP_ADD_TCR_OUTPUT[i]] = {}
            dic_GRP_ADD_TCR_INFO[GRP_ADD_TCR_OUTPUT[i]]['GRP_ADD_TCR_INPUT_RDA'] = GRP_ADD_TCR_INPUT_RDA[i]
            dic_GRP_ADD_TCR_INFO[GRP_ADD_TCR_OUTPUT[i]]['GRP_ADD_TCR_INPUT_CSV_TCR'] = dic_GRP_ADD_TCR_CORRES_SAMPLE_NAMES[GRP_ADD_TCR_SAMPLE_NAME_GRP[i]]['GRP_ADD_TCR_INPUT_CSV_TCR']
    else:
        GRP_ADD_TCR_OUTPUT = [os.path.splitext(x)[0] for x in GRP_ADD_TCR_INPUT_RDA]
        dic_GRP_ADD_TCR_INFO = {}
        for i in range(len(GRP_ADD_TCR_OUTPUT)):
            dic_GRP_ADD_TCR_INFO[GRP_ADD_TCR_OUTPUT[i]] = {}
            dic_GRP_ADD_TCR_INFO[GRP_ADD_TCR_OUTPUT[i]]['GRP_ADD_TCR_INPUT_RDA'] = GRP_ADD_TCR_INPUT_RDA[i]
            dic_GRP_ADD_TCR_INFO[GRP_ADD_TCR_OUTPUT[i]]['GRP_ADD_TCR_INPUT_CSV_TCR'] = GRP_ADD_TCR_INPUT_CSV_TCR[i]

if "Grp_Adding_BCR" in STEPS:
    ### Sample/Project
    if 'Grp_Adding_BCR' in config and 'input.rda' in config['Grp_Adding_BCR'] :
        GRP_ADD_BCR_INPUT_RDA = config['Grp_Adding_BCR']['input.rda']
    elif "Grp_Adding_TCR" in STEPS and not "Grp_Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("Note: No input.rda find in Grp_Adding_BCR section of configfile; input.rda will be determine from Grp_Adding_TCR step for Grp_Adding_BCR step!\n")
        GRP_ADD_BCR_INPUT_RDA = [ x + "_TCR.rda" for x in GRP_ADD_TCR_OUTPUT]
    elif "Grp_Adding_ADT" in STEPS and not "Grp_Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("Note: No input.rda find in Grp_Adding_BCR section of configfile; input.rda will be determine from Grp_Adding_ADT step for Grp_Adding_BCR step!\n")
        GRP_ADD_BCR_INPUT_RDA = [ x + "_ADT.rda" for x in GRP_ADD_ADT_OUTPUT]
    elif "Grp_Clust_Markers_Annot_GE" in STEPS:
        #identification of the end of rda files
        rda_ending=".rda"
        message = True
        if "Grp_Adding_TCR" in STEPS:
            rda_ending = "_TCR" + rda_ending
            sys.stderr.write("Note: No input.rda find in Grp_Adding_BCR section of configfile; input.rda will be determine from Grp_Adding_TCR step for Grp_Adding_BCR step!\n")
            message = False
        if "Grp_Adding_ADT" in STEPS:
            rda_ending = "_ADT" + rda_ending
            if message: sys.stderr.write("Note: No input.rda find in Grp_Adding_BCR section of configfile; input.rda will be determine from Grp_Adding_ADT step for Grp_Adding_BCR step!\n")
            message = False
        if message:
            sys.stderr.write("Note: No input.rda find in Grp_Adding_BCR section of configfile; input.rda will be determine from Clust_Markers_Annot_GE step for Grp_Adding_BCR step!\n")
        #identification of name and rda files.
        GRP_ADD_BCR_INPUT_RDA = []
        GRP_ADD_BCR_SAMPLE_NAME_GRP = []
        for i in range(len(GRP_CMA_INPUT_RDA_GE)):
            for j in range(len(GRP_CMA_CLUST_FOLDERS)):
                GRP_ADD_BCR_INPUT_RDA.append(os.path.normpath(dic_GRP_CMA_INFO[GRP_CMA_INPUT_RDA_GE[i]]['GRP_CMA_OUTPUT_DIR'] + "/" + GRP_CMA_CLUST_FOLDERS[j] + "/" + dic_GRP_CMA_INFO[GRP_CMA_INPUT_RDA_GE[i]]['GRP_CMA_NAME_GRP'] + GRP_CMA_COMPLEMENT[i] + "_" + GRP_CMA_KEEP_DIM_RES[j] + rda_ending))
                GRP_ADD_BCR_SAMPLE_NAME_GRP.append(dic_GRP_CMA_INFO[GRP_CMA_INPUT_RDA_GE[i]]['GRP_CMA_NAME_GRP'])
    else:
        sys.exit("Error: No input.rda in configfile!\n")
    if 'Grp_Adding_BCR' in config and 'vdj.input.files.bcr' in config['Grp_Adding_BCR'] :
        GRP_ADD_BCR_INPUT_CSV_BCR = [ x.replace(", ", ",") for x in config['Grp_Adding_BCR']['vdj.input.files.bcr']]
    else:
        sys.exit("Error: No vdj.input.files.bcr in configfile!\n")
    ### Snakefile parameters
    #Correspondance between ADT sample names, input folder of ADT and GE sample names in yaml
    if "Grp_Norm_DimRed_Eval_GE" in STEPS:
        TMP_SAMPLE_NAME_GRP = copy.deepcopy(GRP_NDRE_NAME_GRP)
    elif "Grp_Clust_Markers_Annot_GE" in STEPS:
        TMP_SAMPLE_NAME_GRP = copy.deepcopy(GRP_CMA_NAME_GRP)
    if "Grp_Norm_DimRed_Eval_GE" in STEPS or "Grp_Clust_Markers_Annot_GE" in STEPS:
        dic_GRP_ADD_BCR_CORRES_SAMPLE_NAMES = {}
        for i in range(len(TMP_SAMPLE_NAME_GRP)):
            dic_GRP_ADD_BCR_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_GRP[i]] = {}
            dic_GRP_ADD_BCR_CORRES_SAMPLE_NAMES[TMP_SAMPLE_NAME_GRP[i]]['GRP_ADD_BCR_INPUT_CSV_BCR'] = GRP_ADD_BCR_INPUT_CSV_BCR[i]
        #Correspondance input/output (output = folder and file name without the extention "_ADT.rda")
        GRP_ADD_BCR_OUTPUT = [os.path.splitext(x)[0] for x in GRP_ADD_BCR_INPUT_RDA]
        dic_GRP_ADD_BCR_INFO = {}
        for i in range(len(GRP_ADD_BCR_OUTPUT)):
            dic_GRP_ADD_BCR_INFO[GRP_ADD_BCR_OUTPUT[i]] = {}
            dic_GRP_ADD_BCR_INFO[GRP_ADD_BCR_OUTPUT[i]]['GRP_ADD_BCR_INPUT_RDA'] = GRP_ADD_BCR_INPUT_RDA[i]
            dic_GRP_ADD_BCR_INFO[GRP_ADD_BCR_OUTPUT[i]]['GRP_ADD_BCR_INPUT_CSV_BCR'] = dic_GRP_ADD_BCR_CORRES_SAMPLE_NAMES[GRP_ADD_BCR_SAMPLE_NAME_GRP[i]]['GRP_ADD_BCR_INPUT_CSV_BCR']
    else:
        GRP_ADD_BCR_OUTPUT = [os.path.splitext(x)[0] for x in GRP_ADD_BCR_INPUT_RDA]
        dic_GRP_ADD_BCR_INFO = {}
        for i in range(len(GRP_ADD_BCR_OUTPUT)):
            dic_GRP_ADD_BCR_INFO[GRP_ADD_BCR_OUTPUT[i]] = {}
            dic_GRP_ADD_BCR_INFO[GRP_ADD_BCR_OUTPUT[i]]['GRP_ADD_BCR_INPUT_RDA'] = GRP_ADD_BCR_INPUT_RDA[i]
            dic_GRP_ADD_BCR_INFO[GRP_ADD_BCR_OUTPUT[i]]['GRP_ADD_BCR_INPUT_CSV_BCR'] = GRP_ADD_BCR_INPUT_CSV_BCR[i]

if "Cerebro" in STEPS:
    ### Sample/Project
    CEREBRO_INPUT_RDA = []
    if 'Cerebro' in config and 'input.rda' in config['Cerebro'] :
        CEREBRO_INPUT_RDA = CEREBRO_INPUT_RDA + config['Cerebro']['input.rda']
    if "Adding_BCR" in STEPS:
        sys.stderr.write("Note: input.rda will be determine from Adding_BCR_GE step for Cerebro step!\n")
        CEREBRO_INPUT_RDA = CEREBRO_INPUT_RDA + [ x + "_BCR.rda" for x in ADD_BCR_OUTPUT]
    elif "Adding_TCR" in STEPS:
        sys.stderr.write("Note: input.rda will be determine from Adding_TCR_GE step for Cerebro step!\n")
        CEREBRO_INPUT_RDA = CEREBRO_INPUT_RDA + [ x + "_TCR.rda" for x in ADD_TCR_OUTPUT]
    elif "Adding_ADT" in STEPS:
        sys.stderr.write("Note: input.rda will be determine from Adding_ADT step for Cerebro step!\n")
        CEREBRO_INPUT_RDA = CEREBRO_INPUT_RDA + [ x + "_ADT.rda" for x in ADD_ADT_OUTPUT]
    elif "Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("Note: input.rda will be determine from Clust_Markers_Annot_GE step for Cerebro step!\n")
        for i in range(len(CMA_INPUT_RDA_GE)):
            for j in range(len(CMA_KEEP_DIM_RES)):
                CEREBRO_INPUT_RDA.append(os.path.normpath(dic_CMA_INFO[CMA_INPUT_RDA_GE[i]]['CMA_OUTPUT_DIR'] + "/" + CMA_CLUST_FOLDERS[j] + "/" + dic_CMA_INFO[CMA_INPUT_RDA_GE[i]]['CMA_SAMPLE_NAME_GE'] + CMA_COMPLEMENT[i] + "_" + CMA_KEEP_DIM_RES[j] + ".rda"))
    if "Int_Adding_BCR" in STEPS:
        sys.stderr.write("Note: input.rda will be determine from Int_Adding_BCR_GE step for Cerebro step!\n")
        CEREBRO_INPUT_RDA = CEREBRO_INPUT_RDA + [ x + "_BCR.rda" for x in INT_ADD_BCR_OUTPUT]
    elif "Int_Adding_TCR" in STEPS:
        sys.stderr.write("Note: input.rda will be determine from Int_Adding_TCR_GE step for Cerebro step!\n")
        CEREBRO_INPUT_RDA = CEREBRO_INPUT_RDA + [ x + "_TCR.rda" for x in INT_ADD_TCR_OUTPUT]
    elif "Int_Adding_ADT" in STEPS:
        sys.stderr.write("Note: input.rda will be determine from Int_Adding_ADT step for Cerebro step!\n")
        CEREBRO_INPUT_RDA = CEREBRO_INPUT_RDA + [ x + "_ADT.rda" for x in INT_ADD_ADT_OUTPUT]
    elif "Int_Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("Note: input.rda will be determine from Int_Clust_Markers_Annot_GE step for Cerebro step!\n")
        for i in range(len(INT_CMA_INPUT_RDA_GE)):
            for j in range(len(INT_CMA_KEEP_DIM_RES)):
                CEREBRO_INPUT_RDA.append(os.path.normpath(dic_INT_CMA_INFO[INT_CMA_INPUT_RDA_GE[i]]['INT_CMA_OUTPUT_DIR'] + "/" + INT_CMA_CLUST_FOLDERS[j] + "/" + dic_INT_CMA_INFO[INT_CMA_INPUT_RDA_GE[i]]['INT_CMA_NAME_INT'] + INT_CMA_COMPLEMENT[i] + "_" + INT_CMA_KEEP_DIM_RES[j] + ".rda"))
    if "Grp_Adding_BCR" in STEPS:
        sys.stderr.write("Note: input.rda will be determine from Grp_Adding_BCR_GE step for Cerebro step!\n")
        CEREBRO_INPUT_RDA = CEREBRO_INPUT_RDA + [ x + "_BCR.rda" for x in GRP_ADD_BCR_OUTPUT]
    elif "Grp_Adding_TCR" in STEPS:
        sys.stderr.write("Note: input.rda will be determine from Grp_Adding_TCR_GE step for Cerebro step!\n")
        CEREBRO_INPUT_RDA = CEREBRO_INPUT_RDA + [ x + "_TCR.rda" for x in GRP_ADD_TCR_OUTPUT]
    elif "Grp_Adding_ADT" in STEPS:
        sys.stderr.write("Note: input.rda will be determine from Grp_Adding_ADT step for Cerebro step!\n")
        CEREBRO_INPUT_RDA = CEREBRO_INPUT_RDA + [ x + "_ADT.rda" for x in GRP_ADD_ADT_OUTPUT]
    elif "Grp_Clust_Markers_Annot_GE" in STEPS:
        sys.stderr.write("Note: input.rda will be determine from Grp_Clust_Markers_Annot_GE step for Cerebro step!\n")
        for i in range(len(GRP_CMA_INPUT_RDA_GE)):
            for j in range(len(GRP_CMA_KEEP_DIM_RES)):
                CEREBRO_INPUT_RDA.append(os.path.normpath(dic_GRP_CMA_INFO[GRP_CMA_INPUT_RDA_GE[i]]['GRP_CMA_OUTPUT_DIR'] + "/" + GRP_CMA_CLUST_FOLDERS[j] + "/" + dic_GRP_CMA_INFO[GRP_CMA_INPUT_RDA_GE[i]]['GRP_CMA_NAME_GRP'] + GRP_CMA_COMPLEMENT[i] + "_" + GRP_CMA_KEEP_DIM_RES[j] + ".rda"))
    if len(CEREBRO_INPUT_RDA) == 0 :
        sys.exit("Error: No input.rda in configfile for Cerebro step!\n")
    ### Analysis Parameters
    CEREBRO_VERSION = getelseDefault('Cerebro', 'version', 'v1.3')
    CEREBRO_GROUPS = getelseDefault('Cerebro', 'groups', 'NULL').replace(", ", ",")
    CEREBRO_REMOVE_OTHER_RED = getelseDefault('Cerebro', 'remove.other.reductions', 'NULL')
    CEREBRO_REMOVE_OTHER_IDENT = getelseDefault('Cerebro', 'remove.other.idents', 'NULL')
    CEREBRO_REMOVE_MT = getelseDefault('Cerebro', 'remove.mt.genes', 'NULL')
    CEREBRO_REMOVE_CRB = getelseDefault('Cerebro', 'remove.crb.genes', 'NULL')
    CEREBRO_REMOVE_STR = getelseDefault('Cerebro', 'remove.str.genes', 'NULL')
    CEREBRO_ONLY_POS_DE = getelseDefault('Cerebro', 'only.pos.DE', 'NULL')
    CEREBRO_REMOVE_CUSTOM_DE = getelseDefault('Cerebro', 'remove.custom.DE', 'NULL')
    CEREBRO_ADD_SURFACE_PROT_INFO = getelseDefault('Cerebro', 'add.surface.prot.info', 'NULL')
    CEREBRO_GMT_FILE = getelseDefault('Cerebro', 'gmt.file', 'NULL')
    CEREBRO_METADATA_FILE = getelseDefault('Cerebro', 'metadata.file', 'NULL').replace(", ", ",")
    ### Snakefile parameters
    #Creation output complement + extention:
    CEREBRO_COMPLEMENT = ""
    if CEREBRO_REMOVE_MT == "TRUE" or CEREBRO_REMOVE_MT == "True": CEREBRO_COMPLEMENT = CEREBRO_COMPLEMENT + "_noMT"
    if CEREBRO_REMOVE_CRB == "TRUE" or CEREBRO_REMOVE_CRB == "True": CEREBRO_COMPLEMENT = CEREBRO_COMPLEMENT + "_noRB"
    if CEREBRO_REMOVE_STR == "TRUE" or CEREBRO_REMOVE_STR == "True": CEREBRO_COMPLEMENT = CEREBRO_COMPLEMENT + "_noSTR"
    if CEREBRO_VERSION == "v1.2":
        CEREBRO_COMPLEMENT_CRB = [CEREBRO_COMPLEMENT + "_v1.2.crb"]
        if CEREBRO_GROUPS != "NULL":
            for group in CEREBRO_GROUPS.split(','):
                CEREBRO_COMPLEMENT_CRB.append(CEREBRO_COMPLEMENT +  '_clusterIs_' + group + "_v1.2.crb")
    else:
        CEREBRO_COMPLEMENT_CRB = [CEREBRO_COMPLEMENT + ".crb"]
    #Correspondance sample/input/output
    CEREBRO_INPUT_RDA_NO_EXTENTION = [os.path.splitext(x)[0] for x in CEREBRO_INPUT_RDA]
    #Singularity environnement
    if CEREBRO_VERSION == "v1.2":
        SINGULARITY_ENV_CEREBRO = PIPELINE_FOLDER + "/envs/singularity/single_cell_oldcerebro.simg"
    elif CEREBRO_VERSION == "v1.3":
        SINGULARITY_ENV_CEREBRO = PIPELINE_FOLDER + "/envs/singularity/single_cell_cerebro.simg"
    else:
        sys.exit("Error: Unknown version of cerebro in configfile!\n")

#singularity containers
SINGULARITY_ENV = PIPELINE_FOLDER + "/envs/singularity/single_cell.simg"
SINGULARITY_ENV_INT = PIPELINE_FOLDER + "/envs/singularity/single_cell_integration.simg"
SINGULARITY_ENV_CR = PIPELINE_FOLDER + "/envs/singularity/single_cell_CellRanger.simg"
SINGULARITY_ENV_TCR_BCR = PIPELINE_FOLDER + "/envs/singularity/single_cell_TCR_BCR.simg"

#conda environments
CONDA_ENV_FASTQC = PIPELINE_FOLDER + "/envs/conda/fastqc.yaml"
CONDA_ENV_FASTQ_SCREEN = PIPELINE_FOLDER + "/envs/conda/fastq-screen.yaml"
CONDA_ENV_MULTIQC = PIPELINE_FOLDER + "/envs/conda/multiqc.yaml"
CONDA_ENV_KALLISTO = PIPELINE_FOLDER + "/envs/conda/kallisto.yaml"
CONDA_ENV_BUSTOOLS = PIPELINE_FOLDER + "/envs/conda/bustools.yaml"
CONDA_ENV_KB_PYTHON = PIPELINE_FOLDER + "/envs/conda/kb-python.yaml"
