#!~/bin/env bash

##################################################
# MANIFEST
##################################################

read -r -d '' MANIFEST <<MANIFEST
`basename ${0}`
*******************************************
`readlink -m $0` ${@}
was called by: `whoami` on `date`
*******************************************
MANIFEST
echo "${MANIFEST}"

read -r -d '' docs <<docs

This script takes in runs the R analysis script

Parameters:
  -i config - experiment config (required)
  -t git dir - path to the SSM code repo 
  -f flowcell - flowcell ID (required)
  -o optimal date - day for analysis
  -h (o) debug - option to print this menu option

Usage:
    $0 -i config -t git_dir -f flowcell -o D14_cont_rat
docs

ECHO="/bin/echo"
#Show help when no parameters are provided
if [ $# -eq 0 ];
then
    ${ECHO} "${docs}" ; exit 1 ;
fi


##################################################
# SET BASH OPTIONS
##################################################

#set -o errexit
#set -o pipefail
#set -o nounset

#################################################
# FIND AND SOURCE PROFILE
##################################################

echo ""
echo "***** Prepping profile properly *****"

if [[ -z ${PROFILE:-} ]]; then
    SCRIPT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    SCRIPT_NAME=$(basename "${SCRIPT}")
    SCRIPT_DIR=$(dirname "${SCRIPT}")
    PROFILE_DIR="${SCRIPT_DIR}/config"
    PROFILE=${PROFILE_DIR}/SSM_config.txt
fi

if [[ ! -f "${PROFILE}" ]]; then
    echo "$(basename ${PROFILE}) was not found. Unable to continue: ${PROFILE}"
    exit 1
fi

echo "Using configuration file at: ${PROFILE}"
source "${PROFILE}"

COMMON_FUNC="${SCRIPT_DIR}/bash/common_functions/commonFunctions.sh"
if [[ ! -f ${COMMON_FUNC} ]]; then
    ${ECHO} -e "\nERROR: Common function UTILITIES not found: ${COMMON_FUNC}\n"
    exit 1
fi
${ECHO} "Using common functions: ${COMMON_FUNC}"
source "${COMMON_FUNC}"

RSCRIPT="~/r/R-4.0.3/bin/Rscript"
RANALYSIS="${SCRIPT_DIR}/R/EventPlots.codonSpliceAI.Rmd"

##################################################
# PARSE PARAMETERS
##################################################=

DEBUG=""

while getopts "i:t:f:o:dh" option
do
    case ${option} in
	i) CONFIG=${OPTARG};;
	t) GIT_DIR=${OPTARG};;
	f) FLOWCELL=${OPTARG};;
	o) OPTIMAL_DAY=${OPTARG};;
        d) DEBUG="-d" ; set -x ;;
        h) ${ECHO} "${docs}" ; exit 1 ;;
        ?) ${ECHO} "${docs}" ; exit 1 ;;
    esac
done


#validateParm "${SAMPLE_PREFIX}"
#validateParm "${AA_BEGIN}"
validateFile "${CONFIG}"
validateParm "${FLOWCELL}"
#validateParm "${ANALYSIS_DAY}"

if [[ -z ${RSTUDIO_PANDOC} ]];then
  RSTUDIO_PANDOC="~/pandoc/2.7.3/bin/"
fi

if [[ -z ${GIT_DIR} ]];then
    GIT_DIR="~/.git"
else
    GIT_DIR=${GIT_DIR}/.git
fi

if [[ -z ${OPTIMAL_DAY} ]];then
    OPTIMAL_DAY="D14_cont_rat"
else
    OPTIMAL_DAY=${OPTIMAL_DAY}
fi


source ${CONFIG}
SAMPLE_PREFIX=${SAMPLE_PREFIX}
OUT_DIR=${DATA_PATH}
AA_begin=${AA_BEGIN}
LEVELS=${LEVELS}
TRUTH_PATH=${TRUTH_DATA}
GENE=${GENE}
validateDirectory "${GIT_DIR}"
OUT_FILE="EventPlots."${GENE}"."${SAMPLE_PREFIX}"codonSpliceAI.Rmd"
function copy_rfile {
    proj_Rmd=${OUT_DIR}/${OUT_FILE}
    cmd="cp ${RANALYSIS} ${OUT_DIR}/${OUT_FILE}"
    if [[ -f ${proj_Rmd} ]];then
	rm -rf ${proj_Rmd}
	logInfo "running command ${cmd}"
	eval ${cmd}
    else
	logInfo "running command ${cmd}"
	eval ${cmd}
    fi
}

function distance_to_start {
    target_array=($TARGET_BOUNDS)
    declare -a distance=()
    for pos in ${POSITIONS};do
	echo ${pos}
	if [[ ${pos} > ${target_array[0]} && ${pos} < ${target_array[1]} ]];then
	    DISTANCE_TO_START="0"
	    distance+=(${DISTANCE_TO_START})
	else
	    DISTANCE_TO_START_LOW=$(expr ${pos} - ${target_array[0]})
	    DISTANCE_TO_START_HIGH=$(expr ${pos} - ${target_array[1]})
	    if [[ ${DISTANCE_TO_START_LOW#-} < ${DISTANCE_TO_START_HIGH#-} ]];then
		distance+=(${DISTANCE_TO_START_HIGH})
	    else
		distance+=(${DISTANCE_TO_START_LOW})
	    fi
	fi
    done
}


function insert_param {
    GIT_COMMIT=$(git --git-dir=${GIT_DIR} log -1 --pretty="%H")
#    eval ${GIT_COMMIT}
    commit_cmd='sed -i "52i COMMIT_ID='"'${GIT_COMMIT}'"' " ${OUT_DIR}/${OUT_FILE}'
    data_path_cmd='sed -i "51i data_path='"'${OUT_DIR}'"' " ${OUT_DIR}/${OUT_FILE}'
    CHROM_CMD='sed -i "53i CHROM='"'${CHROM}'"' " ${OUT_DIR}/${OUT_FILE}'
    SEQUENCE_CMD='sed -i "54i SEQUENCE='"'${SEQUENCE}'"' " ${OUT_DIR}/${OUT_FILE}'
    BASES_CMD='sed -i "55i BASES='"'${BASES}'"' " ${OUT_DIR}/${OUT_FILE}'
    POSITION_CMD='sed -i "56i POSITIONS='"'${POSITIONS}'"' " ${OUT_DIR}/${OUT_FILE}'
    STRAND_CMD='sed -i "57i STRAND='"'${STRAND}'"' " ${OUT_DIR}/${OUT_FILE}'
    TARGET_CMD='sed -i "60i TARGET_BOUNDS='"'${TARGET_BOUNDS}'"' " ${OUT_DIR}/${OUT_FILE}'    
    CODON_DISTANCE_UP_CMD='sed -i "61i CODON_DISTANCE_UP='"'${CODON_DISTANCE_UP}'"' " ${OUT_DIR}/${OUT_FILE}'    
    CUTADAPT_PARAM_CMD='sed -i "62i CUTADAPT_PARAM='"'${CUTADAPT_PARAM}'"' " ${OUT_DIR}/${OUT_FILE}'    
    CODON_DISTANCE_DOWN_CMD='sed -i "63i CODON_DISTANCE_DOWN='"'${CODON_DISTANCE_DOWN}'"' " ${OUT_DIR}/${OUT_FILE}'    
    prefix_cmd='sed -i "65i SAMPLE_PREFIX='"'${SAMPLE_PREFIX}'"' " ${OUT_DIR}/${OUT_FILE}'
    AA_begin_cmd='sed -i "67i AA_begin='"${AA_begin}"' " ${OUT_DIR}/${OUT_FILE}'
    LEVELS_cmd='sed -i "68i LEVELS=c('"'${LEVELS}'"') " ${OUT_DIR}/${OUT_FILE}'
    TRUTH_cmd='sed -i "69i TRUTH=read.csv('"'${TRUTH_PATH}'"', sep='"'\\\\\t'"') " ${OUT_DIR}/${OUT_FILE}'
    START_cmd='sed -i "70i START='"${START}"' " ${OUT_DIR}/${OUT_FILE}'
    SEQPREP_cmd='sed -i "71i SEQPREP_OPTIONS='"'${SEQPREP_OPTIONS}'"' " ${OUT_DIR}/${OUT_FILE}'
    OPTIMAL_DAY_CMD='sed -i "72i Optimal_Date='"'${OPTIMAL_DAY}'"' " ${OUT_DIR}/${OUT_FILE}'
    GENE_cmd='sed -i "73i GENE='"'${GENE}'"' " ${OUT_DIR}/${OUT_FILE}'
    DISTANCE=$(distance_to_start)
    echo "DISTANCE_TO_START" ${DISTANCE}
    DISTANCE_TO_START_cmd='sed -i "72i DISTANCE_TO_START='"'${DISTANCE}'"' " ${OUT_DIR}/${OUT_FILE}'
    FLOWCELL_cmd='sed -i "74i FLOWCELL='"'${FLOWCELL}'"' " ${OUT_DIR}/${OUT_FILE}'
    echo ${commit_cmd}
    echo ${CHROM_CMD}
    echo ${SEQUENCE_CMD}
    echo ${BASES_CMD}
    echo ${POSITION_CMD}
    echo ${STRAND_CMD}
    echo ${TARGET_CMD}
    echo ${CODON_DISTANCE_UP_CMD}
    echo ${CUTADAPT_PARAM_CMD}
    echo ${CODON_DISTANCE_DOWN_CMD}
    echo ${prefix_cmd}
    echo ${data_path_cmd}
    echo ${AA_begin_cmd}
    echo ${LEVELS_cmd}
    echo ${TRUTH_cmd}
    echo ${START_cmd}
    echo ${SEQPREP_cmd}
    echo ${DISTANCE_TO_START_cmd}
    echo ${OPTIMAL_DAY_CMD}
    eval ${commit_cmd}
    eval ${data_path_cmd}
    eval ${CHROM_CMD}
    eval ${SEQUENCE_CMD}
    eval ${BASES_CMD}
    eval ${POSITION_CMD}
    eval ${STRAND_CMD}
    eval ${TARGET_CMD}
    eval ${CODON_DISTANCE_UP_CMD}
    eval ${CUTADAPT_PARAM_CMD}
    eval ${CODON_DISTANCE_DOWN_CMD}
    eval ${prefix_cmd}
    eval ${AA_begin_cmd}
    eval ${LEVELS_cmd}
    eval ${TRUTH_cmd}
    eval ${START_cmd}
    eval ${DISTANCE_TO_START_cmd}
    eval ${SEQPREP_cmd}
    eval ${OPTIMAL_DAY_CMD}
    eval ${GENE_cmd}
    eval ${FLOWCELL_cmd}
}

function print_tracker {
    DISTANCE=$(distance_to_start)
    GIT_COMMIT=$(git --git-dir=${GIT_DIR} log -1 --pretty="%H")
    FILENAME=${GENE}"_"${SAMPLE_PREFIX}${FLOWCELL}".tsv"
    echo -e "CHROM'\t'SEQUENCE'\t'BASES'\t'POSITIONS'\t'STRAND'\t'TARGET'\t'CODON_DISTANCE_UP'\t'CODON_DISTANCE_DOWN'\t'CUTADAPT'\t'AA_BEGIN'\t'START'\t'SEQPREP'\t'DISTANCE_TO_START'\t'FLOWCELL'\t'PROCESSING_DIRECTORY'\t'GIT_COMMIT'\t'Exon" > ${FILENAME}
    echo -e "${CHROM}'\t'${SEQUENCE}'\t'${BASES}'\t'${POSITIONS}'\t'${STRAND}'\t'${TARGET_BOUNDS}'\t'${CODON_DISTANCE_UP}'\t'${CODON_DISTANCE_DOWN}'\t'${CUTADAPT_PARAM}'\t'${AA_BEGIN}'\t'${START}'\t'${SEQPREP_OPTIONS}'\t'${DISTANCE}'\t'${FLOWCELL}'\t'${OUT_DIR}'\t'${GIT_COMMIT}'\t'${SAMPLE_PREFIX}" >> ${FILENAME}
}

function run_markdown {
    pushd ${OUT_DIR}
    cmd="${RSCRIPT} -e 'Sys.setenv(RSTUDIO_PANDOC="'"'${RSTUDIO_PANDOC}'"'");rmarkdown::render("'"'${OUT_FILE}'"'")'"
    echo ${cmd}
    eval ${cmd}
    popd ${OUT_DIR}
}

copy_rfile
insert_param
print_tracker
run_markdown

