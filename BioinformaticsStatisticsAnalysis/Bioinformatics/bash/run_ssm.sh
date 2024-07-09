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

This script takes in run the fastq2bam wdl and count reads code

Parameters:
  -i fastq dir     - path to fastq file
  -o output dir - path to output files
  -c exp_cfg  -  experiment configuration
  -h (o) debug - option to print this menu option

Usage:
    $0 -i fastq dir -o output dir -c config
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

set -o errexit
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

GENERATE_JSON="${SCRIPT_DIR}/python/generate_json.py"
GENERATE_SEQANNO_JSON="${SCRIPT_DIR}/python/generate_seqAnno_json.py"
GENERATE_POSTMAIN_JSON="${SCRIPT_DIR}/python/generate_PostMain_json.py"
FASTQ2BAM_WDL="${SCRIPT_DIR}/wdl/workflows/FastqToSAM.wdl"
SEQANNO_WDL="${SCRIPT_DIR}/wdl/workflows/SeqToAnnoVCF.wdl"
POSTMAIN_WDL="${SCRIPT_DIR}/wdl/workflows/PostMain.wdl"
COUNTREADS="${SCRIPT_DIR}/python/main.py"

##DEFAULT
CIGAR_MIN=60
##################################################
# PARSE PARAMETERS
##################################################=

DEBUG=""

while getopts "i:o:c:dh" option
do
    case ${option} in
        i) FASTQ_DIR=${OPTARG};;
        o) OUT_DIR=${OPTARG};;
	c) EXP_CFG=${OPTARG};;
        d) DEBUG="-d" ; set -x ;;
        h) ${ECHO} "${docs}" ; exit 1 ;;
        ?) ${ECHO} "${docs}" ; exit 1 ;;
    esac
done

validateDirectory "${FASTQ_DIR}"
validateDirectory "${OUT_DIR}"


logInfo "${MANIFEST}"
logInfo "Firing the program: hold on to your beverage"

##################################################
# MAIN
##################################################


function check_job_limit()
{
    job_limit=$JOB_LIMIT
    job_count=$(qstat | sed 1,2d | wc -l)
    while [[ $job_count -ge $job_limit ]]
    do
        echo "User job limit of ${job_limit} reached. Waiting for jobs to finish before submitting more..."
        sleep 30s
        job_count=$(qstat | sed 1,2d | wc -l)
    done
}

function generate_json_cfg {
    echo ${CUTADAPT_PARAM}
    cmd="${PYTHON} ${GENERATE_JSON} -f ${FASTQ_DIR}  -o ${OUT_DIR} -c '${CUTADAPT_PARAM}' -r ${REFERENCE} -j ${SCRIPT_DIR}/wdl/Copy_of_input.json -s '${SEQPREP_OPTIONS}'"
    echo ${cmd}
    eval ${cmd}
    JSON_OUT=`find $OUT_DIR -name "*input.json"`
    echo ${JSON_OUT}
    if [[ -z ${JSON_OUT} ]];then
	logError "No json files created"
    else
	for json in ${JSON_OUT};do
	    logInfo "Parsing json ${json}"
	    SAMPLE_NAME=`basename ${json} | cut -d"." -f1`
	    check_outdir ${SAMPLE_NAME}
	    run_cromwell ${json} ${SAMPLE_NAME} ${OUT_DIR}/${SAMPLE_NAME}
	done
    fi	
}


function run_seqanno {
    cmd="${PYTHON} ${GENERATE_SEQANNO_JSON} -c ${CHROM} -p ${START} -s  ${SEQUENCE} -j ${SCRIPT_DIR}/wdl/SeqToAnnoVCF.json -o ${OUT_DIR}"
    echo ${cmd}
    eval ${cmd}
    JSON_OUT=`ls ${OUT_DIR}/SSM_Allpossible.seqAnno.input.json`
    if [[ -z ${JSON_OUT} ]];then
        logError "No json files created"
    else
	cmd="$QSUB -N SeqAnno ${SEQANNO_MEM} -e ${OUT_DIR} -o ${OUT_DIR} ${CROMWELL} run ${SEQANNO_WDL} --inputs ${JSON_OUT}"
	logInfo "running ${cmd}"
	jobID=$(${cmd})
	logInfo "Run ${jobID} for sample ${SAMPLE_NAME} with json ${JSON}"
    fi
}


function check_outdir {
    SAMPLE_NAME=$1
    SAMPLE_OUT=${OUT_DIR}/${SAMPLE_NAME}
    if [[ -d ${SAMPLE_OUT} ]];then
	logInfo "${SAMPLE_OUT} already exists"
	rm -rf ${SAMPLE_OUT}
	mkdir ${SAMPLE_OUT}
    else
	mkdir ${SAMPLE_OUT}
    fi
}

function check_running_CromwellS {
    VALUE=$1
    IS_STILL_RUNNING=`${QSTAT} | grep "CromwellS" | wc -l`
    echo $IS_STILL_RUNNING "jobs on the grid waiting to complete"
    while [ "$IS_STILL_RUNNING" -gt 0 ]
    do
        sleep 20s
        IS_STILL_RUNNING=`${QSTAT} | grep "CromwellS" | wc -l`
        echo $IS_STILL_RUNNING "jobs on the grid waiting to complete"
    done
    sleep 5s
}


function check_running_SSM {
    VALUE=$1
    IS_STILL_RUNNING=`${QSTAT} | grep "SSM_count" | wc -l`
    echo $IS_STILL_RUNNING "jobs on the grid waiting to complete"
    while [ "$IS_STILL_RUNNING" -gt 0 ]
    do
        sleep 20s
        IS_STILL_RUNNING=`${QSTAT} | grep "SSM_count" | wc -l`
        echo $IS_STILL_RUNNING "jobs on the grid waiting to complete"
    done
    sleep 5s
}


function run_cromwell {
    JSON=$1
    SAMPLE_NAME=$2
    OUTDIR=$3
    cmd="$QSUB -N CromwellSSM.${SAMPLE_NAME} ${FASTQ_MEM} -e ${OUTDIR} -o ${OUTDIR} ${CROMWELL} run ${FASTQ2BAM_WDL} --inputs ${JSON}"
    logInfo "running ${cmd}"
    jobID=$(${cmd})
    logInfo "Run ${jobID} for sample ${SAMPLE_NAME} with json ${JSON}"
}


function count_reads {
    SAM_FILE=`find ${OUT_DIR} -name "*.fastq.gz.sam"`
    echo ${SAM_FILE}
    for sam in ${SAM_FILE};do
	echo ${sam}
	if [[ -f ${sam} ]];then
	    logInfo "${sam} was generated"
	    cmd="${QSUB} -N SSM_count ${COUNTREADS_MEM} -e ${OUT_DIR} -o ${OUT_DIR} ${PYTHON} ${COUNTREADS} -s ${sam} -c ${CHROM} -p ${POSITIONS} -b ${BASES} -r ${SEQUENCE} -z ${STRAND} -d ${CODON_DISTANCE_UP} -P ${TARGET_BOUNDS} -D ${CODON_DISTANCE_DOWN} -m ${CIGAR_MIN} -S ${START}"
	    logInfo "running ${cmd}"
	    standAlone_cmd="${PYTHON} ${COUNTREADS} -s ${sam} -c ${CHROM} -p ${POSITIONS} -b ${BASES} -r ${SEQUENCE} -z ${STRAND} -d ${CODON_DISTANCE_UP} -P ${TARGET_BOUNDS} -D ${CODON_DISTANCE_DOWN} -m ${CIGAR_MIN} -S ${START}"
	    echo ${cmd}
	    eval ${cmd}
	else
	    logError "sam file did not generate"
	fi
    done
}



function run_postmain {
    VCF_LIST=`find ${DATA_PATH} -name "*.fastq.all.vcf" `
    for vcf in ${VCF_LIST};do
	NAME=`basename ${vcf} | cut -d"." -f1`
	cmd="${PYTHON} ${GENERATE_POSTMAIN_JSON} -s ${NAME} -v ${vcf} -j ${SCRIPT_DIR}/wdl/PostMain.json -o ${DATA_PATH}"
	echo ${cmd}
	eval ${cmd}
	JSON_OUT=`ls ${DATA_PATH}/${NAME}.PostMain.input.json`
	if [[ -z ${JSON_OUT} ]];then
           logError "No json files created"
	else
	    cmd="$QSUB -N SeqAnno ${SEQANNO_MEM} -e ${DATA_PATH} -o ${DATA_PATH} ${CROMWELL} run ${POSTMAIN_WDL} --inputs ${JSON_OUT}"
	    logInfo "running ${cmd}"
	    jobID=$(${cmd})
	    logInfo "Run ${jobID} for sample ${SAMPLE_NAME} with json ${JSON}"
	fi
   done
}

if [[ -f ${EXP_CFG} ]];then
    source ${EXP_CFG}
#    validateParm ${DATA_PATH}
    logInfo "Processing data in ${DATA_PATH}"
    run_seqanno
    check_outdir
    generate_json_cfg
    check_running_CromwellS
   count_reads
    check_running_SSM
    run_postmain
else
    logError "No Experiment config supplied"
fi

