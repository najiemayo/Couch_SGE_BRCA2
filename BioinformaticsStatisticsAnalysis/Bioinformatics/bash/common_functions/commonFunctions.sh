#!/bin/bash

PARENT_COMMAND="$(basename ${0})"


##################################################
# Exit Codes
##################################################

INF_OK=0
ERR_GENERAL=1
ERR_MISSINGINPUT=2
ERR_DIRNOTFOUND=5
ERR_FILENOTFOUND=10
ERR_UNEXPECTED=50


##################################################
# Color Codes
##################################################

#main colors
BLACK='\033[0;30m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
GRAY='\033[1;30m'
GREEN='\033[0;32m'
ORANGE='\033[0;33m'
PURPLE='\033[0;35m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
WHITE='\033[1;37m'

#light colors
LCYAN='\033[1;36m'
LGRAY='\033[0;37m'
LGREEN='\033[1;32m'
LRED='\033[1;31m'
LBLUE='\033[1;34m'
LPURPLE='\033[1;35m'

#no color
NC='\033[0m'


##################################################
# DEFAULT VALUES ALL PIPES SHOULD HAVE 
##################################################

if [[ -z ${EMAIL_ADDRESS:-} ]]; then
    EMAIL_ADDRESS="$(~/bin/getent passwd `~/bin/whoami` | ~/bin/cut -d ':' -f 5 | ~/bin/cut -d ';' -f 2)"
fi

LOG_MODE="echo"


##################################################
# Functions
##################################################

function getDate()
{
    echo "$(date +%Y-%m-%d'T'%H:%M:%S%z)"
}

function setupLogging()
{
    LOG_DIR="${1:-${TEMP_LOG_SPACE}}"
    LOG_FILE="${LOG_DIR}/${2:-${SCRIPT_NAME%.*}.log_$$}" #use name, no extension + pid for uniqueness
    LOG_MODE="${3:-echo}"
}
 
# This is "private" function called by the other logging functions, don't call it directly,
# use logError, logWarn, etc.
function _logMsg () {
    if [[ ${LOG_MODE} == "echo" ]] || [[ ${LOG_MODE} == "all" ]]; then
        echo -e "${1}"
    fi

    if [[ ${LOG_MODE} == "file" ]] || [[ ${LOG_MODE} == "all" ]]; then
        echo -e "${1}" | sed -r 's/\\n//'  >> "${LOG_FILE}"
    fi
}
 
function logError()
{
    >&2 _logMsg "[$(getDate)] [ERROR] [${PARENT_COMMAND}] \t${1}"
}

function logDebug()
{
    _logMsg "[$(getDate)] [DEBUG] [${PARENT_COMMAND}] \t${1}"
}
 
function logWarn()
{
    _logMsg "[$(getDate)] [WARN] [${PARENT_COMMAND}] \t${1}"
}
 
function logInfo()
{
    _logMsg "[$(getDate)] [INFO] [${PARENT_COMMAND}] \t${1}"
}

function die() {
    logError "Unexpected error occurred: ${ERR_UNEXPECTED}"
    exit ${ERR_UNEXPECTED};
}

function presskey () {
    read -n 1 -s -p $'Press any key to continue\n'
}

function mailUser () {
	local EMAIL_TO="${1}"
	local SUBJECT="${2}"
	local MESSAGE="${3}"
	local CURRENT_USER="$(~/bin/whoami)"
	local SYSTEM_USER="$(/bin/echo ${CURRENT_USER} | /bin/grep -e "^tu" -e "^wa")"
	local CMD=""

	if [[ ${CURRENT_USER} != ${SYSTEM_USER} ]]; then
		local USER_EMAIL
		USER_EMAIL="$(~/bin/getent passwd `~/bin/whoami` | /bin/cut -d ':' -f 5 | /bin/cut -d ';' -f 2)"
		if [[ ${USER_EMAIL,,} =~ "USEREmail" ]] &&
		   [[ ! ${EMAIL_TO,,} =~ ${USER_EMAIL,,} ]]; then
			EMAIL_TO="${EMAIL_TO},${USER_EMAIL}"
		fi
	fi
	
    CMD="/bin/echo -e \"${MESSAGE}\" | /bin/mail -s \"${SUBJECT}\" \"${EMAIL_TO}\""
    runCmd "${CMD}"
    sleep 1
}

function getIndex() {
    IFS=',' read -r -a array <<< $(head -1 ${1})

    for index in "${!array[@]}"
    do
        if [[ ${array[index],,} == "${2,,}" ]]; then
            break
        fi
    done

    if [[ ${array[index],,} != "${2,,}" ]]; then
        logError "Required index NOT FOUND: ${2}"
        exit ${ERR_GENERAL}
    fi

    index=$(($index+1))
    echo "$index"
}

function runCmd () {
        logInfo "Running: ${1}\n"
        if [[ -f ${LOG_FILE:-} ]]; then
            eval ${1} | tee -a ${LOG_FILE} 2>&1
        else
            eval ${1} 2>&1
        fi
}

function validateParm () {
    local TPARM="${1}"

    if [[ -z "${TPARM:-}" ]]; then
        local TMSG="ERROR: Required parameter has not been defined: ${1}\n${2:-}"
        logError "${TMSG}" ${ERR_GENERAL}
        exit ${ERR_GENERAL}
    fi
}

function validateFile () {
    local TFILE="${1}"

    if [[ ! -f "${TFILE}" ]]; then
        local TMSG="ERROR: File does not exist: ${1}\n${2:-}"
        logError "${TMSG}" ${ERR_FILENOTFOUND} 
        exit ${ERR_FILENOTFOUND}
    fi
}

function validateDirectory () {
    local TDIR="${1}"

    if [[ ! -d "${TDIR}" ]]; then
        local TMSG="ERROR: Directory does not exist: ${1}\n${2:-}"
        logError "${TMSG}" ${ERR_DIRNOTFOUND}
        exit ${ERR_DIRNOTFOUND}
    fi
}

function validateInteger () {
    local TSTR="${1}"

    if [[ ! ${TSTR} =~ ^-?[0-9]+$ ]]; then
        local TMSG="ERROR: Value is NOT an integer: ${1}\n${2:-}"
        logError "${TMSG}" ${ERR_GENERAL}
        exit ${ERR_GENERAL}
    fi
}

function readOption () {
    eval local t="$(awk "/^${2}/{print \$NF; exit;}" "${1}" | sed -r "s/^${2}=(.*)/\1/")"
    echo "${t}"
}

#requires SSH set, optional REMOTE_USER REMOTE_HOST
function runRmtCmd () {
    local RMT_RESULT=""
    RMT_RESULT=$(${SSH} ${REMOTE_USER:-${2}}@${REMOTE_HOST:-${3}} "${1}")
    echo $RMT_RESULT
}

function getCommonPath()
{
  lhs=$1
  rhs=$2
  path=
  OLD_IFS=$IFS; IFS=/
  for w in $rhs; do
    test "$path" = / && try="/$w" || try="$path/$w"
    case $lhs in
      $try*) ;;
      *) break ;;
    esac
    path=$try
  done
  IFS=$OLD_IFS
  echo $path
}

function common_path_all()
{
  local sofar=$1
  shift
  for arg
  do
    sofar=$(getCommonPath "$sofar" "$arg")
  done
  echo ${sofar:-/}
}

function createRelativeLink () {
    if [[ ! -e ${1} ]]; then
        logInfo "Source does not exist: ${1}"
        exit ${ERR_GENERAL}
    fi

    if [[ ! -e ${2} ]]; then
        logInfo "Destination does not exist: ${1}"
        exit ${ERR_GENERAL}
    fi

    local SOURCE_NAME="$(basename ${1})"
    local SOURCE_DIR="$(dirname ${1})"
    local DEST_DIR="${2}"
    local COMMON_PATH=""
    local TPATH=""
    local OFFSHOOT_PATH=""

    COMMON_PATH="$(getCommonPath "${1}" "${2}")"

    #if paths have a common tree
    if [[ "${COMMON_PATH}" != "/" ]]; then
        TPATH="${DEST_DIR#"${COMMON_PATH}"}"
    else
        TPATH="${SOURCE_DIR}"
    fi

    #add branch
    if [[ ${SOURCE_DIR} =~ ^${COMMON_PATH} ]]; then
        OFFSHOOT_PATH="$(echo ${SOURCE_DIR} | sed "s|${COMMON_PATH}||")"
    fi

    #target is in a subdir (source contains dest)
    if [[ ${SOURCE_DIR} =~ ^${DEST_DIR} ]]; then
        OFFSHOOT_PATH=".${OFFSHOOT_PATH}"
    fi     

    #replace path with relative path ../.., append with alternate path /peers
    if [[ ${COMMON_PATH} != "/" ]]; then
        TPATH="$( echo ${TPATH} | sed 's|^/||' | sed -r 's|\w+-?\w+|..|g')${OFFSHOOT_PATH}"
    fi

    cd ${DEST_DIR}
    local CMD="ln -s ${TPATH}/${SOURCE_NAME} ${3:-${SOURCE_NAME}}"
    logInfo "Creating relative link: ${CMD}"
    eval ${CMD}
}

#requires QSTAT var is set
function waitForJob () {
    local WAIT_TIME=0
    local MAX_WAIT=${2:-10800}
    local SLEEP_TIME=${3:-120}

    local JOB_ID="${1}"
    local QSTAT="~/soge/8.1.9b/bin/lx-amd64/qstat"
    local JOB_RUNNING=$(${QSTAT} -u '*' | gawk -F " " -v JID="${JOB_ID}" '$1==JID{print $1}')

    logInfo "Waiting for grid job ${JOB_ID} to complete. MAX_WAIT: ${MAX_WAIT}"

    echo ${JOB_RUNNING}
    while [[ ${JOB_RUNNING} != "" ]]
    do
        if [[ ${WAIT_TIME} -gt ${MAX_WAIT} ]]; then
            logError "Job processing exceeded timeout value. Job ID: ${JOB_ID}" ${ERR_GENERAL}
            exit ${ERR_GENERAL}
        fi

        QSTS=$(${QSTAT} -u '*' | grep ${JOB_ID} || true)
        QSTS=$(echo ${QSTS} | tr -s ' ' | cut -f5 -d " ")
        if [[ ${QSTS} == "Eqw" ]]; then
            logError "QSTAT indicates job ${JOB_ID} failed." ${ERR_GENERAL}
            exit ${ERR_GENERAL}
        fi

        logInfo "Sleeping for ${SLEEP_TIME} seconds"
        sleep ${SLEEP_TIME}
        WAIT_TIME=$(($WAIT_TIME + $SLEEP_TIME))

        JOB_RUNNING=$(${QSTAT} -u '*' | gawk -F " " -v JID="${JOB_ID}" '$1==JID{print $1}')
    done
}

updateWorkbenchFile () {
    local OUTPUT_FILE="${1:-}"
    local TIME="$(date +%Y-%m-%dT%H:%M:%S)"
    local STATE="${2:-}"
    local COMMENT="${3:-}"

    if [[ -z ${OUTPUT_FILE} ]]; then
        logError "OUTPUT_FILE is a required parm"
        exit ${ERR_GENERAL}
    fi

    if [[ -z ${STATE} ]]; then
        logError "STATE is a required parm"
        exit ${ERR_GENERAL}
    fi

    if [[ -z ${COMMENT} ]]; then
        logError "COMMENT is a required parm"
        exit ${ERR_GENERAL}
    fi

    if [[ ! -f ${OUTPUT_FILE} ]]; then
        echo -e "##TIME\tSTATE\tCOMMENT" > "${OUTPUT_FILE}"
    fi
    echo -e "${TIME}\t${STATE}\t${COMMENT}" >> "${OUTPUT_FILE}"
}

# Get current time data for creating run time .metrictsv file
#
# input:	METRIC_NAME - name of metric
#			RUN_ID		- run id
# output:	tab separated data as follows
#           run_time <METRIC_NAME> C <current-time> <RUN_ID> run
function getCurrentTimeForRunTimeMetric() {
	local METRIC_NAME="${1:-}"
	local RUN_ID="${2:-}"

	validateParm "METRIC_NAME"  "METRIC_NAME is a required parm"
	validateParm "RUN_ID" "RUN_ID is a required parm"

	local TIME=$(date +%Y-%m-%d'T'%H:%M:%S%z)
	printf '%s\t%s\t%s\t%s\t%s\t%s\n' "run_time" "${METRIC_NAME}" "C" "${TIME}" "${RUN_ID}" "run"
}

# Get file modify time data for creating run time .metrictsv file
#
# input: 	FILE_NAME 	- fully qualified file name
# 			METRIC_NAME - name of metric
#			RUN_ID		- run id
# output:	tab separated data as follows
#			run_time <METRIC_NAME> C <file-mtime> <RUN_ID> run
function getFileMTimeForRunTimeMetric() {
	local FILE_NAME="${1:-}"
	local METRIC_NAME="${2:-}"
	local RUN_ID="${3:-}"

	validateParm "FILE_NAME" "FILE_NAME is a required parm"
	validateFile "${FILE_NAME}"
	validateParm "METRIC_NAME"  "METRIC_NAME is a required parm"
	validateParm "RUN_ID" "RUN_ID is a required parm"

	local TIME=$(ls -l --time-style=+%Y-%m-%d'T'%H:%M:%S%z ${FILE_NAME} | awk '{print $6}')
	printf '%s\t%s\t%s\t%s\t%s\t%s\n' "run_time" "${METRIC_NAME}" "C" "${TIME}" "${RUN_ID}" "run"
}

# Get file change time data for creating run time .metrictsv
#
# input: 	FILE_NAME 	- fully qualified file name
#       	METRIC_NAME - name of metric
#			RUN_ID		- run id
# output:	tab separated data as follows
#			run_time <METRIC_NAME> C <file-ctime> <RUN_ID> run
function getFileCTimeForRunTimeMetric() {
	local FILE_NAME="${1:-}"
	local METRIC_NAME="${2:-}"
	local RUN_ID="${3:-}"

	validateParm "FILE_NAME" "FILE_NAME is a required parm"
	validateFile "${FILE_NAME}"
	validateParm "METRIC_NAME"  "METRIC_NAME is a required parm"
	validateParm "RUN_ID" "RUN_ID is a required parm"

	local TIME=$(ls -l --time=ctime --time-style=+%Y-%m-%d'T'%H:%M:%S%z ${FILE_NAME} | awk '{print $6}')
	printf '%s\t%s\t%s\t%s\t%s\t%s\n' "run_time" "${METRIC_NAME}" "C" "${TIME}" "${RUN_ID}" "run"
}

# validate chromosomes of a bam file givin a chromose file
#
# input: $1 - input bam filee
#        $2 - input chromosome file
#
# output: 0 - chromosomes of the input bam file is a subset of the input chromosome file (chromosomes of the input bam are valid )
#         1 - at least one of the chromosomes of the input bam file is not found in the input chromosome file (chromosomes of the input bam are invalid )
#
# tools required: (you can source src/cnv.profile to make these tools available)
#                 grep
#                 samtools
#                 cut
#                 sort
#                 uniq
#                 wc
# example calling:
# isValidChromoses=$(validateBamChromosomes /dlmp/qa/runs/IEMCP/IEMCP-NGS69_CAPX-NGS82-VER-0607_CAPX-200407-002_IL2829-P2_HKGYHDRXX_shengbing/samples/NGS82-VE055/mgc/alignment/NGS82-VE055.bam /dlmp/misc-data/pipelinedata/deployments/caps/genomes/GRCh37/hs37d5.chromosomes)
function validateBamChromosomes() {
  # check if each chromose of a bam file is present in a given file
  # $1=bam
  # $2=chromosome file
  local number1
  local number2 # number of chromosomes of the input bam

  number1=$(
    set +e;
    ${GREP} -c -w -f <(
      ${SAMTOOLS} view -H "$1" | ${GREP} '@SQ' | ${CUT} -f2 | ${CUT} -d: -f2 | ${SORT} -V | ${UNIQ}) "$2";
    set -e;
  )

  # number will never be unset given presence of "wc -l" at the end
  number2=$(set +e;
    ${SAMTOOLS} view -H "$1" | ${GREP} '@SQ' | ${CUT} -f2 | ${CUT} -d: -f2 | ${SORT} -V | ${UNIQ} | ${WC} -l;
    set -e;
  )

  if [[  ${number2} -eq 0 ]]; then
    ${ECHO} 1
    exit ${ERR_GENERAL}
  fi

  if [[ ${number1:-} -eq ${number2:-} ]]; then
    ${ECHO} 0
  else
    ${ECHO} 1
  fi
}
