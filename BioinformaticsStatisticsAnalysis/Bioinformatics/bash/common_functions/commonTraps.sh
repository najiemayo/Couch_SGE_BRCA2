#!/bin/bash

# This file can be sourced in any script that needs to trap common signals
# for error reporting.
#
# By default, the traps will echo an error string by calling the reportError
# function which includes the script name, exit code, last line number
# executed, and signal name trapped.
#
# A script may overload any of the traps or reportError function as well as
# define additional traps.  However, you should consider this:
# 1. The DEBUG trap sets LAST_LINENO to the last line executed since
#    bash 4.x does not correctly set BASH_LINENO.
# 2. The reportError function must be the first thing executed in a trap
#    to ensure LAST_LINENO is correct.
# 3. System traps should unset EXIT trap then set the associated exit code
#    (128 + signal number).  The unset of the EXIT trap is to ensure that
#    reportError is not called twice.


##################################################
#Variables
##################################################

LAST_LINENO=""
REAL_LINENO=""


##################################################
# Functions
##################################################

function reportError() {
    local exit_code=${1:-0}
    local trap_string=${2:-""}
    if [[ ${exit_code} -ne "0" ]]; then
        echo "ERROR: ${0}  exit_code=${exit_code}  LAST_LINENO=${LAST_LINENO}  trap_string=${trap_string}"
    fi
}

function enableCommonTraps() {

	# Bash provides these psuedo-signals:
	# - DEBUG   any shell command
	# - EXIT    any exit from shell
	# - ERR     non zero exit from shell
	# - RETURN  any return, any script executed with source or . builtin

	# This will ensure LAST_LINENO is correct when reportError is called.
	trap 'LAST_LINENO=$REAL_LINENO; REAL_LINENO=$LINENO' DEBUG

	# Trap EXIT to catch all bash exits
	trap 'reportError $? EXIT' EXIT

	# Intentionaly do not trap ERR since it only catches non-zero exits
	# from called commands or scripts but not from the current shell.
	# The reportError function checks for non-zero exit codes.


	# System signals can be listed by invoking  'trap -l' or 'kill -l'.
	# For signal details, invoke 'man 7 signal'.

	# Trap system signals associated with processes ending abnormally and exit
	# with associated code (128 + signal number).  These traps unset EXIT so
	# we don't get double reportError calls.
	trap 'reportError 129 SIGHUP;  trap - EXIT; exit 129' SIGHUP
	trap 'reportError 130 SIGINT;  trap - EXIT; exit 130' SIGINT
	trap 'reportError 131 SIGQUIT; trap - EXIT; exit 131' SIGQUIT
	trap 'reportError 134 SIGABRT; trap - EXIT; exit 134' SIGABRT
	trap 'reportError 143 SIGTERM; trap - EXIT; exit 143' SIGTERM

	# Can not trap SIGSTOP or SIGKILL, they will be ignored.

	# jobs executed with qsub -notify recieve SIGUSR1 before SIGSTOP
	trap 'reportError 138 SIGUSR1; trap - EXIT; exit 138' SIGUSR1

	# jobs executed with qsub -notify recieve SIGUSR2 before SIGKILL
	trap 'reportError 140 SIGUSR2; trap - EXIT; exit 140' SIGUSR2

	# qsub jobs that exceed s_cpu, s_vmem receive SIGXCPU
	trap 'reportError 152 SIGXCPU; trap - EXIT; exit 152' SIGXCPU
}

function disableCommonTraps() { 

	# disable bash traps from enableCommonTraps
	trap '' DEBUG
	trap '' EXIT

	# disabled system traps from enableCommonTraps
	trap '' SIGHUP
	trap '' SIGINT
	trap '' SIGQUIT
	trap '' SIGABRT
	trap '' SIGTERM
	trap '' SIGUSR1
	trap '' SIGUSR2
	trap '' SIGXCPU
}
