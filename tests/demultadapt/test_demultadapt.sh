#!/usr/bin/env bash

# Aim: launch a functional test for the demultadapt program
# Copyright (C) 2014 Institut National de la Recherche Agronomique (INRA)
# License: GPL-3+
# Authors: Gautier Sarah, Vincent Maillol

progVersion="1.0"

# Display the help on stdout.
# The format complies with help2man (http://www.gnu.org/s/help2man)
function help () {
    msg="\`${0##*/}' launches a functional test for the demultadapt program.\n"
    msg+="\n"
    msg+="Usage: ${0##*/} [OPTIONS] ...\n"
    msg+="\n"
    msg+="Options:\n"
    msg+=" -h, --help\tdisplay the help and exit\n"
    msg+=" -V, --version\toutput version information and exit\n"
    msg+=" -v, --verbose\tverbosity level (0/default=1/2/3)\n"
    msg+=" -d, --p2d\tabsolute path to the arcad-hts directory\n"
    msg+=" -n, --noclean\tkeep temporary directory with all files\n"
    msg+="\n"
    msg+="Report bugs to <gautier.sarah@supagro.inra.fr>."
    echo -e "$msg"
}

# Display version and license information on stdout.
function version () {
    msg="${0##*/} ${progVersion}\n"
    msg+="\n"
    msg+="Copyright (C) 2014 Institut National de la Recherche Agronomique (INRA).\n"
    msg+="License GPL-3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n"
    msg+="This is free software; see the source for copying conditions. There is NO\n"
    msg+="warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"
    msg+="\n"
    msg+="Written by Gautier Sarah & Vincent Maillol."
    echo -e "$msg"
}

# http://www.linuxjournal.com/content/use-date-command-measure-elapsed-time
function timer () {
    if [[ $# -eq 0 ]]; then
        echo $(date '+%s')
    else
        local startRawTime=$1
        endRawTime=$(date '+%s')
        if [[ -z "$startRawTime" ]]; then startRawTime=$endRawTime; fi
        elapsed=$((endRawTime - startRawTime)) # in sec
        nbDays=$((elapsed / 86400))
        nbHours=$(((elapsed / 3600) % 24))
        nbMins=$(((elapsed / 60) % 60))
        nbSecs=$((elapsed % 60))
        printf "%01dd %01dh %01dm %01ds" $nbDays $nbHours $nbMins $nbSecs
    fi
}

# Parse the command-line arguments.
# http://stackoverflow.com/a/4300224/597069
function parseCmdLine () {
    getopt -T > /dev/null # portability check (say, Linux or Mac OS?)
    if [ $? -eq 4 ]; then # GNU enhanced getopt is available
        TEMP=`getopt -o hVv:d:n -l help,version,verbose:,p2d:,noclean \
-n "$0" -- "$@"`
    else # original getopt is available (no long options, whitespace, sorting)
        TEMP=`getopt hVv:i:n "$@"`
    fi
    if [ $? -ne 0 ]; then
        echo "ERROR: "$(which getopt)" failed"
        getopt -T > /dev/null
        if [ $? -ne 4 ]; then
            echo "did you use long options? they are not handled \
on your system, use -h for help"
        fi
        exit 2
    fi
    eval set -- "$TEMP"
    while [ $# -gt 0 ]; do
        case "$1" in
            -h | --help) help; exit 0; shift;;
            -V | --version) version; exit 0; shift;;
            -v | --verbose) verbose=$2; shift 2;;
            -d | --p2d) pathToArcadHtsDir=$2; shift 2;;
            -n | --noclean) clean=false; shift;;
            --) shift; break;;
            *) echo "ERROR: options parsing failed, use -h for help"; exit 1;;
        esac
    done
    if [ -z "${pathToArcadHtsDir}" ]; then
        echo -e "ERROR: missing compulsory option --p2d\n"
        help
        exit 1
    fi
    if [ ! -d "${pathToArcadHtsDir}" ]; then
        echo -e "ERROR: can't find directory ${pathToArcadHtsDir}\n"
        help
        exit 1
    fi
}

function cmp_or_quit {
    cmp $1 $2
    if [[ "$?" == 0 ]]; then
        rm $2
    else 
        exit 1
    fi
}

function run () {
    cwd=$(pwd)
    cd "${pathToArcadHtsDir}/tests/demultadapt"
    
    # step 1 ------------------------------------------------------------------
    if [ $verbose -gt "0" ]; then
        echo -e "check presence of input data..."
    fi
    if [ ! -f "indi_A_1.fastq" ]; then
        if [ $verbose -gt "0" ]; then
            echo -e "missing indi_A_1.fastq"
            exit 1
        fi
    fi
    
    if [ ! -f "indi_A_2.fastq" ]; then
        if [ $verbose -gt "0" ]; then
            echo -e "missing indi_A_2.fastq"
            exit 1
        fi
    fi
    
    if [ ! -f "adaptator.txt" ]; then
        if [ $verbose -gt "0" ]; then
            echo -e "missing adaptator.txt"
            exit 1
        fi
    fi
    
    # step 2 ------------------------------------------------------------------
    uniqId=$$ # process ID
    testDir=tmp_test_${uniqId}
    rm -rf ${testDir}
    mkdir ${testDir}
    cd ${testDir}
    if [ $verbose -gt "0" ]; then echo "temp dir: "$(pwd); fi
    
    # step 3 ------------------------------------------------------------------
    if [ $verbose -gt "0" ]; then
        echo -e "launch demultadapt paired..."
    fi
    
    export PYTHONPATH=${pathToArcadHtsDir}/lib:$PYTHONPATH
    
    cmd="python"
    cmd+=" ${pathToArcadHtsDir}/sp5_gbs/demultadapt.py"
    cmd+=" -l 1.0"
    cmd+=" -f ${pathToArcadHtsDir}/tests/demultadapt/indi_A_1.fastq"
    cmd+=" -F ${pathToArcadHtsDir}/tests/demultadapt/indi_A_2.fastq"
    cmd+=" -p returned"
    cmd+=" -v ${pathToArcadHtsDir}/tests/demultadapt/adaptator.txt"
    if [ $verbose -le "1" ]; then
      cmd+=" > /dev/null"
    fi
    if [ $verbose -ge "1" ]; then
      echo $cmd
    fi
    eval $cmd
    
    cmp_or_quit ${pathToArcadHtsDir}/tests/demultadapt/expected-indiv_1_1.fastq \
        returned-indiv_1_1.fastq 
    cmp_or_quit ${pathToArcadHtsDir}/tests/demultadapt/expected-indiv_1_2.fastq \
        returned-indiv_1_2.fastq
    cmp_or_quit ${pathToArcadHtsDir}/tests/demultadapt/expected-indiv_2_1.fastq \
        returned-indiv_2_1.fastq
    cmp_or_quit ${pathToArcadHtsDir}/tests/demultadapt/expected-indiv_2_2.fastq \
        returned-indiv_2_2.fastq 
    cmp_or_quit ${pathToArcadHtsDir}/tests/demultadapt/expected-rebut_1.fastq \
        returned-rebut_1.fastq
    cmp_or_quit ${pathToArcadHtsDir}/tests/demultadapt/expected-rebut_2.fastq \
        returned-rebut_2.fastq 
    echo "all tests passed successfully!"
    
    if [ $verbose -gt "0" ]; then
        echo -e "launch demultadapt single..."
    fi
    
    cmd="python"
    cmd+=" ${pathToArcadHtsDir}/sp5_gbs/demultadapt.py"
    cmd+=" -l 1.0"
    cmd+=" -f ${pathToArcadHtsDir}/tests/demultadapt/indi_A_1.fastq"
    cmd+=" -p returned"
    cmd+=" -v ${pathToArcadHtsDir}/tests/demultadapt/adaptator.txt"
    if [ $verbose -le "1" ]; then
      cmd+=" > /dev/null"
    fi
    eval $cmd
    
    cmp_or_quit ${pathToArcadHtsDir}/tests/demultadapt/expected-indiv_1.fastq \
        returned-indiv_1.fastq
    cmp_or_quit ${pathToArcadHtsDir}/tests/demultadapt/expected-indiv_2.fastq \
        returned-indiv_2.fastq
    cmp_or_quit ${pathToArcadHtsDir}/tests/demultadapt/expected-rebut.fastq \
        returned-rebut.fastq
    echo "all tests passed successfully!"
    
    # step 4 ------------------------------------------------------------------
    cd ${cwd}
    if $clean; then rm -rf "${pathToArcadHtsDir}/tests/${testDir}"; fi
}

verbose=1
pathToArcadHtsDir=""
clean=true
parseCmdLine "$@"

if [ $verbose -gt "0" ]; then
    startTime=$(timer)
    msg="START ${0##*/} $(date +"%Y-%m-%d") $(date +"%H:%M:%S")"
    msg+="\ncmd-line: $0 "$@ # comment if an option takes a glob as argument
    msg+="\ncwd: $(pwd)"
    echo -e $msg
fi

run pathToArcadHtsDir clean

if [ $verbose -gt "0" ]; then
    msg="END ${0##*/} $(date +"%Y-%m-%d") $(date +"%H:%M:%S")"
    msg+=" ($(timer startTime))"
    echo $msg
fi
