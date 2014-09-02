#!/bin/bash
# Generate dielectric/kappa/charge files with membrane inserted.
# Requires a customized version 'draw_membrane2a' (can be found in the BornProfiler
# distribution in 'src/drawmembrane/draw_membrane2a.c')
# Written by bornprofiler.membrane instead of MPlaceion.run_drawmembrane()
# :Author: Oliver Beckstein <oliver.beckstein@bioch.ox.ac.uk>
# :Year: 2010-2011
# :Licence: This file is placed in the Public Domain.

# can be set in the environment
: ${APBS:=apbs}
: ${DRAW_MEMBRANE2A:=draw_membrane2a}

# set by the BornProfiler scripts
DUMMY_IN='mem_dummy.in'
DUMMY_OUT='born_dummy.out'
INFICES="_prot_L _prot_M _prot_S _cpx_L _cpx_M _cpx_S"

function die () {
    echo "ERROR: $1";
    exit ${2:-};
}

function checkfile () {
    test -e "$1" || die "File $1 does not exist." 2;
}

function run_apbs () {
    local infile="$1" outfile="$2"
    local targets
    targets=$(ls {dielx,diely,dielz,kappa,charge}*[LMS].dx* 2>/dev/null)
    if [ $(echo $targets | wc -w) -eq 30 ]; then
	echo "APBS-generated diel/kappa/charge dx already exist ... skipping run_apbs('born_dummy')"
	return
    else
	echo "Running ${APBS} ${infile} ... [run_apbs('born_dummy')]"
	${APBS} "${infile}" 2>&1 | tee "${outfile}" \
	    || die "${APBS} ${infile} failed."
    fi
}

function run_drawmembrane () {
    # all draw_membrane interpolations are done here!
    local infix=$1 compression="gz"
    local targets compress_option=""
    targets=$(ls {dielx,diely,dielz,kappa,charge}${infix}m.dx* 2>/dev/null)
    if [ $(echo $targets | wc -w) -eq 5 ]; then
	echo "Files for $infix already exist .. skipping"
	return
    fi
    if [ "$compression" == "gz" ]; then
	compress_option="-Z"
    fi
    ${DRAW_MEMBRANE2A} $compress_option -z 38.000000 -d 45.000000 \
	-p 10.000000 -s 80.000000 -m 2.000000  \
	-R 10.000000 -r 10.000000 -X 201.687906 -Y 64.783443 -c 80.0  \
	-a 0.000000 -i 20.000000 \
	-V 0.000000 -I 0.150000   $infix \
	|| die "${DRAW_MEMBRANE2A} $infix failed."
}    

echo "------------------------------------------------------------"
echo "Generating dx files with membrane"
echo "------------------------------------------------------------"
echo "APBS = ${APBS}"
echo "draw_membrane2a = ${DRAW_MEMBRANE2A}"

checkfile "${DUMMY_IN}"
run_apbs "${DUMMY_IN}" "${DUMMY_OUT}"

echo "Running ${DRAW_MEMBRANE2A} ... [run_drawmembrane()]"
for infix in ${INFICES}; do
    run_drawmembrane $infix
done

echo "Generated all membrane files required for running the window."
exit 0
