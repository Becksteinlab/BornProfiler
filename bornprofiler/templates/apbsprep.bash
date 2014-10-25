#!/bin/bash
# Generate dielectric/kappa/charge files without membrane.
# :Authors: Oliver Beckstein <oliver.beckstein@bioch.ox.ac.uk> Lennard van der Feltz <lvanderf@asu.edu>
# :Year: 2010-2011, 2014
# :Licence: This file is placed in the Public Domain.

# can be set in the environment
: ${APBS:=apbs}

# set by the BornProfiler scripts
DUMMY_IN=%(born_dummy_in)r
DUMMY_OUT=%(born_dummy_out)r
INFICES="%(infices)s"

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

echo "APBS = ${APBS}"

checkfile "${DUMMY_IN}"
run_apbs "${DUMMY_IN}" "${DUMMY_OUT}"


echo "Generated all dx files required for running the window."
exit 0
