#!/bin/bash
# Run one APBS Born window in the local shell.
# This script can be use as the basis for your own submission scripts.

# These binaries need to be on the PATH or supplied with their full path
# (The RUN_DRAWMEMBRANE script can pick them up from the environment.)
# use apbs on PATH or set environment variable APBS
: ${APBS:=apbs}
: ${DRAW_MEMBRANE2A:=draw_membrane2a}
export APBS DRAW_MEMBRANE2A

# If we need to run draw_membrane2a as part of this job
RUN_DRAWMEMBRANE=%(drawmembrane_script)r

# BornProfiler decides if dx files are to be gunzipped: True|1=yes, False|0|*=no
UNPACK_DXGZ=%(unpack_dxgz)s
APBS_VERSION=`${APBS} --version 2>&1 1>/dev/null| awk '/^APBS/ {print $2}'`
if [ -z "${APBS_VERSION}" ]; then
    echo "ERROR: could not find apbs (${APBS}) on PATH. Aborting now."
    exit 1
fi

echo "** %(jobname)s"
echo "-- APBS =${APBS}"
echo "-- Detected apbs version ${APBS_VERSION}"
echo "-- Script was written for version %(apbs_version)s"

if [ -e "${RUN_DRAWMEMBRANE}" ]; then
    echo "-- Running ${RUN_DRAWMEMBRANE}..."
    bash ${RUN_DRAWMEMBRANE}
    # the script checks if it needs to do anything so in principle there
    # is no harm in just running it if we have it
fi

case "${UNPACK_DXGZ}" in
    1|true|True)
	echo "-- unpacking diel/kappa/charge dx files for buggy APBS 1.3..."
	nice gunzip -v {diel,kappa,charge}*m.dx.gz;;
esac

echo "-- APBS Born profile job running on $HOSTNAME"
echo "++ apbs %(infile)s > %(outfile)s"
nice ${APBS} %(infile)s > %(outfile)s
rc=$?
echo "-- job complete: results in %(outfile)s"

case "${UNPACK_DXGZ}" in
    1|true|True)
	echo "-- compressing diel/kappa/charge dx files again to save 98%% of space..."
	nice gzip -v {diel,kappa,charge}*.dx;;
esac

exit $rc
