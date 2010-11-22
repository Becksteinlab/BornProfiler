#$ -N %(jobname)s
#$ -S /bin/bash
#$ -l mem_free=500M,mem_total=500M
#$ -cwd
#$ -j y
#$ -r y

echo "APBS Born profile job running on $HOSTNAME"
ulimit -c 64
#module load abps/32 

# use apbs on PATH or set environment variable APBS
: ${APBS:=apbs}
: ${DRAW_MEMBRANE2A:=draw_membrane2a}
export APBS DRAW_MEMBRANE2A

# If we need to run draw_membrane2a as part of this job
RUN_DRAWMEMBRANE=%(drawmembrane_script)r
# BornProfiler decides if dx files are to be gunzipped: True|1=yes, False|0|*=no
UNPACK_DXGZ=%(unpack_dxgz)s
APBS_VERSION=`${APBS} --version 2>&1 1>/dev/null | awk '/^APBS/ {print $2}'`
if [ -z "${APBS_VERSION}" ]; then
    echo "ERROR: could not find apbs (${APBS}) on PATH. Aborting now."
    exit 1
fi
echo "APBS = ${APBS}"
echo "Detected apbs version ${APBS_VERSION}"
echo "Script was written for version %(apbs_version)s"

if [ -e "${RUN_DRAWMEMBRANE}" ]; then
    echo "Running ${RUN_DRAWMEMBRANE}..."
    nice bash ${RUN_DRAWMEMBRANE}
    # the script checks if it needs to do anything so in principle there
    # is no harm in just running it if we have it
fi

case "${UNPACK_DXGZ}" in
    1|true|True)
	echo "unpacking diel/kappa/charge dx files for buggy APBS 1.3..."
	nice gunzip -v {diel,kappa,charge}*m.dx.gz;;
esac

echo "ensuring single threaded calculation OMP_NUM_THREADS=1"
export OMP_NUM_THREADS=1
nice ${APBS} %(infile)s > %(outfile)s
rc=$?

case "${UNPACK_DXGZ}" in
    1|true|True)
	echo "compressing diel/kappa/charge dx files again to save 98%% of space..."
	nice gzip -v {diel,kappa,charge}*.dx;;
esac

exit $rc
