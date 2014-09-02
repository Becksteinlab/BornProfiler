#$ -N w0044_Na_nAChR
#$ -S /bin/bash
#$ -l mem_free=1024M
#$ -cwd
#$ -j y
#$ -r y

function die () {
    echo "EE ERROR: $1"
    exit ${2:-1}
}

function count_files () {
    local pattern="$*" files
    files=$(ls $pattern 2>/dev/null)
    echo $files | wc -w
}

#------------------------------
# staging
#------------------------------
# requires -cwd flag!
ORIGDIR=$PWD
IS_STAGED=False
WORKDIR=.

function stage () {
# Using local disks to avoid nfs problems
    if [ -n "$TMPDIR" ]; then
	echo "-- Using TMPDIR=$TMPDIR"
        # uses TMPDIR (use BSD mktemp syntax!)
	WORKDIR=$(mktemp -d -t bornprofiler.XXXXXXXX) \
	    || die "Failed to create temp work dir." 2
	echo "-- staging into $WORKDIR..."
	cp -v *.pqr *.in *.bash $WORKDIR
	IS_STAGED=True
	cd $WORKDIR || die "Failed to 'cd $WORKDIR'" 2
    fi
}

function unstage () {
    case ${IS_STAGED} in
	True|1)
	    echo "-- unstaging from $WORKDIR..."
	    cp -v *.dx *.dx.gz *.out $ORIGDIR
	    echo "-- copied files from $HOSTNAME:$WORKDIR --> $ORIGDIR"
	    cd $ORIGDIR || die "Failed to 'cd $ORIGDIR' --- WTF, dude?" 2
	    if [ $(count_files *.out) -lt 2 ] || [ $(count_files *.dx *.dx.gz) -lt 66 ]; then
		die "Missing results from $HOSTNAME:$WORKDIR. Investigate manually!" 2
	    fi
	    rm -rf $WORKDIR
	    echo "-- removed $WORKDIR"
	    ;;
    esac
}

# copy files to TMPDIR
stage

# clean up if we get killed (e.g. qdel): KILL TERM HUP
trap unstage  9 11 1

#
#------------------------------


# These binaries need to be on the PATH or supplied with their full path
# (The RUN_DRAWMEMBRANE script can pick them up from the environment.)
# use apbs on PATH or set environment variable APBS
: ${APBS:=apbs}
: ${DRAW_MEMBRANE2A:=draw_membrane2a}
export APBS DRAW_MEMBRANE2A

# If we need to run draw_membrane2a as part of this job
RUN_DRAWMEMBRANE='run_drawmembrane.bash'

# BornProfiler decides if dx files are to be gunzipped: True|1=yes, False|0|*=no
UNPACK_DXGZ=False
APBS_VERSION=`${APBS} --version 2>&1 1>/dev/null| awk '/^APBS/ {print $2}'`
if [ -z "${APBS_VERSION}" ]; then
    die "ERROR: could not find apbs (${APBS}) on PATH. Aborting now."
fi

echo "** w0044_Na_nAChR"
echo "-- APBS =${APBS}"
echo "-- Detected apbs version ${APBS_VERSION}"
echo "-- Script was written for version 1.4"
echo "-- ensuring single threaded calculation of APBS: OMP_NUM_THREADS=1"
export OMP_NUM_THREADS=1

if [ -e "${RUN_DRAWMEMBRANE}" ]; then
    echo "Running ${RUN_DRAWMEMBRANE}..."
    nice bash ${RUN_DRAWMEMBRANE}
    # the script checks if it needs to do anything so in principle there
    # is no harm in just running it if we have it
fi

case "${UNPACK_DXGZ}" in
    1|true|True)
	echo "-- unpacking diel/kappa/charge dx files for buggy APBS 1.3..."
	nice gunzip -v {diel,kappa,charge}*m.dx.gz;;
esac

echo "-- APBS Born profile job running on $HOSTNAME"
echo "++ apbs job_0044.in 2>&1 | tee job_0044.out"
nice ${APBS} job_0044.in 2>&1 | tee job_0044.out
rc=$?
echo "-- job complete: results in job_0044.out"

case "${UNPACK_DXGZ}" in
    1|true|True)
	echo "-- compressing diel/kappa/charge dx files again to save 98% of space..."
	nice gzip -v {diel,kappa,charge}*.dx;;
esac

# copy files back
unstage

exit $rc
