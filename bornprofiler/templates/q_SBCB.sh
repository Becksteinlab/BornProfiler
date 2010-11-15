#$ -N %(jobname)s
#$ -S /bin/bash
#$ -l mem_free=500M,mem_total=500M
#$ -cwd
#$ -j y
#$ -r y

echo "APBS Born profile job running on $HOSTNAME"
ulimit -c 64
#module load abps/32 

# BornProfiler decides if dx files are to be gunzipped: True|1=yes, False|0|*=no
UNPACK_DXGZ=%(unpack_dxgz)s
APBS_VERSION=`apbs --version 2>&1 1>/dev/null| awk '/^APBS/ {print $2}'`
if [ -z "${APBS_VERSION}" ]; then
    echo "ERROR: could not find apbs on PATH. Aborting now."
    exit 1
fi
echo "Detected apbs version ${APBS_VERSION}"
echo "Script was written for version %(apbs_version)s"

case "${UNPACK_DXGZ}" in
    1|true|True)
	echo "unpacking dx files for buggy APBS 1.3..."
	gunzip -v *m.dx.gz;;
esac

echo "ensuring single threaded calculation OMP_NUM_THREADS=1"
export OMP_NUM_THREADS=1
apbs %(infile)s > %(outfile)s
rc=$?

case "${UNPACK_DXGZ}" in
    1|true|True)
	echo "compressing dx files again to save 98%% of space..."
	gzip -v *.dx;;
esac

exit $rc
