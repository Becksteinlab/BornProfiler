#!/bin/sh
# BornProfiler decides if dx files are to be gunzipped: True|1=yes, False|0|*=no

UNPACK_DXGZ=%(unpack_dxgz)s
APBS_VERSION=`apbs --version 2>&1 1>/dev/null| awk '/^APBS/ {print $2}'`

echo "** %(jobname)s"
echo "-- Detected apbs version ${APBS_VERSION}"
echo "-- Script was written for version %(apbs_version)s"

case "${UNPACK_DXGZ}" in
    1|true|True)
	echo "-- unpacking dx files for buggy APBS 1.3..."
	gunzip -v *m.dx.gz;;
esac

echo "-- APBS Born profile job running on $HOSTNAME"
echo "++ apbs %(infile)s > %(outfile)s"
apbs %(infile)s > %(outfile)s
rc=$?
echo "-- job complete: results in %(outfile)s"

case "${UNPACK_DXGZ}" in
    1|true|True)
	echo "-- compressing dx files again to save 98%% of space..."
	gzip -v *.dx;;
esac

exit $rc
