#!/bin/sh
echo "** %(jobname)s"
echo "-- APBS Born profile job running on $HOSTNAME"
echo "++ apbs %(infile)s > %(outfile)s"
apbs %(infile)s > %(outfile)s
echo "-- job complete: results in %(outfile)s"
