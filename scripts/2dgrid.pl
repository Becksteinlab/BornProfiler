#!/usr/bin/perl -w
# $Id: 2dgrid.pl,v 1.9 2005/10/07 17:41:28 oliver Exp $
# example: how to generate a grid of data points around a given
# point. Use odd dimensions.


use strict;
use FindBin;
use lib $FindBin::Bin; # Messages.pm is in script directory
                       # add other lib dirs in environment var PERL5LIB
use Messages qw($MSG_STDERR $DEBUG msg Progname); 

use Getopt::Long;

my $Version_ID = q($Id: 2dgrid.pl,v 1.9 2005/10/07 17:41:28 oliver Exp $ );
my ($opt_verbose, $opt_version, $opt_help);
my $lengthunit = "Angstrom";
my $lengthfac = 1.0;
my $OUTPUT = "/dev/stdout";

my (@CENTRE,@NDIM,@DELTA);

# msg() to stderr (important, otherwise we pollute the output)
$MSG_STDERR = 1;

sub print_usage {
    my $Prog = Progname;
    print <<EOF;
usage: $Prog [OPTIONS]

Create a (continuous) list of points on a grid, centered on a given point.
The list can be used as input for placeion.sh. All lengths in the output
are in Angstroms (so that it can be used as a pdb).

The OPTIONS can be abbreviated; defaults are given in [].	 
--help              help
--verbose           (alias for DEBUGLEVEL=3)
--version           print version id
--debug=DEBUGLEVEL  0 (almost) quiet, <0 silent,
                    1..3 verbose, 
                    >5 really debugging stuff (20 max) [$DEBUG]
--output=FILE       [$OUTPUT]      
--centre=X --centre=Y --centre=Z
                    centre of the grid; repeat three times 
--ndim=NX [--ndim=NY ]
                    dimension of the grid NX x NY (NX columns, NY rows)
--delta=DX [--delta=DY]
                    spacing in x and y between grid points (if only one
		    --delta is given, the x and y spacing is the same)
--lengthunit=STRING Angstrom|nm  unit of the input [$lengthunit]
EOF
    exit;
}

#----------------------------------------------------------------------
# MAIN
#----------------------------------------------------------------------


# options
&GetOptions( "debug=i"    => \$DEBUG,
             "output=s"   => \$OUTPUT,
             "centre=f@"  => \@CENTRE,
	     "ndim=i@"    => \@NDIM,
	     "delta=f@"   => \@DELTA,
             "lengthunit=s" => \$lengthunit, 
	     "verbose" => \$opt_verbose,
	     "version" => \$opt_version,
             "help" => \$opt_help) or &print_usage;
if ($opt_help) {&print_usage};
if ($opt_version) { msg (0, "%s\n", $Version_ID); exit 0 };
if ($opt_verbose) { $DEBUG = 3 };

# output (Angstrom) = lengthfac * input
$_ = $lengthunit;
LENGTHUNIT: {
    /nm|NM/        && do {$lengthfac = 10.0; last LENGTHUNIT;};
    /ang.*|Ang.*/  && do {$lengthfac = 1.0;  last LENGTHUNIT;};   
    die "Unknown lengthunit '$_', stopped";
}

# grid centre
if (@CENTRE != 3) {
    die "Need (x,y,z) of the centre of the grid. Repeat --centre=XX three times";
}
# convert lengths
@CENTRE = map { $lengthfac*$_ } @CENTRE;
msg(5, "Centre of the grid (in Angstrom): %f %f %f\n",@CENTRE);

# grid dimensions
if (@NDIM == 1) {
    # same ndim in x and y
    push @NDIM,$NDIM[0];
    msg(3,"Using same grid dimension in x and y:\n");
} elsif (@NDIM != 2) {
    die "Need NXxNY number of grid points, repeat --ndim twice.";
}
push @NDIM,1;
msg(5,"Grid dimensions %ix%ix%i\n",@NDIM);

# spacing
if (@DELTA == 1) {
    # same delta in x and y
    push @DELTA,$DELTA[0];
    msg(3,"Using same delta in x and y:\n");
} elsif (@DELTA != 2) {
    die "One or two values for the grid spacing --delta are required.\n";
}
@DELTA = map { $lengthfac*$_ } @DELTA;
push @DELTA,0;
msg(5,"Grid spacing %.2f %.2f %.2f (Angstrom)\n",@DELTA); 
    

my (@i0, @pt);
my ($ix, $iy);

$i0[0] = $NDIM[0]/2.0 + 0.5;
$i0[1] = $NDIM[1]/2.0 + 0.5;

open OUT,"> $OUTPUT" or die "Failed to open $OUTPUT for writing, ";

# header
print OUT "# 2D $NDIM[0] $NDIM[1]\n";

for ($ix=1; $ix<=$NDIM[0]; $ix++) {
    for ($iy=1; $iy<=$NDIM[1]; $iy++) {
	$pt[0] = ($ix-$i0[0]) * $DELTA[0] + $CENTRE[0];
	$pt[1] = ($iy-$i0[1]) * $DELTA[1] + $CENTRE[1];
	$pt[2] = $CENTRE[2];
	printf OUT "%8.4f %8.4f %8.4f\n", $pt[0], $pt[1], $pt[2];
    }
}

msg(5,"Wrote list of points ${OUTPUT}.\n");
close(OUT);
exit 0;
