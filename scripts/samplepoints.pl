#!/usr/bin/perl -w
# $Id: samplepoints.pl,v 1.9 2005/09/13 10:35:53 oliver Exp $

use strict;
use FindBin;
use lib $FindBin::Bin; # Messages.pm is in script directory
                       # add other lib dirs in environment var PERL5LIB
use Messages qw($MSG_STDERR $DEBUG msg Progname); 

use XML::Simple qw(:strict) ;
use Data::Dumper;
use Getopt::Long;

my $Version_ID = q($Id: samplepoints.pl,v 1.9 2005/09/13 10:35:53 oliver Exp $ );
my ($opt_verbose, $opt_version, $opt_help);
my $lengthunit = "nm";
my $lengthfac = 1;

my $OUTPUT = "/dev/stdout";
my $XMLDB = "chl.xml";
my @COM = ();

my $xmlfile;   # name of the xml file
my $xml;       # ref to hash representing the data file

# msg() to stderr (important, otherwise we pollute the output)
$MSG_STDERR = 1;

sub print_usage {
    my $Prog = Progname;
    print <<EOF;
usage: $Prog --com=cX --com=cY [OPTIONS] [LIST ...]

Create a list of points along the axis of a model pore, usable for the
Born profiler scripts (output unit is always Angstrom). The list of 
positions is turned into a list of points with approximately distance
Dz (but Dz is adjusted to fullfill the boundary conditions zmin and zmax).
Typically, zmax(i) and zmin(i+i) are identical. The list can either come 
from standard input or as a file or sequence of files. Each line contains 
the limiting values in space and the spacing:

 zmin(1) zmax(1) Dz(1)  
 zmin(2) zmax(2) Dz(2) 
 ...

The centre of mass is used to determine the x and y coordinates for
the sampling path. It is normally given with the --com= options on the
commandline but when using the protein_modules scripts it can be read 
from an xml data file.

The OPTIONS can be abbreviated; defaults are given in [].	 
--help              help
--verbose           (alias for DEBUGLEVEL=3)
--version           print version id
--debug=DEBUGLEVEL  0 (almost) quiet, <0 silent,
                    1..3 verbose, 
                    >5 really debugging stuff (20 max) [$DEBUG]
--output=FILE       [$OUTPUT]      
--xml=FILE          required unless --com is used [$XMLDB]
--com=X --com=Y     center of mass of protein (in lengthunit)
--lengthunit=STRING Angstrom|nm  unit of the input [$lengthunit]
EOF
    exit;
}



######################################
#
# MAIN
#
######################################


# options
&GetOptions( "debug=i" => \$DEBUG,
             "output=s" => \$OUTPUT,
             "xml=s"    => \$XMLDB,
             "com=f@"   => \@COM,
             "lengthunit=s" => \$lengthunit, 
	     "verbose" => \$opt_verbose,
	     "version" => \$opt_version,
             "help" => \$opt_help) or &print_usage;
if ($opt_help) {&print_usage};
if ($opt_version) { msg (0, "%s\n", $Version_ID); exit 0 };
if ($opt_verbose) { $DEBUG = 3 };

# output (Angstrom) = lengthunit * input
$_ = $lengthunit;
LENGTHUNIT: {
    /nm|NM/        && do {$lengthfac = 10.0; last LENGTHUNIT;};
    /ang.*|Ang.*/  && do {$lengthfac = 1.0;  last LENGTHUNIT;};   
    die "Unknown lengthunit '$_', stopped";
}

if (@COM) {
    if (@COM < 2) {
        die "Need at least x and y coordinates for centre of mass (repeat --com x3)\n";
    }
} else {
    if (! -e $XMLDB) {
        die "XML database file for the channel is required.\n";
    };
}

# get description of z slices (no sanity checks)
# zmin zmax Dz N
my (@zlist, $min, $max, $delta);
while (<>) {
    chomp;
    next if (/^\s*[\#;]/) || (/^$/);
    ($min, $max, $delta) = split;
    push @zlist, { min => $min, max => $max, delta => $delta };
}

msg(9,Dumper(\@zlist));

# now process list
# - add N
# - change Dz
my ($slice, $n, $dz);
my $eps = 1e-6;
for $slice (@zlist) {
    $n = ($slice->{max} - $slice->{min})/$slice->{delta};
    $slice->{N} = int($n+$eps);
    $dz = ($slice->{max} - $slice->{min})/$slice->{N};    
    msg(3,"Adjusted delta for N=%i samples: %g --> %g\n", 
        $slice->{N}, $slice->{delta}, $dz)  unless ($slice->{N} == $n);
    $slice->{delta} = $dz;
}


my @com;
if (@COM) {
    @com = map { $lengthfac*$_ } @COM;
} else {
# For all what follows we want the db in  perl accessible hash
# ($xml is a ref to this hash)
    $xml = XMLin($XMLDB, KeyAttr => ['key'], 
                 ForceArray => 1, 
                 KeepRoot => 0);
    @com = @{ $xml->{ProteinOnly}->{protein}->{center} };
    # XMLDB is in nm (hence factor *10 below)
    @com = map { 10*$_ } @com;
}

msg(5, "centre of mass (in nm): %f %f %f\n",@com);

open OUT, "> $OUTPUT" or die "Failed to open output file $OUTPUT, ";

my ($i, $z, $itot);
for $slice (@zlist) {
    for($i=0;$i<$slice->{N};$i++) {
        $z = $slice->{min} + $i*$slice->{delta};
	$itot++;
        printf OUT "%10.3f %10.3f %10.3f\n", 
          $com[0], $com[1], $lengthfac*$z;
    }
    msg(3,"Generated N=%i points in slice min=%g max=%g with delta=%g\n",
        $slice->{N}, $slice->{min}, $slice->{max}, $slice->{delta});
}
msg(0,"Generated $itot points in total.\n");

close OUT;
exit(0);

