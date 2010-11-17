#!/usr/bin/perl -w
# $Id: xyz2pdb.pl,v 1.5 2005/10/07 17:41:49 oliver Exp $

use File::Basename;
use Getopt::Long;

my $prog = basename($0);
my $Version_ID = q($Id: xyz2pdb.pl,v 1.5 2005/10/07 17:41:49 oliver Exp $ );
my ($opt_help);
my $PDB = "/dev/stdout";

my ($x,$y,$z);
my $atomnr = 0;

my $usage = "usage: $prog [OPTS] --output=samplepoints.pdb  samplepoints.xyz 
This primitive script turns a list of  x y z coordinates (one per line)
into a fake pdb file which can be displayed in eg vmd easily.
--help           help
--output=FILE    write to file instead of stdout [$PDB]
";

sub print_usage {
    print "$usage";
    exit;
}

#------------------------------------------------------------
# main
#------------------------------------------------------------

# options
&GetOptions( "output=s" => \$PDB,
             "help" => \$opt_help) or &print_usage;
if ($opt_help) {&print_usage};

open(PDB,"> ${PDB}") or die "Cannot open ${PDB} for writing, stopped";

my $pdbformat = "ATOM%7i  C   XXX X   1    %8.3f%8.3f%8.3f\n";
while (<>) {
    if (/^#/) {
	print PDB $_;
	next;
    }
    $atomnr++;
    ($x,$y,$z) = split;
    printf PDB $pdbformat,$atomnr,$x,$y,$z;
}

close(PDB);
exit;

#ATOM      1  N   ALA X   1      51.245  51.650  85.060  1.00  0.00
#ATOM      1  C   XXX X   1      63.750  63.750   0.000
