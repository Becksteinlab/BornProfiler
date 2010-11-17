#!/usr/bin/perl -w
# $Id: Messages.pm,v 1.2 2005/12/03 23:30:03 oliver Exp $
# Copyright (c) 2003 Oliver Beckstein <oliver@biop.ox.ac.uk>
# Distributed under the GNU Public License,
# see http://www.gnu.org/copyleft/gpl.html
#
#-------------------------------------------
# Insert into code the use lib path which contains Messages.pm:
# use lib "/sansom/kir/oliver/local/lib/pm/lib/perl";
# use Messages; 
#-------------------------------------------

package      Messages;
require      Exporter;
@ISA       = qw(Exporter);
@EXPORT    = qw(Progname warning errmsg logger log_only msg open_log $DEBUG);
@EXPORT_OK = qw(LOG $MSG_STDERR date xe $SHELL_QUIET);

use POSIX qw(strftime);
use Cwd;
use File::Basename qw(basename);

$DEBUG = 0;
$MSG_STDERR = 0;   # normally print to stdout
$SHELL_QUIET = ">/dev/null 2>&1";
my $use_log = 0;

# somehow main is only happy if the package returns true 
# -- and $DEBUG=0 returns false... weird Perl.
return 1;

END {
    close (LOG);
};

sub Progname () {
    return basename $0;
};

sub date () {
    return POSIX::strftime("%c", localtime(time));
};

sub open_log () {
    my ($error, $LogFn);
    ($LogFn = Progname) =~ s(\.pl$)();  
    $LogFn .= ".log";
    if ($error = open (LOG, "> " . $LogFn)) {
	print LOG Progname . " [Logfile] --- ",
	"opened on ", date(), "\n";
	print LOG "in ", cwd(), "\n";
	print LOG 
	    "-----------------------------------------------------------\n\n";
        $use_log = 1;
    };
    return $error;
};


sub warning {
    my ($format, @args) = @_;
    $format = "Warning: " . $format;
    logger ( 0, $format, @args);
    return;
};

sub errmsg {
    my ($format, @args) = @_;
    $format = "Error: " . $format;
    logger ( 0, $format, @args);
    return;
};


sub logger {
    my ($level, @args) = @_;
    do {
	if ($MSG_STDERR) {
	    printf STDERR "[$level] ";
	    printf STDERR @args;
	} else {
	    printf  "[$level] ";
	    printf  @args;
	};
	printf LOG @args if ($use_log);
    } unless $DEBUG < $level;
    return;
};

sub log_only {
    my ($level, @args) = @_;
    do {
	printf LOG "[$level] ";
	printf LOG @args;
    } unless ($DEBUG < $level or not $use_log);
    return;
};


sub msg {
    my ($level, @args) = @_;
    do {
	if ($MSG_STDERR) {
	    printf STDERR "[$level] ";
	    printf STDERR @args;
	} else {
	    printf  "[$level] ";
	    printf  @args;
	};
    } unless $DEBUG < $level;
    return;
};


# not really logging but...
sub xe ($) {
    my ($CMD) = @_;
    if (LOG) {
	logger (1, ">>> $CMD\n");
    } else {
	print STDERR ">>> $CMD\n";
    };
    system $CMD;
    return $!;
};
