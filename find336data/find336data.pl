#!/usr/bin/env perl

use warnings;
use strict;
$|++;

use Date::Calc;
use POSIX qw(mktime);

# Globals
my $sot_root = "/hao/hinode1/hinode/sot";
my $dest_root = "/home/egeland/Data/hop336";
my $dt_warn = 3600; # [s]; time distance between pointing and observation for which to print a warning
my $dt_skip = 4 * 3600; # [s]; ... to throw an error and skip

sub date2obsdir {
    my @date = @_; # (Y, M, D, h, m, s)
    my $obsdir = sprintf "%04d%02d%02d_%02d%02d%02d", @date;
    return $obsdir;
}

sub obsdir2date {
    my $obsdir = shift;
    my @date = ($obsdir =~ m|(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2})|);
    return @date;
}

sub obsdir_time_diff {
    my ($d1, $d2) = @_;
    my ($Y1, $M1, $D1, $h1, $m1, $s1) = obsdir2date($d1);
    my ($Y2, $M2, $D2, $h2, $m2, $s2) = obsdir2date($d2);
    $s1 = mktime($s1, $m1, $h1, $D1, $M1 - 1, $Y1 - 2000);
    $s2 = mktime($s2, $m2, $h2, $D2, $M2 - 1, $Y2 - 2000);
    return $s2 - $s1;
}

# Prepare data destination
mkdir $dest_root unless -d $dest_root;
mkdir "$dest_root/north" unless -d "$dest_root/north";
mkdir "$dest_root/south" unless -d "$dest_root/south";
mkdir "$dest_root/equator" unless -d "$dest_root/equator";

# Parse the pointing file and put data into a hash keyed on date+time
# Note that this choice eliminates potential duplicates, and prefers the later entry.
#
# File format example:
# 0         1         2         3         4         5         6         7         8         9        10        11
# 012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456
# ORe-point Start     2017/03/07 10:52:00    1      00.0000  00.0000  /*  -172.1/  622.0 : # OP start + 10min + HOP 336 */
# ORe-point Start     2017/03/07 13:52:00    2      00.0000  00.0000  /*   128.0/  622.0 : HOP 336-2 */
# ORe-point Start     2017/03/10 18:40:00    1      00.0000  00.0000  /*  -172.0/ -718.0 : #HOP 336 */

my %obs;
open POINTINGS, '<', 'pointings.txt';
while (<POINTINGS>) {
    # Parse pointing date, e.g. 2017/03/07 10:52:00
    my @datetime = ($_ =~ m|(\d{4}/\d{2}/\d{2}) (\d{2}:\d{2}:\d{2})|);
    my $date = $datetime[0];
    my $time = $datetime[1];
    my $datetime_str = "$date $time";
    my ($Y, $M, $D) = ($date =~ m|(\d{4})/(\d{2})/(\d{2})|);
    my ($h, $m, $s) = ($time =~ m|(\d{2}):(\d{2}):(\d{2})|);

    # Parse pointing coordinates, e.g. -172.1/  622.0
    my ($x, $y) = ($_ =~ m|(\-*\d+\.\d+)/\s+(\-*\d+\.\d+)|);

    # Parse comment
    my $comment = substr $_, 89, -4; # chop off manditory ' */\n'
    $comment =~ s/^\#*\s*(.+)$/$1/;  # chop off optional leading # and spaces

    # Build pointing file data structure
    $obs{$datetime_str} = { DATETIME => $datetime_str,
			    DATE => $date,
			    TIME => $time,
			    Y => $Y, M => $M, D => $D,
			    h => $h, m => $m, s => $s,
			    XCOORD => $x,
			    YCOORD => $y,
			    COMMENT => $comment,
			    LEVEL1 => "$sot_root/level1/$Y/$M/$D/SP3D",
			    LEVEL2 => "$sot_root/level2/$Y/$M/$D/SP3D",
			    MINDIR => date2obsdir($Y, $M, $D, $h, $m, $s)
    };
}
close POINTINGS;

sub find_nearest_obsdir {
    my ($path, $mindir) = @_;
    opendir(my $dh, $path) || die "can't opendir $path: $!";
    my @dirs = sort grep { ! /^\./ && -d "$path/$_" } readdir($dh);
    closedir $dh;

    my $datadir = '';
    for my $d (@dirs) {
	if ($d gt $mindir) {
	    $datadir = $d;
	    last;
	}
    }
    return $datadir;
}

sub tomorrows_path {
    my $obsdir = shift;
    my @date = obsdir2date($obsdir);
    my @day_after = Date::Calc::Add_Delta_Days($date[0], $date[1], $date[2], 1);
    my $path = sprintf "$sot_root/level2/%02d/%02d/%02d/SP3D", @day_after;
    return $path;
}

# Observation sequence is expected to begin a few minutes after the pointing
# Iterate through each observation and check for existence of a nearby 
# following observation in the level2, level1, and level0 directories

my @obs = sort keys %obs;
foreach my $o (@obs) {
    # Get path data
    my $pointdata = $obs{$o};
    my $level2_path = $pointdata->{LEVEL2};
    my $mindir = $pointdata->{MINDIR};

    # Check today's path exists; if not, try tomorrow's path
    if (! -d $level2_path) {
	print "WARN $level2_path does not exist, trying next day\n";
	my $level2_path_orig = $level2_path;
	$level2_path = tomorrows_path($mindir);
	if (! -d $level2_path) {
	    print "ERROR next-day path $level2_path does not exist either, skipping $o\n";
	    next;
	}
    }

    # Get the nearest observation to the pointing time
    # If nothing is found, try tomorrow's path
    my $datadir = find_nearest_obsdir($level2_path, $mindir);
    if (! $datadir ) {
	print "WARN Could not find observation after $o, trying next day\n";
	my $level2_path_orig = $level2_path;
	$level2_path = tomorrows_path($mindir);
	if (! -d $level2_path) {
	    print "ERROR next-day path $level2_path does not exist, skipping $o\n";
	    next;
	}
	$datadir = find_nearest_obsdir($level2_path, $mindir);
	if (! $datadir) {
	    print "ERROR Could not find observation around $o, skipping\n";
	    next;
	}
    }
    
    # Calculate time difference between pointing and observation.
    # Apply thresholds to state of observation lookup.
    my $dt = obsdir_time_diff($mindir, $datadir);
    my $dt_str = '';
    if ($dt < 60) {
	$dt_str = "$dt s";
    } elsif ($dt >= 60 && $dt < 3600) {
	$dt_str = sprintf "%0.2f m", $dt / 60.0;
    } elsif ($dt >= 3600) {
	$dt_str = sprintf "%0.2f h", $dt / 3600.0;
    }
    my $state = "OK";
    my $do_copy = 1;
    if ($dt > $dt_skip) {
	$state = "ERROR (too late)";
	$do_copy = 0;
    } elsif ($dt > $dt_warn) {
	$state = "WARN (late)";
    }
    print "$state $o -> $level2_path/$datadir (dt=$dt_str)\n";
    my $src_dir = "$level2_path/$datadir";
    my $src_fits = "$src_dir/$datadir.fits";
    my $src_sav = "$src_dir/$datadir.sav";
    unless (-f $src_fits) {
	print "ERROR $src_fits not found, skipping $o\n";
	next;
    }
    unless (-f $src_sav) {
	print "ERROR $src_sav not found, skipping $o\n";
	next;
    }

    # If a satisfactory observation has been found, gather it
    if ($do_copy) {
	my $band;
	if ($pointdata->{YCOORD} >= 300) {
	    $band = 'north';
	} elsif ($pointdata->{YCOORD} <= -300) {
	    $band = 'south';
	} else {
	    $band = 'equator';
	}
	my $dest_dir = "$dest_root/$band";
	system "ln -s $src_fits $dest_dir";
	system "ln -s $src_sav $dest_dir";
    }

}

# TODO:
#  - Do the same for level1, level0?
#  - Cross-check with HCR results?
