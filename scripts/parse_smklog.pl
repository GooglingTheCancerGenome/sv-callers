#!/usr/bin/perl
#
# Parse Snakemake log for DELLY jobs.
#

use strict;
use warnings;

$\ = "\n";
$, = ',';
my $read = 0;
my $jobid = 0;
my $sv = '';
my $path = '';

while (<>) {
   chomp;
   if (/rule delly_[sp]/) {
      $read = 1;
   }
   if ($read) {
      $sv = $1 if /sv_type=(\S+)/;
      $jobid = $1 if /jobid '(\d+)'/;
      $path = $1 if /path=([^,]+)/;
      if ($sv && $jobid && $path) {
         print $jobid, $sv, $path;
	 $sv = '';
	 $path = '';
	 $jobid = 0;
	 $read = 0;
      }
   }
}

