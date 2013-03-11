#! /usr/bin/perl

#
# Usage: hpccoutf.pl -a < hpccoputf.txt To show all parameters
#        hpccoutf.pl -w < hpccoputf.txt To show web parameters
#

use strict;
#use diagnostics;
use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION  = 1;
our ( $opt_a, $opt_w, $opt_f, $value, $count, $key ) = 0;
getopts("awf:");

unless ( $opt_a || $opt_w && $opt_f ) {
    print "\n";
    print "Usage: $0 -a -f hpccoutf.txt For all parameters\n";
    print "       $0 -w -f hpccoutf.txt For web parameters\n";
    exit;
}

$/           = undef;
open HPPCOUTF, $opt_f or die $!;
my $file     = <HPPCOUTF>;
my %hpccoutf = $file =~ /^(\w+)=(\d.*)$/mg;
close HPPCOUTF;


my @walkorder = (
  'HPL_Tflops',                      'PTRANS_GBs',
  'MPIRandomAccess_GUPs',            'MPIFFT_Gflops',
  'StarSTREAM_Triad*CommWorldProcs', 'StarSTREAM_Triad',
  'StarDGEMM_Gflops',                'RandomlyOrderedRingBandwidth_GBytes',
  'RandomlyOrderedRingLatency_usec' );

my @walkunits = (
  'Tera Flops per Second',   'Tera Bytes per Second',
  'Giga Updates per Second', 'Tera Flops per Second',
  'Tera Bytes per Second',   'Giga Bytes per Second',
  'Giga Flops per Second',   'Giga Bytes per second',
  'micro-seconds');

my %crosswalk = (
    HPL_Tflops                          => 'G-HPL',
    PTRANS_GBs                          => 'G-PTRANS',
    MPIRandomAccess_GUPs                => 'G-RandomAccess',
    MPIFFT_Gflops                       => 'G-FFT',
    'StarSTREAM_Triad*CommWorldProcs'   => 'EP-STREAM Sys',
    CommWorldProcs                      => 'MPI Processes',
    StarSTREAM_Triad                    => 'EP-STREAM Triad',
    StarDGEMM_Gflops                    => 'EP-DGEMM',
    RandomlyOrderedRingBandwidth_GBytes => 'RandomRing Bandwidth',
    RandomlyOrderedRingLatency_usec     => 'RandomRing Latency'
);

if ( $opt_a ) { show_all(); }
if ( $opt_w ) { show_web(); }

###  Show all parameters from hpcc output file 
sub show_all {
    foreach $key ( sort keys %hpccoutf ) {
        print $key . ': ' . $hpccoutf{"$key"} . "\n";
    }
}

###  Show web parameters from hpcc output file
sub show_web {
    print "\n";
    foreach $key ( @walkorder ) {
        $value = '';
        if ( $key eq 'StarSTREAM_Triad*CommWorldProcs') {
            $value = $hpccoutf{'StarSTREAM_Triad'} * $hpccoutf{'CommWorldProcs'} / 1000;
            write;
        } elsif ( $key eq 'PTRANS_GBs' or $key eq 'MPIFFT_Gflops' ) {
            $value = $hpccoutf{$key} / 10000;
            write;
        } else {
            $value = $hpccoutf{$key};
            write;
        }
        $count++;
    }
    print '--------------------------------------------------------------------------------------------------'."\n";
}

format STDOUT_TOP =
--------------------------------------------------------------------------------------------------
HPCCOUTF NAME                        WEB NAME                      VALUE   UNITS
--------------------------------------------------------------------------------------------------
.

format STDOUT =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  @<<<<<<<<<<<<<<<<<<<  @#######.####   @<<<<<<<<<<<<<<<<<<<<<<
$key,                                $crosswalk{$key},     $value,         $walkunits[$count]
.