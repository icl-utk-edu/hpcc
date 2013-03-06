#!/usr/bin/perl

#
# Usage: hpccoutf.pl < hpccoputf.txt
#

use strict;
use diagnostics;

$/       = undef;
my $file = <STDIN>;
my %hpccout = $file =~ /^(\w+)=(\d.*)$/mg;

my @walkorder = ( 'HPL_Tflops',            'PTRANS_GBs',            'MPIRandomAccess_GUPs',    'MPIFFT_Gflops',        'StarSTREAM_Triad*CommWorldProcs', 'StarSTREAM_Triad',      'StarDGEMM_Gflops',      'RandomlyOrderedRingBandwidth_GBytes', 'RandomlyOrderedRingLatency_usec' );
my @walkunits = ( 'Tera Flops per Second', 'Tera Bytes per Second', 'Giga Updates per Second', 'Tera Flops per Second', 'Tera Bytes per Second',          'Giga Bytes per Second', 'Giga Flops per Second', 'Giga Bytes per second',               'micro-seconds');
my %crosswalk = (
    HPL_Tflops           => 'G-HPL',
    PTRANS_GBs           => 'G-PTRANS',
    MPIRandomAccess_GUPs =>  'G-RandomAccess',
    MPIFFT_Gflops        =>  'G-FFT',
    "StarSTREAM_Triad*CommWorldProcs"   =>  'EP-STREAM Sys',
    CommWorldProcs                      =>  'MPI Processes',
    # StarSTREAM_Triad * CommWorldProcs => EP-STREAM Sys
    StarSTREAM_Triad                    =>  'EP-STREAM Triad',
    StarDGEMM_Gflops                    =>  'EP-DGEMM',
    RandomlyOrderedRingBandwidth_GBytes =>  'RandomRing Bandwidth',
    RandomlyOrderedRingLatency_usec     =>   'RandomRing Latency'
);

my $show_all = 0;
my $show_web = 1;

if ( $show_all ) {
    show_all();
}
if ( $show_web ) {
    show_web(); 
}

sub show_all {
    foreach my $key ( sort keys %hpccout ) {
        print $key . ': ' . $hpccout{"$key"} . "\n";
    }
}

sub show_web {
    my $count = 0;
    foreach my $key ( @walkorder ) {
        if ( $key eq 'StarSTREAM_Triad*CommWorldProcs') {
            print "$key ($crosswalk{$key}) " . $hpccout{'StarSTREAM_Triad'} * $hpccout{'CommWorldProcs'} . " $walkunits[$count]\n";
        } else {
            print "$key ($crosswalk{$key}) $hpccout{$key} $walkunits[$count]\n";
        }
        $count++;
    }
}