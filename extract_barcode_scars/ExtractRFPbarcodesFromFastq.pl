#!/usr/bin/perl -w
use strict;
use List::Util qw(min max);
    
##########################################
#  Counting barcodes for pSMAL vector    #
#  Input: fastq file             #
##########################################
my %Barcodes;
foreach (@ARGV) {
    my $file = $_;
    open FH, $file or die "Ooops";
    open WRITE, ">$file.Barcodes.txt" or die "could not write file\n";
    print WRITE "Barcode\tFrequency\n";
    open WRITE2, ">$file.FrequencySummary.txt" or die "could not write file\n";
    print WRITE2 "Frequency\tnBarcode\n";
    my $Start="CGAATTC"; 
    my $End= "ACCGGTT";
    my $count = 0;
    while (<FH>) {
                chomp;
                my $seq = $_;
                $count ++;
                if (($count == 1) | ($count == 3)) {
                    next;
                } elsif ($count == 2) {
                    if ($seq =~ /($Start)(.*)($End)/) {
                        my $Barcode = $2;
                        if (length $Barcode < 10) {
                            $Barcodes{$Barcode}++;
                        }
                    }
                } elsif ($count == 4) {
                    $count = 0;
                }
    }
	
	my %seen_barcodes;
    for (keys %Barcodes) {
        $seen_barcodes{($Barcodes{$_})}++;
    }
    my $size = scalar keys %Barcodes;
    print "$size\n";
    foreach (keys %Barcodes) {
        print WRITE "$_\t$Barcodes{$_}\n";
    }
    foreach (keys %seen_barcodes) {
        print WRITE2 "$_\t$seen_barcodes{$_}\n";
    }    
}
close FH;
close WRITE;
close WRITE2;
