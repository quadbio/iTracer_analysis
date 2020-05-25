#!/usr/bin/perl -w
use strict;

####################################################
#Identifies and extracts GFP barcodes that can be used for lineage tracing using LINNAEUS
#
#Input: position sorted bam file obtained after running cellranger count against a GFP only reference
#Output: LINNAEUS compatible csv file containing necessary barcode information
#
#Author: Tobias Gerber
####################################################





foreach (@ARGV) {
	open BAM1,"samtools view $_ |" or die "Ooops";
	open BAM2,"samtools view $_ |" or die "Ooops";
    	open WRITEPRE, ">$_.RFPbarcode.csv" or die "could not write file\n";
	print WRITEPRE "Sequence,Barcode,Library,Cell,Reads,UMI,GeneBarcode,plasmid\n";
	open WRITESUMMARY, ">$_.RFPbarcodeFreq.csv" or die "could not write file\n";
    	print WRITESUMMARY "RFPbarcode,freq\n";

	my (%Reads_pSMAL,%Reads_pSBbi,%Count,%Codesum_pSBbi,%Codesum_pSMAL);

	print "Collect nRead\n";
   	while (<BAM1>) {
		chomp;
		my @bamline = split /\t/,$_;
		my ($barcode,$UMI);
		foreach (@bamline) {
			if ($_ =~ /^(CR\:Z\:)(.*)/) {						#extract raw cellular barcode
				$barcode = $2;
			}
			if ($_ =~ /^(CB\:Z\:)(.*)(\-1)/) {					#extract corrected cellular barcode if there
				my $opt_barcode = $2;
					if ($barcode ne $opt_barcode) {
					$barcode = $opt_barcode;
				}
			} 
			if ($_ =~ /^(UR\:Z\:)(.*)/) {						#extract raw UMI
				$UMI = $2;
			}
			if ($_ =~ /^(UB\:Z\:)(.*)/) {						#extract corrected UMI if there
				my $opt_UMI = $2;
				if ($UMI ne $opt_UMI) {
					$UMI = $opt_UMI;
				}
			} 
		}
		my $RFPcode_pSBbi;         
		if ($bamline[9] =~ /(AGTGATCC)([A-Z]{8,})(CCAACC)/ ) {                         #extract barcode for pSBbi barcode that is at least 8 bp long
			$RFPcode_pSBbi = $2;
#			print $2;
		
		}
		if ($RFPcode_pSBbi) {
			my $cell_UMI_barcode = "$barcode\_$UMI\_$RFPcode_pSBbi\n";
			$Reads_pSBbi{$cell_UMI_barcode}++;
			$Codesum_pSBbi{$RFPcode_pSBbi}++;
		}
		my $GFPcode_pSMAL;         
		if ($bamline[9] =~ /(AAACCGGT)([A-Z]{6,})(GAATTCGA)/ ) {                         #extract barcode for pSMAL barcode that is at least 6 bp long
			$GFPcode_pSMAL = $2;
#			print $2;
			
		}
		if ($GFPcode_pSMAL) {
			my $cell_UMI_barcode = "$barcode\_$UMI\_$GFPcode_pSMAL\n";
			$Reads_pSMAL{$cell_UMI_barcode}++;
			$Codesum_pSMAL{$GFPcode_pSMAL}++;
		}
			
     	 
    	}
	foreach (keys %Codesum_pSBbi) {
		my $size = scalar keys %Codesum_pSBbi;
		$Codesum_pSBbi{$_} = $Codesum_pSBbi{$_} / $size;
	}
	foreach (keys %Codesum_pSMAL) {
		my $size = scalar keys %Codesum_pSMAL;
		$Codesum_pSMAL{$_} = $Codesum_pSMAL{$_} / $size;
	}


	my %Code_freq;
	print "Collect Barcode information\n";
    	while (<BAM2>) {
     	   	chomp;
     	 	my @bamline = split /\t/,$_;	
		my ($barcode,$UMI,$start_pos,$cell_seq);								
		foreach (@bamline) {
			if ($_ =~ /^(CR\:Z\:)(.*)/) {						#extract raw cellular barcode
				$barcode = $2;
			}
			if ($_ =~ /^(CB\:Z\:)(.*)(\-1)/) {					#extract corrected cellular barcode if there
				my $opt_barcode = $2;
				if ($barcode ne $opt_barcode) {
					$barcode = $opt_barcode;
				}
			} 
			if ($_ =~ /^(UR\:Z\:)(.*)/) {						#extract raw UMI
				$UMI = $2;
			}
			if ($_ =~ /^(UB\:Z\:)(.*)/) {						#extract corrected UMI if there
				my $opt_UMI = $2;
				if ($UMI ne $opt_UMI) {
					$UMI = $opt_UMI;
				}
			} 
		}
		
		$start_pos = $bamline[3]-1;
		my $cell_UMI_barcode;
		my $RFPcode_pSBbi;
		my $GFPcode_pSMAL;
		if ($bamline[9] =~ /(AGTGATCC)([A-Z]{8,})(CCAACC)/) {
			$RFPcode_pSBbi = $2;
			$cell_UMI_barcode = "$barcode\_$UMI\_$RFPcode_pSBbi\n";
		} elsif ($bamline[9] =~ /(AAACCGGT)([A-Z]{6,})(GAATTCGA)/) {
			$GFPcode_pSMAL = $2;
			$cell_UMI_barcode = "$barcode\_$UMI\_$GFPcode_pSMAL\n";
		} else {
			next;
		}
		
#		print $cell_seq;
		if ($Count{$cell_UMI_barcode}) {						
			next;
		} else {
			if ($Reads_pSBbi{$cell_UMI_barcode} && $Reads_pSBbi{$cell_UMI_barcode} >= 3) {
				print WRITEPRE "$bamline[9],$barcode,na,$barcode,$Reads_pSBbi{$cell_UMI_barcode},$UMI,$RFPcode_pSBbi,pSBbi\n";	
				$Count{$cell_UMI_barcode}++;
				$Code_freq{$bamline[5]}++;
			} elsif ($Reads_pSMAL{$cell_UMI_barcode} && $Reads_pSMAL{$cell_UMI_barcode} >= 3) {
				print WRITEPRE "$bamline[9],$barcode,na,$barcode,$Reads_pSMAL{$cell_UMI_barcode},$UMI,$GFPcode_pSMAL,pSMAL\n";	
				$Count{$cell_UMI_barcode}++;
				$Code_freq{$bamline[5]}++;
				
			}				
		}
		
		
					
    	}
	foreach (keys %Code_freq) {
		print WRITESUMMARY "$_,$Code_freq{$_}\n";
	}
    	close BAM1;
	close BAM2;
	close WRITEPRE;	
	
	
}

#open PREFILE, "possorted_genome_bam.bam.prefilter.Scars.csv" or die "Ooops";
#open WRITEFINAL, ">possorted_genome_bam.bam.final.Scars.csv" or die "could not write file\n";

#my %CountUMI;
#while (<PREFILE>) {
#	chomp;
#	my @line = split /,/,$_;
#	my $cellUMI = "$line[3]\_$line[6]";
#	print "$cellUMI\n";
#	if ($CountUMI{$cellUMI}) {											#check if hash under the cell/UMI id key contains already something if yes skip this line to avoid having the same UMIs for different reads
#		next;
#	} else {
#		my $line = join(',',@line);
#		print WRITEFINAL "$line\n";
#		$CountUMI{$cellUMI}++;
#	}
#	 	
#}	
#close PREFILE;
#close WRITEFINAL;
	

