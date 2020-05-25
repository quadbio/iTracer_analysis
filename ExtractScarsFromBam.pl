#!/usr/bin/perl -w
use strict;

####################################################
#Identifies and extracts scars that can be used for lineage tracing using LINNAEUS
#
#Input: position sorted bam file obtained after running cellranger count against a GFP reference
#Output: LINNAEUS compatible csv file containing necessary scar information
#
#Author: Tobias Gerber
####################################################





foreach (@ARGV) {
	open BAM1,"samtools view $_ |" or die "Ooops";
	open BAM2,"samtools view $_ |" or die "Ooops";
    	open WRITEPRE, ">$_.unique.Scars.csv" or die "could not write file\n";
	print WRITEPRE "Sequence,Barcode,Library,Cell,Reads,CIGAR,UMI,Presence,Experiment,p,Scar\n";
	open WRITESUMMARY, ">$_.CIGARfreq.csv" or die "could not write file\n";
    	print WRITESUMMARY "CIGAR,freq\n";

	my (%Reads,%Count,%CIGARsum);

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
		if ($bamline[5] eq "91M") {										
			if ($bamline[2] eq "Scar_Tomato") {
				my $cell_UMI_scar = "$barcode\_$UMI\_$bamline[5]";
				$Reads{$cell_UMI_scar}++;
				$CIGARsum{$bamline[5]}++;			
			}
		} else {
			if ($bamline[2] eq "Scar_Tomato") {								#select start position
				if (($bamline[5] =~ /^\d+S/) or ($bamline[5] =~ /N$/) or ($bamline[5] =~ /^\d+N/)) {				#filter out bad scars
					next;
				} elsif (($bamline[5] =~ /.*I.*/) or ($bamline[5] =~ /.*D.*/) or ($bamline[5] =~ /.*N.*/) or ($bamline[5] =~ /.*S.*/)) {			#select sacrs with InDels
					my $cell_UMI_scar = "$barcode\_$UMI\_$bamline[5]";
					$Reads{$cell_UMI_scar}++;
					$CIGARsum{$bamline[5]}++;
				}
			}
     	 	}
    	}
	foreach (keys %CIGARsum) {
		my $size = scalar keys %CIGARsum;
		$CIGARsum{$_} = $CIGARsum{$_} / $size;
	}


	my %CIGAR_freq;
	print "Collect Scar information\n";
    	while (<BAM2>) {
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
		if ($bamline[5] eq "91M") {	
			if ($bamline[2] eq "Scar_Tomato") {
				my $cell_UMI_scar = "$barcode\_$UMI\_$bamline[5]";
				#my $cell_UMI = "$barcode\_$UMI";
				my $start_pos = $bamline[3]-1;
				if ($Count{$cell_UMI_scar}) {						
					next;
				} else {
					if ($Reads{$cell_UMI_scar} >= 3) {
						print WRITEPRE "$bamline[9],$barcode,na,$barcode,$Reads{$cell_UMI_scar},$bamline[5],$UMI,na,na,$CIGARsum{$bamline[5]},$start_pos\:$bamline[5]\n";	
						$Count{$cell_UMI_scar}++;
						$CIGAR_freq{$bamline[5]}++;
						#$Count{$cell_UMI} = $Reads{$cell_scar};	
					} else {
						next;
						
					}				
				}
			}

		} else {
			if ($bamline[2] eq "Scar_Tomato") {									#select start position
				if (($bamline[5] =~ /^\d+S/) or ($bamline[5] =~ /S$/) or ($bamline[5] =~ /N$/) or ($bamline[5] =~ /^\d+N/) ) {				#filter out bad scars
					next;
				} elsif (($bamline[5] =~ /.*I.*/) or ($bamline[5] =~ /.*D.*/) or ($bamline[5] =~ /.*N.*/) or ($bamline[5] =~ /.*S.*/)) {				#select sacrs with InDels
					my $cell_UMI_scar = "$barcode\_$UMI\_$bamline[5]";
					#my $cell_UMI = "$barcode\_$UMI";
					my $start_pos = $bamline[3]-1;
					if ($Count{$cell_UMI_scar}) {						
						next;
					} else {
						if ($Reads{$cell_UMI_scar} >= 3) {
							print WRITEPRE "$bamline[9],$barcode,na,$barcode,$Reads{$cell_UMI_scar},$bamline[5],$UMI,na,na,$CIGARsum{$bamline[5]},$start_pos\:$bamline[5]\n";	
							$Count{$cell_UMI_scar}++;
							$CIGAR_freq{$bamline[5]}++;
							#$Count{$cell_UMI} = $Reads{$cell_scar};	
						} else {
							next;
						
						}				
					}
				}
			}	
		}			
    	}
	foreach (keys %CIGAR_freq) {
		print WRITESUMMARY "$_,$CIGAR_freq{$_}\n";
	}
    	close BAM1;
	close BAM2;
	close WRITEPRE;	
	
	
}

#open PREFILE, "possorted_genome_bam.bam.unique.Scars.csv" or die "Ooops";
#open WRITEFINAL, ">possorted_genome_bam.bam.final.Scars.csv" or die "could not write file\n";

#my %CountUMI;
#while (<PREFILE>) {
#	chomp;
#	my @line = split /,/,$_;
#	my $cellUMI = "$line[3]\_$line[6]";
#	
#	if ($CountUMI{$cellUMI}) {
#		print "$cellUMI\t$line[5]\n";									#check if hash under the cell/UMI id key contains already something if yes skip this line to avoid having the same UMIs for different reads
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
	

