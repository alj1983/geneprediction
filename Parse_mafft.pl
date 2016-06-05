#!/usr/bin/perl

use strict;
use warnings;

# check if argument to the script is there.  
if ( $ARGV[0] eq '' ) {
die "A file listing the names of the ClustalW alignment files has to be provided as argument to the
script. The filename must be of the pattern queryname__besthitname__MAFFTalignment.fasta\n"; }


# Load in the alignment file names
open(INFILES, "$ARGV[0]") or die "Could not open file $ARGV[0]\n";
my @filenames;
while (<INFILES>) {
    chomp;
    push @filenames,$_;
}
close INFILES;     

open (OUTLIST,">parsemafft.out");
print OUTLIST "query\thit\tidentity\tsimilarity\n";
foreach (@filenames){
    my $actualfilename=$_;
    my $aln='';
    
    my @lines=(); # this saves the lines that contain the comparison characters of the two alignments
# read the sequence from the input file
    open(IN, $actualfilename) or die "Could not open $actualfilename\n";
    while (<IN>) {
	chomp;

#Remove lines that contain only white space
	if (!/^\s*$/){
	    # Remove lines that have a non-white-space character at the
	    # beginning and thus work only with the lines that show the matching signs '*', ';', and '.'
	    if (!/^\S{1}/) {
		$aln .= $_;
		
		push @lines, $.;
		
	    }
	}
    }
    close IN;

## Getting the two sequences (without their names):
# The lines for the first sequence are
    my @lines1=@lines;
    foreach my $x (@lines1) { $x = $x-1; } #my @lines1 = @lines-1;
# The lines for the second sequence are
    my @lines2=@lines;
    foreach my $x (@lines2) { $x = $x-2; } #my @lines1 = @lines-1;
    
    my $sequence1='';
    my $sequence2='';

    open(IN, "$actualfilename") or die "Could not open file $actualfilename\n";
    while (<IN>) {
	chomp;
# Chose only the lines for sequence 2
	
	if($. ~~ @lines1) {
	    my $s1 = $_;
	    $s1 =~ /^\S*\s*(\S*)$/;
	    my $seq1 = $1;
	    $sequence1 .= $seq1;
	}
	
# Chose only the lines for sequence 2
	if($. ~~ @lines2) {
	    my $s2 = $_;
	    $s2 =~ /^\S*\s*(\S*)$/;
	     my $seq2 = $1;
	    $sequence2 .= $seq2;
	}
    }
    close IN;
    
# Count the number of amino acids in the sequences:
    my $count1= ($sequence1=~s/\w//g);
    my $count2= ($sequence2=~s/\w//g);
    
    my $count='';
    
    if ($count1>$count2){
	$count=$count1;
    }
    else{
	$count=$count2;
    }
    
    my $identical = ( $aln =~ tr/\*// );
    my $similar = ( $aln =~ tr/\.\*\:// );
    
    my $identicalreport=100*($identical/$count);
    my $similarreport=100*($similar/$count);
    my $comparison = $actualfilename =~ /^(.*)\_\_(.*)\_\_.*$/; 
    print OUTLIST "$1\t";
    print OUTLIST "$2\t";
    print OUTLIST "$identicalreport\t";
    print OUTLIST "$similarreport\n";
}
