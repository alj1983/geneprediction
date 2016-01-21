#!/usr/bin/perl -w

# Script to extract Best hit with information such as score and evalue
# and description from a blast text output;
use strict;
use warnings;

if ($ARGV[0] eq '' ) { die "Cfin.compinfo has to be provided as input\n"; }

if ($ARGV[1] eq '' ) { die "Transcriptome fasta file has to be provided as second argument\n"; }

open(IN, $ARGV[0]) or die "Can not open file $ARGV[0]";

# Get the unique comp hit identifiers. 

my $hit=0;
while (<IN>){
    chomp;
    my $line=$_;
    if ($line=~ /^hit:\s*(\S*)/){
	my $identifier=$1;
	if ($hit==0){
	    my $blastdbcmd = `blastdbcmd -db $ARGV[1] -entry $identifier > Comphits.fasta`;    
	}
	else{
	    my $blastdbcmd = `blastdbcmd -db $ARGV[1] -entry $identifier >> Comphits.fasta`;    
	}
	$hit++;
    }
}
close IN;

