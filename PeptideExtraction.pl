#!/usr/bin/perl 

# Extract the longest polypeptide between * (stop
# codons) for each protein sequence in a fasta file

use strict;
use warnings;
use Bio::SeqIO;

if ($ARGV[0] eq '') { die "A fasta input file is needed"};

my $actualfastafile=$ARGV[0];
my %unique;

my $file   = $actualfastafile;
my $seqio  = Bio::SeqIO->new(-file => $file, -format => "fasta");
open(MYOUTFILE, ">LongestPolypeptide.fasta"); #open for write

while(my $seqs = $seqio->next_seq) {
    my $id  = $seqs->display_id;
    print MYOUTFILE "\>$id\n";
    my $seq = $seqs->seq;
    # Split each sequence into an array with * (stop codons) as separators.
    my @seqsplit = split('\*',$seq);
    # Chose for each fasta entry the longest polypeptide sequence, so
    # the longest entry in seqsplit
    my $seqsplitlength=$#seqsplit;
    if ($seqsplitlength>0){
	my @seqlengths;
	# identify the lengths of the polypeptides
	foreach(@seqsplit){
	    push @seqlengths,length($_);
	}
	##	Select the one with the greatest length
	my $idxMax = 0;
	$seqlengths[$idxMax] > $seqlengths[$_] or $idxMax = $_ for 1 .. $#seqlengths;
	my $longest_seqlength=$seqsplit[$idxMax];
	# get the sequence from the starting methionine
	if ($longest_seqlength =~ /.*?(M\w*)\.*.*/){
	    print MYOUTFILE "$1\n\n";
	}
	else{
	    print MYOUTFILE "$longest_seqlength\n\n";
	}
    }
    else{
	my $seqsplit1=$seqsplit[0];
	if ($seqsplit1 =~ /.*?(M\w*)\.*.*/){
	    print MYOUTFILE "$1\n\n";
	}
	else{
	    print MYOUTFILE "$seqsplit1\n\n";
	}

    }
}
close MYOUTFILE;

    
