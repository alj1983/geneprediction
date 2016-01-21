#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::SeqIO;

if ($ARGV[0] eq '' or $ARGV[1] eq '' ) { die "Two input files are needed. First the
output from Parse_DeCypher.pl and second a fasta file with hits in the
C. finmarchicus transcriptome\n"; }

# open the output file from Parse_DeCypherAndConverge.pl

open(IN, $ARGV[0]) or die "Could not open file $ARGV[0]\n";

my @uniquehits;
while(<IN>){
    my $actualline=$_;
    if ($actualline =~ /^hit:\s.*/){
	my $uh = $actualline =~ /^hit:\s(.*)$/;
	push @uniquehits,$1;
    }
}
close IN;


# Search for these uniquehits in the fasta file
open(IN, $ARGV[1]) or die "Could not open the fasta file $ARGV[1]\n";

my $seq_in = Bio::SeqIO->new(
    -file   => "<$ARGV[1]",
    -format => 'Fasta',
    );


my $outfilename="UniqueCompHits.fasta";

my $seq_out = Bio::SeqIO->new(
    -file   => ">$outfilename",
    -format => 'Fasta',
    );

my $fastaid2;
while (my $seq = $seq_in->next_seq) {
    my $fastaid=$seq->display_id;
    if ($fastaid =~ /^.*gb\|(.*)\|\s*(\S*.*)/){
	$fastaid2=$1;
    }
    elsif ($fastaid =~ /^.*emb\|(.*)\|\s*(\S*.*)/){
	$fastaid2=$1;
    }
    elsif ($fastaid =~ /^.*dbj\|(.*)\|\s*(\S*.*)/){
	$fastaid2=$1;
    }
    # NBRF PIR                          pir||entry
    elsif ($fastaid =~ /^.*pir\|\|(\S*)\s*(\S*.*)/){
	$fastaid2=$1;
    }
    # Protein Research Foundation       prf||name
    elsif ($fastaid=~ /^.*prf\|.*\|(.*)\s*(\S*.*)/){
	$fastaid2=$1;
    }
    # SWISS-PROT                        sp|fastaid|name
    elsif ($fastaid=~ /^.*sp\|(.*)\|\s*(\S*.*)/){
	$fastaid2=$1;
    }
    # Brookhaven Protein Data Bank (1)  pdb|entry|chain
    elsif ($fastaid=~ /^.*pdb\|(.*)\|\s*(\S*.*)/){
	$fastaid2=$1;
    }
    # GenInfo Backbone Id               bbs|number
    elsif ($fastaid=~ /^.*bbs\|(.*)/){
	$fastaid2=$1;
    }
    # General database identifier       gnl|database|identifier
    elsif ($fastaid=~ /^.*gnl\|(.*)\|\s*(\S*.*)/){
	$fastaid2=$1;
    }
    # NCBI Reference Sequence           ref|fastaid|locus
    elsif ($fastaid=~ /^.*ref\|(.*)\|\s*(\S*.*)/){
	$fastaid2=$1;
    }
    # tpg|DAA34093.1| 
    elsif ($fastaid=~ /^.*tpg\|(.*)\|\s*(\S*.*)/){
	$fastaid2=$1;
    }
    # Local Sequence identifier         lcl|identifier
    elsif ($fastaid=~ /^.*lcl\|(.*)\|\s*(\S*.*)/){
	$fastaid2=$1;
    }
    elsif ($fastaid=~ /^.*lcl\|(.*).*/){
	$fastaid2=$1;
    }
    elsif ($fastaid){
	$fastaid2="";
    }
    foreach (@uniquehits){
	my $actual_uniquehit=$_;
	if ($fastaid2 eq $actual_uniquehit){
#		    print "$fastaid2\n";
	    $seq_out->write_seq($seq);
	}
	#print $fastaid,"\n";
#	    print $seq->seq()."\n";
    }
}
