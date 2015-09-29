#!/usr/bin/perl
use strict;
use warnings;

use Bio::SeqIO;

if ( $ARGV[0] eq '' ) {
    die "fasta file with protein sequence(s) has to be provided as argument to the
script\n"; }

if ( $ARGV[1] eq '' ) {
    die "The protein database (nr) of taxa related to your target species has to be provided as second argument to the script\n"; }

my $seqio  = Bio::SeqIO->new(-file => $ARGV[0], -format => "fasta");


open (MYOUTFILE, ">Blastp\.outfiles"); #open for write - this will contain
				   #the filenames of the output files
				   #this script produces
print ("waiting for blastp against the protein database (nr)...");

my $comp;
while (my $seqs = $seqio->next_seq) {
    my $seqout = Bio::SeqIO->new(-file => '>temporarysequence.fasta', -format => 'fasta' );
    $seqout->write_seq($seqs);
    my $id  = $seqs->display_id;
    if ($id =~ /^.*gb\|(.*)\|\s*(\S*.*)/){
	$comp=$1;
    }
    elsif ($id =~ /^.*emb\|(.*)\|\s*(\S*.*)/){
	$comp=$1;
    }
    elsif ($id =~ /^.*dbj\|(.*)\|\s*(\S*.*)/){
	$comp=$1;
    }
    # NBRF PIR                          pir||entry
    elsif ($id =~ /^.*pir\|\|(\S*)\s*(\S*.*)/){
	$comp=$1;
    }
    # Protein Research Foundation       prf||name
    elsif ($id=~ /^.*prf\|.*\|(.*)\s*(\S*.*)/){
	$comp=$1;
    }
    # SWISS-PROT                        sp|fasta$id|name
    elsif ($id=~ /^.*sp\|(.*)\|\s*(\S*.*)/){
	$comp=$1;
    }
    # Brookhaven Protein Data Bank (1)  pdb|entry|chain
    elsif ($id=~ /^.*pdb\|(.*)\|\s*(\S*.*)/){
	$comp=$1;
    }
    # GenInfo Backbone $Id               bbs|number
    elsif ($id=~ /^.*bbs\|(.*)/){
	$comp=$1;
    }
    # General database $identifier       gnl|database|$identifier
    elsif ($id=~ /^.*gnl\|(.*)\|\s*(\S*.*)/){
	$comp=$1;
    }
    # NCBI Reference Sequence           ref|fasta$id|locus
    elsif ($id=~ /^.*ref\|(.*)\|\s*(\S*.*)/){
	$comp=$1;
    }
    # tpg|DAA34093.1| 
    elsif ($id=~ /^.*tpg\|(.*)\|\s*(\S*.*)/){
	$comp=$1;
    }
    # Local Sequence $identifier         lcl|$identifier
    elsif ($id=~ /^.*lcl\|(.*)\|\s*(\S*.*)/){
	$comp=$1;
    }
    elsif ($id=~ /^.*lcl\|(\S*)\s*\S*/){
	$comp=$1;
    }
    elsif ($id){
	$comp="";
    }

    my $filename=$comp."\_blastpresult\.out";
    my $blastprun = `blastp -query temporarysequence.fasta -db $ARGV[1] -out $filename`;
    print (".") ;
    print MYOUTFILE "$filename\n";  
}
print ("\n");

close MYOUTFILE;
      
