#!/usr/bin/perl
use strict;
use warnings;

use Bio::SeqIO;

if ( $ARGV[0] eq '' ) {
    die "fasta file with protein sequence(s) has to be provided as argument to the
script\n"; }

if ( $ARGV[1] eq '' ) {
    die "The fasta file of the local transcriptome database needs to be provided as second argument to the script\n"; }


my $seqio = Bio::SeqIO->new(-file=>$ARGV[0] , '-format' => 'fasta' );

open (MYOUTFILE, ">tblastn\.LocalDatabase"); #open for write - this will contain
				   #the filenames of the output files
				   #this script produces

print ("waiting for tblastn against the local transcriptome database...");
while (my $seqs = $seqio->next_seq) {
    my $seqout = Bio::SeqIO->new(-file => '>temporarysequence.fasta', -format => 'fasta' );
    $seqout->write_seq($seqs);
    my $id  = $seqs->display_id;
    my $comp;
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
    # SWISS-PROT                        sp|accession|name
    elsif ($id=~ /^.*sp\|(.*)\|\s*(\S*.*)/){
	$comp=$1;
    }
    # Brookhaven Protein Data Bank (1)  pdb|entry|chain
    elsif ($id=~ /^.*pdb\|(.*)\|\s*(\S*.*)/){
	$comp=$1;
    }
    # GenInfo Backbone Id               bbs|number
    elsif ($id=~ /^.*bbs\|(.*)/){
	$comp=$1;
    }
    # General database identifier       gnl|database|identifier
    elsif ($id=~ /^.*gnl\|(.*)\|\s*(\S*.*)/){
	$comp=$1;
    }
    # NCBI Reference Sequence           ref|accession|locus
    elsif ($id=~ /^.*ref\|(.*)\|\s*(\S*.*)/){
	$comp=$1;
    }
    # tpg|DAA34093.1| 
    elsif ($id=~ /^.*tpg\|(.*)\|\s*(\S*.*)/){
	$comp=$1;
    }
    # Local Sequence identifier         lcl|identifier
    elsif ($id=~ /^.*lcl\|(.*)\|\s*(\S*.*)/){
	$comp=$1;
    }
    elsif ($id=~ /^.*lcl\|(\S*)\s*(\S*.*)/){
	$comp=$1;
    }
    elsif ($id=~ /^(\S*)/){
	$comp=$1;
    }
    elsif ($id){
	$comp=$1;
    }
    my $filename=$comp."\_tblastnLocalDatabase\.out";
    my $tblastnrun = `tblastn -query temporarysequence.fasta -db $ARGV[1] -out $filename`;
    print (".");
    print MYOUTFILE "$filename\n";
}
print ("\n");

close MYOUTFILE;
      
