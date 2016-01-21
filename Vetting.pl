#!/usr/bin/perl
use strict;
use warnings;

use Bio::SeqIO;

if ( $ARGV[0] eq '' ) {
    die "fasta file with protein sequence(s) has to be provided as first argument to the
script\n"; }

if ( $ARGV[1] eq '' ) {
    die "fasta file with EST sequence(s) of your target species has to be provided as second argument to the
script\n"; }

#my $prog = 'tblastn';
#my $db   = 'est';
#my $e_val= '1e-10';

#$Bio::Tools::Run::RemoteBlast::HEADER{'ENTREZ_QUERY'} = 'Calanus finmarchicus [ORGN]';
#$Bio::Tools::Run::RemoteBlast::RETRIEVALHEADER{'FORMAT_TYPE'} = 'Text';

my $seqio = Bio::SeqIO->new(-file=>$ARGV[0] , '-format' => 'fasta' );

open (MYOUTFILE, ">tblastn\.ESTout"); #open for write - this will contain
				   #the filenames of the output files
				   #this script produces

print ("waiting for tblastn against the EST database...");
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
    # GenInfo Backbone $Id               bbs|number
    elsif ($id=~ /^.*bbs\|(.*)/){
	$comp=$1;
    }
    # General database $identifier       gnl|database|$identifier
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
    # Local Sequence $identifier         lcl|$identifier
    elsif ($id=~ /^.*lcl\|(.*)\|\s*(\S*.*)/){
	$comp=$1;
    }
    elsif ($id=~ /^.*lcl\|(.*)\s*\S*/){
	$comp=$1;
    }
    elsif ($id){
	$comp="";
    }
    #my $q = $id =~ /^.*\|(.*)\|/;
    #my $comp=$1;
    my $filename=$comp."\_tblastnresult\.out";
    my $blastprun = `tblastn -query temporarysequence.fasta -db $ARGV[1] -out $filename`;
    print (".");
    print MYOUTFILE "$filename\n";
}
print ("\n");

close MYOUTFILE;
      