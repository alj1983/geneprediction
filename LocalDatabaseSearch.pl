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
    my $q = $id =~ /^.*\|(.*)\|/;
    my $comp=$1;
    my $filename=$comp."\_tblastnLocalDatabase\.out";
    my $tblastnrun = `tblastn -query temporarysequence.fasta -db $ARGV[1] -out $filename`;
    print (".");
    print MYOUTFILE "$filename\n";
}
print ("\n");

close MYOUTFILE;
      
