#!/usr/bin/perl -w

package ProteinScreen;
use match::smart; ## Replaces the smartmatch operator ~~
use strict;
use warnings;
use Bio::SeqIO;
use LWP::Simple;
use List::MoreUtils qw(uniq); # This allows to find all unique values in an array


### Function definitions

## local_database_search
sub local_database_search {

    if ( $_[0] eq '' ) {
	die "fasta file with protein sequence(s) has to be provided as argument to the
script\n"; }

    if ( $_[1] eq '' ) {
	die "The fasta file of the local transcriptome database needs to be provided as second argument to the script\n"; }


    my $seqio = Bio::SeqIO->new(-file=>$_[0] , '-format' => 'fasta' );

    open (MYOUTFILE, ">tblastn\.LocalDatabase"); #open for write - this will contain
    #the filenames of the output files
    #this script produces

    print ("waiting for tblastn against the local transcriptome database...");
    while (my $seqs	=  $seqio->next_seq) {
	my $seqout	=  Bio::SeqIO->new(-file => '>temporarysequence.fasta', -format => 'fasta' );
	$seqout->write_seq($seqs);
	my $id	=  $seqs->display_id;
	my $comp;
	if ($id	=~ /^.*gb\|(.*)\|\s*(\S*.*)/){
	    $comp		=  $1;
	}
	elsif ($id							       =~ /^.*emb\|(.*)\|\s*(\S*.*)/){
	    $comp							       =  $1;
	}
	elsif ($id							       =~ /^.*dbj\|(.*)\|\s*(\S*.*)/){
	    $comp							       =  $1;
	}
	# NBRF PIR                          pir||entry
	elsif ($id							       =~ /^.*pir\|\|(\S*)\s*(\S*.*)/){
	    $comp							       =  $1;
	}
	# Protein Research Foundation       prf||name
	elsif ($id							       =~ /^.*prf\|.*\|(.*)\s*(\S*.*)/){
	    $comp							       =  $1;
	}
	# SWISS-PROT                        sp|accession|name
	elsif ($id							       =~ /^.*sp\|(.*)\|\s*(\S*.*)/){
	    $comp							       =  $1;
	}
	# Brookhaven Protein Data Bank (1)  pdb|entry|chain
	elsif ($id							       =~ /^.*pdb\|(.*)\|\s*(\S*.*)/){
	    $comp							       =  $1;
	}
	# GenInfo Backbone Id               bbs|number
	elsif ($id							       =~ /^.*bbs\|(.*)/){
	    $comp							       =  $1;
	}
	# General database identifier       gnl|database|identifier
	elsif ($id							       =~ /^.*gnl\|(.*)\|\s*(\S*.*)/){
	    $comp							       =  $1;
	}
	# NCBI Reference Sequence           ref|accession|locus
	elsif ($id							       =~ /^.*ref\|(.*)\|\s*(\S*.*)/){
	    $comp							       =  $1;
	}
	# tpg|DAA34093.1| 
	elsif ($id							       =~ /^.*tpg\|(.*)\|\s*(\S*.*)/){
	    $comp							       =  $1;
	}
	# Local Sequence identifier         lcl|identifier
	elsif ($id							       =~ /^.*lcl\|(.*)\|\s*(\S*.*)/){
	    $comp							       =  $1;
	}
	elsif ($id							       =~ /^.*lcl\|(\S*)\s*(\S*.*)/){
	    $comp							       =  $1;
	}
	elsif ($id							       =~ /^(\S*)/){
	    $comp							       =  $1;
	}
	elsif ($id){
	    $comp							       =  $1;
	}
	my $filename						       =  $comp."\_tblastnLocalDatabase\.out";
	my $tblastnrun						       =  `tblastn -query temporarysequence.fasta -db $_[1] -out $filename`;
	print (".");
	print MYOUTFILE "$filename\n";
    }
    print ("\n");

    close MYOUTFILE;
    
}      


## parse_local_database_search
# Script to extract Best hit in the local transcriptome/genome with
# information such as evalue and reading frame from tblastn output	and
# then selecting only those with an e-value <1e-05
sub parse_local_database_search {

    if ($_[0] eq '' ) { die "text file has to be provided that lists
the names of tblastn report files which shall be treated in this
script\n"; }

    open(IN, $_[0]) or die "Can not open file $_[0]";
# Get an array containing the filenames that shall be parsed here
    my @files;
    while (<IN>){
	chomp;
	push @files,$_;
	
    }
    close IN;



    open(MYOUTFILE1, ">Cfin\.queryinfo"); #open for write




    my @querynames1; ## This will hold the query names repeated for each hit.
    my @hitnames1; ## This will hold the hit names
    my @hitdescriptions1; ## This will hold the hit descriptions
    my @querydescriptions1; ## This will hold the query descriptions
    my @hitlength1; ## This will hold the hit length
    my @evals1; ## This will hold the e-values
    my @ranks1; ## This will hold the ranks
    my @queryframe1;  # this will hold the reading frame of the query
    my @targetframe1; # this will hold the reading frame of the target

#    open(IN, $_[0]) or die "Can not open file $_[0]";

    foreach (@files){
	my $f=$_;
	
	open (IN,$f) or die "Can not open file $f";
	
	my $hitnumber=0;
	my $hitnumber2=0;
	my $hitnumber3=0;
	my $querynumber=1;
	
	my $queryname1; # This will hold the single query name for one file
	my $querydescription1; # This will hold the single query description
	while(<IN>){
	    chomp;
	    my $line=$_;
	    if ($line =~ /^Query=\s(.*)/){
		my $accession = $1;
		### Test what kind of Fasta definition line format this is:
		# GenBank                           gi|gi-number|gb|accession|locus
		# EMBL Data Library                 gi|gi-number|emb|accession|locus
		# DDBJ, DNA Database of Japan       gi|gi-number|dbj|accession|locus
		if ($accession =~ /^.*gb\|(.*)\|\s*(\S*.*)/){
		    $queryname1=$1;
		    $querydescription1=$2;
		}
		elsif ($accession =~ /^.*emb\|(.*)\|\s*(\S*.*)/){
		    $queryname1=$1;
		    $querydescription1=$2;
		}
		elsif ($accession =~ /^.*dbj\|(.*)\|\s*(\S*.*)/){
		    $queryname1=$1;
		    $querydescription1=$2;
		}
		# NBRF PIR                          pir||entry
		elsif ($accession =~ /^.*pir\|\|(\S*)\s*(\S*.*)/){
		    $queryname1=$1;
		    $querydescription1=$2;
		}
		# Protein Research Foundation       prf||name
		elsif ($accession=~ /^.*prf\|.*\|(.*)\s*(\S*.*)/){
		    $queryname1=$1;
		    $querydescription1=$2;
		}
		# SWISS-PROT                        sp|accession|name
		elsif ($accession=~ /^.*sp\|(.*)\|\s*(\S*.*)/){
		    $queryname1=$1;
		    $querydescription1=$2;
		}
		# Brookhaven Protein Data Bank (1)  pdb|entry|chain
		elsif ($accession=~ /^.*pdb\|(.*)\|\s*(\S*.*)/){
		    $queryname1=$1;
		    $querydescription1=$2;
		}
		# GenInfo Backbone Id               bbs|number
		elsif ($accession=~ /^.*bbs\|(.*)/){
		    $queryname1=$1;
		    $querydescription1=$2;
		}
		# General database identifier       gnl|database|identifier
		elsif ($accession=~ /^.*gnl\|(.*)\|\s*(\S*.*)/){
		    $queryname1=$1;
		    $querydescription1=$2;
		}
		# NCBI Reference Sequence           ref|accession|locus
		elsif ($accession=~ /^.*ref\|(.*)\|\s*(\S*.*)/){
		    $queryname1=$1;
		    $querydescription1=$2;
		}
		# tpg|DAA34093.1| 
		elsif ($accession=~ /^.*tpg\|(.*)\|\s*(\S*.*)/){
		    $queryname1=$1;
		    $querydescription1=$2;
		}
		# Local Sequence identifier         lcl|identifier
		elsif ($accession=~ /^.*lcl\|(.*)\|\s*(\S*.*)/){
		    $queryname1=$1;
		    $querydescription1=$2;
		}
		elsif ($accession=~ /^(\S*)/){
		    $queryname1=$1;
		    $querydescription1="No description";
		}
		elsif ($accession){
		    $queryname1="";
		    $querydescription1="No description";
		}
	    }
	    if ($line =~ /^>/){
		push @querynames1, $queryname1;
		push @querydescriptions1, $querydescription1;
		### Test what kind of Fasta definition line format this is:
		# GenBank                           gi|gi-number|gb|accession|locus
		# EMBL Data Library                 gi|gi-number|emb|accession|locus
		# DDBJ, DNA Database of Japan       gi|gi-number|dbj|accession|locus
		if ($line =~ /^>.*gb\|(.*)\|\s*(\S*.*)/){
		    push @hitnames1, $1 ;
		    push @hitdescriptions1,$2;
		}
		elsif ($line =~ /^>.*emb\|(.*)\|\s*(\S*.*)/){
		    push @hitnames1,$1;
		    push @hitdescriptions1,$2;
		}
		elsif ($line =~ /^.*dbj\|(.*)\|\s*(\S*.*)/){
		    push @hitnames1,$1;
		    push @hitdescriptions1,$2;
		}
		# NBRF PIR                          pir||entry
		elsif ($line =~ /^.*pir\|\|(\S*)\s*(\S*.*)/){
		    push @hitnames1,$1;
		    push @hitdescriptions1,$2;
		}
		# Protein Research Foundation       prf||name
		elsif ($line=~ /^.*prf\|.*\|(.*)\s*(\S*.*)/){
		    push @hitnames1,$1;
		    push @hitdescriptions1,$2;
		}
		# SWISS-PROT                        sp|accession|name
		elsif ($line=~ /^.*sp\|(.*)\|\s*(\S*.*)/){
		    push @hitnames1,$1;
		    push @hitdescriptions1,$2;
		}
		# Brookhaven Protein Data Bank (1)  pdb|entry|chain
		elsif ($line=~ /^.*pdb\|(.*)\|\s*(\S*.*)/){
		    push @hitnames1,$1;
		    push @hitdescriptions1,$2;
		}
		# GenInfo Backbone Id               bbs|number
		elsif ($line=~ /^.*bbs\|(.*)/){
		    push @hitnames1,$1;
		    push @hitdescriptions1,$1;
		}
		# General database identifier       gnl|database|identifier
		elsif ($line=~ /^.*gnl\|(.*)\|\s*(\S*.*)/){
		    push @hitnames1,$1;
		    push @hitdescriptions1,$2;
		}
		# NCBI Reference Sequence           ref|accession|locus
		elsif ($line=~ /^.*ref\|(.*)\|\s*(\S*.*)/){
		    push @hitnames1,$1;
		    push @hitdescriptions1,$2;
		}
		# tpg|DAA34093.1| 
		elsif ($line=~ /^.*tpg\|(.*)\|\s*(\S*.*)/){
		    push @hitnames1,$1;
		    push @hitdescriptions1,$2;
		}
		# Local Sequence identifier         lcl|identifier
		elsif ($line=~ /^.*lcl\|(.*)\|\s*(\S*.*)/){
		    push @hitnames1,$1;
		    push @hitdescriptions1,$2;
		}
		elsif ($line=~ /^.*lcl\|(\S*)\s*(\S*.*)/){
		    push @hitnames1,$1;
		    push @hitdescriptions1,$2;
		}
		elsif ($line=~ /^(\S*)/){
		    push @hitnames1, $1;
		    push @hitdescriptions1,"";
		}
		elsif ($line){
		    push @hitnames1,"";
		    push @hitdescriptions1,"";
		}
		$hitnumber++;
		$hitnumber2=0;
		$hitnumber3=0;
		$querynumber=$hitnumber;
		if ($hitnumber>0){
		    push @ranks1, $hitnumber;
		}
	    }
	    

	    if ($line =~ /Expect.*=\s*(\S*),/g){
		my $e=$1;
		if ($hitnumber2<2){
		    #push @estscores,$1;
		    push @evals1,$e;
		    $hitnumber2=2;
		}
	    }
	    if ($hitnumber>0){
		#push @estscores,$1;
		if ($line =~ /Length=(.*)$/g){
		    push @hitlength1,$1; 
		}
	    }
	    
	    if ($line =~ /Frame = (.*)$/g){
		if ($hitnumber3<2){
		    push @targetframe1,$1; 
		    push @queryframe1,"1";
		    $hitnumber3=2;
		}
	    }
	}
	close IN;
    }


##### Remove those hits with an e-value >=1e-05

# Exclude those records without evalues
    my @excluding;
    @excluding = grep { $evals1[$_] eq "" } 0..$#evals1;

    my $ex;
    my $evalsl1 = @evals1;
    my @evals;
    my @querynames; ## This will hold the query names
    my @querydescriptions; ## This will hold the query descriptions
    my @hitnames; ## This will hold the hit names
    my @hitlength; ## This will hold the hit length
    my @ranks; ## This will hold the ranks
    my @queryframe;
    my @targetframe;


    for ($ex=0;$ex<$evalsl1;$ex++){
	if ($ex ~~ @excluding){}
	else{
	    push @evals,$evals1[$ex];
	    push @querynames,$querynames1[$ex];
	    push @hitnames,$hitnames1[$ex];
	    push @hitlength,$hitlength1[$ex];
	    push @ranks,$ranks1[$ex];
	    push @queryframe,$queryframe1[$ex];
	    push @targetframe,$targetframe1[$ex];
	    push @querydescriptions,$querydescriptions1[$ex];
	}
    }


    my $evalsl = @evals;
    my @querynames_sig;
    my @querydescriptions_sig;
    my @hitnames_sig;
    my @hitlength_sig;
    my @evals_sig;
    my @ranks_sig;
    my @queryframe_sig;
    my @targetframe_sig;

    my $o;
    for ($o=0;$o<$evalsl;$o++){
	if ($evals[$o]<1e-05){
	    push @querynames_sig, $querynames[$o];
	    push @querydescriptions_sig, $querydescriptions[$o];
	    push @hitnames_sig, $hitnames[$o];
	    push @hitlength_sig, $hitlength[$o];
	    push @evals_sig, $evals[$o];
	    push @ranks_sig, $ranks[$o];
	    push @queryframe_sig,$queryframe[$o];
	    push @targetframe_sig,$targetframe[$o];
	}
	
	
    }

#####
# Identify the unique comps for each query #####
# First get an array of the unique queries that were used: 

    my @unique_queries;
    @unique_queries = uniq @querynames_sig;
    my $unique_queriesl=@unique_queries;
    print MYOUTFILE1 "\nQueries used:\n";
    print MYOUTFILE1 "-------------\n";

    my $r;
    for ($r=0;$r<$unique_queriesl;$r++){
	print MYOUTFILE1 "Query: $unique_queries[$r]";
	my @querydescriptions_sig_index = grep { $querynames_sig[$_] eq $unique_queries[$r] } 0..$#querynames_sig;
	my $actualquerydescription=$querydescriptions_sig[$querydescriptions_sig_index[0]];
	print MYOUTFILE1 "\t$actualquerydescription\n";    
    }

    my $q;
    my @all_unique_comp_hits; # this will collect the unique comp hits of
    # all queries in a single array
    for ($q=0;$q<$unique_queriesl;$q++){
	print MYOUTFILE1 "\nUsedQuery: $unique_queries[$q]\n";
	print MYOUTFILE1 "--------------------\n";
	
	
	my $querynames_sigl = @querynames_sig;
	
	my $i;
	my @comp_hits;
	my @queryindices;
	for ($i=0;$i<$querynames_sigl;$i++){
	    if ("$querynames_sig[$i]" eq "$unique_queries[$q]"){
# get the comp hits for that query
		push @comp_hits, $hitnames_sig[$i];
		push @queryindices, $i;
	    }
	}
	my @compnames;
	foreach(@comp_hits)
	{
	    if ($_ =~ /^(.*)_.*_.*$/){
		push @compnames, $1;
	    }
	    else {
		push @compnames, $_;
	    }
	}
	
	print MYOUTFILE1 "\nUnique comp hits \n";
	print MYOUTFILE1 "---------------- \n";
#get the unique comp hits for that query
	my @unique_comp_hits = uniq @compnames;
	
	my @index;
	foreach(@unique_comp_hits){
	    
# identify the number of hits for each unique comp. If it is more than
# one, select only the one with the longest length, otherwise select
# the only one
	    my $comp_hit=$_;
	    @index = grep { $compnames[$_] eq $comp_hit } 0..$#compnames;
#    print $#index;
	    if ($#index>0){
		my @lengths=@hitlength_sig[@queryindices];
		my @complengths=@lengths[@index];
#	Select the one with the greatest length
		my $idxMax = 0;
		$complengths[$idxMax] > $complengths[$_] or $idxMax = $_ for 1 .. $#complengths;
		
		my @queryhits=@hitnames_sig[@queryindices];
		my @comphits=@queryhits[@index];
		my $longest_comphit=$comphits[$idxMax];
		print MYOUTFILE1 "$unique_queries[$q] unique: $longest_comphit\n";
		push @all_unique_comp_hits,$longest_comphit;
	    }
	    else {
		my @queryhits=@hitnames_sig[@queryindices];
		my @comphits=@queryhits[@index];
		print MYOUTFILE1 "$unique_queries[$q] unique: @comphits\n";
		push @all_unique_comp_hits,@comphits;
	    }
	    
	}
	

	
#########################################################
######## Identify the comps with multiple hits ##########
#########################################################
	

	print MYOUTFILE1 "\nMultiple comp hits \n";
	print MYOUTFILE1 "---------------- \n";
# get the unique comp hits for that query
	
	
	my @index2;
	my $n=1;
	foreach(@unique_comp_hits){
	    
# identify the number of hits for each unique comp. If it is more than
# one, select only the one with the longest length, otherwise select
# the only one
	    my $comp_hit=$_;
	    @index2 = grep { $compnames[$_] eq $comp_hit } 0..$#compnames;
#    print $#index2;
	    if ($#index2>0){
		my @lengths=@hitlength_sig[@queryindices];
		my @complengths=@lengths[@index2];
		
		my @queryhits=@hitnames_sig[@queryindices];
		my @comphits=@queryhits[@index2];
		my $x=1;
		foreach(@comphits)
		{
		    print MYOUTFILE1 "$unique_queries[$q] multiple $n.$x: $_\n";
		    $x++;
		}
		$n++;	
	    }
	    
	}
	
    }


    close MYOUTFILE1;

########################################################################
###### converge unique comp hits between the different queries #########
########################################################################
    open(MYOUTFILE2, ">Cfin\.compinfo"); #open for write
    my @overall_unique_comp_hits = uniq @all_unique_comp_hits;

# scan through all the table and print the queries, ranks and evals
# for all the unique comp hits:
    print MYOUTFILE2 "Unique comp hits in all queries\n";
    print MYOUTFILE2 "-------------------------------\n\n";
    foreach (@overall_unique_comp_hits){
	print MYOUTFILE2 "hit: $_\n";
    }
    print MYOUTFILE2 "-------------------------------\n\n";
    foreach (@overall_unique_comp_hits){
	my $actual_comp_hit=$_;
	print MYOUTFILE2 "Unique comp hit: $actual_comp_hit\n";
	print MYOUTFILE2 "---------------------------------\n";
	my @index;
	@index = grep { $hitnames_sig[$_] eq $actual_comp_hit } 0..$#hitnames_sig;
	print MYOUTFILE2 "Comp\tQuery\tRank\tE-value\tQueryframe\tTargetframe\n";
	foreach(@index){
	    print MYOUTFILE2 "$actual_comp_hit\t$querynames_sig[$_]\t$ranks_sig[$_]\t$evals_sig[$_]\t$queryframe_sig[$_]\t$targetframe_sig[$_]\n";
	}
	print MYOUTFILE2 "\n";
    }

    close MYOUTFILE2;

}



## best_hits_to_fasta

sub best_hits_to_fasta {


    if ($_[0] eq '' ) { die "Cfin.compinfo has to be provided as input\n"; }

    if ($_[1] eq '' ) { die "Transcriptome fasta file has to be provided as second argument\n"; }

    open(IN, $_[0]) or die "Can not open file $_[0]";

# Get the unique comp hit identifiers. 

    my $hit=0;
    while (<IN>){
	chomp;
	my $line=$_;
	if ($line=~ /^hit:\s*(\S*)/){
	    my $identifier=$1;
	    if ($hit==0){
		my $blastdbcmd = `blastdbcmd -db $_[1] -entry $identifier > Comphits.fasta`;    
	    }
	    else{
		my $blastdbcmd = `blastdbcmd -db $_[1] -entry $identifier >> Comphits.fasta`;    
	    }
	    $hit++;
	}
    }
    close IN;

}

## translation
sub translation {



# The genetic code is represented as a hash. Each codon (hash key)
# is associated with an amino acid (hash value).
# A stop codon is shown as '*'
    my %code = (
	'UUU' => 'F', 'UUC' => 'F', 'UUA' => 'L', 'UUG' => 'L',
	'CUU' => 'L', 'CUC' => 'L', 'CUA' => 'L', 'CUG' => 'L',
	'AUU' => 'I', 'AUC' => 'I', 'AUA' => 'I', 'AUG' => 'M',
	'GUU' => 'V', 'GUC' => 'V', 'GUA' => 'V', 'GUG' => 'V',
	'UCU' => 'S', 'UCC' => 'S', 'UCA' => 'S', 'UCG' => 'S',
	'CCU' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
	'ACU' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
	'GCU' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',
	'UAU' => 'Y', 'UAC' => 'Y', 'UAA' => '*', 'UAG' => '*',
	'CAU' => 'H', 'CAC' => 'H', 'CAA' => 'Q', 'CAG' => 'Q',
	'AAU' => 'N', 'AAC' => 'N', 'AAA' => 'K', 'AAG' => 'K',
	'GAU' => 'D', 'GAC' => 'D', 'GAA' => 'E', 'GAG' => 'E',
	'UGU' => 'C', 'UGC' => 'C', 'UGA' => '*', 'UGG' => 'W',
	'CGU' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',
	'AGU' => 'S', 'AGC' => 'S', 'AGA' => 'R', 'AGG' => 'R',
	'GGU' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G' 
	);

    if ($_[0] eq '') { die "A fasta input file is needed"};
    if ($_[1] eq '') { die "The file Cfin.compinfo is needed as second input file"};


    open(COMPINFO, $_[1]) or die "Can not open file $_[1]";

# Get the comp ids and in which target frame they shall be translated
    my @comps;
    my @frames;
    while(<COMPINFO>) {
	chomp;
	my $actualline=$_;
	if ($actualline =~ /\S*\t\S*\t\d*\t\S*\t\S*\t\S*/){
	    my @columns=split(/\t/);
	    push @comps,$columns[0];
	    push @frames,$columns[5];
	}
    }
    close COMPINFO;

    my $compnumber=$#comps;

# Split these lines and get only compID and targetframe

    my $fastafile   = $_[0];
    my $seqio  = Bio::SeqIO->new(-file => $fastafile, -format => "fasta");
    open(MYOUTFILE, ">TranslatedCompHits.fasta"); #open for write

# Get only the unique sequence * frame combinations
    my %unique;      
    my @uniquecomps;
    my @uniqueframes;
    for (my $a = 0 ; $a < $compnumber+1 ; $a++){
	my $hashkey ="$comps[$a]"."$frames[$a]";
	unless(exists($unique{$hashkey})) {
	    push @uniquecomps,$comps[$a];
	    push @uniqueframes,$frames[$a];
	    $unique{$hashkey} +=1;	  
	}
    }
    my $uniquecompnumber=$#uniquecomps;

    while(my $seqs = $seqio->next_seq) {
	my $id  = $seqs->display_id;
	my $seq = $seqs->seq;
	my $translation;
	for (my $b = 0 ; $b < $uniquecompnumber+1 ; $b++){
	    #print "$uniquecomps[$b]\n";
	    if ($id =~ /$uniquecomps[$b]/){
		#print "ok\n";
		print MYOUTFILE ">$id\n";    
		# translastion for positive readframes
		if ($uniqueframes[$b]>0){
		    my $rnaseq = $seq;
		    $rnaseq =~ tr/T/U/;       # the RNA sequence
		    my $startpoint = (abs $uniqueframes[$b]) - 1;
		    for ( my $i = $startpoint ; $i < length($rnaseq) - 2 ; $i = $i + 3 ) {
			my $codon = substr( $rnaseq, $i, 3 );
			my $amino_acid = $code{$codon};
			$translation=$translation.$amino_acid;
		    }
		    print MYOUTFILE "$translation\n";    
		}
		
		
		# translation for negative read frames
		if ($uniqueframes[$b]<0){
		    my $revseq=reverse($seq);
		    $revseq =~ tr/natgcbdkrvhmyxwsNATGCBDKRVHMYXWS/ntacgvhmybdkrxswNTACGVHMYBDKRXSW/;
		    my $rnaseq = $revseq;
		    $rnaseq =~ tr/T/U/;       # the RNA sequence
		    my $startpoint = (abs $uniqueframes[$b]) - 1;
		    for ( my $i = $startpoint ; $i < length($rnaseq) - 2 ; $i = $i + 3 ) {
			my $codon = substr( $rnaseq, $i, 3 );
			my $amino_acid = $code{$codon};
			$translation=$translation.$amino_acid;
			
		    }
		    print MYOUTFILE "$translation\n";    
		}
		
		
		
		
	    }
	    
	}
	
    }

    close MYOUTFILE;

}

## deduplicate_fasta
sub deduplicate_fasta {
    my %unique;

    my $file   = $_[0];
    my $filenamecode = $file =~/^(.*).fasta/;
    my $seqio  = Bio::SeqIO->new(-file => $file, -format => "fasta");
    my $outseq = Bio::SeqIO->new(-file => ">deduplicated.fasta", -format => "fasta");

    while(my $seqs = $seqio->next_seq) {
	my $id  = $seqs->display_id;
	my $seq = $seqs->seq;
	unless(exists($unique{$seq})) {
	    $outseq->write_seq($seqs);
	    $unique{$seq} +=1;
	}
    }
}

## parse_decypher2
sub parse_decypher2 {
    if ($_[0] eq '' or $_[1] eq '' ) { die "Two input files are needed. First the
output from Parse_DeCypher.pl and second a fasta file with hits in the
C. finmarchicus transcriptome\n"; }

# open the output file from Parse_DeCypherAndConverge.pl

    open(IN, $_[0]) or die "Could not open file $_[0]\n";

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
    open(IN, $_[1]) or die "Could not open the fasta file $_[1]\n";

    my $seq_in = Bio::SeqIO->new(
	-file   => "<$_[1]",
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
}

## peptide_extraction
sub peptide_extraction {

    if ($_[0] eq '') { die "A fasta input file is needed"};

    my $actualfastafile=$_[0];
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

}

## reciprocal_blast
sub reciprocal_blast {

    if ( $_[0] eq '' ) {
	die "fasta file with protein sequence(s) has to be provided as argument to the
script\n"; }

    if ( $_[1] eq '' ) {
	die "The protein database (nr) of taxa related to your target species has to be provided as second argument to the script\n"; }

    my $seqio  = Bio::SeqIO->new(-file => $_[0], -format => "fasta");


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
	my $blastprun = `blastp -query temporarysequence.fasta -db $_[1] -out $filename`;
	print (".") ;
	print MYOUTFILE "$filename\n";  
    }
    print ("\n");

    close MYOUTFILE;
}

## parse_blast
sub parse_blastp {

    if ($_[0] eq '' ) { die "Text file has to be provided that lists
the names of blast report files which shall be treated in this
script\n"; }

    if ($_[1] eq '' ) {die "The protein database (nr) of taxa related to your target species has to be provided as fourth argument to the script\n"; }

    open(IN, $_[0]) or die "Can not open file $_[0]";

# Get an array containing the filenames that shall be parsed here
    my @files;
    while (<IN>){
	chomp;
	push @files,$_;
	
    }
    close IN;


    my @besthits;
    my @bestscores;
    my @bestevalues;
    my @bestdescriptions;
    my @queries;


    foreach (@files){
	my $f=$_;

	open (IN,$f) or die "Can not open file $f";



	
	my @hits;
	my @scores;
	my @evalues;
	my @descriptions;
	my $query;
	
	while (<IN>){
	    chomp;
	    my $line=$_;
	    if ($line =~ /^Query=\s(.*)/){
		my $accession = $1;
		### Test what kind of Fasta definition line format this is:
		# GenBank                           gi|gi-number|gb|accession|locus
		# EMBL Data Library                 gi|gi-number|emb|accession|locus
		# DDBJ, DNA Database of Japan       gi|gi-number|dbj|accession|locus
		if ($accession =~ /^.*gb\|(.*)\|\s*(\S*.*)/){
		    $query=$1;
		}
		elsif ($accession =~ /^.*emb\|(.*)\|\s*(\S*.*)/){
		    $query=$1;
		}
		elsif ($accession =~ /^.*dbj\|(.*)\|\s*(\S*.*)/){
		    $query=$1;
		}
		# NBRF PIR                          pir||entry
		elsif ($accession =~ /^.*pir\|\|(\S*)\s*(\S*.*)/){
		    $query=$1;
		}
		# Protein Research Foundation       prf||name
		elsif ($accession=~ /^.*prf\|.*\|(.*)\s*(\S*.*)/){
		    $query=$1;
		}
		# SWISS-PROT                        sp|accession|name
		elsif ($accession=~ /^.*sp\|(.*)\|\s*(\S*.*)/){
		    $query=$1;
		}
		# Brookhaven Protein Data Bank (1)  pdb|entry|chain
		elsif ($accession=~ /^.*pdb\|(.*)\|\s*(\S*.*)/){
		    $query=$1;
		}
		# GenInfo Backbone Id               bbs|number
		elsif ($accession=~ /^.*bbs\|(.*)/){
		    $query=$1;
		}
		# General database identifier       gnl|database|identifier
		elsif ($accession=~ /^.*gnl\|(.*)\|\s*(\S*.*)/){
		    $query=$1;
		}
		# NCBI Reference Sequence           ref|accession|locus
		elsif ($accession=~ /^.*ref\|(.*)\|\s*(\S*.*)/){
		    $query=$1;
		}
		# tpg|DAA34093.1| 
		elsif ($accession=~ /^.*tpg\|(.*)\|\s*(\S*.*)/){
		    $query=$1;
		}
		# Local Sequence identifier         lcl|identifier
		elsif ($accession=~ /^.*lcl\|(.*)\|\s*(\S*.*)/){
		    $query=$1;
		}
		elsif ($accession=~ /^.*lcl\|(.*)\s*\S*/){
		    $query=$1;
		}
		elsif ($accession){
		    $query="";
		}
	    }
	    
	    
	    if ($line =~ /^>/){
		### Test what kind of Fasta definition line format this is:
		# GenBank                           gi|gi-number|gb|accession|locus
		# EMBL Data Library                 gi|gi-number|emb|accession|locus
		# DDBJ, DNA Database of Japan       gi|gi-number|dbj|accession|locus
		if ($line =~ /^>.*gb\|(.*)\|\s*(\S*.*)/){
		    push @hits,$1;
		    push @descriptions,$2;
		}
		if ($line =~ /^>.*emb\|(.*)\|\s*(\S*.*)/){
		    push @hits,$1;
		    push @descriptions,$2;
		}
		if ($line =~ /^>.*dbj\|(.*)\|\s*(\S*.*)/){
		    push @hits,$1;
		    push @descriptions,$2;
		}
		# NBRF PIR                          pir||entry
		if ($line =~ /^>.*pir\|\|(\S*)\s*(\S*.*)/){
		    push @hits,$1;
		    push @descriptions,$2;
		}
		# Protein Research Foundation       prf||name
		if ($line=~ /^>.*prf\|.*\|(.*)\s*(\S*.*)/){
		    push @hits,$1;
		    push @descriptions,$2;
		}
		# SWISS-PROT                        sp|accession|name
		if ($line=~ /^>.*sp\|(.*)\|\s*(\S*.*)/){
		    push @hits,$1;
		    push @descriptions,$2;
		}
		# Brookhaven Protein Data Bank (1)  pdb|entry|chain
		if ($line=~ /^>.*pdb\|(.*)\|\s*(\S*.*)/){
		    push @hits,$1;
		    push @descriptions,$2;
		}
		# GenInfo Backbone Id               bbs|number
		if ($line=~ /^>.*bbs\|(.*)/){
		    push @hits,$1;
		    push @descriptions,$1;
		}
		# General database identifier       gnl|database|identifier
		if ($line=~ /^>.*gnl\|(.*)\|\s*(\S*.*)/){
		    push @hits,$1;
		    push @descriptions,$2;
		}
		# NCBI Reference Sequence           ref|accession|locus
		if ($line=~ /^>.*ref\|(.*)\|\s*(\S*.*)/){
		    push @hits,$1;
		    push @descriptions,$2;
		}
		# tpg|DAA34093.1| 
		if ($line=~ /^>.*tpg\|(.*)\|\s*(\S*.*)/){
		    push @hits,$1;
		    push @descriptions,$2;
		}
		# Local Sequence identifier         lcl|identifier
		if ($line=~ /^>.*lcl\|(.*)\|\s*(\S*.*)/){
		    push @hits,$1;
		    push @descriptions,$2;
		}
		if ($line=~ /^>.*lcl\|(.*)\s*\S*/){
		    push @hits,$1;
		    push @descriptions,$2;
		}
	    }
	    if ($line =~ /Score\s*=\s*(.*)\s*bits.*Expect\s*=\s*(.*),.*$/g){
		push @scores,$1;
		push @evalues,$2;
	    }
	}
	close IN;
### Now I have to exclude those hits that are putative or hypthetical

# I remove everything that contains "predict", "hypothetic",
# "putativ", "provisional". If I use //i, this will match the string
# independent of wether it is upper- or lowercase.

# the word partial was often only indicated when I retrieve the
# sequence from ncbi (in the web version of this program), but here
# the sequences were directly from the local database so the
# identifiers in the blastpresult.out file is all I got and trust
# on. 
	print "... finding best hit for query $query\n";
	my @efetch_results=@descriptions;
	# Take only non-hypothetical proteins
	my @index = grep { $efetch_results[$_] !~ /predict|hypothetic|putativ|provisional|partial|unknown/i } 0..$#efetch_results;
	if ($#index<0) # if no best hit is found then take a hit that contains "predict, hypothetic..."
	{
	    @index = grep { $efetch_results[$_] =~ /predict|hypothetic|putativ|provisional|partial|unknown/i } 0..$#efetch_results;
	}
	

	# remove the hits that are putative or hypothetical
	my @hits_selected = @hits[@index];
	my @scores_selected = @scores[@index];
	my @descriptions_selected = @descriptions[@index];
	my @evalues_selected = @evalues[@index];
	push @queries,$query;
	push @besthits,$hits_selected[0];
	push @bestscores,$scores_selected[0];
	push @bestdescriptions,$descriptions_selected[0];
	push @bestevalues,$evalues_selected[0];
    }

    open (MYOUTFILE,">blastp\.besthits");
    print MYOUTFILE "query\tbesthit\tdescription\tscore\tevalue\n";
    my $n;
    my $querynumber = $#queries;
    for ($n=0;$n<$querynumber+1;$n++){
	print MYOUTFILE "$queries[$n]\t$besthits[$n]\t$bestdescriptions[$n]\t$bestscores[$n]\t$bestevalues[$n]\n";
    }
    close MYOUTFILE;

    open (MYOUTFILE2,">besthits\.fasta");
############# Now get the sequences of the best hits and write them to a fastafile


    foreach (@besthits){

	my $db     = "protein";
	my $query  = $_;
	print "$query\n";
	my $blastdbcmd = `blastdbcmd -db $_[1] -entry $query >> besthits.fasta`;

    }
    close MYOUTFILE2;
}

## hmmscan
sub hmmscan {



    if ( $_[0] eq '' ) {
	die "Two fasta files with protein sequence(s) have to be provided as argument to the
script. One contains the contigs of your target species the other contains the sequences of the best blastp hits\n"; }

    if ( $_[1] eq '' ) {
	die "Two fasta files with protein sequence(s) have to be provided as argument to the
script. One contains the contigs of your target species the other contains the sequences of the best blastp hits\n"; }

    if ( $_[2] eq '' ) {
	die "Provide the Pfam path (including the file) as third argument\n",}


    open my $fh, '<', $_[0] or die "error opening $_[0]: $!";
    my $data1 = do { local $/; <$fh> };

    open my $fh2, '<', $_[1] or die "error opening $_[1]: $!";
    my $data2 = do { local $/; <$fh2> };

    open (MYOUTFILE, ">alldata.fasta");
    my $alldata=$data1.$data2;
    print MYOUTFILE $alldata;
    close MYOUTFILE;

# running hmmscan from the command line
    my $output = `hmmscan --domtblout hmmscanoutput.tab $_[2] alldata.fasta`;
    open (MYOUTFILE2, ">hmmscan.out");
    print MYOUTFILE2 "queryname\tlen\tdomname\tbegin\tend\tdomaccession\tdomlength\tquerylength\tEval\tscore\tiEval\tdomscore\n";


    open(IN, 'hmmscanoutput.tab') or die "could not open file hmmscanoutput.tab\n";
# Since one domain can be reported more than once, I will select only the unique ones
    my @uniquelines;
    while (<IN>) {
	unless (/^\#/) {    # avoid all lines beginning
	    # with the '#' character
	    my @columns  = split(/ +/);
	    my $domname  = $columns[0];
	    my $queryname = $columns[3];
#	$queryname1 =~ /.*\|.*$/g;
	    #$queryname    =~ s/.*\|//;
	    my $len      = $columns[5];
	    my $begin    = $columns[19]; # This is the domain envelope that is also reported in SMART.
	    my $end      = $columns[20];
	    my $domaccession      = $columns[1];
	    my $domlength      = $columns[2];
	    my $querylength      = $columns[5];
	    my $Eval      = $columns[6];
	    my $score      = $columns[7];
	    my $iEval     = $columns[12];
	    my $domscore    = $columns[13];
	    
	    if ( $iEval < 1e-10 ) {
		push @uniquelines, "$queryname\t$len\t$domname\t$begin\t$end\t$domaccession\t$domlength\t$querylength\t$Eval\t$score\t$iEval\t$domscore\n";
	    }
	}
    }
    close IN;

    my @unique = uniq @uniquelines;
    foreach ( @unique ) {
	print MYOUTFILE2 $_;
    }

    close MYOUTFILE2;
}

## mafft
sub mafft {

    if ( $_[0] eq '' ) {
	die "Two fasta files with protein sequence(s) have to be provided as argument to the
script. The first contains the contigs of your target species and the second contains the sequences of the best blastp hits\n"; }

    if ( $_[1] eq '' ) {
	die "Two fasta files with protein sequence(s) have to be provided as argument to the
script. The first contains the contigs of your target species and the second contains the sequences of the best blastp hits\n"; }

    if ( $_[2] eq '' ) {
	die "Provide the file ending with .besthits as third argument\n",}


# First loading in the information which fasta sequences form a pair 
    open(IN, $_[2]) or die "could not open file $_[2]\n";

    my @queries;
    my @besthits;
    while (<IN>) {
	chomp;
	if ($.!=1){# skip the header line
	    my @row = split(/\t/);
	    push @queries,$row[0];
	    push @besthits,$row[1];
	}
    }
    close IN;

#### Extract the protein sequences for the query-besthit pairs and
#### submit them for alignment to mafft

    open(OUT, ">mafftalignments.outfilenames"); # this will list the outputfilenames in one file.

# Creating a temporary output file that stores the sequence pairs in fasta format
    my $t;
    my $tlength=@queries;
    for ($t=0;$t<$tlength;$t++){
	open (TEMPORARY, ">temporary_fastapair.fasta");
	print "$t\n";
## Getting the query sequence
	my $query_in = Bio::SeqIO->new(
	    -file   => "<$_[0]",
	    -format => 'Fasta',
	    );
	
	my $actualquery;
	while (my $seq = $query_in->next_seq) {
	    my $fastaid=$seq->display_id;
	    if ($fastaid =~ /.*$queries[$t].*/){
		$actualquery = $seq->seq();
	    }
	}
	
	
	print TEMPORARY ">$queries[$t]\n";
	print TEMPORARY "$actualquery\n\n";
	
## getting the besthit sequence
	my $besthit_in = Bio::SeqIO->new(
	    -file   => "<$_[1]",
	    -format => 'Fasta',
	    );
	
	my $actualbesthit;
	while (my $seq = $besthit_in->next_seq) {
	    my $fastaid=$seq->display_id;
	    if ($fastaid =~ /.*$besthits[$t].*/){
		$actualbesthit = $seq->seq();
	    }
	}
	print TEMPORARY ">$besthits[$t]\n";
	print TEMPORARY "$actualbesthit\n\n";

	close TEMPORARY;

	my $fastapairname1="$queries[$t]\_\_$besthits[$t]\_\_MAFFTalignmentClustalW\.txt";
	print OUT "$fastapairname1\n";
	my $maffrun1 = `mafft --clustalout temporary_fastapair.fasta  > $fastapairname1`;

	my $fastapairname2="$queries[$t]\_\_$besthits[$t]\_\_MAFFTalignment\.fasta";
	my $maffrun2 = `mafft temporary_fastapair.fasta  > $fastapairname2`;
	
	


    }
    close OUT;
}

## parse_mafft
sub parse_mafft {


# check if argument to the script is there.  
    if ( $_[0] eq '' ) {
	die "A file listing the names of the ClustalW alignment files has to be provided as argument to the
script. The filename must be of the pattern queryname__besthitname__MAFFTalignment.fasta\n"; }


# Load in the alignment file names
    open(INFILES, "$_[0]") or die "Could not open file $_[0]\n";
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

}

## vetting
sub vetting {


    if ( $_[0] eq '' ) {
	die "fasta file with protein sequence(s) has to be provided as first argument to the
script\n"; }

    if ( $_[1] eq '' ) {
	die "fasta file with EST sequence(s) of your target species has to be provided as second argument to the
script\n"; }

#my $prog = 'tblastn';
#my $db   = 'est';
#my $e_val= '1e-10';

#$Bio::Tools::Run::RemoteBlast::HEADER{'ENTREZ_QUERY'} = 'Calanus finmarchicus [ORGN]';
#$Bio::Tools::Run::RemoteBlast::RETRIEVALHEADER{'FORMAT_TYPE'} = 'Text';

    my $seqio = Bio::SeqIO->new(-file=>$_[0] , '-format' => 'fasta' );

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
	my $blastprun = `tblastn -query temporarysequence.fasta -db $_[1] -out $filename`;
	print (".");
	print MYOUTFILE "$filename\n";
    }
    print ("\n");

    close MYOUTFILE;
    
}

## parse_vetting
sub parse_vetting {

    if ($_[0] eq '' ) { die "text file has to be provided that lists
the names of tblastn report files which shall be treated in this
script\n"; }

    open(IN, $_[0]) or die "Can not open file $_[0]";
    open(OUT, ">Vetting.out");

# Get an array containint the filenames that shall be parsed here
    my @files;
    while (<IN>){
	chomp;
	push @files,$_;
	
    }
    close IN;


    print OUT ("query\tquerylength\test\tlength\tscore\tidentity\tevalue\tqueryalignmentstart\tdescription\tframe\n");

    foreach (@files){
	my $f=$_;

	open (IN,$f) or die "Can not open file $f";

	my $query;
	my @querylengths;
	my @esthits;
	my @estscores;
	my @estidentities;
	my @estevalues;
	my @estdescriptions;
	my @estlengths;
	my @queryalignmentstart;
	my @frame;

	my $estnumber=0;
	my $querynumber=1;
	while (<IN>){
	    chomp;
	    my $line=$_;
	    if ($line =~ /^Query=\s(.*)/){
		my $accession = $1;
		### Test what kind of Fasta definition line format this is:
		# GenBank                           gi|gi-number|gb|accession|locus
		# EMBL Data Library                 gi|gi-number|emb|accession|locus
		# DDBJ, DNA Database of Japan       gi|gi-number|dbj|accession|locus
		if ($accession =~ /^.*gb\|(.*)\|\s*(\S*.*)/){
		    $query=$1;
		}
		elsif ($accession =~ /^.*emb\|(.*)\|\s*(\S*.*)/){
		    $query=$1;
		}
		elsif ($accession =~ /^.*dbj\|(.*)\|\s*(\S*.*)/){
		    $query=$1;
		}
		# NBRF PIR                          pir||entry
		elsif ($accession =~ /^.*pir\|\|(\S*)\s*(\S*.*)/){
		    $query=$1;
		}
		# Protein Research Foundation       prf||name
		elsif ($accession=~ /^.*prf\|.*\|(.*)\s*(\S*.*)/){
		    $query=$1;
		}
		# SWISS-PROT                        sp|accession|name
		elsif ($accession=~ /^.*sp\|(.*)\|\s*(\S*.*)/){
		    $query=$1;
		}
		# Brookhaven Protein Data Bank (1)  pdb|entry|chain
		elsif ($accession=~ /^.*pdb\|(.*)\|\s*(\S*.*)/){
		    $query=$1;
		}
		# GenInfo Backbone Id               bbs|number
		elsif ($accession=~ /^.*bbs\|(.*)/){
		    $query=$1;
		}
		# General database identifier       gnl|database|identifier
		elsif ($accession=~ /^.*gnl\|(.*)\|\s*(\S*.*)/){
		    $query=$1;
		}
		# NCBI Reference Sequence           ref|accession|locus
		elsif ($accession=~ /^.*ref\|(.*)\|\s*(\S*.*)/){
		    $query=$1;
		}
		# tpg|DAA34093.1| 
		elsif ($accession=~ /^.*tpg\|(.*)\|\s*(\S*.*)/){
		    $query=$1;
		}
		# Local Sequence identifier         lcl|identifier
		elsif ($accession=~ /^.*lcl\|(.*)\|\s*(\S*.*)/){
		    $query=$1;
		}
		elsif ($accession=~ /^.*lcl\|(.*)\s*\S*/){
		    $query=$1;
		}
		elsif ($accession){
		    $query="";
		}
	    }
	    if ($line =~ /^>(.*)/){
		### Test what kind of Fasta definition line format this is:
		# GenBank                           gi|gi-number|gb|accession|locus
		# EMBL Data Library                 gi|gi-number|emb|accession|locus
		# DDBJ, DNA Database of Japan       gi|gi-number|dbj|accession|locus
		my $line2=$1;
		
		if ($line2 =~ /^.*gb\|(.*)\|\s*(\S*.*)/){
		    push @esthits, $1 ;
		    push @estdescriptions,$2;
		}
		elsif ($line2 =~ /^.*emb\|(.*)\|\s*(\S*.*)/){
		    push @esthits,$1;
		    push @estdescriptions,$2;
		}
		elsif ($line2 =~ /^.*dbj\|(.*)\|\s*(\S*.*)/){
		    push @esthits,$1;
		    push @estdescriptions,$2;
		}
		# NBRF PIR                          pir||entry
		elsif ($line2 =~ /^.*pir\|\|(\S*)\s*(\S*.*)/){
		    push @esthits,$1;
		    push @estdescriptions,$2;
		}
		# Protein Research Foundation       prf||name
		elsif ($line2=~ /^.*prf\|.*\|(.*)\s*(\S*.*)/){
		    push @esthits,$1;
		    push @estdescriptions,$2;
		}
		# SWISS-PROT                        sp|accession|name
		elsif ($line2=~ /^.*sp\|(.*)\|\s*(\S*.*)/){
		    push @esthits,$1;
		    push @estdescriptions,$2;
		}
		# Brookhaven Protein Data Bank (1)  pdb|entry|chain
		elsif ($line2=~ /^.*pdb\|(.*)\|\s*(\S*.*)/){
		    push @esthits,$1;
		    push @estdescriptions,$2;
		}
		# GenInfo Backbone Id               bbs|number
		elsif ($line2=~ /^.*bbs\|(.*)/){
		    push @esthits,$1;
		    push @estdescriptions,$1;
		}
		# General database identifier       gnl|database|identifier
		elsif ($line2=~ /^.*gnl\|(.*)\|\s*(\S*.*)/){
		    push @esthits,$1;
		    push @estdescriptions,$2;
		}
		# NCBI Reference Sequence           ref|accession|locus
		elsif ($line2=~ /^.*ref\|(.*)\|\s*(\S*.*)/){
		    push @esthits,$1;
		    push @estdescriptions,$2;
		}
		# tpg|DAA34093.1| 
		elsif ($line2=~ /^.*tpg\|(.*)\|\s*(\S*.*)/){
		    push @esthits,$1;
		    push @estdescriptions,$2;
		}
		# Local Sequence identifier         lcl|identifier
		elsif ($line2=~ /^.*lcl\|(.*)\|\s*(\S*.*)/){
		    push @esthits,$1;
		    push @estdescriptions,$2;
		}
		elsif ($line2=~ /^>.*lcl\|(.*)\s*\S*/){
		    push @esthits,$1;
		    push @estdescriptions,$2;
		}
		elsif ($line2){
		    push @esthits,"";
		    push @estdescriptions,"";
		}
		
		$estnumber++;
		$querynumber=$estnumber;
	    }	
	    if ($line =~ /Score.*\s*=\s*(\S*)\s*bits.*Expect.*\s*=\s*(\S*)\,/){
		my $a=$1;
		my $b=$2;
		unless  ($line =~ /.*mRNA.*/){
		    push @estscores,$a;
		    push @estevalues,$b;
		}

	    }
	    if ($line =~ /Identities\s*=.*\((.*)\).*\(.*\).*\(.*\)$/){
		push @estidentities,$1;
	    }
	    if ($line =~ /Length=(.*)$/g){
		push @estlengths,$1; #the first hit here will be the querylength the rest will be the est lengths
	    }
	    if ($line =~ /Frame = (.*)$/g){
		push @frame,$1; #the first hit here will be the reading frame.
	    }
	    if ($line =~ /^Query\s*(\d*)\s*\D*\s*\d*/g){
		my $as=$1;
		if ($querynumber == $estnumber){
		    push @queryalignmentstart,$as; #this collects the start of
		    #the alignment of the query
		    #file to identify if the hit
		    #aligns at the N- or
		    #c-terminus or in between
		    $querynumber++
		}
	    }
	}
	my $querylength=$estlengths[0];
	my @estlengths2;
	if ($#estlengths>=1){
	    @estlengths2=@estlengths[1..$#estlengths];
	}
	else {
	    @estlengths2=0;
	}
	close IN;
	# Print one line for each est hit
	my $lines=@esthits;
	for (my $l=0;$l<$lines;$l++){

	    print OUT "$query\t$querylength\t$esthits[$l]\t$estlengths2[$l]\t$estscores[$l]\t$estidentities[$l]\t$estevalues[$l]\t$queryalignmentstart[$l]\t$estdescriptions[$l]\t$frame[$l]\n";
	}
    }

    close OUT;

}


## create_resport
sub create_report {

    open (MYOUTFILE,">ResultReport.tex");
#if ( $_[0] eq '' ) {
#    die "A fasta file from DeCypher has to be provided as first input\n";}

#if ( $_[1] eq '' ) {
#    die "A text file from DeCypher hast to be provided as second input\n";}

#if ( $_[2] eq '' ) {
#    die "The filepath of the Pfam-A database has to be provided as third input\n";}
    print MYOUTFILE "\\documentclass\[a4paper\]\{article\}\n";
    print MYOUTFILE "\\usepackage\{graphicx\}\n";
    print MYOUTFILE "\\usepackage\{hyperref\}\n";
    print MYOUTFILE "\\hypersetup\{\n";
    print MYOUTFILE "colorlinks   = true,\n";
    print MYOUTFILE "urlcolor     = blue,\n";
    print MYOUTFILE "linkcolor    = blue,\n";
    print MYOUTFILE "citecolor   = red\n";
    print MYOUTFILE "\}\n";
    print MYOUTFILE "\\usepackage\[top=2cm, bottom=2cm, left=2.5cm, right=2.5cm\]\{geometry\}\n";
    print MYOUTFILE "\\usepackage\{longtable\}\n";

    print MYOUTFILE "\\begin\{document\}\n";
    print MYOUTFILE "\\title\{Report\}\n";
    print MYOUTFILE "\\author\{Protein prediction\}\n";
    print MYOUTFILE "\\maketitle\n";
    print MYOUTFILE "\\section\{Queries\}\n";
    print MYOUTFILE "Accession numbers and descriptions of the queries used\:\\\\\\\\\n";
    open(QUERYINFO,"Cfin.queryinfo");

    while (<QUERYINFO>){
	chomp;
	my $line=$_;
	if ($line =~ /^Query\:\s(.*)\t(.*)$/){
	    my $query=$1;
	    my $querydescription=$2;
	    $querydescription =~ s/\_/\\_/g;
	    $querydescription =~ s/\=/\\=/g;
	    my $searchquery=$query;
	    $query =~ s/\_/\\_/g;
	    print MYOUTFILE "\\href\{http://www.ncbi.nlm.nih.gov/protein/$searchquery\}{$query} - $querydescription\\\\\n";
	}
    }
    close QUERYINFO;

    print MYOUTFILE "\\section\{Unique contig hits in the transcriptome of your target species}\n";
    print MYOUTFILE "The following sections will summarize the results for each of these contigs\:\\\\\\\\\n";
    open(COMPINFO,"Cfin.compinfo");
    my @comphits;
    my @comphitslatex;
    while (<COMPINFO>){
	chomp;
	my $line=$_;
	if ($line =~ /^hit\:\s(.*)$/){
	    my $hit=$1;
	    push @comphits,$hit;
	    $hit =~ s/\_/\\_/g;
	    push @comphitslatex,$hit;
	    print MYOUTFILE "$hit\\\\\n";
	}
    }
    close COMPINFO;

    print MYOUTFILE "\\\\The fasta file \\href\{UniqueComphits.fasta\}\{UniqueComphits.fasta\} contains the polypeptide sequences of these contigs.\\\\\n";        
    print MYOUTFILE "The fasta file \\href\{LongestPolypeptide.fasta\}\{LongestPolypeptide.fasta\} contains the polypeptide sequences of the longest open reading frames for these contigs.\\\\\n";        

    foreach (@comphitslatex){
	my $actualcomphit=$_;
	my $comphit= $actualcomphit;
	$comphit =~ s/\\_/\_/g;
	print MYOUTFILE "\\clearpage\n";        
	print MYOUTFILE "\\section\{$actualcomphit\}\n";

############### Include a table on the information to which querys
############### this comphit showed up as a hit and in which quality
############### print MYOUTFILE "\\begin\{center\}\n";
	print MYOUTFILE "\\textbf\{Contig $actualcomphit showed up as a hit to the following queries (sorted by e-values). The targetframe refers to the  contig of your target species.\}\\\\\\\\\n";
	print MYOUTFILE "\\begin\{longtable\}\{c c c r r\}\n";
	print MYOUTFILE "\\hline\n";
	print MYOUTFILE "Query \& Rank \& E\-value \& Queryframe \& Targetframe\\\\\n";
	print MYOUTFILE "\\hline\n";
	print MYOUTFILE "\\endfirsthead\n";
	print MYOUTFILE "\\multicolumn\{5\}\{c\}\n";
	print MYOUTFILE "\{\\tablename\ -- \\textit\{Continued from previous page\}\}\\\\\n";
	print MYOUTFILE "\\hline\n";
	print MYOUTFILE "Query \& Rank \& E\-value \& Queryframe \& Targetframe\\\\\n";
	print MYOUTFILE "\\hline\n";
	print MYOUTFILE "\\endhead\n";
	print MYOUTFILE "\\multicolumn\{5\}\{r\}\n";
	print MYOUTFILE "\{\\textit\{Continued on next page\}\} \\\\\n";
	print MYOUTFILE "\\endfoot\n";
	print MYOUTFILE "\\hline\n";
	print MYOUTFILE "\\endlastfoot\n";

	open(COMPINFO2,"Cfin.compinfo");
	my @queries;
	my @ranks;
	my @evalues;
	my @queryframes;
	my @targetframes;
	while (<COMPINFO2>){
	    chomp;
	    my $line=$_;
	    if ($line =~ /^$comphit\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)$/){
		push @queries,$1;
		push @ranks,$2;
		push @evalues,$3;
		push @queryframes,$4;
		push @targetframes,$5;
	    }
	}
	close COMPINFO2;
	# rankig
	my @list_order = sort { $evalues[$a] <=> $evalues[$b] } 0 .. $#evalues;
	@queries=@queries[@list_order];
	@ranks=@ranks[@list_order];
	@evalues=@evalues[@list_order];
	@queryframes=@queryframes[@list_order];
	@targetframes=@targetframes[@list_order];
	my $e;
	my $elength=@evalues;
	for ($e=0;$e<$elength;$e++){
	    my $q=$queries[$e];
	    $q =~ s/\_/\\_/g;
	    print MYOUTFILE "\\href\{http://www.ncbi.nlm.nih.gov/protein/$queries[$e]\}{$q} \& $ranks[$e] \& $evalues[$e] \& $queryframes[$e] \& $targetframes[$e] \\\\\n";
	}
	print MYOUTFILE "\\end\{longtable\}\\\\\n";
#    print MYOUTFILE "\\end\{center\}\n";    



################# proteins on NCBI
	
	print MYOUTFILE "\\\\ \\\\ \\textbf\{Best blastp hit in the protein database:\}\\\\\\\\\n";
	my $besthits;
	open(BESTHIT,"blastp.besthits");
	while (<BESTHIT>){
	    chomp;
	    my $line=$_;
	    if ($line =~ /^$comphit\t(.*)\t(.*)\t(.*)\t(.*)/g){
		my $besthit=$1;
		$besthits=$besthit;
		my $desc=$2;	    
		my $s=$3;
		my $ev=$4;
		if ($besthit =~ /\_/){
		    $besthit =~ s/\_/\\_/;
		}
		if ($desc =~ /\_/){
		    $desc =~ s/\_/\\_/;
		}
		if ($besthits){
		    print MYOUTFILE "Accession: \\href\{http://www.ncbi.nlm.nih.gov/protein/$besthits\}\{$besthit}\\\\\n";
		    print MYOUTFILE "Description: $desc\\\\\n";
		    print MYOUTFILE "Score: $s\\\\\n";
		    print MYOUTFILE "E-value: $ev\\\\\n";
		}
		else {
		    print MYOUTFILE "No best hit found !\\\\\n";
		}
	    }
	}
	close BESTHIT;
	



############## Include information on the identity and similarity from the mafft alignment
	if ($besthits){    
	    open(MAFFTOUT,"parsemafft.out");
	    while (<MAFFTOUT>){
		chomp;
		my $line=$_;
		if ($line =~ /^$comphit\t.*\t(.*)\t(.*)/g){
		    my $identity=sprintf "%.2f", $1;
		    my $similarity=sprintf "%.2f", $2;
		    print MYOUTFILE "\\% Identity (based on alignment): $identity\\\\\n";
		    print MYOUTFILE "\\% Similarity (based on alignment): $similarity\\\\\n";
		}
	    }
	    close MAFFTOUT;
	    
########### Make a link to the clustalw alignment file
	    
	    my $printcomphit=$comphit;
	    $printcomphit =~ s/\_/\\_/g;
	    my $printbesthit=$besthits;
	    $printbesthit =~ s/\_/\\_/g;
	    print MYOUTFILE "Link to alignment file:\\\\\n";
	    print MYOUTFILE "\\href\{$comphit\_\_$besthits\_\_MAFFTalignmentClustalW\.txt\}\{$printcomphit\\_\\_$printbesthit\\_\\_MAFFTalignmentClustalW\.txt\}\\\\\n";
	    

############# Include information on the protein domains in the Cfin contig and the best hit.
	    print MYOUTFILE "\\\\\\textbf\{Domain predictions. Start and End refer to the non-aligned protein-sequences. Slength: length of the query sequence; Dlength: Length of the domain.\}\\\\\\\\\n";    
	    
	    print MYOUTFILE "\\begin\{longtable\}\{l r l l r r r l\}\n";
	    print MYOUTFILE "\\hline\n";
	    print MYOUTFILE "Sequence \& Slength \& Domain \& Accession \& Dlength \& Start \& End \& i\-E\-value \\\\\n";
	    print MYOUTFILE "\\hline\n";
	    print MYOUTFILE "\\endfirsthead\n";
	    print MYOUTFILE "\\multicolumn\{8\}\{c\}\n";
	    print MYOUTFILE "\{\\tablename\ -- \\textit\{Continued from previous page\}\}\\\\\n";
	    print MYOUTFILE "\\hline\n";
	    print MYOUTFILE "Sequence \& Slength \& Domain \& Accession \& Dlength \& Start \& End \& i\-E\-value \\\\\n";
	    print MYOUTFILE "\\hline\n";
	    print MYOUTFILE "\\endhead\n";
	    print MYOUTFILE "\\multicolumn\{8\}\{r\}\n";
	    print MYOUTFILE "\{\\textit\{Continued on next page\}\} \\\\\n";
	    print MYOUTFILE "\\endfoot\n";
	    print MYOUTFILE "\\hline\n";
	    print MYOUTFILE "\\endlastfoot\n";

	    
	    open(PROTEINDOMAINS,"hmmscan.out");
	    while (<PROTEINDOMAINS>){
		chomp;
		my $line=$_;
		if ($line =~ /$comphit/g){
		    $line =~ /.*\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)/g;
		    my $seqlength=$1;
		    my $domain=$2;
		    my $accession=$5;
		    my $domlength=$6;
		    my $start=$3;
		    my $end=$4;
		    my $evalue=$10;
		    $domain =~ s/\_/\\_/g;
		    print MYOUTFILE "$actualcomphit \& $seqlength \& $domain \& \\href\{http://pfam.xfam.org/family/$accession\}\{$accession\} \& $domlength \& $start \& $end \& $evalue \\\\\n";
		}
		if ($line =~ /$besthits/g){
		    $line =~ /.*\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)/g;
		    my $seqlength=$1;
		    my $domain=$2;
		    my $accession=$5;
		    my $domlength=$6;
		    my $start=$3;
		    my $end=$4;
		    my $evalue=$10;
		    $domain =~ s/\_/\\_/g;
		    print MYOUTFILE "$printbesthit \& $seqlength \& $domain \& \\href\{http://pfam.xfam.org/family/$accession\}\{$accession\} \& $domlength \& $start \& $end \& $evalue \\\\\n";
		}
	    }
	    close PROTEINDOMAINS;
	    print MYOUTFILE "\\hline\n";
	    print MYOUTFILE "\\end\{longtable\}\\\\\n";


############### Include the figure of alignment and protein domain
	    print MYOUTFILE "\\begin\{figure\}\[h\]\n";
	    print MYOUTFILE "\\centering\n";
	    print MYOUTFILE "\\includegraphics\[width\=\.6\\textwidth\]\{\{$comphit\_Pfam\}.png\}\n";
	    print MYOUTFILE "\\caption\{Protein domains predicted along the protein sequence of your target species derived from $actualcomphit and its aligned best blastp hit $printbesthit.\}\n";
	    print MYOUTFILE "\\end\{figure\}\n";
	}
	
	
############## Include information on the EST hits in the Cfin database
	print MYOUTFILE "\\\\\\textbf\{EST hits\:\}\\\\\\\\\n";    
	
	print MYOUTFILE "\\begin\{longtable\}\{l r r r l r r\}\n";
	print MYOUTFILE "\\hline\n";
	print MYOUTFILE "EST \& Length \& Score \& Identity \& E\-value \& Alignmentstart \& Frame\\\\\n";
	print MYOUTFILE "\\hline\n";
	print MYOUTFILE "\\endfirsthead\n";
	print MYOUTFILE "\\multicolumn\{7\}\{c\}\n";
	print MYOUTFILE "\{\\tablename\ -- \\textit\{Continued from previous page\}\}\\\\\n";
	print MYOUTFILE "\\hline\n";
	print MYOUTFILE "EST \& Length \& Score \& Identity \& E\-value \& Alignmentstart \& Frame\\\\\n";
	print MYOUTFILE "\\hline\n";
	print MYOUTFILE "\\endhead\n";
	print MYOUTFILE "\\multicolumn\{7\}\{r\}\n";
	print MYOUTFILE "\{\\textit\{Continued on next page\}\} \\\\\n";
	print MYOUTFILE "\\endfoot\n";
	print MYOUTFILE "\\hline\n";
	print MYOUTFILE "\\endlastfoot\n";
	
	
	
	open(VETTINGOUT,"Vetting.out");
	while (<VETTINGOUT>){
	    chomp;
	    my $line=$_;
	    if ($line =~ /^$comphit\t.*\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t(.*)\t.*\t(.*)/g){
		my $est=$1;
		my $searchest=$est;
		my $length=$2;
		my $score=$3;
		my $identity=$4;
		my $evalue=$5;
		my $start=$6;
		my $readingframe=$7;
		if ($est =~ /\_/){
		    $est =~ s/\_/\\_/;
		}
		$identity =~ s/%/\\%/;
		print MYOUTFILE "\\href\{http://www.ncbi.nlm.nih.gov/nucest/$searchest\}\{$est\} \& $length \& $score \& $identity \& $evalue \& $start \& $readingframe \\\\\n";
	    }
	}
	close VETTINGOUT;
	print MYOUTFILE "\\hline\n";
	print MYOUTFILE "\\end\{longtable\}\\\\\n";
	

    }



    print MYOUTFILE "\\end\{document\}\n";


#\section{Conclusion}
#Write your conclusion here.




    close MYOUTFILE;

}


### Function calls
sub ProteinScreen {
# argument 0: A fasta file with query sequences
# argument 1: The filepath of the Pfam-A database
# argument 2: The fasta file of the local transcriptome database
# argument 3: The protein database (nr) of taxa related to your target species
# argument 4: fasta file with EST sequence(s) of your target species
    local_database_search($_[0], $_[2]);
    parse_local_database_search("tblastn.LocalDatabase");
    best_hits_to_fasta("Cfin.compinfo");
    translation("Comphits.fasta", "Cfin.compinfo");
    deduplicate_fasta("TranslatedCompHits.fasta");
    parse_decypher2("Cfin.compinfo", "deduplicated.fasta");
    peptide_extraction("UniqueCompHits.fasta");
    reciprocal_blast("LongestPolypeptide.fasta", $_[3]);
    parse_blastp("Blastp.outfiles");
    hmmscan("LongestPolypeptide.fasta", "besthits.fasta", $_[1]);
    mafft("LongestPolypeptide.fasta", "besthits.fasta", "blastp.besthits");
    parse_mafft("mafftalignments.outfilenames");
    vetting("LongestPolypeptide.fasta", $_[4]);
    parse_vetting("tblastn.ESTout");
#my $retVal12							       =  `R CMD BATCH Pfam.r`;
    create_report();
#my $retVal14							       =  `pdflatex --interaction=nonstopmode ResultReport.tex`;
#my $retVal15							       =  `pdflatex --interaction=nonstopmode ResultReport.tex`;
}

#ProteinScreen::ProteinScreen("/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/Fishes_reproduction_embl.fasta", "/home/alj/Dropbox.personal/Dropbox/Pfam/Pfam-A.hmm", "/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/nonredundant_dbv1.fasta", "/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/nr_bonyfishes", "/home/alj/Dropbox.personal/Dropbox/Programming/2014CopepodHSPs/201408ProgramAdjustable/TestFiles/GuppyESTs.fasta");

1;

# Need to fix the uninitialized value warning?!

#XXX How to include Pfam.r execution??


# Running R script from within perl: http://www.perlmonks.org/?node_id =  1009021
# XX Create reports in markdown and then convert it to pdf and html!

#XX Create usage of model in separate script file and use subroutines with NAM::function(args);
## XX Make an argument for an optional output folder
# XX Export only the ProteinScreen function and make the others inaccessible

# http://www.perlmonks.org/?node_id				       =  304000
