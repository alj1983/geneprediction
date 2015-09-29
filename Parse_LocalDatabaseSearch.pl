#!/usr/bin/perl -w

# Script to extract Best hit in the local transcriptome/genome with
# information such as evalue and reading frame from tblastn output and
# then selecting only those with an e-value <1e-05
use strict;
use warnings;
use LWP::Simple;
use List::MoreUtils qw(uniq); # This allows to find all unique values in an array


if ($ARGV[0] eq '' ) { die "text file has to be provided that lists
the names of tblastn report files which shall be treated in this
script\n"; }

open(IN, $ARGV[0]) or die "Can not open file $ARGV[0]";
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
    my @hitlength1; ## This will hold the hit length
    my @evals1; ## This will hold the e-values
    my @ranks1; ## This will hold the ranks
    my @queryframe1;  # this will hold the reading frame of the query
    my @targetframe1; # this will hold the reading frame of the target
    my @querydescriptions1;
    
#    open(IN, $ARGV[0]) or die "Can not open file $ARGV[0]";
foreach (@files){
    my $f=$_;

    open (IN,$f) or die "Can not open file $f";
    
    my $hitnumber=0;
    my $querynumber=1;

my $queryname1; # This will hold the single query name for one file
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
		push @querydescriptions1,$2;
	    }
	    elsif ($accession =~ /^.*emb\|(.*)\|\s*(\S*.*)/){
		$queryname1=$1;
		push @querydescriptions1,$2;
	    }
	    elsif ($accession =~ /^.*dbj\|(.*)\|\s*(\S*.*)/){
		$queryname1=$1;
		push @querydescriptions1,$2;
	    }
	    # NBRF PIR                          pir||entry
	    elsif ($accession =~ /^.*pir\|\|(\S*)\s*(\S*.*)/){
		$queryname1=$1;
		push @querydescriptions1,$2;
	    }
	    # Protein Research Foundation       prf||name
	    elsif ($accession=~ /^.*prf\|.*\|(.*)\s*(\S*.*)/){
		$queryname1=$1;
		push @querydescriptions1,$2;
	    }
	    # SWISS-PROT                        sp|accession|name
	    elsif ($accession=~ /^.*sp\|(.*)\|\s*(\S*.*)/){
		$queryname1=$1;
		push @querydescriptions1,$2;
	    }
	    # Brookhaven Protein Data Bank (1)  pdb|entry|chain
	    elsif ($accession=~ /^.*pdb\|(.*)\|\s*(\S*.*)/){
		$queryname1=$1;
		push @querydescriptions1,$2;
	    }
	    # GenInfo Backbone Id               bbs|number
	    elsif ($accession=~ /^.*bbs\|(.*)/){
		$queryname1=$1;
		push @querydescriptions1,$1;
	    }
	    # General database identifier       gnl|database|identifier
	    elsif ($accession=~ /^.*gnl\|(.*)\|\s*(\S*.*)/){
		$queryname1=$1;
		push @querydescriptions1,$2;
	    }
	    # NCBI Reference Sequence           ref|accession|locus
	    elsif ($accession=~ /^.*ref\|(.*)\|\s*(\S*.*)/){
		$queryname1=$1;
		push @querydescriptions1,$2;
	    }
	    # tpg|DAA34093.1| 
	    elsif ($accession=~ /^.*tpg\|(.*)\|\s*(\S*.*)/){
		$queryname1=$1;
		push @querydescriptions1,$2;
	    }
	    # Local Sequence identifier         lcl|identifier
	    elsif ($accession=~ /^.*lcl\|(.*)\|\s*(\S*.*)/){
		$queryname1=$1;
		push @querydescriptions1,$2;
	    }
	    elsif ($accession){
		$queryname1="";
		push @querydescriptions1,"";
	    }
	}
	if ($line =~ /^>/){
	    push @querynames1, $queryname1;
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
	    elsif ($line){
		push @hitnames1,"";
		push @hitdescriptions1,"";
	    }
	    $hitnumber++;
	    $querynumber=$hitnumber;
	    if ($hitnumber>0){
		push @ranks1, $hitnumber;
	    }
	}
	
	if ($line =~ /Score\s*=\s*(.*)\s*bits.*Expect\s*=\s*(\S*),.*$/g){
	    #push @estscores,$1;
	    push @evals1,$2;
	}
	if ($hitnumber>0){
	    #push @estscores,$1;
	    if ($line =~ /Length=(.*)$/g){
		push @hitlength1,$1; 
	    }
	}
	
	if ($line =~ /Frame = (.*)$/g){
	    push @targetframe1,$1; 
	    push @queryframe1,"1";
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
my @hitnames; ## This will hold the hit names
my @hitlength; ## This will hold the hit length
my @ranks; ## This will hold the ranks
my @queryframe;
my @targetframe;
my @querydescriptions;

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
my @hitnames_sig;
my @hitlength_sig;
my @evals_sig;
my @ranks_sig;
my @queryframe_sig;
my @targetframe_sig;
my @querydescriptions_sig;

my $o;
for ($o=0;$o<$evalsl;$o++){
    if ($evals[$o]<1e-05){
       push @querynames_sig, $querynames[$o];
       push @hitnames_sig, $hitnames[$o];
       push @hitlength_sig, $hitlength[$o];
       push @evals_sig, $evals[$o];
       push @ranks_sig, $ranks[$o];
       push @queryframe_sig,$queryframe[$o];
       push @targetframe_sig,$targetframe[$o];
       push @querydescriptions_sig,$querydescriptions[$o];
    }


}


############################################################
############# Identify the unique comps for each query #####
############################################################
# First get an array of the unique queries that were used: 
my @unique_queries = uniq @querynames_sig;
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
	$_ =~ /^(.*)_.*_.*$/;
	push @compnames, $1;
    }
    
    print MYOUTFILE1 "\nUnique comp hits \n";
    print MYOUTFILE1 "---------------- \n";
# get the unique comp hits for that query
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


