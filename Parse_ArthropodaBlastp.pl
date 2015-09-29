#!/usr/bin/perl -w

# Script to extract Best hit with information such as score and evalue
# and description from a blast text output;
use strict;
use warnings;
use LWP::Simple;

if ($ARGV[0] eq '' ) { die "text file has to be provided that lists
the names of blast report files which shall be treated in this
script\n"; }

open(IN, $ARGV[0]) or die "Can not open file $ARGV[0]";

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
    my $blastdbcmd = `blastdbcmd -db nr_arthropoda -entry $query >> besthits.fasta`;

}
close MYOUTFILE2;
