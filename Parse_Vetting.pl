#!/usr/bin/perl -w

# Script to extract Best hit with information such as score and evalue
# and description from a blast text output;
use strict;
use warnings;
use LWP::Simple;

if ($ARGV[0] eq '' ) { die "text file has to be provided that lists
the names of tblastn report files which shall be treated in this
script\n"; }

open(IN, $ARGV[0]) or die "Can not open file $ARGV[0]";
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
