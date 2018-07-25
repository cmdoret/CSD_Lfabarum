#!/usr/bin/env perl

#This script takes a sam file as input (plain or compressed) and splits it into three sam files. 1. unique hits; 2. reads that hit more than one location; 3. reads that hit nothing

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use Cwd qw(abs_path);
use POSIX qw(strftime);

my $target_dir=abs_path;
my ($in_sam, $prefix, $line,$hitcount);
my (@headers,@reads,@nohits,@multi,@hits,@unique);
my ($nohits, $multi, $unique, $limit);
my ($c_nohits, $c_multi, $c_unique)=(0,0,0);

GetOptions (
	"in:s"	=> \$in_sam,
	"out:s"	=> \$prefix,
	"mem:i"	=> \$limit
#	"fq!"	=> \$fq
);

my $USAGE="\nThis script takes as input a sam file (plain or compressed) and splits it (if applicable) into three sam files containing: 1. unique hits; 2. reads with multiple hits; 3. reads without hits. 
		Usage: split_sam.pl -i infile.sam -o outprefix 

			optional:
			-m <int>	maximum memory usage in Gb\n\n";
################

if ((!$in_sam) || (!$prefix)){
	print "$USAGE\n";
	exit;
}

my $sam_fh = &read_fh ($in_sam);
print "\nsplit_sam.pl is processing $in_sam .. ";
open (UNIQUE, ">$prefix-uniq.sam") or die $!;
open (MULTI, ">$prefix-multi.sam") or die $!;
open (NOHIT, ">$prefix-nohit.sam") or die $!;
################

#chdir $base_dir or die $!;
#$base_dir=abs_path;

while (<$sam_fh>) {
	chomp;
	if ($_ =~ /^@/){
		push (@headers,$_);
		print UNIQUE "$_\n";
		print MULTI "$_\n";
		print NOHIT "$_\n";
	}else{
		$line=$_;
		$hitcount=0;
#		print "testing $line\n";
		@reads = split ("\t");
		if ($reads[2] =~ /\*/){
#			push (@nohits,$_);
			$nohits.="$_\n";
			$c_nohits++;
		}else {
			@hits=grep(/X[0-1]:i/, @reads);
#			print "@hits\n";
			for (@hits){
				$_=~s/X[0-1]:i://;
#				print "$_\n";
				$hitcount+=$_;
			}
#			print "hitcount is $hitcount\n";		
			if ($hitcount >= 2){
#				push (@multi,$line);
				$multi.="$line\n";
				$c_multi++;
			}else{
				$unique.="$line\n";
				$c_unique++;
#				push (@unique,$line);
			}
		}
	}
	if ($limit){
		if ($. % 1000000 == 0){
			if (&check_memory($limit, $prefix)){
				print NOHIT "$nohits";
				print MULTI "$multi";
				print UNIQUE "$unique";
				undef $nohits;
				undef $multi;
				undef $unique;
				($nohits, $multi, $unique)=("","","");
			}
		}
	}
}

#writing to files
print NOHIT "$nohits";
print MULTI "$multi";
print UNIQUE "$unique";

#for (@nohits){
#	print NOHIT "$_\n";
#}
#for (@multi){
#	print MULTI "$_\n";		
#}
#for (@unique){
#	print UNIQUE "$_\n";
#}

my $total=$c_nohits+$c_multi+$c_unique;
#my $total=(scalar @nohits) + (scalar @multi) + (scalar @unique);
my $unique_perc = sprintf "%.2f", ($c_unique/$total*100);
#my $unique_perc = sprintf "%.2f", (scalar @unique/$total*100);
my $multi_perc = sprintf "%.2f", ($c_multi/$total*100);
#my $multi_perc = sprintf "%.2f", (scalar @multi/$total*100);
my $nohit_perc = sprintf "%.2f", ($c_nohits/$total*100);
#my $nohit_perc = sprintf "%.2f", (scalar @nohits/$total*100);

my $summary="
Summary for $in_sam:

total number of reads: $total
unique hit reads: $c_unique ($unique_perc %)
multi hit reads: $c_multi ($multi_perc %)
no hit reads: $c_nohits ($nohit_perc %)\n";

print "done!\n";

#open (SUM, ">$prefix.log") or die $!;
#print SUM "$summary\n";
print "$summary\n";
#close SUM;

######################
sub read_fh {
	my $filename = shift @_;
	my $filehandle;
	if ($filename =~ /gz$/) {
		open $filehandle, "gunzip -dc $filename |" or die $!;
	}
	else {
		open $filehandle, "<$filename" or die $!;
	}
	return $filehandle;
}

sub check_memory{
        my $return=0;
        my $limit=shift;
        my $pattern=shift;
        my (@out, @var);
        my $mem;
#       print "$0\n";
        print "checking memory usage\t". strftime("%b %e %H:%M:%S", localtime)."\n";
        my @ps_aux=`ps aux | grep $0`;
        print "@ps_aux\n";
#       print "looking for $pattern\n";
        for (@ps_aux){
                if ($_ =~ /$pattern/){
                        print "ps aux output: \n$_\n";
                        @out=split("\n",$_);
                        @var=split('\s+');
#                       print "values: @var\n";
#                       print "percentage memory used: $var[3]\n";
#                       print "physical memory used: $var[5]\n";
                        $mem=$var[5]/1048576;
                        print "$mem Gb\n";
                        if ($mem > $limit){
                                $return=1;
                                print "memory usage exceeds maximum ($limit G) -> will be emptied\n";
                        }
        
                }
        }
        return $return;
}

