#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my($input, $output);
GetOptions(

    "input=s" => \$input,
    "output=s"=> \$output,
);

open IN, $input or die $!;

my $i=0;
while(<IN>){

    chomp;
    $i++;
}

close IN;


my $dir = dirname($input);
my $sample = basename($dir);

if($i > 0){

    print "$sample\t$i\n";

}














