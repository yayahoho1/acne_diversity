#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
my($input, $threshold, $output);
GetOptions(

    "input=s" => \$input,
    "threshold=s"=> \$threshold,
    "output=s" => \$output,
);

open IN, $threshold or die $!;
my %hash;
while(<IN>){

    chomp;
    my @arry=split /\t/,$_;
    #print Dumper(@arry);
    $hash{$arry[0]}=$arry[1];
}
close IN;

#print "$hash{"Threshold:"}"


open IN, $input or die $!;

while(<IN>){

    chomp;
    next if $_=~/^variant/;
    my @tmp=split /\t/,$_;

    if($tmp[3] < $hash{"Threshold:"}){

        print "$_\n";
    }


}

close IN;











