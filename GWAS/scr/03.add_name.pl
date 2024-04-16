#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my($ref, $input, $output);
GetOptions(

    "ref=s" => \$ref,
    "input=s" => \$input,
    "output=s" => \$output,
);

open IN, $ref or die $!;
my %hash;
while(<IN>){

    chomp;
    next if $_=~/^ID/;

    my @arry=split /\t/,$_;
    $hash{$arry[0]}= $arry[1];
}
close IN;

open IN, $input or die $!;

while(<IN>){

    chomp;
    my @tmp=split /\t/,$_;

    if(exists $hash{$tmp[0]}){

         print "$tmp[0]\t$hash{$tmp[0]}\t$tmp[1]\n";
    }else{


         print "$tmp[0]\tNA\t$tmp[1]\n";
    }

}


close IN;























