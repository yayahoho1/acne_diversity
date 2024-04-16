#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my($input, $output);
GetOptions(

    "input=s" => \$input,
    "output=s" => \$output,
);

open IN, $input or die $!;

my %hash;
while(<IN>){

    chomp;
    my @arry=split /\t/,$_;
    if(exists $hash{$arry[1]}){

        $hash{$arry[1]}= $hash{$arry[1]} + $arry[2];
    }else{

        $hash{$arry[1]}=$arry[2];
    }

}

close IN;

foreach my $key (sort keys %hash){

    print "$key\t$hash{$key}\n";
}











