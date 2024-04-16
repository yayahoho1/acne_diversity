#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my($ref, $output);
GetOptions(

    "ref=s" => \$ref,
    "output=s" => \$output,
);

open IN, $ref or die $!;

my %hash;
while(<IN>){

    chomp;
    my @arry=split /\t/,$_;

    my $pos = $arry[1] . "-" . $arry[2];
    $hash{$pos} = $arry[4];
}

close IN;


my @group=@ARGV;
my $length=@ARGV;

my %hash_sig;
my %hash_phenotype;
my @suffixlist = qw(_significant_SNP.txt);

for(my $i=0; $i < $length;$i++){

    open IN, $group[$i] or die $!;
    my $phenotype= basename($group[$i],@suffixlist);
    while(<IN>){

        chomp;
        next if $_=~/^variant/;

        my @tmp=split /\t/,$_;
        my $snppos=(split/\_/,$tmp[0])[1];
        my $j=0;
        foreach my $key (sort keys %hash){

            my($start, $end)=split /\-/,$key;
            if($snppos >= $start and $snppos <= $end){
                $j++;
                $hash_sig{$hash{$key}}++;
                $hash_phenotype{$phenotype}{$hash{$key}}++;
            }
        }

        if($j==0){

            print "$phenotype\t$_\n";
        }
        #else{
        #    print "This is $j\n";
        #}
    }

    close IN;


}


my $gene_snp= "gene2snpnumber_stat.txt";

my $pheno_gene_snp = "phenotype_gene2snpnumber_stat.txt";

open OUT1,">>$gene_snp" or die $!;
open OUT2,">>$pheno_gene_snp" or die $!;

foreach my $key (sort keys %hash_sig){

    print OUT1 "$key\t$hash_sig{$key}\n";

}

foreach my $item (sort keys %hash_phenotype){

    foreach my $tt (sort keys %{$hash_phenotype{$item}}){

        print OUT2 "$item\t$tt\t$hash_phenotype{$item}{$tt}\n";
    }

}

close OUT1;
close OUT2;







