#!/usr/bin/perl
#Created by Duan Zhongqu, 2018-10-22
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h);
getopts("h");
my $usage="\nUsage: $0 [options] <sequence_file> <blast_file> <accession_file> <output_directory> <min_length> <min_identity>

$0 is used to detect and discard the potentail contaminated sequence.

Necessary input description:

  sequence_file        <string>      The sequence file.
  
  blast_file           <string>      The blast result of sequence file.

  accession_file       <string>      The accession file.

  output_directory     <string>      Both final output files and intermediate results 
                                     will be found in this directory. To avoid 
                                     overwriting of existing files. We kindly request
                                     that the output_directory should not exist. It is
                                     to say, this directory will be created by the 
                                     script itself.

  min_length           <int>         The minimum alignment.

  min_identity         <int>         The minimun identity.

Options:

     -h                            Print this usage page.

";

die $usage if @ARGV!=6;
die $usage if defined($opt_h); 
my ($seq_file,$blast_file,$accession_file,$out_dir,$min_length,$min_identity)=@ARGV;

#adjust directory names and create output directory
$out_dir.="/" unless($out_dir=~/\/$/);
mkdir($out_dir);
my $log_file=$out_dir."RemoveContaminate.log";
    
#Classify the accessions of target sequences according to the orangism information
open(IN,$accession_file)||die("Error9: cannot read the file: $accession_file!\n");
my %accession_Human;
my %accession_OtherPrimate;
my %accession_Microbiology;
my %accession_Other;
my $count_Human=0;
my $count_OtherPrimate=0;
my $count_Microbiology=0;
my $count_Other=0;
while(my $line=<IN>){
    chomp $line;
    my @string=split "\t", $line;
    if($string[3]=~/genus:Homo/){
        $accession_Human{$string[1]}=$string[2];
        $count_Human+=1;
    }
    elsif($string[3]=~/order:Primates/){
        $accession_OtherPrimate{$string[1]}=$string[2];
        $count_OtherPrimate+=1;
    }
    elsif($string[3]=~/superkingdom:Viruses/||$string[3]=~/superkingdom:Bacteria/||$string[3]=~/superkingdom:Archaea/){
        $accession_Microbiology{$string[1]}=$string[2];
        $count_Microbiology+=1;
    }
    else{
        $accession_Other{$string[1]}=$string[2];
        $count_Other+=1;
    } 
}    
close IN;
my $count_total=$count_Human+$count_OtherPrimate+$count_Microbiology+$count_Other;
open(LOG,">$log_file")||die("Error10: cannot write the file: $log_file.\n");
print LOG "Accession information:\n";
print LOG "In total, $count_total accessiones in the '$accession_file' file, of these\n\t$count_Human are belong to human;\n\t$count_OtherPrimate are belong to other primates;\n\t$count_Microbiology are belong to microbiology (including bacteria, virses and archaes);\n\tand $count_Other are belong to other orangism (including non-primate animals and plants).\n\n";
close LOG;
#Classify the contigs according to the blast result
my %contig_Human;
my %contig_OtherPrimate;
my %contig_Microbiology;
my %contig_Other;
open(IN,$blast_file)||die("Error9: cannot read the file: $blast_file!\n");
while(my $line=<IN>){
    chomp $line;
    next if $line=~/^#/;
    my @string=split "\t",$line;
    if((int($string[2])>=$min_identity)&&($string[3]>=$min_length)){
        if(exists($accession_Human{$string[1]})){
            $contig_Human{$string[0]}=$string[1];
        }
        elsif(exists($accession_OtherPrimate{$string[1]})){
            $contig_OtherPrimate{$string[0]}=$string[1];
        }
        elsif(exists($accession_Microbiology{$string[1]})){
            $contig_Microbiology{$string[0]}=$string[1];
        }
        elsif(exists($accession_Other{$string[1]})){
            $contig_Other{$string[0]}=$string[1];
        }
        else{
            print STDERR ("Warning: $string[1] have not in $accession_file\n");
        }
    }
}   
close IN;
    
#Output the classified contigs 
open(IN,$seq_file)||die("Error9: cannot read the file: $seq_file!\n");
my $seq_orangism_file=$out_dir."sequence.orangism.txt";
my $Human_seq_file=$out_dir."human.fa";
my $OtherPrimate_seq_file=$out_dir."other_primate.fa";
my $Microbiology_seq_file=$out_dir."microbiology.fa";
my $Other_seq_file=$out_dir."other.fa";
my $NotInNt_seq_file=$out_dir."not_in_nt.fa";
my $Novel_seq_file=$out_dir."novel_sequence.fa";
open(OUT,">$seq_orangism_file")||die("Error10: cannot write the file: $seq_orangism_file.\n");
open(OUT1,">$Human_seq_file")||die("Error10: cannot write the file: $Human_seq_file.\n");
open(OUT2,">$OtherPrimate_seq_file")||die("Error10: cannot write the file: $OtherPrimate_seq_file.\n");
open(OUT3,">$Microbiology_seq_file")||die("Error10: cannot write the file: $Microbiology_seq_file.\n");
open(OUT4,">$Other_seq_file")||die("Error10: cannot write the file: $Other_seq_file.\n");
open(OUT5,">$NotInNt_seq_file")||die("Error10: cannot write the file: $NotInNt_seq_file.\n");
open(NOVEL,">$Novel_seq_file")||die("Error10: cannot write the file: $Novel_seq_file.\n");
my $flag=0;
my $count_contig_Human=0;
my $count_contig_OtherPrimate=0;
my $count_contig_Microbiology=0;
my $count_contig_Other=0;
my $count_contig_NotInNt=0;
my $length_contig_Human=0;
my $length_contig_OtherPrimate=0;
my $length_contig_Microbiology=0;
my $length_contig_Other=0;
my $length_contig_NotInNt=0;
while(my $line=<IN>){
    chomp $line;
    if($line=~/^>/){
        my @string=split /\s+/,$line;
        my $name=substr($string[0],1,length($string[0])-1);
        if(exists($contig_Human{$name})){
            print OUT $name."\t".$contig_Human{$name}."\t"."human\n";
            print OUT1 $line."\n";
            print NOVEL $line."\n";
            $count_contig_Human+=1;
            $flag=1;
        }
        elsif(exists($contig_OtherPrimate{$name})){
            print OUT $name."\t".$contig_OtherPrimate{$name}."\t"."other primate\n";
            print OUT2 $line."\n";
            print NOVEL $line."\n";
            $count_contig_OtherPrimate+=1;
            $flag=2;
        }
        elsif(exists($contig_Microbiology{$name})){
            print OUT $name."\t".$contig_Microbiology{$name}."\t"."microbiology\n";
            print OUT3 $line."\n";
            $count_contig_Microbiology+=1;
            $flag=3;
        }
        elsif(exists($contig_Other{$name})){
            print OUT $name."\t".$contig_Other{$name}."\t"."other\n";
            print OUT4 $line."\n";
            $count_contig_Other+=1;
            $flag=4;
        }
        else{
            print OUT $name."\t--\t"."not in nt database\n";
            print OUT5 $line."\n";
            print NOVEL $line."\n";
            $count_contig_NotInNt+=1;
            $flag=5;
        }
    }
    else{
        if($flag==1){
            print OUT1 $line."\n";
            print NOVEL $line."\n";
            $length_contig_Human+=length($line);
        }
        elsif($flag==2){
            print OUT2 $line."\n";
            print NOVEL $line."\n";
            $length_contig_OtherPrimate+=length($line);
        }
        elsif($flag==3){
            print OUT3 $line."\n";
            $length_contig_Microbiology+=length($line);
        }
        elsif($flag==4){
            print OUT4 $line."\n";
            $length_contig_Other+=length($line);
        }
        elsif($flag==5){
            print OUT5 $line."\n";
            print NOVEL $line."\n";
            $length_contig_NotInNt+=length($line);
        }
    }
}
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;
close NOVEL;
my $count_contig_total=$count_contig_Human+$count_contig_OtherPrimate+$count_contig_Microbiology+$count_contig_Other+$count_contig_NotInNt;
my $length_contig_total=$length_contig_Human+$length_contig_OtherPrimate+$length_contig_Microbiology+$length_contig_Other+$length_contig_NotInNt;
open(LOG,">>$log_file")||die("Error10: cannot write the file: $log_file.\n");
print LOG "Result summary:\n";
print LOG "Type\tNumber\tLength\n";
print LOG "Human\t$count_contig_Human\t$length_contig_Human\n";
print LOG "Other Primate\t$count_contig_OtherPrimate\t$length_contig_OtherPrimate\n";
print LOG "Microbiology\t$count_contig_Microbiology\t$length_contig_Microbiology\n";
print LOG "Other\t$count_contig_Other\t$length_contig_Other\n";
print LOG "Not in nt database\t$count_contig_NotInNt\t$length_contig_NotInNt\n";
print LOG "Total\t$count_contig_total\t$length_contig_total\n";
close LOG;
