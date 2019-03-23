#!/usr/bin/perl
#Created by Duan Zhongqu at 2018-11-5
use strict;
use warnings;
use Getopt::Std;

my $usage="Usage: $0 seq_file gff_file tra_file pro_file genome_ref tran_ref cdhit_exe blast_exe mkblastdb_exe repeat_exe gene_identity ref_identity thread_num out_dir\n";

die $usage if @ARGV!=14;

my ($seq_file,$gff_file,$tra_file,$pro_file,$genome_ref,$tran_ref,$cdhit_exe,$blast_exe,$mkblastdb_exe,$repeat_exe,$gene_identity,$ref_identity,$thread_num,$out_dir)=@ARGV;

my $log_file=$out_dir."FilterNovelGene.log";
open(LOG,">$log_file")||die("Cannot write the file.\n");
print LOG "\nCommand: $0 $seq_file $gff_file $tra_file $pro_file $genome_ref $tran_ref $cdhit_exe $blast_exe $mkblastdb_exe $repeat_exe $gene_identity $ref_identity $thread_num $out_dir\n\n";
print LOG "\tInput file:\n\t$seq_file\n\t$gff_file\n\t$tra_file\n\t$pro_file\n\n";
print LOG "\tReference information:\n\t$genome_ref\n\t$tran_ref\n\n";
print LOG "\tOutput directory:\n\t$out_dir\n\n";
print LOG "\tBinary executable file:\n\t$cdhit_exe\n\t$blast_exe\n\t$mkblastdb_exe\n\t$repeat_exe\n\n";
print LOG "\tOptions:\n\tGene identity: $gene_identity\n\tReference identity: $ref_identity\n\tThread number: $thread_num\n\n";

#read the gff file and get the number of novel predicted genes.
my $raw_out_dir=$out_dir."raw/";
mkdir($raw_out_dir);
my $filter_100bp_dir=$out_dir."100bps/";
mkdir($filter_100bp_dir);
open(GFF,$gff_file)||die("Error9: cannot read the gff file.\n");
my $new_gff_file=$raw_out_dir."maker.raw.gff";
open(OUT1,">$new_gff_file")||die("Error10: cannot write the file: $new_gff_file.\n");
my $raw_gene_map=$raw_out_dir."maker.gene.map";
open(OUT,">$raw_gene_map")||die("Error10: cannot write the file: $raw_gene_map.\n");
print OUT "sequence name\tstart\tend\tgene name\n";
my $gene_map=$filter_100bp_dir."gene_100bps.map";
open(OUT2,">$gene_map")||die("Error10: cannot write the file: $gene_map.\n");
print OUT2 "sequence name\tstart\tend\tgene name\n";
my $count_raw_gene=0;
my $count_gene=0;
print LOG "--------------------\nRead $gff_file file...\n\tall the annotation informations from maker will be selected and writtend into $new_gff_file.\n\tthe positions of each predicted gene were recorded in $raw_gene_map.\n";
while(my $line=<GFF>){
    chomp $line;
    if($line=~/maker/){
        print OUT1 $line."\n";
        if($line=~/maker\tgene/){
            $count_raw_gene+=1;
            my @string=split "\t",$line;
            my $gene_length=$string[4]-$string[3]+1;
            my @tmp=split ";",$string[8];
            my $gene_id=substr($tmp[0],3);
            print OUT $string[0]."\t".$string[3]."\t".$string[4]."\t".$gene_id."\n";
            if($gene_length>=100){
                $count_gene+=1;
                print OUT2 $string[0]."\t".$string[3]."\t".$string[4]."\t".$gene_id."\n";
            }
        }
    }        
}
print LOG "\tthere are $count_raw_gene novel gene were predicted by the maker tools.\n";
print LOG "\tof these, $count_gene genes are longer than 100 bp.\n\tthe positions of each predicted gene were recorded in $gene_map.\n";
close GFF;
close OUT1;
close OUT2;
close OUT;

#read the novel, transcript and protein sequences onto hash
my %novel_seq=readSeq($seq_file);
my %tran_seq=readSeq($tra_file);
my %prot_seq=readSeq($pro_file);
#obtain the gene sequences from novel sequences based on the gff file.
open(IN,$gene_map)||die("Error9: cannot read the file: $gene_map.\n");
my $gene_file=$filter_100bp_dir."gene.100bps.fasta";
open(OUT,">$gene_file")||die("Error10: cannot write the file: $gene_file.\n");
readline IN;
my %gene_seq;
my %gene_length;
print LOG "\twrite these genes into $gene_file.\n";
while(my $line=<IN>){
    chomp $line;
    my @string=split "\t",$line;
    my $name=$string[0];
    my $start=$string[1];
    my $end=$string[2];
    my $length=$end-$start+1;
    my $gene_id=$string[3];
    my $all_seq=$novel_seq{$name};
    my $seq=substr($all_seq,$start-1,$length);
    $gene_seq{$gene_id}=$seq;
    $gene_length{$gene_id}=length($seq);
    print OUT ">".$gene_id."\t".$name."\t".$start."\t".$end."\n".$seq."\n";
}
close IN;
close OUT;
my $iden=$gene_identity*100;
print LOG "Done.\n--------------------\nRemoving the redundancy genes with sequence identity of ".$iden."%...\n";
#remove the redundant gene sequence
my $rm_redundancy_dir=$out_dir."Non-redundancy/";
mkdir ($rm_redundancy_dir);
my $cdhit_out=$rm_redundancy_dir."non.redundant.gene.fasta";
my $cdhit_log=$rm_redundancy_dir."cdhit.log";
my $com="$cdhit_exe -i $gene_file -o $cdhit_out -T $thread_num -c $gene_identity >$cdhit_log";
system($com);
#Count the number of rest novel gene
my @non_redun_gene=extract_gene_name($cdhit_out);
my $non_redun_gene_count=@non_redun_gene;
print LOG "\tthere are $non_redun_gene_count novel genes after removing homology gene sequences in $cdhit_out.\n";
   
#align the novel sequence to reference genome seqeuces
print LOG "Done.\n-------------------\nAlign the gene sequences to reference genome...\n";
my $genome_ref_dir=$out_dir."ref_genome/";
mkdir($genome_ref_dir);
my $genome_db=$genome_ref_dir."blast_db_genome";
my $genome_blast=$genome_ref_dir."non.redundant.gene.blast2genome";
runblast($blast_exe,$mkblastdb_exe,$genome_ref,$genome_db,$cdhit_out,$genome_blast,$thread_num);
 
#discard the gene could align to reference genome at the global identity cutoff 50%
my %gene_similar_genome=filterblast($genome_blast,\%gene_length);
my @gene_non_similar_genome=grep {!$gene_similar_genome{$_}} @non_redun_gene;
my $gene_non_similar_genome_count=@gene_non_similar_genome;
print LOG "\tthere are $gene_non_similar_genome_count novel genes after removing genes similar with the reference genome.\n";
my $gene_non_similar_genome_file=$genome_ref_dir."Non-similar2genome.gene.fasta";
my @genes=extractSub($cdhit_out,$gene_non_similar_genome_file,"gene",\@gene_non_similar_genome);
my $tran_non_similar_genome_file=$genome_ref_dir."Non-similar2genome.transcripts.fasta";
my @tran_non_similar_genome=extractSub($tra_file,$tran_non_similar_genome_file,"transcripts",\@gene_non_similar_genome);
print LOG "Done.\n--------------------\nAlign the gene sequences to reference transcript...\n";
#align the transcripts sequence to reference transcript sequences
my $tran_ref_dir=$out_dir."ref_tran/";
mkdir($tran_ref_dir);
my $tran_db=$tran_ref_dir."blast_db_transcripts";
my $tran_blast=$tran_ref_dir."non.similar2genome.transcripts.blast2transcript";
runblast($blast_exe,$mkblastdb_exe,$tran_ref,$tran_db,$tran_non_similar_genome_file,$tran_blast,$thread_num);

#discard the gene could align to transcript genome at the global identity cutoff 50%
my %tran_length;
while((my $key,my $value)=each %tran_seq){
    $tran_length{$key}=length($value);
}
my %tran_similar_genome=filterblast($tran_blast,\%tran_length);
my @tran_non_similar_transcript=grep {!$tran_similar_genome{$_}} @tran_non_similar_genome;
my %hash1;
foreach my $sub (@tran_non_similar_transcript){
    my $tmp=substr($sub,0,length($sub)-3);
    if(grep /^$tmp$/,@gene_non_similar_genome){
        $hash1{$tmp}=1;
    }
}
my @gene_non_similar_transcript=keys %hash1;
my $gene_non_similar_transcript_count=@gene_non_similar_transcript;
print LOG "\tthere are $gene_non_similar_transcript_count novel genes after removing genes similar with the reference transcript.\n";
my $gene_non_similar_transcript_file=$tran_ref_dir."Non-similar2transcript.gene.fasta";
@genes=extractSub($gene_non_similar_genome_file,$gene_non_similar_transcript_file,"gene",\@gene_non_similar_transcript);
my $tran_non_similar_transcript_file=$tran_ref_dir."Non-similar2transcript.transcripts.fasta";
my @trans=extractSub($tran_non_similar_genome_file,$tran_non_similar_transcript_file,"transcript",\@gene_non_similar_transcript);
print LOG "Done.\n--------------------\nDetermine the full-length genes...\n";
#define whether the novel gene is complete according to start codon and stop codon.
my %filtered_tran_seq=readSeq($tran_non_similar_transcript_file);
my %complete_gene_transcript_seq;
while((my $key,my $value)=each %filtered_tran_seq){
    my $start=substr($value,0,3);
    my $end=substr($value,length($value)-3);
    if(($start eq "ATG")&&($end eq "TAG"||$end eq "TAA"||$end eq "TGA")){
        $complete_gene_transcript_seq{$key}=$value;
    }
}
my %hash2;
while((my $key,my $value)=each %complete_gene_transcript_seq){
    my $tmp=substr($key,0,length($key)-3);
    $hash2{$tmp}=1;
}
my @complete_gene_name=keys %hash2;
my $complete_gene_dir=$out_dir."Full-length/";
mkdir($complete_gene_dir);
my $complete_gene_file=$complete_gene_dir."Full-length.gene.fasta";
@genes=extractSub($gene_non_similar_genome_file,$complete_gene_file,"gene",\@complete_gene_name);
my $complete_tran_file=$complete_gene_dir."Full-length.transcript.fasta";
@trans=extractSub($tran_non_similar_transcript_file,$complete_tran_file,"transcript",\@complete_gene_name);
my $complete_gene_count=@complete_gene_name;
print LOG "\tthere are $complete_gene_count full-length novel genes.\n";
print LOG "Done.\n--------------------\nRepeatmask the sequences of rest genes...\n";
#repeat component of full-length novel gene.
my $repeat_out_dir=$out_dir."Repeat/";
mkdir($repeat_out_dir);
$com="$repeat_exe -parallel $thread_num -species human -html -gff -dir $repeat_out_dir $complete_gene_file";
system($com);
my $masked_gene_file=$repeat_out_dir."Full-length.gene.fasta.masked";
my %masked_gene_seq=readSeq($masked_gene_file);
my %final_gene;
while((my $key,my $value)=each %masked_gene_seq){
    my $length=length($value);
    my $count = $value =~ tr/N/N/;
    if($count/$length<=0.5){
        $final_gene{$key}=1;   
    }
}
print LOG "Done.\n--------------------\nProcude final data set of novel genes...\n";
my $final_gene_count=keys %final_gene;
print LOG "\tthere are $final_gene_count full-length novel genes.\n";
my @final_gene_name=keys %final_gene;
my $final_out_dir=$out_dir."Final/";
mkdir($final_out_dir);
my $final_gene_file=$final_out_dir."Final.gene.fasta";
@genes=extractSub($complete_gene_file,$final_gene_file,"gene",\@final_gene_name);
my $final_tran_file=$final_out_dir."Final.transcript.fasta";
@trans=extractSub($complete_tran_file,$final_tran_file,"transcript",\@final_gene_name);
my $final_prot_file=$final_out_dir."Final.protein.fasta";
my @proteins=extractSub($pro_file,$final_prot_file,"protein",\@final_gene_name);
my $final_gff_file=$final_out_dir."Final.gff";
open(IN,$new_gff_file)||die("Error09: cannot read the file: $new_gff_file");
open(OUT,">$final_gff_file")||die("Errpr10: cannot write the file: $final_gff_file");
while(my $line=<IN>){
    chomp $line;
    my @string=split "\t", $line;
    my @tmp=split ";", $string[8];
    my $name=substr($tmp[0],3,length($final_gene_name[0]));
    if(grep /^$name$/, @final_gene_name){
        print OUT $line."\n";
    }
}
print LOG "Done.\n--------------------\n";
print LOG "The final data set includes four files:\n\tgene sequences: $final_gene_file\t\ntranscript sequences: $final_tran_file\t\nprotein sequences: $final_prot_file\t\nannotation information: $final_gff_file.\n";
close IN;
close OUT;
close LOG;

sub readSeq{
    use strict;
    use warnings;
    my ($seq_file)=@_;
    my %sequence;
    open(IN,$seq_file)||die("Error09: Could not read file: $seq_file.\n");
    my $count=0;
    my $seq;
    my $name;
    while(my $line=<IN>){
        chomp $line;
        if($line=~/^>/){
            my @string=split /\s+/, $line;
            if($count==0){
                $count+=1;
                $name=substr($string[0],1);
            }
            else{
                $sequence{$name}=$seq;
                $name=substr($string[0],1);
                $seq="";
            }
        }
        else{
            $seq.=$line;
            if(eof){
                $sequence{$name}=$seq;
            }
        }
    }
    close IN;
    close OUT;
    return %sequence;
}

sub runblast{
    use strict;
    use warnings;
    my($blast_exe,$mkblastdb_exe,$ref,$db,$input,$output,$thread_num)=@_;
    my $com="$mkblastdb_exe -in $ref -dbtype nucl -out $db";
    $com.="\n$blast_exe -query $input -db $db -out $output -outfmt 7 -max_target_seqs 1 -num_threads $thread_num\n";
    #print $com."\n";
    system("$com");
}
1;

sub extract_gene_name{
    use strict;
    use warnings;
    my ($input_file)=@_;
    open(IN,$input_file)||die("Error09: cannot read the file: $input_file.\n");
    my @names;
    while(my $line=<IN>){
        chomp $line;
        if($line=~/^>/){
            my @string=split "\t",$line;
            push @names, substr($string[0],1);
        }
    }
    return @names;
}

sub filterblast{
    use strict;
    use warnings;
    my ($blast_file,$gene_length)=@_;
    my %gene_length=%$gene_length;
    open(IN,$blast_file)||die("Error09: cannot read the file: $blast_file");
    my %names;
    while(my $line=<IN>){
        chomp $line;
        next if $line=~/^#/;
        my @string=split "\t", $line;
        my $name=$string[0];
        my $identity=$string[2];
        my $align_rate=$string[3]/$gene_length{$name};
        if($align_rate*$identity>=50){
        #if(($align_rate>=0.5)&&($identity>=80)){
            $names{$name}=1;
        }
    }
    return %names;
}
1;
sub extractSub{
    use strict;
    use warnings;
    my ($input,$output,$type,$subsets)=@_;
    my @subsets=@$subsets;
    open(IN,$input)||die("Error09: cannot read the file: $input.\n");
    open(OUT,">$output")||die("Error10: cannot write the file: $output.\n");
    my $flag=0;
    my %names;
    while(my $line=<IN>){
        chomp $line;
        if($line=~/^>/){
            my @string=split /\s+/,$line;
            my $name;
            my $tmp=substr($string[0],1);
            $names{$tmp}=1;
            if($type=~/gene/){
                $name=substr($string[0],1);
            }
            elsif($type=~/transcript/||$type=~/protein/){
                $name=substr($string[0],1,length($string[0])-4);   
            }
            else{
                die("Please define the sequence type: 'gene', 'transcript' or 'protein'.\n");
            }
            if(grep /^$name$/, @subsets){
                print OUT $line."\n";
                $flag=1;
            }
            else{
                $flag=0;
            }
        }
        else{
            if($flag==1){
                print OUT $line,"\n";
            }
        }
    }
    close IN;
    close OUT;
    return %names;
}
