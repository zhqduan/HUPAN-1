use strict;
use warnings;
package fastaSta;

sub sta{
    my $usage="hupan fastaSta <fasta>
";
    die $usage if @ARGV!=1;

    open(IN,$ARGV[0]);
    my $i=0;
    my $j=0;
    while(<IN>){
	chomp;
	if(/^>(.+)$/){
	    $j++;
	}
	else{
	    $i+=length($_);
	}
    }

    print STDOUT "Total sequence number: ",$j,"\n";
    print STDOUT "Total base number: ",$i,"\n";
}
1;
