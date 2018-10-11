
use strict;
use warnings;
use lib "/home/local/ARCS/nz2274/Pipeline/NA_script/";
use INFO2map;
sub cell2array {
	my @arr=($_[0]);
	return \@arr;
}

print "This script is used to correct the genotype and info filed for each individuals when several ALT occured in the same site in the joint calling.\n ";
print "Input format: raw_file_tab_delimited  header  OUTPUT_FOLDER\n";
print "please use ## to comment out unwantted lines, the first line (not start with ##) is taken as the colnames.\n";
my $file=$ARGV[0];
my @bases=split(/\//,$file);
my $fb=$bases[@bases-1];
print $fb.".extract.txt\n";
#exit;;
my %samples=();
my $dout="./";
if( @ARGV >2 ) {$dout=$ARGV[2];}
#my $ouf_folder=
open IN,"less $file|" or die "Error: $file cannot find!\n";

my $proband=uc ($ARGV[3]);
my $fheader=$ARGV[1];
my $freq=0.01;
my $laf=0.1;
my $ped=$ARGV[4];
if(@ARGV>5){$freq=$ARGV[5]}
my %header=();
my @keys=();
my @parents=();
#print $proband."\t".$ped."\n";
open P, "egrep $proband -i $ped|";
print "search $proband in $ped \n"; 
my $relation=<P>;
my @pedigree=split(/\s+/,$relation);
for(my $i=1;$i<4;$i++){
	
	if(length($pedigree[$i])>2){
		push(@parents,$pedigree[$i]);

	}
}
close P;

open HEAD ,$fheader or die "exit! cannot find $fheader, please specify the header file\n";
my $c=0;
while(my $line=<HEAD>){
   chomp($line);
    foreach my $s (split("\t",$line)){
     $s=uc($s);
     $header{$s}=$c;
     push(@keys,$s);
    
	 #print $s."\t".$c."\n";
    # print OUT $s."\t";
	 $c=$c+1;
 }
}
#print OUT "\n";
close HEAD;
foreach my $id(@parents){
	$id=uc($id);
	$header{$id}=$c;
	push(@keys,$id);
	
	$c=$c+1;
}


open OUT, ">$dout/$fb.extract.txt";
print "Please find out in $dout/$fb.extract.txt\n";

my $lc=0;
#push(@header,"FILTER");
while(my $line=<IN>){
	chomp($line);
	#print $line."\n";
	if ($line =~ /^##/){
		next;
	}elsif ($lc==0){
		my %check=();
		my @sets=split("\t",$line);
		print "Total columns: ".@sets."\n";
		my @keeps=();
		for(my $i=0;$i<@sets;$i++){
			$sets[$i]=uc($sets[$i]);
			if(exists($check{$sets[$i]})){next;}
			$check{$sets[$i]}=1;
			#	print "keep\t".$sets[$i]."\t".exists($header{$sets[$i]})."\n";
		    if(exists($header{$sets[$i]})){
		     	push(@keeps,$header{$sets[$i]});
			
				 $header{$sets[$i]}=$i;
				 
			}else{
			 #print $sets[$i]."\t".$proband."\t".(index($sets[$i],$proband))."\n";
		#	 if(index($sets[$i],$proband) !=-1){
			# print "Found ".$sets[$i]."\n";
			 my @sub_h=split(/\./,$sets[$i]);
			 my $subk=join(".",@sub_h[0..(@sub_h-2)]);
 			#if(exists($check{$subk})){next;}
 			#$check{$subk}=1;
			#  print "keep2 $sets[$i]\t".$subk."\t".(@sub_h)."\t"."\n";
			 if(index($sets[$i],$proband) !=-1 && exists($header{$subk})){
                	   push(@keeps,$header{$subk});
			  #print "keep ".$subk."\t".(@sub_h)."\t".$header{$subk}."\n";
	                   $header{$subk}=$i;
	  			   

        	   }
		      #}
		    }	
			#$samples{$i}=$sets[$i];
		}
		@keeps=sort { $a <=> $b } @keeps;
	    @keys=@keys[@keeps];
	 	foreach my $item(@keys){
	 #		   print "find\t".$item."\t".$header{$item}."\n";	
	 	    print OUT $item."\t";
	 		
	 	}
	 	print OUT "\n";
        $lc=$lc+1;
		next;
	}
		
		my @sets=split("\t",$line);
		my $pc=0;
		if(exists($header{"POPSCORE"})){
			if($sets[$header{"POPSCORE"}] >$freq){next;}
		}
		if(exists($header{"AF"})){
                        if($sets[$header{"AF"}] >$laf){next;}
                }
		foreach my $item(@keys){
		#   print "find\t$pc\t".$item."\t".$header{$item}."\n";	
		    print OUT $sets[$header{$item}]."\t";
			$pc=$pc+1;
		}
		print OUT "\n";
  #     exit;
	
}


close IN;
close OUT;
