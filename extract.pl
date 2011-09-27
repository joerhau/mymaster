#!/usr/bin/perl
use Math::Round;

$nr=$ARGV[0];

$dir="./2_1_15taxa/";
$file1="prot.phy";
$file2="partitionFile.txt";
open DAT1, $dir.$file1 or die "Could not open file!";
open DAT2, $dir.$file2 or die "Could not open file!";
my @lines1=<DAT1>;
my @lines2=<DAT2>;

if(@lines2 < $nr){
exit "less than ".$nr." partitions...";
}
$size=@lines2;
print "nr partitions: ".$size."\n";
#print @lines2;
print "\n";

my @part = ();
my @parte = ();
my @parts = ();

#if($nr<$size){
for($i=0;$i<$nr;$i++){

$r=round(rand()*($size));
$part[$i]=$lines2[$r];
if($lines[nr]=~ m/\d*/){
print "hit: ".$1;
}
#print $part[$i]."\n";
}
#} else {

#}
print "\n";
#use Data::Dumper;
#print Dumper $lines;
