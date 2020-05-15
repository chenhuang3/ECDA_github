@data=`ls libpaw/src/*F90`;
$nlines=@data;

$n=0;
#print "You have  $nlines files in this directory \n";

#$file=$ARGV[0];
foreach(@data){
 chop($_);
 #print "$_\n";
 $file=$_;
 FORMATfile(); 

}

sub FORMATfile {
  open TEMP, "<$file" or die "Opening '$file': $!\n";
  $fileout=$file.".tmp";
  open OUT, ">$fileout" or die "Opening '$fileout': $!\n";

  $iline=0;
  while ($line = <TEMP>){
   $iline=$iline+1;
   if($line=~/\s*use\s*interfaces/){
    if( ($line=~/\s*use\s*interfaces_\d+_hidewrite/) ||
        ($line=~/\s*use\s*interfaces_\d+_hideleave/)){
     print(OUT "$line");
    }
    else{
     print(OUT "!$line");
    }
   }
   else{
    print(OUT $line);
   }#if
  }
 
  close OUT; 
  close TEMP;
  rename($fileout,$file);
 
}

###############
sub ARGVerror {
###############
    chop ( $case = qx /  echo \$PWD \| awk -F \/ '{print\$NF} '  /  );   #directory at which the user is working
    system("clear");
    print "\e[0;31m*******************************************************\e[00m \n";
    print "The scripts needs to be run as follows:\n";
    print "\n";
    print "format-bib.pl \e[0;34m file.bib \e[0;34m  \e[0m\n";
    print "\n";
    print "file.bib is a .bib file to be formated. \n";
    print "\e[0;31m******************************************************\e[00m \n";
    exit(1);
}

# Declare the subroutines


# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
# Left trim function to remove leading whitespace
sub ltrim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	return $string;
}
# Right trim function to remove trailing whitespace
sub rtrim($)
{
	my $string = shift;
	$string =~ s/\s+$//;
	return $string;
}
