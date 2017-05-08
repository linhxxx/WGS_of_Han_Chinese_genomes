#!/usr/bin/perl -w
use strict;
my $t="/ifs5/BC_MEDEA/PROJECT/linhx/bin/bin/jdk1.7.0_60/bin/java";
while (<>){
   if (/Pending: (.*)/){
      my $buff=$1;
      $buff=~s/java/$t/;
      $buff=~s/'//g;
      print $buff,"\n";
   }
}
