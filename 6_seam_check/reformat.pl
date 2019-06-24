#!/usr/bin/perl
use strict;

my $notprinted = 1;
my @head;
while(<STDIN>) {
    chomp;
    my @a = split;
    if (  /_rln(\w+)/ ) {
        push @head , $1
    }
    elsif ( $#a > 3 && /\w+/ ) {
        if ( $notprinted ) {
            print join( "\t", @head ) . "\n" ;
            $notprinted = 1;
        }
        print $_ . "\n";
    }
}
