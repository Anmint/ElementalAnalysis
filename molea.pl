#!/usr/bin/perl

#"read_atom" was cited from http://www.shylph.info/blog/programming/perl/257

#use strict;
use warnings;
use utf8;
use Encode;
use Data::Dumper;
#use Getopt::std;
binmode STDOUT, ":utf8";

#constant declaration
my $atomfile = "molea.ini";

#--------------------------#
#molea.ini have .csv style.
#[Atom]
#H,1.008
#C,12.01
#atom,weight
#--------------------------#

#variable declaration
my $peritable = &read_atom;
my $molecule = &set_mol;
&calc_mol($molecule);

#sort a argument to elemental analysis style -> return as a reference of 2Dlist
#list = ([atom,quantity],[atom,quantity]...)
sub set_mol{
    my @mol = ();
    eval{
  my $str = $ARGV[0] or die $!;
#separate between atoms.
	@mol = split /(?=[A-Z])/, $str;
    };
    if($@){
	print "ERROR!: $@";
	exit 1;
    }
    my @mol2D = ();
    foreach(@mol){
	my @tmp = ();
	if($_ !~ /\d{1,}/){
	    push @mol2D, [$_, 1];
	}else{
	@tmp = $_ =~ m/[A-Z][a-z]*|\d{1,}/g;
	push @mol2D, \@tmp;
	}
    }
#arrange overlapped atoms to unique
    my %uniq;
    foreach my $tmp (@mol2D){
	if(defined $uniq{$tmp->[0]}){
	    $uniq{$tmp->[0]}=$uniq{$tmp->[0]} + $tmp->[1];
	}else{
	    $uniq{$tmp->[0]}=$tmp->[1];
	}
    }
    foreach my $key (keys %uniq){
	my @tmp = ($key, $uniq{$key});
	push @moluniq, \@tmp;
    }
    return \@moluniq;
}

#calculation of molecular weight and ratio -> return as a ref of list
sub calc_mol{
    my $refmol2D = shift;
    my @moltmpwght =();
    my $Mwght = 0;
    foreach my $tmp (@$refmol2D){
	my @tmpmol = ();
	push @tmpmol, @$tmp;
	my $Awght = ($peritable->{$tmpmol[0]})->[1];
	$Awght = $Awght * $tmpmol[1];
	push @tmpmol, $Awght;
	$Mwght = $Mwght + $Awght;
	push @moltmpwght, \@tmpmol;
    }
    foreach my $tmp (@moltmpwght){
	my $wpercent = ($tmp->[2] / $Mwght) * 100;
	$wpercent = sprintf("%.4g", $wpercent);
	push @$tmp, $wpercent;
    }
    $" = "\t";
    print "Atom\t","N\t", "sub MW\t", "subtotal\n";
    print "\t\t", "g\/mol\t", "\%\n";
    foreach(@moltmpwght){
	print "@$_\n";
    }
    print "------------------------------------------\n";
    print "Fomula weight\t", $Mwght, "\n";
}

#read atomic average weight file -> return as a ref. of hash
#hash={{atom}=[Atomic Number,Atomic weight]...}
sub read_atom{
    my $filename = \$atomfile;
    my %hash;

    eval{
	open my $fh, '<:encoding(utf8)', $$filename or die $!;
	my $head = <$fh>;
	#catch a header code
	if($head =~ /\[Atom\]/){
	    while(<$fh>){
		chomp;
       		my ($Anum, $Aname, $weight) = split(/,/, $_, 3);
		$hash{$Aname} = [$Anum, $weight];
	    }
	}else{
	    die $!;
	}
	close $fh;
    };
    if($@){
	print "ERROR!: $@";
	exit 1;
    }
#    print scalar(keys(%hash)), "\n";
    return \%hash;
}
