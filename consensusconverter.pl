use strict;
use warnings;

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
package IdSequence;

sub new {
    my ($class, $id, $sequence) = @_;
    my $self = {
        _id => $id,
        _sequence => $sequence,
    };
    return bless($self, $class);
};
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
package main;
my $file_path = shift @ARGV;

my $file = $file_path =~ s/.*\///r;
my ($file_name, $file_extension) = $file =~ /^(.+)\.([^.]+)$/;

open my $input, "<:utf8", $file_path or die;
my (@IdSequence_array, $id);
while (<$input>) {
    chomp;
    if (m/\A>/) {
        ($id) = m/\A>(\S+)/;
    } else {
        my $a_IdSequence = IdSequence->new($id, uc $_);
        push @IdSequence_array, $a_IdSequence;
    };
    last if eof $input;
};
close $input;

sub delineate {
    my $sequence = $_[0];
    my %hash = (
        "A" => ["A"],
        "C" => ["C"],
        "G" => ["G"],
        "T" => ["T"],
        "U" => ["U"],
        "R" => ["A","G"],
        "Y" => ["C","T"],
        "S" => ["G","C"],
        "W" => ["A","T"],
        "K" => ["G","T"],
        "M" => ["A","C"],
        "B" => ["C","G","T"],
        "D" => ["A","G","T"],
        "H" => ["A","C","T"],
        "V" => ["A","C","G"],
        "N" => ["A","T","G","C"],
        "." => [""],
        "-" => [""],
    );
    
# A.................Adenine
# C.................Cytosine
# G.................Guanine
# T (or U)..........Thymine (or Uracil)
# R.................A or G
# Y.................C or T
# S.................G or C
# W.................A or T
# K.................G or T
# M.................A or C
# B.................C or G or T
# D.................A or G or T
# H.................A or C or T
# V.................A or C or G
# N.................any base
# . or -............gap
    
    my $solutions = 1;
    $solutions *= $_ for my @entries = map { 0 + @{$hash{$_}} } my @keys = sort keys %hash; # put reverse before sort to have the keys alphabetically
    
    my @pairs = map {
        my $i = $_;
        [ map { # put reverse before map to produce list alphabetically
            [$keys[$_], $hash{$keys[$_]}[$i % $entries[$_]]],
            ($i = int($i / $entries[$_])) x 0
        } 0 .. $#entries];
    } 0 .. $solutions - 1;
    
    my @results;
    my $i = 0;
    for my $pair (@pairs) {
        my $temporary = $sequence;
        for my $replacement (@{$pair}) {
            if (@{$replacement}[0] eq "." || @{$replacement}[0] eq "-") {
                @{$replacement}[0] = "\\" . @{$replacement}[0];
            };
            $temporary =~ s/@{$replacement}[0]/@{$replacement}[1]/g;
        };
        push @results, $temporary;
    };
    return uniq(@results);
};

open my $output, ">:utf8", "delineated_${file_name}.${file_extension}" or die;
for my $a_IdSequence (@IdSequence_array) {
    my @delineated = delineate($a_IdSequence->{_sequence});
    for (my $i = 0; $i < scalar @delineated; $i++) {
        print $output ">" . $a_IdSequence->{_id} . "_${i}" . "\n";
        print $output $delineated[$i] . "\n";
    };
};

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
#     grep !$seen{$_->[0]}{$_->[1]}++, @_;
#     grep !$seen{"$_->[0],$_->[1]"}++, @_;
};
