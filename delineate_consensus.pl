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
my $path = shift @ARGV;
my $opt = shift @ARGV;

my $file = $path =~ s/.*\///r;
my ($file_name, $file_extension) = $file =~ /^(.+)\.([^.]+)$/;

open my $input, "<:utf8", $path or die;
my (@IdSequence_array, $id);
while (<$input>) {
    chomp;
    if (m/\A>/) {
        ($id) = m/\A>(.+)/;
    } else {
        my $a_IdSequence = IdSequence->new($id, uc $_);
        push @IdSequence_array, $a_IdSequence;
    };
    last if eof $input;
};
close $input;

sub na_delineate {
    my $seq = $_[0];
    my %hash = (
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
    );

    my $sols = 1;
    $sols *= $_ for my @entrs = map { 0 + @{$hash{$_}} } my @keys = sort keys %hash;

    my @pairs = map {
        my $i = $_;
        [ map {
            [$keys[$_], $hash{$keys[$_]}[$i % $entrs[$_]]],
            ($i = int($i / $entrs[$_])) x 0
        } 0 .. $#entrs];
    } 0 .. $sols - 1;

    my @res;
    for my $pair (@pairs) {
        my $temp = $seq;
        for my $repl (@{$pair}) {
            $temp =~ s/@{$repl}[0]/@{$repl}[1]/g;
        };
        push @res, $temp;
    };
    return uniq(@res);
};

sub aa_delineate {
    my $seq = $_[0];
#     my @subs = $seq =~ /(?<=\[)\b\w+\b(?=[\]])/g; # match content inside the square brackets
#     my @not_subs = $seq =~ /(?<!\[)\b\w+\b(?![\]])/g; # match content outside the square brackets
    my @subs = $seq =~ /\[[^]]*\]/g;
    my %hash;
    for my $sub (@subs) {
        my @chars = split("", $sub =~ s/[][]//gr);
        $hash{$sub} = [@chars] unless $hash{$sub};
    };

    my $sols = 1;
    $sols *= $_ for my @entrs = map { 0 + @{$hash{$_}} } my @keys = sort keys %hash;

    my @pairs = map {
        my $i = $_;
        [ map {
            [$keys[$_], $hash{$keys[$_]}[$i % $entrs[$_]]],
            ($i = int($i / $entrs[$_])) x 0
        } 0 .. $#entrs];
    } 0 .. $sols - 1;

    my @res;
    for my $pair (@pairs) {
        my $temp = $seq;
        for my $repl (@{$pair}) {
            @{$repl}[0] =~ s/\[/\\\[/g;
            @{$repl}[0] =~ s/\]/\\\]/g;
            $temp =~ s/@{$repl}[0]/@{$repl}[1]/g;
        };
        push @res, $temp;
    };
    return uniq(@res);
};

open my $output, ">:utf8", "delineated_${file_name}.${file_extension}" or die;
for my $a_IdSequence (@IdSequence_array) {
    my @delin;
    if ($opt eq "na") {
        @delin = na_delineate($a_IdSequence->{_sequence});
    } elsif ($opt eq "aa") {
        @delin = aa_delineate($a_IdSequence->{_sequence});
    };
    for (my $i = 0; $i < scalar @delin; $i++) {
        print $output ">" . $a_IdSequence->{_id} . "_" . ($i + 1) . "\n";
        print $output $delin[$i] . "\n";
    };
};

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
#     grep !$seen{$_->[0]}{$_->[1]}++, @_;
#     grep !$seen{"$_->[0],$_->[1]"}++, @_;
};
