use strict;
use warnings;

my %codons = (
    "TTT" => "F", "TTC" => "F",
    "TTA" => "L", "TTG" => "L", "CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
    "ATT" => "I", "ATC" => "I", "ATA" => "I",
    "ATG" => "M",
    "GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
    "TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S", "AGT" => "S", "AGC" => "S",
    "CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
    "ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
    "GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
    "TAT" => "Y", "TAC" => "Y",
    "TAA" => "*", "TAG" => "*", "TGA" => "*",
    "CAT" => "H", "CAC" => "H",
    "CAA" => "Q", "CAG" => "Q",
    "AAT" => "N", "AAC" => "N",
    "AAA" => "K", "AAG" => "K",
    "GAT" => "D", "GAC" => "D",
    "GAA" => "E", "GAG" => "E",
    "TGT" => "C", "TGC" => "C",
    "TGG" => "W",
    "CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R", "AGA" => "R", "AGG" => "R",
    "GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G",
);

my $path = shift @ARGV;
my $file = $path =~ s/.*\///r;
my ($file_name, $file_extension) = $file =~ /^(.+)\.([^.]+)$/;

open my $input, "<:utf8", $path or die;
my (%id_sequence, $id);
while (<$input>) {
    chomp;
    if (m/\A>/) {
        ($id) = m/\A>(\S+)/;
    } else {
        $id_sequence{$id} = uc $_;
    };
    last if eof $input;
};
close $input;

open my $output, ">:utf8", "translated_${file_name}.faa";
for my $id (sort keys %id_sequence) {
    for my $i (0 .. 2) {
        my $translated = translate($id_sequence{$id}, $i);
        my $reverse_translated = translate(reverse_complement($id_sequence{$id}), $i);
        print $output ">${id} F" . ($i + 1) . "\n" . $translated . "\n" if $translated ne "";
        print $output ">${id} R" . ($i + 1) . "\n" . $translated . "\n" if $translated ne "";
    };
};
close $output;

sub reverse_complement {
    my $seq = $_[0];
    $seq =~ tr/ATGC/TACG/;
    return reverse($seq);
};

sub translate {
    my $seq = $_[0];
    my $step = $_[1];
    my @seq_arr;
    my $idx = 0;
    for (my $i = 0; $i + $step <= length($seq); $i += 3) {
        my $codon = substr($seq, $i + $step, 3);
        push @seq_arr, $codons{$codon} if $codons{$codon};
    };
    $idx++ if (grep {/\*/} @seq_arr);
    if ($idx == 0 || $idx == scalar @seq_arr) {
        return join("", @seq_arr);
    } else {
        return "";
    };
};
