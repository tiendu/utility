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

my @files = @ARGV;
for my $file (@files) {
    open my $input, "<:utf8", $file or die;
    my (%id_sequence, $id);
    while (<$input>) {
        chomp;
        if (m/\A>/) {
            ($id) = m/\A>(.+)/;
        } else {
            $id_sequence{$id} = uc $_;
        };
        last if eof $input;
    };
    close $input;
    $file =~ s/.*\///; 
    my ($name, $extension) = $file =~ /^(.+)\.([^.]+)$/;
    open my $output, ">:utf8", "translated_${name}.faa";
    for my $id (sort keys %id_sequence) {
        for my $i (0 .. 2) {
            my $fw_trans = translate($id_sequence{$id}, $i);
            my $rv_trans = translate(reverse_complement($id_sequence{$id}), $i);
            print $output ">${id}_F" . ($i + 1) . "\n" . $fw_trans . "\n" if $fw_trans ne "";
            print $output ">${id}_R" . ($i + 1) . "\n" . $rv_trans . "\n" if $rv_trans ne "";
        };
    };
    close $output;
};

sub reverse_complement {
    my $seq = $_[0];
    $seq =~ tr/ATGC/TACG/;
    return reverse($seq);
};

sub translate {
    my $seq = $_[0];
    my $step = $_[1];
    my @arr;
    for (my $i = 0; $i + $step <= length($seq); $i += 3) {
        my $codon = substr($seq, $i + $step, 3);
        push @arr, $codons{$codon} if $codons{$codon};
    };
    my ($stop_idx) = grep {$arr[$_] eq "*"} 0 .. $#arr;
    if (!defined($stop_idx) || $stop_idx + 1 == scalar(@arr)) {
        return join("", @arr);
    } else {
        return "";
    };
};
