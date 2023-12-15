use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

# GetOptions(
#     "input|i=s" => \my $file,
# );
# 
# open my $input, "<:utf8", $file or die;
# my (%id_sequence, $id);
# while (<$input>) {
# 	chomp;
# 	if (m/\A>/) {
# 		($id) = m/\A>(.+)/;
# 	} else {
# 		$id_sequence{$id} = uc $_;
# 	};
# 	last if eof $input;
# };
# close $input;

our $m = 2; # match
our $mm = -1; # mismatch 
our $g = -2; # gap penalty
our @mat; # alignment matrix filled by similarity

# my @keys = sort keys %id_sequence;
# for my $i (0 .. $#keys) {
# 	for my $j ($i + 1 .. $#keys) {
# 		my $merged = get_merged($id_sequence{$keys[$i]}, $id_sequence{$keys[$j]});
# 		if ($merged ne "") {
# 			print ">" . $keys[$i] . "_" . $keys[$j] . "\n" . $merged . "\n";
# 		};
# 	};
# };

sub max_by {
	my $max = shift;
	foreach (@_) {
		$max = $_ if @{$_}[0] > @{$max}[0];
	};
	return $max;
};

sub get_merged {
	my ($seq1, $seq2) = @_;
	if (length($seq1) >= length($seq2)) {
		$seq1 = $seq1;
		$seq2 = $seq2;
	} else {
		my $swap = $seq1;
		$seq1 = $seq2;
		$seq2 = $swap;
	};
	my ($left, $diag, $up) = (1, 2, 3);
	my @lst1 = split("", $seq1);
	my @lst2 = split("", $seq2);
	my $seq1len = length($seq1);
	my $seq2len = length($seq2);
	my $max_score = [0, 0];
	my $pos = [0, 0];
	for my $i (0 .. $seq1len) {
		for my $j (0 .. $seq2len) {
			$mat[$i][$j] = [0, 0];
		};
	};
	for my $i (1 .. $seq1len) {
		for my $j (1 .. $seq2len) {
			my $p = (substr($seq1, $i - 1, 1) eq substr($seq2, $j - 1, 1)) ? $::m : $::mm;
			# scoring function for Smith-Waterman
			my $score = max_by(
				[0, 0],
				[@{$mat[$i - 1][$j]}[0] + $::g, $up],
				[@{$mat[$i][$j - 1]}[0] + $::g, $left],
				[@{$mat[$i - 1][$j - 1]}[0] + $p, $diag]
			);
			$mat[$i][$j] = $score;
			if (@{$score}[0] >= @{$max_score}[0]) {
				$max_score = $score;
				$pos = [$i, $j];
			};
		};
	};
	
	my ($x, $y) = @{$pos};
	my $aln1 = "";
	my $aln2 = "";
	my $cons = "";
	my $i = $seq1len - 1;
	my $j = $seq2len - 1;
	for (my $i = $seq1len - 1; $i >= $x; $i--) {
		$cons .= $lst1[$i];
	};
	my ($m_, $mm_, $g_) = (0, 0, 0);
	while (1) {
		if ($mat[$x][$y][1] == 1) {
			$y -= 1;
			$aln1 .= "-";
			$aln2 .= $lst2[$y];
			$g_ += 1;
		} elsif ($mat[$x][$y][1] == 2) {
			$x -= 1;
			$y -= 1;
			$aln1 .= $lst1[$x];
			$aln2 .= $lst2[$y];
			if ($lst1[$x] eq $lst2[$y]) {
				$m_ += 1;
			} else {
				$mm_ += 1;
			};
		} elsif ($mat[$x][$y][1] == 3) {
			$x -= 1;
			$aln1 .= $lst1[$x];
			$aln2 .= "-";
			$g_ += 1;
		} elsif ($mat[$x][$y][1] == 0) {
			last;
		};
	};
	my $test1 = reverse $aln1;
	my $test2 = reverse $aln2;
	print $test1 . "\n";
	print $test2 . "\n";
	if ($aln1 eq $aln2) {
		$cons .= $aln1;
	} else {
		return "";
	};
	for (my $i = $y; $i >= 0; $i--) {
		$cons .= $lst2[$i];
	};
	return $cons;
};

my $seq2 = "AGGTTG";
my $seq1 = "GGTTGAT";

print get_merged($seq1, $seq2);
