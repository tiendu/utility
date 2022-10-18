use strict;
use warnings;

# Inspired by Abhinav Mishra https://gist.github.com/bibymaths

my $file_path = shift @ARGV;
open my $input, "<:utf8", $file_path or die;
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

my $global = 1;
our $m = 1; # match
our $mm = -1; # mismatch 
our $g = 0; # gap penalty
our @M; # alignment matrix filled by similarity [NW]
our @N; # alignment matrix filled by similarity [SW]

my @keys = sort keys %id_sequence;
for my $i (0 .. $#keys) {
	for my $j ($i + 1 .. $#keys) {
		print "=" x 30 . "\n";
		print $keys[$i] . " and " . $keys[$j] . "\n";
		if ($global) {
			foreach (get_global_alignment($id_sequence{$keys[$i]}, $id_sequence{$keys[$j]})) {
				print $_ . "\n";
			};
		} else {
			foreach (get_local_alignment($id_sequence{$keys[$i]}, $id_sequence{$keys[$j]})) {
				print $_ . "\n";
			};
		};
	};
};

sub p {  
	my ($loc_1, $loc_2) = @_;
	return ($loc_1 eq $loc_2) ? $::m : $::mm; 
};

sub max {
	my $max = shift;
	foreach (@_) {
		$max = $_ if $_ > $max;
	};
	return $max;
};

sub local_similarity {
	my ($s, $t) = @_;
	for my $i (0 .. length($s)) { $N[$i][0] = $::g * $i; }
	for my $j (0 .. length($t)) { $N[0][$j] = $::g * $j; }
	for my $i (1 .. length($s)) {
		for my $j (1 .. length($t)) {
			my $p = p(substr($s, $i - 1, 1), substr($t, $j - 1, 1)); 
			# scoring function for [SW]
			$N[$i][$j] = max(
				0, 
				$N[$i - 1][$j] + $::g,
				$N[$i][$j - 1] + $::g,
				$N[$i - 1][$j - 1] + $p
				);
			};
	};
	return ($N[length($s)][length($t)]);
};

sub global_similarity {
	my ($s, $t) = @_;
	for my $i (0 .. length($s)) { $M[$i][0] = $::g * $i; }
	for my $j (0 .. length($t)) { $M[0][$j] = $::g * $j; }
	for my $i (1 .. length($s)) {
		for my $j (1 .. length($t)) {
			my $p =  p(substr($s, $i - 1, 1), substr($t, $j - 1, 1)); 
			# scoring function for [NW]
			$M[$i][$j] = max(
				$M[$i - 1][$j] + $::g,
				$M[$i][$j - 1] + $::g,
				$M[$i - 1][$j - 1] + $p
				);
		};
	};
	return ($M[length($s)][length($t)]);
};

# Recursively reconstructs best alignment of strings $s and $t  
# using information stored in alignment matrix similarity.
# Returns a list of two strings representing best alignments,
# i.e. $s and $t with gap symbols inserted.
sub get_local_alignment {
	my ($s, $t) = @_;
	my ($i, $j) = (length($s), length($t));
	my ($N) = local_similarity($s, $t);
	return ("-" x $j, $t) if ($i == 0);
	return ($s, "-" x $i) if ($j == 0);
	my ($s_last, $t_last) = (substr($s, -1), substr($t, -1));
	if ($N[$i][$j] == $N[$i - 1][$j - 1] + p($s_last, $t_last)) {
	# case 1: last letters are paired in the best alignment
		my ($sa, $ta) = get_local_alignment(substr($s, 0, -1), substr($t, 0, -1));
		return ($sa . $s_last, $ta . $t_last);
	} elsif ($N[$i][$j] == $N[$i - 1][$j] + $g) {
	# case 2: last letter of the first string is paired with a gap
		my ($sa, $ta) = get_local_alignment(substr($s, 0, -1), $t);
		return ($sa . $s_last, $ta . "-");
	} elsif ($N[$i][$j] == $N[$i][$j - 1] + $g) {
	# case 3: last letter of the second string is paired with a gap
		my ($sa, $ta) = get_local_alignment($s, substr($t, 0, -1));
		return ($sa . "-", $ta . $t_last);
	};
};

sub get_global_alignment {
	my ($s, $t) = @_;
	my ($i, $j) = (length($s), length($t));
	my ($M) = global_similarity($s, $t);
	return ("-" x $j, $t) if ($i == 0);
	return ($s, "-" x $i) if ($j == 0);
	my ($s_last, $t_last) = (substr($s, -1), substr($t, -1));
	if ($M[$i][$j] == $M[$i - 1][$j - 1] + p($s_last, $t_last)) {
	# case 1: last letters are paired in the best alignment
		my ($sa, $ta) = get_global_alignment(substr($s, 0, -1), substr($t, 0, -1));
		return ($sa . $s_last, $ta . $t_last);
	} elsif ($M[$i][$j] == $M[$i - 1][$j] + $g) {
	# case 2: last letter of the first string is paired with a gap
		my ($sa, $ta) = get_global_alignment(substr($s, 0, -1), $t);
		return ($sa . $s_last, $ta . "-");
	} elsif ($M[$i][$j] == $M[$i][$j - 1] + $g) {
	# case 3: last letter of the second string is paired with a gap
		my ($sa, $ta) = get_global_alignment($s, substr($t, 0, -1));
		return ($sa . "-", $ta . $t_last);
	};
};

