use strict;

my $file_path = shift @ARGV;
my $threshold = shift @ARGV;
my $file_path_cp = $file_path;
$file_path_cp =~ s/.*\///;
my ($file_name, $file_extension) = $file_path_cp =~ /^(.+)\.([^.]+)$/;

open my $input, "<:utf8", $file_path or die;
my (%id_sequence_hash, $stored_id, $flag);
while (<$input>) {
    chomp;
    if ($_ =~ m/\A>(.+)/) {
        my @fields = split " ";
        my $id = $fields[0];
        $id =~ s/>//;
        $stored_id = $id;
        $flag = 1;
    } elsif ($flag) {
        $id_sequence_hash{$stored_id} = uc $_;
        $flag = 0;
    };
    last if eof $input;
};
close $input;

our @tetranucleotides;
foreach (kmer_generator(4)) {
    push @tetranucleotides, $_ if palindrome(\$_);
};

my %id_kmer_normusage;
for my $id (sort keys %id_sequence_hash) {
    my $no_repeat_sequence = repeat_modifier($id_sequence_hash{$id}, 4, 0);
    my $kmer_counts = kmer_count_generator($no_repeat_sequence, 4);
    my $base_frequencies = base_frequency_generator($no_repeat_sequence);
    foreach (@tetranucleotides) {
        my $expected_count = 1;
        my $normalised_usage;
        my $kmer_base_counts = base_count_generator($_);
        for my $base ("A", "T", "G", "C") {
            $expected_count *= ($base_frequencies->{$base} ** $kmer_base_counts->{$base});
        };
        $expected_count *= (length($no_repeat_sequence) - 4 + 1);
        $normalised_usage = $kmer_counts->{$_} / $expected_count;
        $id_kmer_normusage{$id}{$_} = $normalised_usage;
    };
};

my %result = agglomerative_clustering(\%id_kmer_normusage, $threshold);

my $i = 1;
foreach (sort keys %result) {
    my @ids = split ",";
    open my $output, ">:utf8", "cluster_${i}_${file_name}.${file_extension}" or die;
    for my $id (@ids) {
        print $output ">${id}\n$id_sequence_hash{$id}\n";
    };
    $i++;
    close $output;
};

sub agglomerative_clustering {
    my %data = %{$_[0]};
    my $threshold = $_[1];
    my $size = keys %data;
    my %clusters;
    for (my $i = 1; $i < $size; $i++) {
        my (%distances, $find, @keys);
        @keys = sort keys %data;
        for my $index_1 (0 .. $#keys) {
            for my $index_2 (1 + $index_1 .. $#keys) {
                my ($distance, $key_1, $key_2) = (0, $keys[$index_1], $keys[$index_2]);
                $distance += ($data{$key_1}{$_} - $data{$key_2}{$_}) ** 2 foreach @::tetranucleotides;
                $distance = sqrt($distance);
                $distances{$key_1}{$key_2} = $distance;
                $find->{min} = $distance unless $find->{min};
                $find->{key} = [$key_1, $key_2] unless $find->{key};
                if ($find->{min} > $distance) {
                    $find->{min} = $distance;
                    $find->{key} = [$key_1, $key_2];
                };
            };
        };
        my ($key_1, $key_2) = $find->{key}->@*;
        $data{"$key_1,$key_2"}{$_} = mean($data{$key_1}{$_}, $data{$key_2}{$_}) foreach @::tetranucleotides;
        delete @data{($key_1, $key_2)};
        last if $find->{min} >= $threshold;
        %clusters = %data;
    };
    return %clusters;
};

sub mean {
    my @data = @_;
    my $sum;
    $sum += $_ for @data;    
    return $sum / @data;
};

sub palindrome {
    my $word = $_;
    my $reversed_word = reverse $word;
    if ($word eq $reversed_word) {
        return 1;
    } else {
        return 0;
    };
};

sub repeat_modifier {
    my ($string, $threshold, $change) = @_;
    $string =~ s/(.)\1{$threshold,}/$1 x $change/ge;
    return $string;
};

sub kmer_generator {
    my ($k) = @_;
    my @bases_1 = ("A", "T", "G", "C");
    my @bases_2 = @bases_1;
    for (my $i = 1; $i < $k; $i++) {
        my @temporary;
        foreach my $base_1 (@bases_1) {
            foreach my $base_2 (@bases_2) {
                push @temporary, "$base_1" . "$base_2";
            };
        };
        undef @bases_2;
        @bases_2 = @temporary;
    };
    return @bases_2;
};

sub base_count_generator {
    my $sequence = $_[0];
    my @base_list = ("A", "T", "G", "C");
    my $base;
    map {$base->{$_} = 0} @base_list;
    for (my $i = 0; $i <= length($sequence); $i++) {
        $base->{substr($sequence, $i, 1)}++ if exists $base->{substr($sequence, $i, 1)};
    };
    return $base;
};

sub base_frequency_generator {
    my $sequence = $_[0];
    my @base_list = ("A", "T", "G", "C");
    my $base;
    map {$base->{$_} = 0} @base_list;
    for (my $i = 0; $i <= length($sequence); $i++) {
        $base->{substr($sequence, $i, 1)}++ if exists $base->{substr($sequence, $i, 1)};
    };
    my $total = 0;
    $total += $base->{$_} for keys %{$base};
    $base->{$_} /= $total for keys %{$base};
    return $base;
};

sub kmer_count_generator {
    my $sequence = $_[0];
    my $k = $_[1];
    my @kmers_list = kmer_generator($k);
    my $kmers;
    map {$kmers->{$_} = 0} @kmers_list;
    for (my $i = 0; $i <= length($sequence) - $k; $i++) {
        $kmers->{substr($sequence, $i, $k)}++ if exists $kmers->{substr($sequence, $i, $k)};
    };
    return $kmers;
};

sub kmer_frequency_generator {
    my $sequence = $_[0];
    my $k = $_[1];
    my @kmers_list = kmer_generator($k);
    my $kmers;
    map {$kmers->{$_} = 0} @kmers_list;
    for (my $i = 0; $i <= length($sequence) - $k; $i++) {
        $kmers->{substr($sequence, $i, $k)}++ if exists $kmers->{substr($sequence, $i, $k)};
    };
    my $total = 0;
    $total += $kmers->{$_} for keys %{$kmers};
    $kmers->{$_} /= $total for keys %{$kmers};
    return $kmers;
};
