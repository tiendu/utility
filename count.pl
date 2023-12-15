use strict;

my $file_path = shift @ARGV;
my $file = $file_path =~ s/.*\///r;
my ($file_name, $file_extension) = $file =~ /^(.+)\.([^.]+)$/;

open my $file, "<:utf8", $file_path or die;
my %counter;
while (<$file>) {
    chomp;
    $counter{$_}++;
};
close $file;

open my $count, ">:utf8", "count_${file_name}.${file_extension}";
for my $key (sort keys %counter) {
        print $count join("\t", $key, $counter{$key}, "\n");
};
close $count;
