use strict;
use warnings;

my $file_path = shift @ARGV;
my $file = $file_path =~ s/.*\///r;
my ($file_name, $file_extension) = $file =~ /^(.+)\.([^.]+)$/;

open my $input, "<:utf8", $file_path or die;
open my $output, ">:utf8", "singleline_${file_name}.${file_extension}";
while (<$input>) {
    chomp unless m/\A>(.+)/;
    if ($. > 1 && m/\A>(.+)/) {
        $_ =~ s/\A/\n/;
    };
    print $output "$_";
    print $output "\n" if eof;
};
close $input;
close $output;
