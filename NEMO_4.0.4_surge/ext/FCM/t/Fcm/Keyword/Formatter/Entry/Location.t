#!/usr/bin/perl

use strict;
use warnings;

use Fcm::Keyword::Entry::Location;
use Test::More qw{no_plan};

main();

sub main {
    my $class = 'Fcm::Keyword::Formatter::Entry::Location';
    use_ok($class);
    test_normal($class);
}

################################################################################
# Tests normal usage
sub test_normal {
    my ($class) = @_;
    my $prefix = 'normal';
    my $formatter = $class->new();
    isa_ok($formatter, $class, $prefix);
    my $entry = Fcm::Keyword::Entry::Location->new({key => 'k', value => 'v'});
    like($formatter->format($entry), qr{k \s = \s v}xms, "$prefix: format");
}

__END__
