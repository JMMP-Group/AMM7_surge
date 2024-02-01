#!/usr/bin/perl

use strict;
use warnings;

use Fcm::Keyword::Entry;
use Test::More qw{no_plan};

main();

sub main {
    my $class = 'Fcm::Keyword::Formatter::Entry';
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
    my $entry = Fcm::Keyword::Entry->new({key => 'k', value => 'v'});
    is($formatter->format($entry), "k = v\n", "$prefix: format");
}

__END__
