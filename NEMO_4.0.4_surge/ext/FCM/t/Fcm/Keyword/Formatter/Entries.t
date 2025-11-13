#!/usr/bin/perl

use strict;
use warnings;

use Fcm::Keyword::Entries;
use Test::More qw{no_plan};

main();

sub main {
    my $class = 'Fcm::Keyword::Formatter::Entries';
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
    my $entries = Fcm::Keyword::Entries->new();
    $entries->add_entry('foo', 'food');
    $entries->add_entry('bar', 'barley');
    is($formatter->format($entries), "BAR = barley\nFOO = food\n",
        "$prefix: format");
}

__END__
