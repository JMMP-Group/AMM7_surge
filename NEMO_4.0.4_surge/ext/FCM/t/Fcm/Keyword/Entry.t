#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw{no_plan};

main();

sub main {
    my $class = 'Fcm::Keyword::Entry';
    use_ok($class);
    test_normal($class);
}

################################################################################
# Tests normal usage
sub test_normal {
    my ($class) = 'Fcm::Keyword::Entry';
    my $entry = $class->new({key => 'key', value => 'value'});
    isa_ok($entry, $class);
    is($entry->get_key(), 'key', "normal: get_key()");
    is($entry->get_value(), 'value', "normal: get_value()");
}

__END__
