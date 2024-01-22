#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw{no_plan};

main();

sub main {
    my $class = 'Fcm::Keyword::Exception';
    use_ok($class);
    test_normal($class);
}

################################################################################
# Tests normal usage
sub test_normal {
    my ($class) = @_;
    my $prefix = 'normal';
    my $e = $class->new({message => 'message'});
    isa_ok($e, $class, $prefix);
    is("$e", "$class: message\n", "$prefix: as_string()");
    is($e->get_message(), 'message', "$prefix: get_message()");
}

__END__
