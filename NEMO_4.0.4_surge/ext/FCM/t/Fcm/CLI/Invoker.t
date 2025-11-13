#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw{no_plan};

main();

sub main {
    my $class = 'Fcm::CLI::Invoker';
    use_ok($class);
    test_normal($class);
}

################################################################################
# Tests normal usage
sub test_normal {
    my ($class) = @_;
    my $prefix = "normal";
    my %OPTIONS   = (option1 => 1, option2 => 2, option3 => 3);
    my @ARGUMENTS = ('argument 1', 'argument 2');
    my $invoker = $class->new({
        command   => 'command',
        options   => \%OPTIONS,
        arguments => \@ARGUMENTS,
    });
    isa_ok($invoker, $class, $prefix);
    is($invoker->get_command(), 'command', "$prefix: command");
    is_deeply({$invoker->get_options()}, \%OPTIONS, "$prefix: options");
    is_deeply([$invoker->get_arguments()], \@ARGUMENTS, "$prefix: arguments");
    eval {
        $invoker->invoke();
    };
    isa_ok($@, 'Fcm::CLI::Exception', "$prefix: invoke");
}

__END__
