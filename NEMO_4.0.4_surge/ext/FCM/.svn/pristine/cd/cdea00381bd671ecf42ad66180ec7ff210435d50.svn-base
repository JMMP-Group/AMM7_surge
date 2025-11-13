#!/usr/bin/perl

use strict;
use warnings;


################################################################################
# A sub-class of Fcm::CLI::Invoker for testing
{
    package TestInvoker;
    use base qw{Fcm::CLI::Invoker};

    our $LATEST_INSTANCE;

    ############################################################################
    # Returns a test attrib
    sub get_test_attrib {
        my ($self) = @_;
        return $self->{test_attrib};
    }

    ############################################################################
    # Invokes the sub-system
    sub invoke {
        my ($self) = @_;
        $LATEST_INSTANCE = $self;
    }
}

use Fcm::CLI::Config;
use Fcm::CLI::Subcommand;
use Test::More (tests => 25);

main();

sub main {
    use_ok('Fcm::CLI');
    test_invalid_subcommand();
    test_invoker_not_implemented();
    test_normal_invoke();
    test_help_invoke();
    test_get_invoker_normal();
    test_load_invoker_class();
}

################################################################################
# Tests to ensure that an invalid subcommand results in an exception
sub test_invalid_subcommand {
    Fcm::CLI::Config->instance({core_subcommands => [], vc_subcommands => []});
    eval {
        local(@ARGV) = ('foo');
        Fcm::CLI::invoke();
    };
    like($@, qr{foo: unknown command}, 'invalid subcommand');
}

################################################################################
# Tests to ensure that an unimplemented invoker results in an exception
sub test_invoker_not_implemented {
    Fcm::CLI::Config->instance({
        core_subcommands => [
            Fcm::CLI::Subcommand->new({names => ['foo']}),
            Fcm::CLI::Subcommand->new({
                names => ['bar'], invoker_class => 'barley',
            }),
        ],
        vc_subcommands   => [],
    });
    eval {
        local(@ARGV) = ('foo');
        Fcm::CLI::invoke();
    };
    like($@, qr{foo: \s command \s not \s implemented}xms, 'not implemented');
    eval {
        local(@ARGV) = ('bar');
        Fcm::CLI::invoke();
    };
    like($@, qr{barley: \s class \s loading \s failed}xms, 'not implemented');
}

################################################################################
# Tests normal usage of invoke
sub test_normal_invoke {
    my $prefix = "normal invoke";
    Fcm::CLI::Config->instance({
        core_subcommands => [
            Fcm::CLI::Subcommand->new({
                names          => ['foo'],
                invoker_class  => 'TestInvoker',
                invoker_config => {test_attrib => 'test_attrib value'},
            }),
        ],
        vc_subcommands   => [],
    });
    ok(!$TestInvoker::LATEST_INSTANCE, "$prefix: invoker not called");
    local(@ARGV) = ('foo', 'bar', 'baz');
    Fcm::CLI::invoke();
    my $invoker = $TestInvoker::LATEST_INSTANCE;
    if (!$invoker) {
        fail($prefix);
    }
    else {
        is($invoker->get_command(), 'foo', "$prefix: invoker command");
        is_deeply({$invoker->get_options()}, {}, "$prefix: invoker options");
        is_deeply([$invoker->get_arguments()], ['bar', 'baz'],
            "$prefix: invoker arguments");
        is($invoker->get_test_attrib(), 'test_attrib value',
            "$prefix: invoker test attrib");
    }
    $TestInvoker::LATEST_INSTANCE = undef;
}

################################################################################
# Tests help usage of invoke
sub test_help_invoke {
    my $prefix = "help invoke";
    Fcm::CLI::Config->instance({
        core_subcommands => [
            Fcm::CLI::Subcommand->new({
                names          => ['foo'],
                invoker_class  => 'TestInvoker',
                invoker_config => {test_attrib => 'test_attrib value normal'},
                options        => [
                    Fcm::CLI::Option->new({name => 'foo', is_help => 1}),
                ],
            }),
            Fcm::CLI::Subcommand->new({
                names          => [q{}],
                invoker_class  => 'TestInvoker',
            }),
        ],
        vc_subcommands   => [],
    });
    ok(!$TestInvoker::LATEST_INSTANCE, "$prefix: invoker not called");
    local(@ARGV) = ('foo', '--foo');
    Fcm::CLI::invoke();
    my $invoker = $TestInvoker::LATEST_INSTANCE;
    if (!$invoker) {
        fail($prefix);
    }
    else {
        is_deeply([$invoker->get_arguments()], ['foo'],
            "$prefix: invoker argument");
    }
    $TestInvoker::LATEST_INSTANCE = undef;
}

################################################################################
# Tests getting an invoker
sub test_get_invoker_normal {
    my $prefix = 'get invoker normal';
    my @options = (
        Fcm::CLI::Option->new({name => 'foo'}),
        Fcm::CLI::Option->new({name => 'bar'}),
        Fcm::CLI::Option->new({name => 'baz'}),
        Fcm::CLI::Option->new({
            name      => q{pork},
            delimiter => q{,},
            has_arg   => Fcm::CLI::Option->ARRAY_ARG,
        }),
    );
    my $subcommand = Fcm::CLI::Subcommand->new({options => \@options});
    my %TEST = (
        test1 => {
            argv      => ['--foo', '--bar', 'egg', 'ham', 'sausage'],
            command   => 'command',
            options   => {foo => 1, bar => 1},
            arguments => ['egg', 'ham', 'sausage'],
        },
        test2 => {
            argv      => ['--baz', '--foo', '--bar'],
            command   => 'test',
            options   => {foo => 1, bar => 1, baz => 1},
            arguments => [],
        },
        test3 => {
            argv      => ['egg', 'ham', 'sausage'],
            command   => 'meal',
            options   => {},
            arguments => ['egg', 'ham', 'sausage'],
        },
        test4 => {
            argv      => ['--pork', 'ham', '--pork', 'sausage'],
            command   => 'pig',
            options   => {pork => ['ham', 'sausage']},
            arguments => [],
        },
        test5 => {
            argv      => ['--pork', 'ham,sausage', '--pork', 'bacon', 'liver'],
            command   => 'pig',
            options   => {pork => ['ham', 'sausage', 'bacon']},
            arguments => ['liver'],
        },
    );
    for my $key (keys(%TEST)) {
        local(@ARGV) = @{$TEST{$key}{argv}};
        my ($opts_ref, $args_ref) = Fcm::CLI::_parse_argv_using($subcommand);
        is_deeply($opts_ref, $TEST{$key}{options},
            "$prefix $key: get options");
        is_deeply($args_ref, $TEST{$key}{arguments},
            "$prefix $key: get arguments");
    }
    my %BAD_TEST = (
        test1 => {
            argv => ['--egg', '--bar', 'foo', 'ham', 'sausage'],
        },
        test2 => {
            argv => ['--foo=egg'],
        },
    );
    for my $key (keys(%BAD_TEST)) {
        local(@ARGV) = @{$BAD_TEST{$key}{argv}};
        eval {
            Fcm::CLI::_parse_argv_using($subcommand);
        };
        isa_ok($@, 'Fcm::CLI::Exception', "$prefix $key");
    }
}

################################################################################
# Tests loading an invoker with a different class
sub test_load_invoker_class {
    my $prefix = 'get invoker class';
    eval {
        my $subcommand = Fcm::CLI::Subcommand->new({invoker_class => 'foo'});
        Fcm::CLI::_load_invoker_class_of($subcommand);
    };
    isa_ok($@, 'Fcm::Exception', "$prefix");

    my $invoker_class = 'Fcm::CLI::Invoker::ConfigSystem';
    my $subcommand
        = Fcm::CLI::Subcommand->new({invoker_class => $invoker_class});
    my $class = Fcm::CLI::_load_invoker_class_of($subcommand);
    is($class, $invoker_class, "$prefix: $invoker_class");
}

__END__
