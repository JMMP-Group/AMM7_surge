#!/usr/bin/perl

use strict;
use warnings;

################################################################################
# A mock Fcm::ConfigSystem object
{
    package MockConfigSystem;
    use base qw{Fcm::ConfigSystem};

    our $LATEST_INVOKED_INSTANCE;

    ############################################################################
    # Returns the arguments to the last invoke() call
    sub get_invoke_args {
        my ($self) = @_;
        return $self->{invoke_args};
    }

    ############################################################################
    # Does nothing but captures the arguments
    sub invoke {
        my ($self, %args) = @_;
        $LATEST_INVOKED_INSTANCE = $self;
        $self->{invoke_args} = \%args;
        return 1;
    }
}

use Cwd;
use Test::More qw{no_plan};

main();

sub main {
    my $class = 'Fcm::CLI::Invoker::ConfigSystem';
    use_ok($class);
    test_invoke($class);
}

################################################################################
# Tests normal usage of invoke()
sub test_invoke {
    my ($class) = @_;
    my $prefix = "invoke";
    my %TEST = (
        test1 => {
            command            => 'pig',
            options            => {'egg' => 1},
            arguments          => ['bacon'],
            expected_options   => {FOO => undef, BAR_BAZ => undef, EGGS => 1},
            expected_arguments => 'bacon',
        },
        test2 => {
            command            => 'pig',
            options            => {'foo' => 1, 'bar-baz' => 1},
            arguments          => [],
            expected_options   => {FOO => 1, BAR_BAZ => 1, EGGS => undef},
            expected_arguments => cwd(),
        }
    );
    for my $key (keys(%TEST)) {
        my $invoker = $class->new({
            command            => $TEST{$key}{command},
            options            => $TEST{$key}{options},
            arguments          => $TEST{$key}{arguments},
            impl_class         => 'MockConfigSystem',
            cli2invoke_key_map => {
                'foo' => 'FOO', 'bar-baz' => 'BAR_BAZ', 'egg' => 'EGGS',
            },
        });
        isa_ok($invoker, 'Fcm::CLI::Invoker::ConfigSystem', "$prefix: $key");
        $invoker->invoke();
        my $config_system_instance = $MockConfigSystem::LATEST_INVOKED_INSTANCE;
        isa_ok(
            $config_system_instance,
            'Fcm::ConfigSystem',
            "$prefix: $key: Fcm::ConfigSystem",
        );
        is(
            $config_system_instance->cfg()->src(),
            $TEST{$key}{expected_arguments},
            "$prefix: $key: cfg()->src()",
        );
        is_deeply(
            $config_system_instance->get_invoke_args(),
            $TEST{$key}{expected_options},
            "$prefix: $key: invoke args",
        );
    }
}

__END__
