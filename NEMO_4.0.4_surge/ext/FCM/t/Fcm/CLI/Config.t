#!/usr/bin/perl

use strict;
use warnings;

use Fcm::CLI::Config::Default;
use Fcm::CLI::Subcommand;
use Test::More qw{no_plan};

main();

sub main {
    my $class = 'Fcm::CLI::Config';
    use_ok($class);
    test_get_instance($class);
    test_get_subcommand_of_string($class);
}

################################################################################
# Tests normal usage of getting an instance
sub test_get_instance {
    my ($class) = @_;
    my $prefix = 'constructor';
    my $cli_config = $class->instance();
    isa_ok($cli_config, $class, "$prefix");
    is_deeply(
        [$cli_config->get_core_subcommands()],
        \@Fcm::CLI::Config::Default::CORE_SUBCOMMANDS,
        "$prefix: default core",
    );
    is_deeply(
        [$cli_config->get_vc_subcommands()],
        \@Fcm::CLI::Config::Default::VC_SUBCOMMANDS,
        "$prefix: default vc",
    );
    is_deeply(
        [$cli_config->get_subcommands()],
        [$cli_config->get_core_subcommands(), $cli_config->get_vc_subcommands()],
        "$prefix: default",
    );
    is($class->instance(), $cli_config, "$prefix: same instance");
    isnt($class->instance({}), $cli_config, "$prefix: not the same instance");
    my $empty_cli_config = $class->instance({
        core_subcommands => [],
        vc_subcommands   => [],
    });
    is_deeply(
        [$empty_cli_config->get_core_subcommands()],
        [],
        "$prefix: empty core",
    );
    is_deeply(
        [$empty_cli_config->get_vc_subcommands()],
        [],
        "$prefix: empty vc",
    );
    is_deeply(
        [$empty_cli_config->get_subcommands()],
        [],
        "$prefix: empty",
    );
}

################################################################################
# Tests getting a subcommand of a matching string
sub test_get_subcommand_of_string {
    my ($class) = @_;
    my $prefix = 'get_subcommand_of';
    my $foo_subcommand = Fcm::CLI::Subcommand->new({names => ['food', 'foo']});
    my $bar_subcommand = Fcm::CLI::Subcommand->new({names => ['barley', 'bar']});
    my $cli_config = $class->instance({
        core_subcommands => [$foo_subcommand, $bar_subcommand],
        vc_subcommands   => [],
    });
    for my $key ('food', 'foo') {
        is($cli_config->get_subcommand_of($key), $foo_subcommand,
            "$prefix: $key");
    }
    for my $key ('barley', 'bar') {
        is($cli_config->get_subcommand_of($key), $bar_subcommand,
            "$prefix: $key");
    }
    is($cli_config->get_subcommand_of('baz'), undef, "$prefix: baz");
}

__END__
