# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::CLI;

use Carp qw{croak};
use Fcm::CLI::Config;
use Fcm::CLI::Exception;
use Fcm::Util::ClassLoader;
use File::Basename qw{basename};
use Getopt::Long qw{GetOptions};
use Scalar::Util qw{blessed};

################################################################################
# Invokes the FCM command line interface
sub invoke {
    local(@ARGV) = @ARGV;
    my $config          = Fcm::CLI::Config->instance();
    my $subcommand_name = @ARGV ? shift(@ARGV) : q{};
    my $subcommand      = $config->get_subcommand_of($subcommand_name);
    eval {
        if (!$subcommand) {
            croak(Fcm::CLI::Exception->new({message => 'unknown command'}));
        }
        my ($opts_ref, $args_ref, $is_help) = _parse_argv_using($subcommand);
        my ($invoker_class, $invoker);
        if ($is_help) {
            $invoker_class
                = _load_invoker_class_of($config->get_subcommand_of(q{}));
            $invoker = $invoker_class->new({
                command   => $subcommand_name,
                arguments => [$subcommand_name],
            });
        }
        else {
            $invoker_class = _load_invoker_class_of($subcommand);
            $invoker = $invoker_class->new({
                command   => $subcommand_name,
                options   => $opts_ref,
                arguments => $args_ref,
                (
                    $subcommand->get_invoker_config()
                    ? %{$subcommand->get_invoker_config()}
                    : ()
                ),
            });
        }
        $invoker->invoke();
    };
    if ($@) {
        if (Fcm::CLI::Exception->caught($@)) {
            die(sprintf(
                qq{%s%s: %s\nType "%s help%s" for usage\n},
                basename($0),
                ($subcommand_name ? qq{ $subcommand_name} : q{}),
                $@->get_message(),
                basename($0),
                defined($subcommand) ? qq{ $subcommand_name} : q{},
            ));
        }
        else {
            die($@);
        }
    }
}

################################################################################
# Parses options in @ARGV using the options settings of a subcommand
sub _parse_argv_using {
    my ($subcommand) = @_;
    my %options = ();
    my $is_help = undef;
    if (($subcommand->get_options())) {
        my $problem = q{};
        local($SIG{__WARN__}) = sub {
            ($problem) = @_;
        };
        my $success = GetOptions(
            \%options,
            (map {$_->get_arg_for_getopt_long()} ($subcommand->get_options())),
        );
        if (!$success) {
            croak(Fcm::CLI::Exception->new({message => sprintf(
                "option parse failed: %s", $problem,
            )}));
        }

        OPTION:
        for my $option ($subcommand->get_options()) {
            if (!exists($options{$option->get_name()})) {
                next OPTION;
            }
            if ($option->is_help()) {
                $is_help = 1;
            }
            if (
                $option->has_arg() == $option->ARRAY_ARG
                && $option->get_delimiter()
            ) {
                $options{$option->get_name()} = [split(
                    $option->get_delimiter(),
                    join(
                        $option->get_delimiter(),
                        @{$options{$option->get_name()}},
                    ),
                )];
            }
        }
    }
    return (\%options, [@ARGV], $is_help);
}

################################################################################
# Loads and returns the invoker class of a subcommand
sub _load_invoker_class_of {
    my ($subcommand) = @_;
    my $invoker_class
        = $subcommand->get_invoker_class() ? $subcommand->get_invoker_class()
        :                                    'Fcm::CLI::Invoker'
        ;
    return Fcm::Util::ClassLoader::load($invoker_class);
}

1;
__END__

=head1 NAME

Fcm::CLI

=head1 SYNOPSIS

    use Fcm::CLI
    Fcm::CLI::invoke();

=head1 DESCRIPTION

Invokes the FCM command line interface.

=head1 FUNCTIONS

=over 4

=item invoke()

Invokes the FCM command line interface.

=back

=head1 TO DO

Move option/argument parsing to L<Fcm::CLI::Invoker|Fcm::CLI::Invoker>?

Use an OO interface?

=head1 SEE ALSO

L<Fcm::CLI::Config|Fcm::CLI::Config>,
L<Fcm::CLI::Invoker|Fcm::CLI::Invoker>,
L<Fcm::CLI::Subcommand|Fcm::CLI::Subcommand>,
L<Fcm::CLI::Option|Fcm::CLI::Option>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
