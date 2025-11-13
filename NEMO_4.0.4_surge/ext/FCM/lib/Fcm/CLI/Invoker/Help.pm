# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::CLI::Invoker::Help;
use base qw{Fcm::CLI::Invoker};

use Carp qw{croak};
use Fcm::CLI::Exception;
use Fcm::CLI::Config;
use Fcm::Config;
use Fcm::Util qw{run_command};
use IO::File;

################################################################################
# Invokes the sub-system
sub invoke {
    my ($self) = @_;
    my @subcommand_names = $self->get_arguments();
    if (@subcommand_names) {
        for my $subcommand_name (@subcommand_names) {
            my $help_string = $self->_get_help_for($subcommand_name);
            if (!defined($help_string)) {
                croak(Fcm::CLI::Exception->new({message => sprintf(
                    "%s: unknown command", $subcommand_name,
                )}));
            }
            print($help_string, "\n");
        }
    }
    else {
        print($self->_get_help());
    }
}

################################################################################
# Returns the help string for a subcommand matching $subcommand_name
sub _get_help_for {
    my ($self, $subcommand_name) = @_;
    my $subcommand
        = Fcm::CLI::Config->instance()->get_subcommand_of($subcommand_name);
    if (!$subcommand) {
        return;
    }
    if ($subcommand->is_vc()) {
        my $invoker = $subcommand->get_invoker($subcommand_name);
        local(@ARGV) = '--help';
        $invoker->invoke();
        return q{};
    }
    my $prog = Fcm::Config->instance()->setting('FCM_COMMAND');
    # FIXME: can do with using Text::Template or Perl6::Form
    my $help = sprintf(
        "%s %s: %s\n",
        $prog,
        $subcommand->as_string(),
        $subcommand->get_synopsis(),
    );
    $help .= sprintf(
        "usage: %s %s %s\n",
        $prog, $subcommand->get_names()->[0], $subcommand->get_usage(),
    );
    if ($subcommand->get_description()) {
        my @lines = (q{}, split("\n", $subcommand->get_description()), q{});
        $help .= join(qq{\n  }, @lines) . "\n";
    }
    if ($subcommand->get_options()) {
        $help .= "Valid options:\n";
        my $max_length_of_name = 0;
        my @option_names;
        for my $option ($subcommand->get_options()) {
            if (length($option->get_name()) > $max_length_of_name) {
                $max_length_of_name = length($option->get_name());
            }
        }
        for my $option ($subcommand->get_options()) {
            $help .= sprintf(
                "  --%s%s%s%s : %s\n",
                $option->get_name(),
                q{ } x ($max_length_of_name - length($option->get_name())),
                (
                    $option->get_letter()
                    ? q{ [-} . $option->get_letter() . q{]} : q{     }
                ),
                ($option->has_arg() ? q{ arg} : q{ } x 4),
                $option->get_description(),
            );
        }
    }
    return $help;
}

################################################################################
# Returns the general help string
sub _get_help {
    my ($self) = @_;
    my $release = $self->_get_release();

    # FIXME: can do with using Text::Template or Perl6::Form
    my $prog = Fcm::Config->instance()->setting('FCM_COMMAND');
    my $return = sprintf(
          qq{usage: %s <subcommand> [options] [args]\n}
        . qq{Flexible configuration management system, release %s.\n}
        . qq{Type "%s help <subcommand>" for help on a specific subcommand\n}
        . qq{\n}
        . qq{Available subcommands:\n}
        ,
        $prog, $release, $prog,
    );
    for my $subcommand (Fcm::CLI::Config->instance()->get_core_subcommands()) {
        $return .= sprintf(qq{  %s\n}, $subcommand->as_string());
    }

    my @lines = run_command(
        [qw/svn help/], DEVNULL => 1, METHOD => 'qx', ERROR => 'ignore',
    );
    if (@lines) {
        for my $subcommand (Fcm::CLI::Config->instance()->get_vc_subcommands()) {
            if (defined($subcommand->get_synopsis())) {
                $return .= sprintf(qq{  %s\n}, $subcommand->as_string());
            }
            else {
                $return .= qq{  <version control system commands, see below>\n};
            }
        }
        $return .= "\n=> svn help\n". join(q{}, @lines);
    }
    return $return;
}

################################################################################
# Returns the release number of the current program
sub _get_release {
    my ($self) = @_;
    my $release  = Fcm::Config->instance()->setting('FCM_RELEASE');
    my $rev_file = Fcm::Config->instance()->setting('FCM_REV_FILE');
    if (-r $rev_file) {
        my $handle = IO::File->new($rev_file, 'r');
        if ($handle) {
            my $rev = $handle->getline();
            $handle->close();
            chomp($rev);
            if ($rev) {
                $release .= qq{ (r$rev)};
            }
        }
    }
    return $release;
}

1;
__END__

=head1 NAME

Fcm::CLI::Invoker::Help

=head1 SYNOPSIS

    use Fcm::CLI::Invoker::Help;
    $invoker = Fcm::CLI::Invoker::Help->new({
        command   => $command,
        options   => \%options,
        arguments => $arguments,
    });
    $invoker->invoke();

=head1 DESCRIPTION

This class extends L<Fcm::CLI::Invoker|Fcm::CLI::Invoker> an inherits all its
methods. An object of this class is used to provide help on the command line
interface.

=head1 METHODS

See L<Fcm::CLI::Invoker|Fcm::CLI::Invoker> for a list of inherited methods.

=over 4

=item invoke()

Provides help. If a subcommand name is specified in the argument, provides help
for the specified subcommand. If a subcommand name is not specified, provides
general CLI help.

=back

=head1 DIAGNOSTICS

=over 4

=item L<Fcm::CLI::Exception|Fcm::CLI::Exception>

The invoke() method can croak() with this exception if the specified subcommand
cannot be identified.

=back

=head1 TO DO

Unit tests.

Separate logic in this module with that of L<Fcm::CLI::Config|Fcm::CLI::Config>.

Decouples help formatter with this invoker.

=head1 SEE ALSO

L<Fcm::CLI::Exception|Fcm::CLI::Exception>,
L<Fcm::CLI::Subcommand|Fcm::CLI::Subcommand>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
