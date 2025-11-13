# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::CLI::Invoker::GUI;
use base qw{Fcm::CLI::Invoker};

use Fcm::Util qw{run_command};

################################################################################
# Invokes the sub-system
sub invoke {
    my ($self) = @_;
    my ($target) = $self->get_arguments();
    run_command(['fcm_gui', ($target ? $target : ())], METHOD => 'exec');
}

1;
__END__

=head1 NAME

Fcm::CLI::GUIInvoker

=head1 SYNOPSIS

    use Fcm::CLI::Invoker::GUI;
    $invoker = Fcm::CLI::Invoker::GUI->new({
        command   => $command,
        options   => \%options,
        arguments => $arguments,
    });
    $invoker->invoke($command, \%options, $target);

=head1 DESCRIPTION

This class extends L<Fcm::CLI::Invoker|Fcm::CLI::Invoker> an inherits all its
methods. An object of this class is used to invoke the FCM GUI.

=head1 METHODS

See L<Fcm::CLI::Invoker|Fcm::CLI::Invoker> for a list of inherited methods.

=over 4

=item invoke()

Invokes the FCM GUI. If a target is specified as argument, it is the initial
working directory of the GUI.

=back

=head1 TO DO

Unit tests.

=head1 SEE ALSO

L<Fcm::CLI::Invoker|Fcm::CLI::Invoker>,
L<Fcm::CLI::Subcommand|Fcm::CLI::Subcommand>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
