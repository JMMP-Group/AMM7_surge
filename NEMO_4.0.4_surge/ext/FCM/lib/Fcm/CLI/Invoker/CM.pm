# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::CLI::Invoker::CM;
use base qw{Fcm::CLI::Invoker};

use Fcm::Cm qw{cli};

################################################################################
# Invokes the sub-system
sub invoke {
    my ($self) = @_;
    return cli($self->get_command(), @ARGV);
}

1;
__END__

=head1 NAME

Fcm::CLI::Invoker::CM

=head1 SYNOPSIS

    use Fcm::CLI::Invoker::CM;
    $invoker = Fcm::CLI::Invoker::CM->new();
    $invoker->invoke();

=head1 DESCRIPTION

This class extends L<Fcm::CLI::Invoker|Fcm::CLI::Invoker> an inherits all its
methods. An object of this class is used to invoke a command in the CM
sub-system.

It is worth noting that this is not yet a full implementation.

=head1 METHODS

See L<Fcm::CLI::Invoker|Fcm::CLI::Invoker> for a list of inherited methods.

=over 4

=item invoke()

Invokes a command in the CM sub-system.

=back

=head1 TO DO

Bring the CM system into this framework.

Unit tests.

=head1 SEE ALSO

L<Fcm::Cm|Fcm::Cm>,
L<Fcm::CLI::Subcommand|Fcm::CLI::Subcommand>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
