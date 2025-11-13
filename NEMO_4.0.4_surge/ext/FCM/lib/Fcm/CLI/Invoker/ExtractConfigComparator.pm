# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::CLI::Invoker::ExtractConfigComparator;
use base qw{Fcm::CLI::Invoker};

use Cwd qw{cwd};
use Fcm::ExtractConfigComparator;
use Fcm::Config;

################################################################################
# Invokes the sub-system
sub invoke {
    my ($self) = @_;
    my ($cfg_file1, $cfg_file2) = $self->get_arguments();
    if (exists($self->get_options()->{verbose})) {
        Fcm::Config->instance()->verbose($self->get_options()->{verbose});
    }

    my $system = Fcm::ExtractConfigComparator->new({
        files => [$cfg_file1, $cfg_file2], wiki => $self->get_options()->{wiki},
    });
    $system->invoke();
}

1;
__END__

=head1 NAME

Fcm::CLI::ExtractInvoker

=head1 SYNOPSIS

    use Fcm::CLI::Invoker::ExtractConfigComparator;
    $invoker = Fcm::CLI::Invoker::ExtractConfigComparator->new({
        command   => $command,
        options   => \%options,
        arguments => $arguments,
    });
    $invoker->invoke();

=head1 DESCRIPTION

This class extends L<Fcm::CLI::Invoker|Fcm::CLI::Invoker> an inherits all its
methods. An object of this class is used to invoke the extract configuration
comparator.

=head1 METHODS

See L<Fcm::CLI::Invoker|Fcm::CLI::Invoker> for a list of inherited methods.

=over 4

=item invoke()

Invokes the extract configuration comparator.

The I<wiki> option is mapped directly to that of the constructor of
L<Fcm::ExtractConfigComparator|Fcm::ExtractConfigComparator> object.

=back

=head1 TO DO

Unit tests.

=head1 SEE ALSO

L<Fcm::ExtractConfigComparator|Fcm::ExtractConfigComparator>,
L<Fcm::CLI::Invoker|Fcm::CLI::Invoker>,
L<Fcm::CLI::Subcommand|Fcm::CLI::Subcommand>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
