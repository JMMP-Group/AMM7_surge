# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::CLI::Invoker::CfgPrinter;
use base qw{Fcm::CLI::Invoker};

use Carp qw{croak};
use Fcm::Exception;
use Fcm::CfgFile;
use Fcm::Config;

################################################################################
# Invokes the sub-system
sub invoke {
    my ($self) = @_;
    my ($cfg_file) = $self->get_arguments();
    if (!$cfg_file) {
        croak(Fcm::CLI::Exception->new({message => 'no CFGFILE specified'}));
    }
    my $cfg = Fcm::CfgFile->new(SRC => $cfg_file);
    Fcm::Config->instance()->verbose(0); # suppress message printing to STDOUT
    my $read = $cfg->read_cfg();
    if (!$read) {
        croak(Fcm::Exception->new({message => sprintf(
            "% :cannot read", $cfg_file,
        )}));
    }
    $cfg->print_cfg($self->get_options()->{output});
 }

1;
__END__

=head1 NAME

Fcm::CLI::Invoker::CfgPrinter

=head1 SYNOPSIS

    use Fcm::CLI::Invoker::CfgPrinter;
    $invoker = Fcm::CLI::Invoker::CfgPrinter->new({
        command   => $command,
        options   => \%options,
        arguments => $arguments,
    });
    $invoker->invoke();

=head1 DESCRIPTION

This class extends L<Fcm::CLI::Invoker|Fcm::CLI::Invoker> an inherits all its
methods. An object of this class is used to invoke the pretty printer for FCM
configuration files.

=head1 METHODS

See L<Fcm::CLI::Invoker|Fcm::CLI::Invoker> for a list of inherited methods.

=over 4

=item invoke()

Invokes the pretty printer for a FCM configuration file.

If the I<output> option is set, output goes to the location specified by this
value.  Otherwise, it prints to STDOUT.

=back

=head1 DIAGNOSTICS

=over 4

=item L<Fcm::Exception|Fcm::Exception>

The invoke() method can croak() with this exception if the configuration file
cannot be read.

=item L<Fcm::CLI::Exception|Fcm::CLI::Exception>

The invoke() method can croak() with this exception if no configuration file is
specified.

=back

=head1 TO DO

Unit tests.

=head1 SEE ALSO

L<Fcm::CfgFile|Fcm::CfgFile>,
L<Fcm::CLI::Exception|Fcm::CLI::Exception>,
L<Fcm::CLI::Invoker|Fcm::CLI::Invoker>,
L<Fcm::CLI::Subcommand|Fcm::CLI::Subcommand>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
