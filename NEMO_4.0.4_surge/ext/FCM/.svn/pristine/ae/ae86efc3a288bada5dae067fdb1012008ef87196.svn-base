# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::CLI::Invoker::ConfigSystem;
use base qw{Fcm::CLI::Invoker};

use Cwd qw{cwd};
use Fcm::CLI::Exception;
use Fcm::Config;
use Fcm::Util::ClassLoader;

################################################################################
# Returns a hash map to convert CLI options to system invoke options.
sub get_cli2invoke_key_map {
    my ($self) = @_;
    return (
          wantarray() ? %{$self->{cli2invoke_key_map}}
        :               $self->{cli2invoke_key_map}
    );
}

################################################################################
# Returns the Fcm::ConfigSystem class for invoking the sub-system.
sub get_impl_class {
    my ($self) = @_;
    return $self->{impl_class};
}

################################################################################
# Invokes the sub-system
sub invoke {
    my ($self) = @_;
    my $options_ref = $self->get_options();
    if (exists($options_ref->{verbose})) {
        Fcm::Config->instance()->verbose($options_ref->{verbose});
    }

    Fcm::Util::ClassLoader::load($self->get_impl_class());
    my $system = $self->get_impl_class()->new();
    my ($cfg_file) = $self->get_arguments();
    $system->cfg()->src($cfg_file ? $cfg_file : cwd());

    my %map = $self->get_cli2invoke_key_map();
    my %invoke_args = map {($map{$_}, $options_ref->{$_})} keys(%map);
    $system->invoke(%invoke_args);
}

1;
__END__

=head1 NAME

Fcm::CLI::Invoke::ConfigSystem

=head1 SYNOPSIS

    use Fcm::CLI::Invoker::ConfigSystem;
    $invoker = Fcm::CLI::Invoker::ConfigSystem->new({
        command            => $command,
        options            => \%options,
        arguments          => $arguments,
        impl_class         => $class_name,
        cli2invoke_key_map => {
            option => 'OPTION',
            # ... more keys
        },
    });
    $invoker->invoke();

=head1 DESCRIPTION

This class extends L<Fcm::CLI::Invoker|Fcm::CLI::Invoker> and inherits all its
methods. An object of this class is used to invoke a
L<Fcm::ConfigSystem|Fcm::ConfigSystem>, e.g. extract and build.

=head1 METHODS

See L<Fcm::CLI::Invoker|Fcm::CLI::Invoker> for a list of inherited methods.

=over 4

=item get_cli2invoke_key_map()

Returns a hash containing a mapping table from the names of the relevant command
line options to the names to be given to the invoke() method of the implementing
L<Fcm::ConfigSystem|Fcm::ConfigSystem> object.

=item get_impl_class()

Returns the actual class that implements L<Fcm::ConfigSystem|Fcm::ConfigSystem>.
An object of this implementation will be created and used by invoke().

=item invoke()

Invokes the L<Fcm::ConfigSystem|Fcm::ConfigSystem> sub-system. If a
configuration file is not specified in the argument, it uses the current working
directory.

=back

=head1 SEE ALSO

L<Fcm::ConfigSystem|Fcm::ConfigSystem>,
L<Fcm::CLI::Exception|Fcm::CLI::Exception>,
L<Fcm::CLI::Invoker|Fcm::CLI::Invoker>,
L<Fcm::CLI::Subcommand|Fcm::CLI::Subcommand>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
