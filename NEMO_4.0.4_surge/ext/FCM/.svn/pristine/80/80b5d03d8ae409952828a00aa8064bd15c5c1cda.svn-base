# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::CLI::Config;

use Fcm::CLI::Config::Default;
use List::Util qw{first};
use Scalar::Util qw{blessed};

my $INSTANCE;

################################################################################
# Class method: returns an instance of this class
sub instance {
    my ($class, $args_ref) = @_;
    if ($args_ref || !$INSTANCE) {
        $INSTANCE = bless({
            core_subcommands => [@Fcm::CLI::Config::Default::CORE_SUBCOMMANDS],
            vc_subcommands   => [@Fcm::CLI::Config::Default::VC_SUBCOMMANDS],
            (defined($args_ref) ? %{$args_ref} : ()), 
        }, $class);
    }
    return $INSTANCE;
}

################################################################################
# Returns a subcommand matching $key
sub get_subcommand_of {
    my ($self, $key) = @_;
    if (blessed($key) && $key->isa('Fcm::CLI::Subcommand')) {
        return first {"$_" eq "$key"} ($self->get_subcommands());
    }
    else {
        return first {$_->has_a_name($key)} ($self->get_subcommands());
    }
}

################################################################################
# Returns the subcommands
sub get_subcommands {
    my ($self) = @_;
    my @return = ($self->get_core_subcommands(), $self->get_vc_subcommands());
    return (wantarray() ? @return : \@return);
}

################################################################################
# Returns the core subcommands
sub get_core_subcommands {
    my ($self) = @_;
    return (
        wantarray() ? @{$self->{core_subcommands}} : $self->{core_subcommands}
    );
}

################################################################################
# Returns the subcommands that are relevant only with a VC system
sub get_vc_subcommands {
    my ($self) = @_;
    return (wantarray() ? @{$self->{vc_subcommands}} : $self->{vc_subcommands});
}

1;
__END__

=head1 NAME

Fcm::CLI::Config

=head1 SYNOPSIS

    use Fcm::CLI::Config;
    $cli_config = Fcm::CLI::Config->instance();
    $subcommand = $cli_config->get_subcommand_of($key);
    @subcommands = $cli_config->get_subcommands();
    @core_subcommands = $cli_config->get_core_subcommands();
    @vc_subcommands = $cli_config->get_vc_subcommands();

=head1 DESCRIPTION

This class provides the configuration of the FCM command line interface.

=head1 METHODS

=over 4

=item instance($arg_ref)

Returns an instance of this class.

Creates the instance on first call, or replaces it with a new one if $args_ref
is defined in subsequent call. $args_ref should be a reference to a hash. The
hash can contain I<core_subcommands> and I<vc_subcommands>. Each of these
settings should point to an array reference containing L<Fcm::CLI::Subcommand>
objects. If the setting is unspecified, it uses the default from
L<Fcm::CLI::Config::Default|Fcm::CLI::Config::Default>.

=item get_subcommand_of($key)

Returns a L<Fcm::CLI::Subcommand|Fcm::CLI::Subcommand> object matching the
search $key. Returns undef if there is no match.

=item get_subcommands()

Short-hand for:
    ($self->get_core_subcommands(), $self->get_vc_subcommands())

=item get_core_subcommands()

Returns the core subcommands.

=item get_vc_subcommands()

Returns the subcommands that are relevant only in the presence of a VC system.

=back

=head1 SEE ALSO

L<Fcm::CLI::Config::Default|Fcm::CLI::Config::Default>,
L<Fcm::CLI::Invoker|Fcm::CLI::Invoker>,
L<Fcm::CLI::Subcommand|Fcm::CLI::Subcommand>,
L<Fcm::CLI::Option|Fcm::CLI::Option>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
