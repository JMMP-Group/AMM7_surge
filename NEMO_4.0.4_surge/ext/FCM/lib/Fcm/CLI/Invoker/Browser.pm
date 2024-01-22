# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::CLI::Invoker::Browser;
use base qw{Fcm::CLI::Invoker};

use Carp qw{croak};
use Fcm::CLI::Exception;
use Fcm::Config;
use Fcm::Keyword;
use Fcm::Util qw{expand_tilde get_url_of_wc is_wc run_command};

################################################################################
# Invokes the sub-system
sub invoke {
    my ($self) = @_;
    my $config = Fcm::Config->instance();
    my $browser
        = $self->get_options()->{browser} ? $self->get_options()->{browser}
        :                                   $config->setting(qw/WEB_BROWSER/)
        ;
    my ($target) = $self->get_arguments();
    if (!$target) {
        if (is_wc()) {
            $target = q{.};
        }
        else {
            croak(Fcm::CLI::Exception->new({
                message => 'no TARGET specified and . not a working copy',
            }));
        }
    }
    $target = expand_tilde($target);
    if (-e $target) {
        $target = get_url_of_wc($target);
    }

    my $browser_url = Fcm::Keyword::get_browser_url($target);
    my @command = (split(qr{\s+}xms, $browser), $browser_url);
    run_command(\@command, METHOD => 'exec', PRINT => 1);
}

1;
__END__

=head1 NAME

Fcm::CLI::Invoker::Browser

=head1 SYNOPSIS

    use Fcm::CLI::Invoker::Browser;
    $invoker = Fcm::CLI::Invoker::Browser->new({
        command   => $command,
        options   => \%options,
        arguments => $arguments,
    });
    $invoker->invoke();

=head1 DESCRIPTION

This class extends L<Fcm::CLI::Invoker|Fcm::CLI::Invoker> an inherits all its
methods. An object of this class is used to invoke a web browser of a VC
location.

=head1 METHODS

See L<Fcm::CLI::Invoker|Fcm::CLI::Invoker> for a list of inherited methods.

=over 4

=item invoke()

Invokes a web browser for a VC target, if it can be mapped to a browser URL. If
a target is not specified in arguments, it uses the current working directory
as the target.

If the browser option is set, it is used as the browser command. Otherwise, the
default browser is used.

=back

=head1 DIAGNOSTICS

=over 4

=item L<Fcm::CLI::Exception|Fcm::CLI::Exception>

The invoke() method can croak() with this exception if no target is specified
and a target cannot be deduced from the current working directory.

=item L<Fcm::Keyword::Exception|Fcm::Keyword::Exception>

The invoke() method can croak() with this exception if the target cannot be
mapped to a browser URL.

=back

=head1 TO DO

Unit tests.

=head1 SEE ALSO

L<Fcm::CLI::Exception|Fcm::CLI::Exception>,
L<Fcm::CLI::Invoker|Fcm::CLI::Invoker>,
L<Fcm::CLI::Subcommand|Fcm::CLI::Subcommand>,
L<Fcm::Keyword|Fcm::Keyword>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
