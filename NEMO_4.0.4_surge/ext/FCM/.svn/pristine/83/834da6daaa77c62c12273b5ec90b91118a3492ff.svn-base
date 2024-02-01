# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::CLI::Invoker::KeywordPrinter;
use base qw{Fcm::CLI::Invoker};

use Carp qw{croak};
use Fcm::CLI::Exception;
use Fcm::Keyword;
use Fcm::Keyword::Formatter::Entries;
use Fcm::Keyword::Formatter::Entry::Location;
use Fcm::Keyword::Exception;
use Fcm::Util qw{get_url_of_wc};

################################################################################
# Invokes the sub-system
sub invoke {
    my ($self) = @_;
    my @targets = $self->get_arguments();
    if (@targets) {
        for my $target (@targets) {
            my $entry_list = Fcm::Keyword::get_location_entries_for($target);
            my $loc = $target;
            if (-e $target) {
                $loc = get_url_of_wc($target);
                if (!$loc) {
                    croak(Fcm::Keyword::Exception->new({message => sprintf(
                        "%s: unrecognised version control resource", $target,
                    )}));
                }
            }
            my @entry_list = Fcm::Keyword::get_location_entries_for($loc);
            if (!@entry_list) {
                croak(Fcm::Keyword::Exception->new({message => sprintf(
                    "%s: no FCM location keyword found for this target", $target,
                )}));
            }
            my $formatter = Fcm::Keyword::Formatter::Entry::Location->new();
            for my $entry (
                sort {$a->get_key() cmp $b->get_key()}
                grep {!$_->is_implied()}
                @entry_list
            ) {
                print($formatter->format($entry), "\n");
            }
        }
    }
    else {
        my $formatter = Fcm::Keyword::Formatter::Entries->new();
        print($formatter->format(Fcm::Keyword::get_entries()));
    }
}

1;
__END__

=head1 NAME

Fcm::CLI::Invoker::KeywordPrinter

=head1 SYNOPSIS

    use Fcm::CLI::Invoker::KeywordPrinter;
    $invoker = Fcm::CLI::Invoker::KeywordPrinter->new({
        command   => $command,
        options   => \%options,
        arguments => $arguments,
    });
    $invoker->invoke();

=head1 DESCRIPTION

This class extends L<Fcm::CLI::Invoker|Fcm::CLI::Invoker> an inherits all its
methods. An object of this class is used to invoke the location keyword printer.

=head1 METHODS

See L<Fcm::CLI::Invoker|Fcm::CLI::Invoker> for a list of inherited methods.

=over 4

=item invoke()

Invokes the location keyword printer. If a namespace is specified in the
argument, prints revision keywords and browser mapping templates for the
specified namespace. If a namespace is not specified, prints all registered
location keywords.

=back

=head1 DIAGNOSTICS

=over 4

=item L<Fcm::Keyword::Exception|Fcm::Keyword::Exception>

The invoke() method can croak() with this exception if there is no matching
namespace matching that of the specified.

=back

=head1 TO DO

Unit tests.

=head1 SEE ALSO

L<Fcm::CLI::Exception|Fcm::CLI::Exception>,
L<Fcm::CLI::Invoker|Fcm::CLI::Invoker>,
L<Fcm::CLI::Subcommand|Fcm::CLI::Subcommand>,
L<Fcm::Keyword|Fcm::Keyword>,
L<Fcm::Keyword::Formatter::Entries|Fcm::Keyword::Formatter::Entries>,
L<Fcm::Keyword::Formatter::Entry::Location|Fcm::Keyword::Formatter::Entry::Location>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
