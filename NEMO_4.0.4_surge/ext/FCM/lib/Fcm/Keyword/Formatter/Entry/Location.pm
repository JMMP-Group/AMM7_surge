# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::Keyword::Formatter::Entry::Location;
use base qw{Fcm::Keyword::Formatter::Entry};

use Fcm::Config;
use Fcm::Keyword::Formatter::Entries;

################################################################################
# Formats a keyword entry
sub format {
    my ($self, $entry) = @_;
    my $return = $self->SUPER::format($entry);
    for my $implied_entry (@{$entry->get_implied_entry_list()}) {
        $return .= $self->SUPER::format($implied_entry);
    }
    if (@{$entry->get_revision_entries()->get_all_entries()}) {
        my $formatter = Fcm::Keyword::Formatter::Entries->new();
        $return .= "\n[revision keyword]\n";
        $return .= $formatter->format($entry->get_revision_entries());
    }
    return $return;
}

1;
__END__

=head1 NAME

Fcm::Keyword::Formatter::Entry::Location

=head1 SYNOPSIS

    use Fcm::Keyword::Formatter::Entry::Location;
    $formatter = Fcm::Keyword::Formatter::Entry::Location->new();
    print($formatter->format($entry));

=head1 DESCRIPTION

An object of this class is used to format the detail in a location keyword entry.

=head1 METHODS

=over 4

=item new()

Constructor.

=item format($entry)

Returns a string representation of $entry.

=back

=head1 SEE ALSO

L<Fcm::Keyword::Entry|Fcm::Keyword::Entry>,
L<Fcm::Keyword::Formatter::Entry|Fcm::Keyword::Formatter::Entry>,
L<Fcm::Keyword::Formatter::Entries|Fcm::Keyword::Formatter::Entries>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
