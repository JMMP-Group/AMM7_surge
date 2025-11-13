# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::Keyword::Formatter::Entries;

use Fcm::Keyword::Formatter::Entry;

################################################################################
# Constructor
sub new {
    my ($class) = @_;
    return bless(\do{my $annon_scalar}, $class);
}

################################################################################
# Formats a keyword entry
sub format {
    my ($self, $entries) = @_;
    my $formatter = Fcm::Keyword::Formatter::Entry->new();
    my $return = q{};
    for my $entry (
        sort {$a->get_key() cmp $b->get_key()}
        grep {!$_->can('is_implied') || !$_->is_implied()}
        $entries->get_all_entries()
    ) {
        $return .= $formatter->format($entry);
    }
    return $return;
}

1;
__END__

=head1 NAME

Fcm::Keyword::Formatter::Entries

=head1 SYNOPSIS

    use Fcm::Keyword::Formatter::Entries;
    $formatter = Fcm::Keyword::Formatter::Entries->new();
    print($formatter->format($entries));

=head1 DESCRIPTION

An object of this class is used to format a keyword entries object.

=head1 METHODS

=over 4

=item new()

Constructor.

=item format($entries)

Returns a simple string representation of $entries.

=back

=head1 SEE ALSO

L<Fcm::Keyword::Entries|Fcm::Keyword::Entries>,
L<Fcm::Keyword::Entry|Fcm::Keyword::Entry>,
L<Fcm::Keyword::Formatter::Entry|Fcm::Keyword::Formatter::Entry>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
