# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::Keyword::Formatter::Entry;

################################################################################
# Constructor
sub new {
    my ($class) = @_;
    return bless(\do{my $annon_scalar}, $class);
}

################################################################################
# Formats a keyword entry
sub format {
    my ($self, $entry) = @_;
    return sprintf("%s = %s\n", $entry->get_key(), $entry->get_value());
}

1;
__END__

=head1 NAME

Fcm::Keyword::Formatter::Entry

=head1 SYNOPSIS

    use Fcm::Keyword::Formatter::Entry;
    $formatter = Fcm::Keyword::Formatter::Entry->new();
    print($formatter->format($entry));

=head1 DESCRIPTION

An object of this class is used to format a keyword entry.

=head1 METHODS

=over 4

=item new()

Constructor.

=item format($entry)

Returns a simple string representation of $entry.

=back

=head1 SEE ALSO

L<Fcm::Keyword::Entry|Fcm::Keyword::Entry>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
