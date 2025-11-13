# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::Keyword::Entry;

sub new {
    my ($class, $args_ref) = @_;
    if (!$args_ref) {
        $args_ref = {};
    }
    return bless({%{$args_ref}}, $class);
}

################################################################################
### Methods: get_*
for my $name (
    # Returns the key of this entry
    'key',
    # Returns the value of this entry
    'value',
) {
    no strict qw{refs};
    my $getter = "get_$name";
    *$getter = sub {
        my ($self) = @_;
        return $self->{$name};
    }
}

1;
__END__

=head1 NAME

Fcm::Keyword::Entry

=head1 SYNOPSIS

    use Fcm::Keyword::Entry;

    $entry = Fcm::Keyword::Entry->new({key => $key, value => $value});
    $key = $entry->get_key();
    $value = $entry->get_value();

=head1 DESCRIPTION

An object of this class represents a FCM keyword entry.

=head1 METHODS

=over 4

=item C<new({key =E<gt> $key, value =E<gt> $value})>

Constructor.

=item get_key()

Returns the key of this keyword entry.

=item get_value()

Returns the value of this keyword entry.

=back

Simple formatter for displaying an entry.

=head1 SEE ALSO

L<Fcm::Keyword|Fcm::Keyword>,
L<Fcm::Keyword::Entries|Fcm::Keyword::Entries>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
