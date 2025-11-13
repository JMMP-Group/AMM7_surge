# ------------------------------------------------------------------------------
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# ------------------------------------------------------------------------------
use strict;
use warnings;

package Fcm::Keyword::Entries;

use Carp qw{croak};
use Fcm::Util::ClassLoader;

sub new {
    my ($class, $args_ref) = @_;
    return bless(
        {
            entry_class             => 'Fcm::Keyword::Entry',
            entry_by                => {key => {}, value => {}},
            has_loaded_entries_from => {},
            loaders                 => [],
            ($args_ref && ref($args_ref) eq 'HASH' ? %{$args_ref} : ()),
        },
        $class,
    );
}

################################################################################
# Returns the class of entries stored by this entries list
sub get_entry_class {
    my ($self) = @_;
    return $self->{entry_class};
}

################################################################################
# Returns all entries
sub get_all_entries {
    my ($self) = @_;
    if (!%{$self->{entry_by}{key}}) {
        # Nothing set, attempt to load entries
        $self->load_entries();
    }
    if (wantarray()) {
        return values(%{$self->{entry_by}{key}});
    }
    else {
        return [values(%{$self->{entry_by}{key}})];
    }
}

################################################################################
# Methods: get_entry_by_*
for my $name (
    ### Returns an entry with a matching key
    'key',
    ### Returns an entry with a matching value
    'value'
) {
    no strict qw{refs};
    my $method = "get_entry_by_$name";
    *$method = sub {
        my ($self, $search_key) = @_;
        if (!defined($search_key)) {
            return;
        }
        my $sk = ($name eq 'key') ? uc($search_key) : $search_key;
        if (!exists($self->{entry_by}{$name}{$sk})) {
            $self->load_entries($name, $sk);
        }
        if (exists($self->{entry_by}{$name}{$sk})) {
            return $self->{entry_by}{$name}{$sk};
        }
        else {
            return;
        }
    }
}

################################################################################
# Adds an entry
sub add_entry {
    my ($self, $key, $value, $args_ref) = @_;
    Fcm::Util::ClassLoader::load($self->get_entry_class());
    my $entry = $self->get_entry_class()->new({
        key   => uc($key),
        value => $value,
        ($args_ref && ref($args_ref) eq 'HASH' ? %{$args_ref} : ()),
    });
    $self->{entry_by}{key}{uc($key)} = $entry;
    $self->{entry_by}{value}{$value} = $entry;
    return $entry;
}

################################################################################
# Returns the loaders for this entries list
sub get_loaders {
    my ($self) = @_;
    return (wantarray() ? @{$self->{loaders}} : $self->{loaders});
}

################################################################################
# Loads entries using its loaders
sub load_entries {
    my ($self, $name, $search_key) = @_;
    LOADER:
    for my $loader ($self->get_loaders()) {
        if ($self->{has_loaded_entries_from}{$loader->get_source()}) {
            next LOADER;
        }
        $self->{has_loaded_entries_from}{$loader->get_source()}
            = $loader->load_to($self);
        if ($name && exists($self->{entry_by}{$name}{$search_key})) {
            last LOADER;
        }
    }
}

1;
__END__

=head1 NAME

Fcm::Keyword::Entries

=head1 SYNOPSIS

    use Fcm::Keyword::Entries;

    my $entries = Fcm::Keyword::Entries->new({
        entry_class => $entry_class,
        loaders     => \@loaders,
    });
    $entry = $entries->get_entry_by_key($key);
    $entry = $entries->get_entry_by_value($value);

    for my $entry ($entries->get_entries()) {
        # ...
    }

    $entries->add_entry($key, $value);

=head1 DESCRIPTION

This module is used to manipulate FCM keyword entries. It is used by
L<Fcm::Keyword|Fcm::Keyword> to store keyword entries, which are
L<Fcm::Keyword::Entry|Fcm::Keyword::Entry> objects.

=head1 METHODS

=over 4

=item C<new({entry_class =E<gt> $entry_class, loaders =E<gt> \@loaders})>

Constructor. The argument should be a reference to hash, where:

I<entry_class> is a string representing the class name of entries in this
object. The class must be a sub-class of
L<Fcm::Keyword::Entry|Fcm::Keyword::Entry>. The default is
"L<Fcm::Keyword::Entry|Fcm::Keyword::Entry>".

I<loaders> is a reference to an array of
L<Fcm::Keyword::Loader|Fcm::Keyword::Loader> objects, which will be used to
load entries into this object. The default is an empty array.

=item add_entry($key,$value)

Adds an entry. Returns the added entry. (Keys are converted to uppercases
automatically.)

=item get_all_entries()

Returns all entries that are currently loaded.

=item get_entry_by_key($key)

Return an entry, whose key matches $key. (Search is case-insensitive.) Returns
undef if there is no matching entry.

=item get_entry_by_value($value)

Return an entry, whose value matches $value. (Search is case-sensitive.)
Returns undef if there is no matching entry.

=item get_loaders()

Returns the loaders for loading entries.

=item load_entries()

Loads entries from its loaders, as returned by get_loaders(). This method can
also be triggered by get_all_entries(), if the entry list is empty, or by
get_entry_by_key($key) and get_entry_by_value($value) methods, if there is no
matching entry in the current lookup lists.

=back

=head1 TO DO

Handle duplicated entries in add_entry($key,$value).

=head1 SEE ALSO

L<Fcm::Keyword|Fcm::Keyword>,
L<Fcm::Keyword::Entry|Fcm::Keyword::Entry>,
L<Fcm::Keyword::Loader|Fcm::Keyword::Loader>

=head1 COPYRIGHT

E<169> Crown copyright Met Office. All rights reserved.

=cut
