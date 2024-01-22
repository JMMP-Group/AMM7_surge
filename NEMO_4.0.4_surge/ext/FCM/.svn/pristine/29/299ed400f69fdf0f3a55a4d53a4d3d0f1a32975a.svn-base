# ------------------------------------------------------------------------------
# NAME
#   Fcm::Base
#
# DESCRIPTION
#   This is base class for all FCM OO packages.
#
# COPYRIGHT
#   (C) Crown copyright Met Office. All rights reserved.
#   For further details please refer to the file COPYRIGHT.txt
#   which you should have received as part of this distribution.
# ------------------------------------------------------------------------------

package Fcm::Base;

# Standard pragma
use strict;
use warnings;

use Fcm::Config;

my @scalar_properties = (
  'config', # instance of Fcm::Config, configuration setting
);

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $obj = Fcm::Base->new;
#
# DESCRIPTION
#   This method constructs a new instance of the Fcm::Base class.
# ------------------------------------------------------------------------------

sub new {
  my $this  = shift;
  my %args  = @_;
  my $class = ref $this || $this;

  my $self  = {};
  for (@scalar_properties) {
    $self->{$_} = exists $args{uc ($_)} ? $args{uc ($_)} : undef;
  }

  bless $self, $class;
  return $self;
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $value = $obj->X;
#   $obj->X ($value);
#
# DESCRIPTION
#   Details of these properties are explained in @scalar_properties.
# ------------------------------------------------------------------------------

for my $name (@scalar_properties) {
  no strict 'refs';

  *$name = sub {
    my $self = shift;

    # Argument specified, set property to specified argument
    if (@_) {
      $self->{$name} = $_[0];
    }

    # Default value for property
    if (not defined $self->{$name}) {
      if ($name eq 'config') {
        # Configuration setting of the main program
        $self->{$name} = Fcm::Config->instance();
      }
    }

    return $self->{$name};
  }
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $value = $self->setting (@args); # $self->config->setting
#   $value = $self->verbose (@args); # $self->config->verbose
# ------------------------------------------------------------------------------

for my $name (qw/setting verbose/) {
  no strict 'refs';

  *$name = sub {
    my $self = shift;
    return $self->config->$name (@_);
  }
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $value = $self->cfglabel (@args);
#
# DESCRIPTION
#   This is an alias to $self->config->setting ('CFG_LABEL', @args);
# ------------------------------------------------------------------------------

sub cfglabel {
  my $self = shift;
  return $self->setting ('CFG_LABEL', @_);
}

# ------------------------------------------------------------------------------

1;

__END__
