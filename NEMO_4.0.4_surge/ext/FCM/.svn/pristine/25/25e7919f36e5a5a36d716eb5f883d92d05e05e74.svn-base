# ------------------------------------------------------------------------------
# NAME
#   Fcm::CmCommitMessage
#
# DESCRIPTION
#   This class contains methods to read, write and edit the commit message file
#   in a working copy.
#
# COPYRIGHT
#   (C) Crown copyright Met Office. All rights reserved.
#   For further details please refer to the file COPYRIGHT.txt
#   which you should have received as part of this distribution.
# ------------------------------------------------------------------------------

package Fcm::CmCommitMessage;
@ISA = qw(Fcm::Base);

# Standard pragma
use warnings;
use strict;

# Standard modules
use Carp;
use Cwd;
use File::Spec;
use File::Temp qw/tempfile/;

# FCM component modules
use Fcm::Base;
use Fcm::Util qw/e_report run_command/;

# List of property methods for this class
my @scalar_properties = (
  'auto_mesg',   # the automatically inserted part of a commit message
  'base',        # the base name of the commit message file
  'dir',         # the directory container of the commit message file
  'ignore_mesg', # the ignored part of a commit message
  'user_mesg',   # the user defined part of a commit message
);

# Commit log delimiter messages
my $log_delimiter        = '--Add your commit message ABOVE - ' .
                           'do not alter this line or those below--';
my $auto_delimiter       = '--FCM message (will be inserted automatically)--';
my $auto_delimiter_old   = '--This line will be ignored and those below ' .
                           'will be inserted automatically--';
my $status_delimiter     = '--Change summary ' .
                           '(not part of commit message)--';
my $status_delimiter_old = '--This line, and those below, will be ignored--';

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $obj = Fcm::CmCommitMessage->new ();
#
# DESCRIPTION
#   This method constructs a new instance of the Fcm::CmCommitMessage class.
# ------------------------------------------------------------------------------

sub new {
  my $this  = shift;
  my %args  = @_;
  my $class = ref $this || $this;

  my $self = Fcm::Base->new (%args);

  $self->{$_} = undef for (@scalar_properties);

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
      if ($name eq 'base') {
        # Reference to an array
        $self->{$name} = '#commit_message#';

      } elsif ($name eq 'dir') {
        # Current working directory
        $self->{$name} = &cwd ();

      } elsif ($name =~ /_mesg$/) {
        # Reference to an array
        $self->{$name} = [];
      }
    }

    return $self->{$name};
  }
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $file = $obj->file;
#   $obj->file ($file);
#
# DESCRIPTION
#   This method returns the full name of the commit message file. If an
#   argument is specified, the file is reset using the value of the argument.
# ------------------------------------------------------------------------------

sub file {
  my ($self, $file) = @_;

  if ($file) {
    $self->dir  (dirname  ($file));
    $self->base (basename ($file));
  }

  return File::Spec->catfile ($self->dir, $self->base);
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   ($user, $auto) = $obj->read_file ();
#
# DESCRIPTION
#   This function reads from the commit log message file. It resets the user
#   and the automatic messages after reading the file. It returns the message
#   back in two array references.
# ------------------------------------------------------------------------------

sub read_file {
  my $self = shift;

  my @user = ();
  my @auto = ();
  my $file = $self->file;

  if (-r $file) {
    open FILE, '<', $file or croak 'Cannot open ', $file, '(', $!, '), abort';

    my $in_auto = 0;
    while (<FILE>) {

      next if (index ($_, $log_delimiter) == 0); 

      if (index ($_, $status_delimiter)     == 0 ||
          index ($_, $status_delimiter_old) == 0) {
        # Ignore after the ignore delimiter
        last;
      }

      if (index ($_, $auto_delimiter)     == 0 ||
          index ($_, $auto_delimiter_old) == 0) {
        # Beginning of the automatically inserted message
        $in_auto = 1;
        next;
      }

      if ($in_auto) {
        push @auto, $_;

      } else {
        push @user, $_;
      }
    }

    close FILE;

    $self->user_mesg (\@user);
    $self->auto_mesg (\@auto);
  }

  return (\@user, \@auto);
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $obj->write_file ();
#
# DESCRIPTION
#   This function writes to the commit log message file based on the content of
#   the user defined message, and the automatically inserted message.
# ------------------------------------------------------------------------------

sub write_file {
  my $self = shift;
  my %args = @_;

  my @user = @{ $self->user_mesg };
  my @auto = @{ $self->auto_mesg };
  my $file = $self->file;

  open FILE, '>', $file or die 'Cannot open ', $file, '(', $!, '), abort';
  print FILE @user;
  print FILE $log_delimiter, "\n", $auto_delimiter, "\n", @auto if @auto;
  close FILE or croak 'Cannot close ', $file, '(', $!, '), abort';

  return;
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $file = $obj->edit_file ([TEMP => 1,] [BATCH => 1,]);
#
# DESCRIPTION
#   This function normally triggers an editor for editing the commit message.
#   If TEMP is set, it edits a temporary file. Otherwise, it edits the current
#   commit message file. It resets the user defined message on success. Returns
#   the name of the commit log file. Do not start the editor if BATCH is set.
# ------------------------------------------------------------------------------

sub edit_file {
  my $self  = shift;
  my %args  = @_;
  my $temp  = exists $args{TEMP}  ? $args{TEMP}  : 0;
  my $batch = exists $args{BATCH} ? $args{BATCH} : 0;

  my @user   = @{ $self->user_mesg };
  my @auto   = @{ $self->auto_mesg };
  my @ignore = @{ $self->ignore_mesg };
  my $file   = $self->file;

  if ($temp) {
    my $fh;
    ($fh, $file) = tempfile (SUFFIX => ".fcm", UNLINK => 1);
    close $fh;
  }

  # Add original or code driven message and status information to the file
  my $select = select;
  open FILE, '>', $file or croak 'Cannot open ', $file, ' (', $!, '), abort';
  select FILE;

  print @user;
  print (@auto || @user ? '' : "\n");
  print $log_delimiter, "\n";
  print $auto_delimiter, "\n", @auto, "\n" if @auto;
  print $status_delimiter, "\n\n";
  print @ignore if @ignore;

  close FILE or die 'Cannot close ', $file, ' (', $!, '), abort';
  select $select;

  if (not $batch) {
    # Select editor
    my $editor = 'nedit';
    
    if ($ENV{'SVN_EDITOR'}) {
      $editor = $ENV{'SVN_EDITOR'};

    } elsif ($ENV{'VISUAL'}) {
      $editor = $ENV{'VISUAL'};

    } elsif ($ENV{'EDITOR'}) {
      $editor = $ENV{'EDITOR'};
    }

    # Execute command to start the editor
    print 'Starting ', $editor, ' to edit commit message ...', "\n";
    &run_command ([split (/\s+/, $editor), $file]);
  }

  # Read the edited file, and extract user log message from it
  open FILE, '<', $file or croak 'Cannot open ', $file, ' (', $!, '), abort';
  my (@log);
  my $delimiter_found = 0;

  while (<FILE>) {
    if (index ($_, $log_delimiter) == 0) {
      $delimiter_found = 1;
      last;
    }
    push @log, $_;
  }

  close FILE;

  # Ensure log delimiter line was not altered
  e_report 'Error: the line "', $log_delimiter, '" has been altered, abort.'
    if not $delimiter_found;

  # Check for empty commit log
  e_report 'Error: log message unchanged or not specified, abort.'
    if join (' ', (@log, @auto)) =~ /^\s*$/;

  # Echo the commit message to standard output
  my $separator = '-' x 80 . "\n";
  print 'Change summary:', "\n";
  print $separator, @ignore, $separator;
  print 'Commit message is as follows:', "\n";
  print $separator, @log, @auto, $separator;

  open FILE, '>', $file or croak 'Cannot open ', $file, ' (', $!, '), abort';
  print FILE @log, @auto;
  close FILE or croak 'Cannot close ', $file, ' (', $!, '), abort';

  # Reset the array for the user specified log message
  $self->user_mesg (\@log);

  return $file;
}

# ------------------------------------------------------------------------------

1;

__END__
