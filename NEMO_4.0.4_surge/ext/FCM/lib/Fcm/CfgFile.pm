# ------------------------------------------------------------------------------
# NAME
#   Fcm::CfgFile
#
# DESCRIPTION
#   This class is used for reading and writing FCM config files. A FCM config
#   file is a line-based text file that provides information on how to perform
#   a particular task using the FCM system.
#
# COPYRIGHT
#   (C) Crown copyright Met Office. All rights reserved.
#   For further details please refer to the file COPYRIGHT.txt
#   which you should have received as part of this distribution.
# ------------------------------------------------------------------------------

package Fcm::CfgFile;
@ISA = qw(Fcm::Base);

# Standard pragma
use warnings;
use strict;

# Standard modules
use Carp;
use File::Basename;
use File::Path;
use File::Spec;

# FCM component modules
use Fcm::Base;
use Fcm::CfgLine;
use Fcm::Config;
use Fcm::Keyword;
use Fcm::Util;

# List of property methods for this class
my @scalar_properties = (
  'actual_src', # actual source of configuration file
  'lines',      # list of lines, Fcm::CfgLine objects
  'pegrev',     # peg revision of configuration file
  'src',        # source of configuration file
  'type',       # type of configuration file
  'version',    # version of configuration file
);

# Local module variables
my $expand_type   = 'bld|ext'; # config file type that needs variable expansions

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $obj = Fcm::CfgFile->new (%args);
#
# DESCRIPTION
#   This method constructs a new instance of the Fcm::CfgFile class. See above
#   for allowed list of properties. (KEYS should be in uppercase.)
# ------------------------------------------------------------------------------

sub new {
  my $this  = shift;
  my %args  = @_;
  my $class = ref $this || $this;

  my $self = Fcm::Base->new (%args);

  bless $self, $class;

  for (@scalar_properties) {
    $self->{$_} = exists $args{uc ($_)} ? $args{uc ($_)} : undef;
  }

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

    if (@_) {
      $self->{$name} = $_[0];
    }

    if (not defined $self->{$name}) {
      if ($name eq 'lines') {
        $self->{$name} = [];
      }
    }

    return $self->{$name};
  }
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $mtime = $obj->mtime ();
#
# DESCRIPTION
#   This method returns the modified time of the configuration file source.
# ------------------------------------------------------------------------------

sub mtime {
  my $self  = shift;
  my $mtime = undef;

  if (-f $self->src) {
    $mtime = (stat $self->src)[9];
  }

  return $mtime;
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $read = $obj->read_cfg ();
#
# DESCRIPTION
#   This method reads the current configuration file. It returns the number of
#   lines read from the config file, or "undef" if it fails. The result is
#   placed in the LINES array of the current instance, and can be accessed via
#   the "lines" method.
# ------------------------------------------------------------------------------

sub read_cfg {
  my $self = shift;

  my @lines = $self->_get_cfg_lines;

  # List of CFG types that need INC declarations expansion
  my %exp_inc    = ();
  for (split (/$Fcm::Config::DELIMITER_LIST/, $self->setting ('CFG_EXP_INC'))) {
    $exp_inc{uc ($_)} = 1;
  }

  # List of CFG labels that are reserved keywords
  my %cfg_keywords = ();
  for (split (/$Fcm::Config::DELIMITER_LIST/, $self->setting ('CFG_KEYWORD'))) {
    $cfg_keywords{$self->cfglabel ($_)} = 1;
  }

  # Loop each line, to separate lines into label : value pairs
  my $cont = undef;
  my $here = undef;
  for my $line_num (1 .. @lines) {
    my $line = $lines[$line_num - 1];
    chomp $line;

    my $label   = '';
    my $value   = '';
    my $comment = '';

    # If this line is a continuation, set $start to point to the line that
    # starts this continuation. Otherwise, set $start to undef
    my $start = defined ($cont) ? $self->lines->[$cont] : undef;
    my $warning = undef;

    if ($line =~ /^(\s*#.*)$/) { # comment line
      $comment = $1;

    } elsif ($line =~ /\S/) {    # non-blank line
      if (defined $cont) {
        # Previous line has a continuation mark
        $value = $line;

        # Separate value and comment
        if ($value =~ s/((?:\s+|^)#\s+.*)$//) {
          $comment = $1;
        }

        # Remove leading spaces
        $value =~ s/^\s*\\?//;

        # Expand environment variables
        my $warn;
        ($value, $warn) = $self->_expand_variable ($value, 1) if $value;
        $warning .= ($warning ? ', ' : '') . $warn if $warn;

        # Expand internal variables
        ($value, $warn) = $self->_expand_variable ($value, 0) if $value;
        $warning .= ($warning ? ', ' : '') . $warn if $warn;

        # Get "line" that begins the current continuation
        my $v = $start->value . $value;
        $v =~ s/\\$//;
        $start->value ($v);

      } else {
        # Previous line does not have a continuation mark
        if ($line =~ /^\s*(\S+)(?:\s+(.*))?$/) {
          # Check line contains a valid label:value pair
          $label = $1;
          $value = defined ($2) ? $2 : '';

          # Separate value and comment
          if ($value =~ s/((?:\s+|^)#\s+.*)$//) {
            $comment = $1;
          }

          # Remove trailing spaces
          $value =~ s/\s+$//;

          # Value begins with $HERE?
          $here  = ($value =~ /\$\{?HERE\}?(?:[^A-Z_]|$)/);

          # Expand environment variables
          my $warn;
          ($value, $warn) = $self->_expand_variable ($value, 1) if $value;
          $warning .= ($warning ? ', ' : '') . $warn if $warn;

          # Expand internal variables
          ($value, $warn) = $self->_expand_variable ($value, 0) if $value;
          $warning .= ($warning ? ', ' : '') . $warn if $warn;
        }
      }

      # Determine whether current line ends with a continuation mark
      if ($value =~ s/\\$//) {
        $cont = scalar (@{ $self->lines }) unless defined $cont;

      } else {
        $cont = undef;
      }
    }

    if (    defined($self->type())
        &&  exists($exp_inc{uc($self->type())})
        &&  uc($start ? $start->label() : $label) eq $self->cfglabel('INC')
        &&  !defined($cont)
    ) {
      # Current configuration file requires expansion of INC declarations
      # The start/current line is an INC declaration
      # The current line is not a continuation or is the end of the continuation

      # Get lines from an "include" configuration file
      my $src = ($start ? $start->value : $value);
      $src   .= '@' . $self->pegrev if $here and $self->pegrev;

      if ($src) {
        # Invoke a new instance to read the source
        my $cfg = Fcm::CfgFile->new (
          SRC => expand_tilde ($src), TYPE => $self->type,
        );

        $cfg->read_cfg;

        # Add lines to the lines array in the current configuration file
        $comment = 'INC ' . $src . ' ';
        push @{$self->lines}, Fcm::CfgLine->new (
          comment => $comment . '# Start',
          number  => ($start ? $start->number : $line_num),
          src     => $self->actual_src,
          warning => $warning,
        );
        push @{ $self->lines }, @{ $cfg->lines };
        push @{$self->lines}, Fcm::CfgLine->new (
          comment => $comment . '# End',
          src     => $self->actual_src,
        );

      } else {
        push @{$self->lines}, Fcm::CfgLine->new (
          number  => $line_num,
          src     => $self->actual_src,
          warning => 'empty INC declaration.'
        );
      }

    } else {
      # Push label:value pair into lines array
      push @{$self->lines}, Fcm::CfgLine->new (
        label   => $label,
        value   => ($label ? $value : ''),
        comment => $comment,
        number  => $line_num,
        src     => $self->actual_src,
        warning => $warning,
      );
    }

    next if defined $cont; # current line not a continuation

    my $slabel = ($start ? $start->label : $label);
    my $svalue = ($start ? $start->value : $value);
    next unless $slabel;

    # Check config file type and version
    if (index (uc ($slabel), $self->cfglabel ('CFGFILE')) == 0) {
      my @words = split /$Fcm::Config::DELIMITER_PATTERN/, $slabel;
      shift @words;

      my $name = @words ? lc ($words[0]) : 'type';

      if ($self->can ($name)) {
        $self->$name ($svalue);
      }
    }

    # Set internal variable
    $slabel =~ s/^\%//; # Remove leading "%" from label

    $self->config->variable ($slabel, $svalue)
      unless exists $cfg_keywords{$slabel};
  }

  # Report and reset warnings
  # ----------------------------------------------------------------------------
  for my $line (@{ $self->lines }) {
    w_report $line->format_warning if $line->warning;
    $line->warning (undef);
  }

  return @{ $self->lines };

}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $rc = $obj->print_cfg ($file, [$force]);
#
# DESCRIPTION
#   This method prints the content of current configuration file. If no
#   argument is specified, it prints output to the standard output. If $file is
#   specified, and is a writable file name, the output is sent to the file.  If
#   the file already exists, its content is compared to the current output.
#   Nothing will be written if the content is unchanged unless $force is
#   specified. Otherwise, for typed configuration files, the existing file is
#   renamed using a prefix that contains its last modified time. The method
#   returns 1 if there is no error.
# ------------------------------------------------------------------------------

sub print_cfg {
  my ($self, $file, $force) = @_;

  # Count maximum number of characters in the labels, (for pretty printing)
  my $max_label_len = 0;
  for my $line (@{ $self->lines }) {
    next unless $line->label;
    my $label_len  = length $line->label;
    $max_label_len = $label_len if $label_len > $max_label_len;
  }

  # Output string
  my $out = '';

  # Append each line of the config file to the output string
  for my $line (@{ $self->lines }) {
    $out .= $line->print_line ($max_label_len - length ($line->label) + 1);
    $out .= "\n";
  }

  if ($out) {
    my $old_select = select;

    # Open file if necessary
    if ($file) {
      # Make sure the host directory exists and is writable
      my $dirname = dirname $file;
      if (not -d $dirname) {
        print 'Make directory: ', $dirname, "\n" if $self->verbose;
        mkpath $dirname;
      }
      croak $dirname, ': cannot write to config file directory, abort'
        unless -d $dirname and -w $dirname;

      if (-f $file and not $force) {
        if (-r $file) {
          # Read old config file to see if content has changed
          open IN, '<', $file or croak $file, ': cannot open (', $!, '), abort';
          my $in_lines = '';
          while (my $line = <IN>) {
            $in_lines .= $line;
          }
          close IN or croak $file, ': cannot close (', $!, '), abort';

          # Return if content is up-to-date
          if ($in_lines eq $out) {
            print 'No change in ', lc ($self->type), ' cfg: ', $file, "\n"
              if $self->verbose > 1 and $self->type;
            return 1;
          }
        }

        # If config file already exists, make sure it is writable
        if (-w $file) {
          if ($self->type) {
            # Existing config file writable, rename it using its time stamp
            my $mtime = (stat $file)[9];
            my ($sec, $min, $hour, $mday, $mon, $year) = (gmtime $mtime)[0 .. 5];
            my $timestamp = sprintf '%4d%2.2d%2.2d_%2.2d%2.2d%2.2d_',
                            $year + 1900, $mon + 1, $mday, $hour, $min, $sec;
            my $oldfile   = File::Spec->catfile (
              $dirname, $timestamp . basename ($file)
            );
            rename $file, $oldfile;
            print 'Rename existing ', lc ($self->type), ' cfg: ',
                  $oldfile, "\n" if $self->verbose > 1;
          }

        } else {
          # Existing config file not writable, throw an error
          croak $file, ': config file not writable, abort';
        }
      }

      # Open file and select file handle
      open OUT, '>', $file
        or croak $file, ': cannot open config file (', $!, '), abort';
      select OUT;
    }

    # Print output
    print $out;

    # Close file if necessary
    if ($file) {
      select $old_select;
      close OUT or croak $file, ': cannot close config file (', $!, '), abort';

      if ($self->type and $self->verbose > 1) {
        print 'Generated ', lc ($self->type), ' cfg: ', $file, "\n";

      } elsif ($self->verbose > 2) {
        print 'Generated cfg: ', $file, "\n";
      }
    }

  } else {
    # Warn if nothing to print
    my $warning = 'Empty configuration';
    $warning   .= ' - nothing written to file: ' . $file if $file;
    carp $warning if $self->type;
  }

  return 1;
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   @lines = $self->_get_cfg_lines ();
#
# DESCRIPTION
#   This internal method reads from a configuration file residing in a
#   Subversion repository or in the normal file system.
# ------------------------------------------------------------------------------

sub _get_cfg_lines {
  my $self  = shift;
  my @lines = ();

  my $verbose = $self->verbose;

  my ($src) = $self->src();
  if ($src =~ qr{\A([A-Za-z][\w\+-\.]*):}xms) { # is a URI
    $src = Fcm::Keyword::expand($src);
    # Config file resides in a SVN repository
    # --------------------------------------------------------------------------
    # Set URL source and version
    my $rev = 'HEAD';

    # Extract version from source if it exists
    if ($src =~ s{\@ ([^\@]+) \z}{}xms) {
      $rev = $1;
    }

    $src = Fcm::Util::tidy_url($src);

    # Check whether URL is a config file
    my $rc;
    my @cmd = (qw/svn cat/, $src . '@' . $rev);
    @lines = &run_command (
      \@cmd, METHOD => 'qx', DEVNULL => 1, RC => \$rc, ERROR => 'ignore',
    );

    # Error in "svn cat" command
    if ($rc) {
      # See whether specified config file is a known type
      my %cfgname = %{ $self->setting ('CFG_NAME') };
      my $key     = uc $self->type;
      my $file    = exists $cfgname{$key} ? $cfgname{$key} : '';

      # If config file is a known type, specified URL may be a directory
      if ($file) {
        # Check whether a config file with a default name exists in the URL
        my $path = $src . '/' . $file;
        my @cmd  = (qw/svn cat/, $path . '@' . $rev);

        @lines = &run_command (
          \@cmd, METHOD => 'qx', DEVNULL => 1, RC => \$rc, ERROR => 'ignore',
        );

        # Check whether a config file with a default name exists under the "cfg"
        # sub-directory of the URL
        if ($rc) {
          my $cfgdir = $self->setting (qw/DIR CFG/);
          $path   = $src . '/' . $cfgdir . '/' . $file;
          my @cmd = (qw/svn cat/, $path . '@' . $rev);

          @lines  = &run_command (
            \@cmd, METHOD => 'qx', DEVNULL => 1, RC => \$rc, ERROR => 'ignore',
          );
        }

        $src = $path unless $rc;
      }
    }

    if ($rc) {
      # Error in "svn cat"
      croak 'Unable to locate config file from "', $self->src, '", abort';

    } else {
      # Print diagnostic, if necessary
      if ($verbose and $self->type and $self->type =~ /$expand_type/) {
        print 'Config file (', $self->type, '): ', $src;
        print '@', $rev if $rev;
        print "\n";
      }
    }

    # Record the actual source location
    $self->pegrev ($rev);
    $self->actual_src ($src);

  } else {
    # Config file resides in the normal file system
    # --------------------------------------------------------------------------
    my $src = $self->src;

    if (-d $src) { # Source is a directory
      croak 'Config file "', $src, '" is a directory, abort' if not $self->type;

      # Get name of the config file by looking at the type
      my %cfgname = %{ $self->setting ('CFG_NAME') };
      my $key     = uc $self->type;
      my $file    = exists $cfgname{$key} ? $cfgname{$key} : '';

      if ($file) {
        my $cfgdir = $self->setting (qw/DIR CFG/);

        # Check whether a config file with a default name exists in the
        # specified path, then check whether a config file with a default name
        # exists under the "cfg" sub-directory of the specified path
        if (-f File::Spec->catfile ($self->src, $file)) {
          $src = File::Spec->catfile ($self->src, $file);

        } elsif (-f File::Spec->catfile ($self->src, $cfgdir, $file)) {
          $src = File::Spec->catfile ($self->src, $cfgdir, $file);

        } else {
          croak 'Unable to locate config file from "', $self->src, '", abort';
        }

      } else {
        croak 'Unknown config file type "', $self->type, '", abort';
      }
    }

    if (-r $src) {
      open FILE, '<', $src;
      print 'Config file (', $self->type, '): ', $src, "\n"
        if $verbose and $self->type and $self->type =~ /$expand_type/;

      @lines = readline 'FILE';
      close FILE;

    } else {
      croak 'Unable to read config file "', $src, '", abort';
    }

    # Record the actual source location
    $self->actual_src ($src);
  }

  return @lines;
}

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $string = $self->_expand_variable ($string, $env[, \%recursive]);
#
# DESCRIPTION
#   This internal method expands variables in $string. If $env is true, it
#   expands environment variables. Otherwise, it expands local variables. If
#   %recursive is set, it indicates that this method is being called
#   recursively. In which case, it must not attempt to expand a variable that
#   exists in the keys of %recursive.
# ------------------------------------------------------------------------------

sub _expand_variable {
  my ($self, $string, $env, $recursive) = @_;

  # Pattern for environment/local variable
  my @patterns = $env
    ? (qr#\$([A-Z][A-Z0-9_]+)#, qr#\$\{([A-Z][A-Z0-9_]+)\}#)
    : (qr#%(\w+(?:::[\w\.-]+)*)#, qr#%\{(\w+(?:(?:::|/)[\w\.-]+)*)\}#);

  my $ret = '';
  my $warning = undef;
  while ($string) {
    # Find the first match in $string
    my ($prematch, $match, $postmatch, $var_label);
    for my $pattern (@patterns) {
      next unless $string =~ /$pattern/;
      if ((not defined $prematch) or length ($`) < length ($prematch)) {
        $prematch = $`;
        $match = $&;
        $var_label = $1;
        $postmatch = $';
      }
    }

    if ($match) {
      $ret .= $prematch;
      $string = $postmatch;

      # Get variable value from environment or local configuration
      my $variable = $env
                     ? (exists $ENV{$var_label} ? $ENV{$var_label} : undef)
                     : $self->config->variable ($var_label);

      if ($env and $var_label eq 'HERE' and not defined $variable) {
        $variable = dirname ($self->actual_src);
        $variable = File::Spec->rel2abs ($variable) if not &is_url ($variable);
      }

      # Substitute match with value of variable
      if (defined $variable) {
        my $cyclic = 0;
        if ($recursive) {
          if (exists $recursive->{$var_label}) {
            $cyclic = 1;

          } else {
            $recursive->{$var_label} = 1;
          }

        } else {
          $recursive = {$var_label => 1};
        }

        if ($cyclic) {
          $warning .= ', ' if $warning;
          $warning .= $match . ': cyclic dependency, variable not expanded';
          $ret .= $variable;

        } else {
          my ($r, $w) = $self->_expand_variable ($variable, $env, $recursive);
          $ret .= $r;
          if ($w) {
            $warning .= ', ' if $warning;
            $warning .= $w;
          }
        }

      } else {
        $warning .= ', ' if $warning;
        $warning .= $match . ': variable not expanded';
        $ret .= $match;
      }

    } else {
      $ret .= $string;
      $string = "";
    }
  }

  return ($ret, $warning);
}

1;

__END__
