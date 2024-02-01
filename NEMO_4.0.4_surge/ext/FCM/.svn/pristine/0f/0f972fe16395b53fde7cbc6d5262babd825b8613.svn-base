# ------------------------------------------------------------------------------
# NAME
#   Fcm::Timer
#
# DESCRIPTION
#   This is a package of timer utility used by the FCM command.
#
# COPYRIGHT
#   (C) Crown copyright Met Office. All rights reserved.
#   For further details please refer to the file COPYRIGHT.txt
#   which you should have received as part of this distribution.
# ------------------------------------------------------------------------------

package Fcm::Timer;

# Standard pragma
use warnings;
use strict;

# Exports
our (@ISA, @EXPORT, @EXPORT_OK);

sub timestamp_command;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(timestamp_command);

# ------------------------------------------------------------------------------

# Module level variables
my %cmd_start_time = (); # Command start time, (key = command, value = time)

# ------------------------------------------------------------------------------
# SYNOPSIS
#   $string = &Fcm::Timer::timestamp_command ($command[, $status]);
#
# DESCRIPTION
#   This function returns a string adding to $command a prefix according the
#   value of $status. If $status is not specified or does not match the word
#   "end", the status is assumed to be "start". At "start", the prefix will
#   contain the current timestamp. If $status is the word "end", the prefix
#   will contain the total time taken since this function was called with the
#   same $command at the "start" status.
# ------------------------------------------------------------------------------

sub timestamp_command {
  (my $command, my $status) = @_;

  my $prefix;
  if ($status and $status =~ /end/i) {
    # Status is "end", insert time taken
    my $lapse = time () - $cmd_start_time{$command};
    $prefix = sprintf "# Time taken: %12d s=> ", $lapse;

  } else {
    # Status is "start", insert time stamp
    $cmd_start_time{$command} = time;

    (my $sec, my $min, my $hour, my $mday, my $mon, my $year) = localtime;
    $prefix = sprintf "# Start: %04d-%02d-%02d %02d:%02d:%02d=> ",
              $year + 1900, $mon + 1, $mday, $hour, $min, $sec;
  }

  return $prefix . $command . "\n";
}

# ------------------------------------------------------------------------------

1;

__END__
