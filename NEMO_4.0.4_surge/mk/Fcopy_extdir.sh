#!/bin/bash
#set -x
set -o posix
#set -u
#set -e
#+
# ============
# Fcopy_extdir.sh
# ============
# --------------------------
# Copy a reference directory
# --------------------------
# SYNOPSIS
# ========
# ::
#  $ Fcopy_extdir.sh
# DESCRIPTION
# ===========
# When an unsupported configuration is requested, 
# Prepare sources for the NEMO sub-directories needed (OCE, TOP ...)
# EXAMPLES
# ========
# ::
#  $ ./Fcopy_extdir.sh ORCA2_LIM
# TODO
# ====
# option debug
# EVOLUTIONS
# ==========
# $Id: Fcopy_extdir.sh 3294 2012-01-28 16:44:18Z rblod $
#   * creation
#-
grep "$1 " ${CONFIG_DIR}/uspcfg.txt > ${CONFIG_DIR}/cfg.tmp
#
LOCAL_REF=$(cat cfg.tmp | awk 'BEGIN {FS = "#" }{print $2}')
TAB=$(cat cfg.tmp | awk 'BEGIN {FS = "#" }{print $3}')
REMOTE_CTL=$(cat cfg.tmp | awk 'BEGIN {FS = "#" }{print $4}')
#
\rm ${CONFIG_DIR}/cfg.tmp
