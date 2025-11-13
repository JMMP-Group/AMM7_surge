#!/bin/bash
#set -x
set -o posix
#set -u
#set -e
#+
# ===============
# Ffetch_extdir.sh
# ===============
# ---------------
# Make the config 
# ---------------
# SYNOPSIS
# ========
# ::
#  $ Ffetch_extdir.sh
# DESCRIPTION
# ===========
# - Make the config directory 
# - Create repositories needed :
#  
#  - EXP00 for namelist
#  - MY_SRC for user sources
#  - BLD for compilation 
# EXAMPLES
# ========
# ::
#  $ ./Ffetch_extdir.sh CONFIG_NAME REMOTE_CTL 
# TODO
# ====
# option debug
# EVOLUTIONS
# ==========
# $Id: Ffetch_extdir.sh 3715 2012-11-28 16:06:02Z acc $
#   * creation
#-
basedir=$(pwd)
cd ${1}
wget ${2} -O remote_file.list
#
if [ -f remote_file.list ] ; then
 cat remote_file.list | grep -v '^#' |
  while
   read remfile locfile
   do
    if [ $remfile == 'create_directory' ] ;then
      mkdir $locfile
    else
      wget $remfile -O $locfile
    fi
   done
fi
cd $basedir
