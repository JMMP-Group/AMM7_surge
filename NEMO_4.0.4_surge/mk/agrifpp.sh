#!/bin/bash
#set -x
set -o posix
#set -u
#set -e
#+
#
# ==========
# agrifpp.sh
# ==========
#
# ----------------------------
# Preform AGrif pre-processing
# ----------------------------
#
# SYNOPSIS
# ========
#
# ::
#
#  $ agrifpp.sh
#
#
# DESCRIPTION
# ===========
#
#
# Preprocess file using the conv in NEMOFILES directory
# Standard preprocessed files are stored in NEMOFILES/ppsrc/nemo
# Source files are stored under NEMOFILES/obj
# Include filess  in NEMOFILES/inc
# Note that agrif2model.F90 should not be preprocess (standard one) 
#
# EXAMPLES
# ========
#
# ::
#
#  $ ./agrifpp.sh FILE_TO_PROCESS
# 
# TODO
# ====
#
# option debug
#
#
# EVOLUTIONS
# ==========
#
# $Id: agrifpp.sh 2143 2010-10-04 12:49:55Z rblod $
#
#
#
#   * creation
#
#-
MYFILE=$(basename "$1")
if [ "$MYFILE" == "agrif2model.f90" ];then
   \cp ${NEMO_TDIR}/${NEW_CONF}/WORK/${MYFILE/.f90/.F90} ${NEMO_TDIR}/${NEW_CONF}/NEMOFILES/obj/$MYFILE
else
cd ${NEMO_TDIR}/${NEW_CONF}/NEMOFILES/ppsrc/nemo ; ${NEMO_TDIR}/${NEW_CONF}/NEMOFILES/conv ${NEMO_TDIR}/${NEW_CONF}/NEMOFILES/agrif_oce.in -rm -incdir ${NEMO_TDIR}/${NEW_CONF}/NEMOFILES/inc -comdirout ${NEMO_TDIR}/${NEW_CONF}/NEMOFILES/obj -convfile ${MYFILE} > /dev/null 
fi