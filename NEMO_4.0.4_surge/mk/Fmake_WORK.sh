#!/bin/bash
######################################################
# Author : Rachid Benshila for NEMO
# Contact : rblod@locean-ipsl.upmc.fr
#
# Some functions called from makenemo
# Fmake_WORK      : create links in the WORK
######################################################
#set -vx
set -o posix
#set -u
#set -e
#+
#
# =============
# Fmake_WORK.sh
# =============
#
# -----------------------
# Make the WORK directory
# -----------------------
#
# SYNOPSIS
# ========
#
# ::
#
#  $ Fmake_WORK.sh
#
#
# DESCRIPTION
# ===========
#
#
# Make the WORK directory:
#
# - Create line in NEW_CONF/WORK
# - Use specified sub-directories previously
# - OPA has to be done first !!!
#
#
# EXAMPLES
# ========
#
# ::
#
#  $ ./Fmake_WORK.sh ORCA2_LIM OCE ICE
#
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
# $Id: Fmake_WORK.sh 9651 2018-05-28 06:47:14Z nicolasmartin $
#
#
#
#   * creation
#
#-
declare ZSRC=$1 ; shift 
declare ZCONF=$1 ; shift
ZTAB=( $@ )
declare i=0 ; declare NDIR=${#ZTAB[@]}

echo 'Creating '${ZCONF}'/WORK = '${ZTAB[*]}' for '${ZCONF}

[ ! -d ${ZCONF}/MY_SRC ] && \mkdir ${ZCONF}/MY_SRC
[   -d ${ZCONF}/WORK   ] || \mkdir ${ZCONF}/WORK

if [ "${ZSRC}" != 'none' ] ; then 

	if [ -d ${ZSRC} ] ; then 
		ln -sf ${ZSRC}/*.[Ffh]90 ${ZCONF}/MY_SRC/. 
		echo 'MY_SRC content is linked to '${ZSRC}
	else
		echo 'External directory for MY_SRC does not exist. Using default.'
	fi

else 
	echo 'MY_SRC directory is : '${ZCONF}'/MY_SRC'
fi

#\rm -f ../${1}/WORK/*

for comp in ${ZTAB[*]}; do
	find ${NEMO_DIR}/$comp -name *.[Ffh]90 -exec ln -sf {} ${ZCONF}/WORK \;
done

for i in `(cd ${ZCONF}/MY_SRC ; \ls *.[Ffh]90 2>/dev/null ) `; do
	[ -f ${ZCONF}/MY_SRC/$i ] &&  ln -sf $PWD/${ZCONF}/MY_SRC/${i} ${ZCONF}/WORK/.
done

unset -v ZCONF ZTAB i NDIR
