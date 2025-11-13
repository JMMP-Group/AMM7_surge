MODULE defprec 
!-
! $Id: defprec.f90 2281 2010-10-15 14:21:13Z smasson $
!-
! This software is governed by the CeCILL license
! See IOIPSL/IOIPSL_License_CeCILL.txt
!!--------------------------------------------------------------------
!! The module "defprec" set default precision for computation
!!
!! This module should be used by every modules
!! to keep the right precision for every variable
!!--------------------------------------------------------------------
!?INTEGERS of KIND 1 are not supported on all computers
!?INTEGER,PARAMETER :: i_1=SELECTED_INT_KIND(2)
  INTEGER,PARAMETER :: i_2=SELECTED_INT_KIND(4)
  INTEGER,PARAMETER :: i_4=SELECTED_INT_KIND(9)
  INTEGER,PARAMETER :: i_8=SELECTED_INT_KIND(13)
  INTEGER,PARAMETER :: r_4=SELECTED_REAL_KIND(6,37)
  INTEGER,PARAMETER :: r_8=SELECTED_REAL_KIND(15,307)
  INTEGER,PARAMETER :: i_std=i_4, r_std=r_8
!-----------------
END MODULE defprec
