MODULE julian
   !!======================================================================
   !!                       ***  MODULE julian   ***
   !! Ocean          : Julian data utilities
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   jul2greg        : Convert relative time to date
   !!   greg2jul        : Convert date to relative time
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind, ONLY : &       ! Precision variables
      & wp, &
      & dp  
   !USE in_out_manager           ! I/O manager
   USE lib_mpp,  ONLY : &
      & ctl_warn, ctl_stop

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE
   PUBLIC jul2greg,        &  ! Convert relative time to date
      &   greg2jul            ! Convert date to relative time 
  
   !! $Id$
CONTAINS
 
#include "jul2greg.h90"

#include "greg2jul.h90"
   
END MODULE julian
