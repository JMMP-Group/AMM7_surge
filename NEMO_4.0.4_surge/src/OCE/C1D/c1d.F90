MODULE c1d
   !!======================================================================
   !!                     ***  MODULE  c1d  ***
   !! Ocean domain  :  1D configuration
   !!=====================================================================
   !! History :  2.0  !  2004-09 (C. Ethe)     Original code
   !!            3.0  !  2008-04 (G. Madec)    adaptation to SBC
   !!            3.5  !  2013-10 (D. Calvert)  add namelist
   !!----------------------------------------------------------------------
#if defined key_c1d
   !!----------------------------------------------------------------------
   !!   'key_c1d'                                   1D column configuration
   !!----------------------------------------------------------------------
   !!   c1d_init      : read in the C1D namelist
   !!----------------------------------------------------------------------
   USE par_kind       ! kind parameters
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   c1d_init   ! called by nemogcm.F90

   LOGICAL , PUBLIC, PARAMETER ::  lk_c1d = .TRUE.   ! 1D config. flag

   REAL(wp), PUBLIC ::  rn_lat1d     !: Column latitude
   REAL(wp), PUBLIC ::  rn_lon1d     !: Column longitude
   LOGICAL , PUBLIC ::  ln_c1d_locpt !: Localization (or not) of 1D column in a grid

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!======================================================================
CONTAINS

   SUBROUTINE c1d_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE c1d_init  ***
      !! 
      !! ** Purpose :   Initialization of C1D options
      !!
      !! ** Method  :   Read namelist namc1d 
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !!
      NAMELIST/namc1d/ rn_lat1d, rn_lon1d , ln_c1d_locpt
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_ref )              ! Namelist namc1d in reference namelist : Tracer advection scheme
      READ  ( numnam_ref, namc1d, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namc1d in reference namelist' )
      !
      REWIND( numnam_cfg )              ! Namelist namtra_adv in configuration namelist : Tracer advection scheme
      READ  ( numnam_cfg, namc1d, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namc1d in configuration namelist' )
      IF(lwm) WRITE ( numond, namc1d )
      !
      IF(lwp) THEN                    ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'c1d_init : Initialize 1D model configuration options'
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) '   Namelist namc1d : set options for the C1D model'
         WRITE(numout,*) '      column latitude                 rn_lat1d     = ', rn_lat1d
         WRITE(numout,*) '      column longitude                rn_lon1d     = ', rn_lon1d
         WRITE(numout,*) '      column localization in a grid   ln_c1d_locpt = ', ln_c1d_locpt
      ENDIF
      !
   END SUBROUTINE c1d_init

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                           No use of 1D configuration
   !!----------------------------------------------------------------------
   USE par_kind         ! kind parameters
   LOGICAL, PUBLIC, PARAMETER ::   lk_c1d = .FALSE.   !: 1D config. flag de-activated
   REAL(wp)                   ::   rn_lat1d, rn_lon1d
   LOGICAL , PUBLIC           ::   ln_c1d_locpt = .FALSE. 
CONTAINS
   SUBROUTINE c1d_init               ! Dummy routine
   END SUBROUTINE c1d_init
#endif

   !!======================================================================
END MODULE c1d
