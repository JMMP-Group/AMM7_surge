MODULE iscplini
   !!======================================================================
   !!                       ***  MODULE  sbciscpl  ***
   !! Ocean forcing:  ?????
   !!=====================================================================
   !! History :  NEMO  ! 2015-01 P. Mathiot: original 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   iscpl_init     : initialisation routine (namelist)
   !!   iscpl_alloc    : allocation of correction variables
   !!----------------------------------------------------------------------
   USE oce             ! global tra/dyn variable
   USE dom_oce         ! ocean space and time domain
   !
   USE lib_mpp         ! MPP library
   USE lib_fortran     ! MPP library
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE
   
   PUBLIC   iscpl_init      
   PUBLIC   iscpl_alloc 
   
   !                                 !!* namsbc_iscpl namelist *
   LOGICAL , PUBLIC ::   ln_hsb       !:
   INTEGER , PUBLIC ::   nn_fiscpl    !:
   INTEGER , PUBLIC ::   nn_drown     !:
   
   INTEGER , PUBLIC ::   nstp_iscpl   !:
   REAL(wp), PUBLIC ::   rdt_iscpl    !: 
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   hdiv_iscpl   !:
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   htsc_iscpl   !:

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION iscpl_alloc()
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE sbc_iscpl_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( htsc_iscpl(jpi,jpj,jpk,jpts) , hdiv_iscpl(jpi,jpj,jpk) , STAT=iscpl_alloc )
         !
      CALL mpp_sum ( 'iscplini', iscpl_alloc )
      IF( iscpl_alloc > 0 )   CALL ctl_warn('iscpl_alloc: allocation of arrays failed')
   END FUNCTION iscpl_alloc


   SUBROUTINE iscpl_init()
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      INTEGER ::   ios           ! Local integer output status for namelist read
      NAMELIST/namsbc_iscpl/ nn_fiscpl, ln_hsb, nn_drown
      !!----------------------------------------------------------------------
      !
      nn_fiscpl = 0
      ln_hsb    = .FALSE.
      REWIND( numnam_ref )              ! Namelist namsbc_iscpl in reference namelist : Ice sheet coupling
      READ  ( numnam_ref, namsbc_iscpl, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_iscpl in reference namelist' )
      REWIND( numnam_cfg )              ! Namelist namsbc_iscpl in configuration namelist : Ice Sheet coupling
      READ  ( numnam_cfg, namsbc_iscpl, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_iscpl in configuration namelist' )
      IF(lwm) WRITE ( numond, namsbc_iscpl )
      !
      nstp_iscpl=MIN( nn_fiscpl, nitend-nit000+1 ) ! the coupling period have to be less or egal than the total number of time step
      rdt_iscpl = nstp_iscpl * rn_rdt
      !
      IF (lwp) THEN
         WRITE(numout,*) 'iscpl_rst:'
         WRITE(numout,*) '~~~~~~~~~'
         WRITE(numout,*) ' coupling     flag (ln_iscpl )            = ', ln_iscpl
         WRITE(numout,*) ' conservation flag (ln_hsb   )            = ', ln_hsb
         WRITE(numout,*) ' nb of stp for cons (rn_fiscpl)           = ', nstp_iscpl
         IF (nstp_iscpl .NE. nn_fiscpl) WRITE(numout,*) 'W A R N I N G: nb of stp for cons has been modified &
            &                                           (larger than run length)'
         WRITE(numout,*) ' coupling time step                       = ', rdt_iscpl
         WRITE(numout,*) ' number of call of the extrapolation loop = ', nn_drown
      ENDIF
      !
   END SUBROUTINE iscpl_init

   !!======================================================================
END MODULE iscplini
