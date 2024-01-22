MODULE usrdef_nam
   !!======================================================================
   !!                     ***  MODULE usrdef_nam   ***
   !!
   !!                     ===  GYRE configuration  ===
   !!
   !! User defined : set the domain characteristics of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-03  (S. Flavoni, G. Madec)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_nam   : read user defined namelist and set global domain size
   !!   usr_def_hgr   : initialize the horizontal mesh 
   !!----------------------------------------------------------------------
   USE dom_oce  , ONLY: nimpp, njmpp       ! ocean space and time domain
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_nam   ! called in nemogcm.F90 module

   !                              !!* namusr_def namelist *!!
   LOGICAL, PUBLIC ::   ln_bench   ! =T benchmark test with gyre: the gridsize is constant (no need to adjust timestep or viscosity)
   INTEGER, PUBLIC ::   nn_GYRE    ! 1/nn_GYRE = the resolution chosen in degrees and thus defining the horizontal domain size

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_nam( cd_cfg, kk_cfg, kpi, kpj, kpk, kperio )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_nam  ***
      !!                    
      !! ** Purpose :   read user defined namelist and define the domain size
      !!
      !! ** Method  :   read in namusr_def containing all the user specific namelist parameter
      !!
      !!                Here GYRE configuration
      !!
      !! ** input   : - namusr_def namelist found in namelist_cfg
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(out) ::   cd_cfg          ! configuration name
      INTEGER         , INTENT(out) ::   kk_cfg          ! configuration resolution
      INTEGER         , INTENT(out) ::   kpi, kpj, kpk   ! global domain sizes 
      INTEGER         , INTENT(out) ::   kperio          ! lateral global domain b.c. 
      !
      INTEGER ::   ios   ! Local integer
      !!
      NAMELIST/namusr_def/ nn_GYRE, ln_bench, jpkglo
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_cfg )          ! Namelist namusr_def (exist in namelist_cfg only)
      READ  ( numnam_cfg, namusr_def, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namusr_def in configuration namelist' )
      !
      IF(lwm)   WRITE( numond, namusr_def )
      !
      cd_cfg = 'GYRE'               ! name & resolution (not used)
#if defined key_agrif
      IF (.NOT.Agrif_root()) nn_GYRE = Agrif_parent(nn_GYRE) * Agrif_irhox()
#endif
      kk_cfg = nn_GYRE
      !
      kpi = 30 * nn_GYRE + 2        ! Global Domain size
      kpj = 20 * nn_GYRE + 2
#if defined key_agrif
      IF( .NOT. Agrif_Root() ) THEN
         kpi  = nbcellsx + 2 + 2*nbghostcells
         kpj  = nbcellsy + 2 + 2*nbghostcells
      ENDIF
#endif
      kpk = jpkglo
      !                             ! Set the lateral boundary condition of the global domain
      kperio = 0                    ! GYRE configuration : closed domain
      !
      !                             ! control print
      IF(lwp) THEN
         WRITE(numout,*) '   '
         WRITE(numout,*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'
         WRITE(numout,*) '~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namusr_def : GYRE case'
         WRITE(numout,*) '      GYRE used as Benchmark (=T)                      ln_bench  = ', ln_bench
         WRITE(numout,*) '      inverse resolution & implied domain size         nn_GYRE   = ', nn_GYRE
#if defined key_agrif
         IF( Agrif_Root() ) THEN
#endif
         WRITE(numout,*) '         jpiglo = 30*nn_GYRE+2                            jpiglo = ', kpi
         WRITE(numout,*) '         jpjglo = 20*nn_GYRE+2                            jpjglo = ', kpj
#if defined key_agrif
         ENDIF
#endif
         WRITE(numout,*) '      number of model levels                              jpkglo = ', kpk
         WRITE(numout,*) '   '
         WRITE(numout,*) '   Lateral b.c. of the global domain set to closed        jperio = ', kperio
      ENDIF
      !
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam
