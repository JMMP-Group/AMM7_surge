MODULE cpl_oasis3
   !!======================================================================
   !!                    ***  MODULE cpl_oasis  ***
   !! Coupled O/A : coupled ocean-atmosphere case using OASIS3-MCT
   !!=====================================================================
   !! History :  1.0  !  2004-06  (R. Redler, NEC Laboratories Europe, Germany) Original code
   !!             -   !  2004-11  (R. Redler, NEC Laboratories Europe; N. Keenlyside, W. Park, IFM-GEOMAR, Germany) revision
   !!             -   !  2004-11  (V. Gayler, MPI M&D) Grid writing
   !!            2.0  !  2005-08  (R. Redler, W. Park) frld initialization, paral(2) revision
   !!             -   !  2005-09  (R. Redler) extended to allow for communication over root only
   !!             -   !  2006-01  (W. Park) modification of physical part
   !!             -   !  2006-02  (R. Redler, W. Park) buffer array fix for root exchange
   !!            3.4  !  2011-11  (C. Harris) Changes to allow mutiple category fields
   !!            3.6  !  2014-11  (S. Masson) OASIS3-MCT
   !!----------------------------------------------------------------------
   
   !!----------------------------------------------------------------------
   !!   'key_oasis3'                    coupled Ocean/Atmosphere via OASIS3-MCT
   !!   'key_oa3mct_v3'                 to be added for OASIS3-MCT version 3
   !!----------------------------------------------------------------------
   !!   cpl_init     : initialization of coupled mode communication
   !!   cpl_define   : definition of grid and fields
   !!   cpl_snd      : snd out fields in coupled mode
   !!   cpl_rcv      : receive fields in coupled mode
   !!   cpl_finalize : finalize the coupled mode communication
   !!----------------------------------------------------------------------
#if defined key_oasis3
   USE mod_oasis                    ! OASIS3-MCT module
#endif 
#if defined key_iomput && defined key_xios25  
   USE xios, ONLY:xios_oasis_enddef
#endif
   USE par_oce                      ! ocean parameters
   USE dom_oce                      ! ocean space and time domain
   USE in_out_manager               ! I/O manager
   USE lbclnk                       ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   cpl_init
   PUBLIC   cpl_define
   PUBLIC   cpl_snd
   PUBLIC   cpl_rcv
   PUBLIC   cpl_freq
   PUBLIC   cpl_finalize

   INTEGER, PUBLIC            ::   OASIS_Rcv  = 1    !: return code if received field
   INTEGER, PUBLIC            ::   OASIS_idle = 0    !: return code if nothing done by oasis
   INTEGER                    ::   ncomp_id          ! id returned by oasis_init_comp
   INTEGER                    ::   nerror            ! return error code
#if ! defined key_oasis3
   ! OASIS Variables not used. defined only for compilation purpose
   INTEGER                    ::   OASIS_Out         = -1
   INTEGER                    ::   OASIS_REAL        = -1
   INTEGER                    ::   OASIS_Ok          = -1
   INTEGER                    ::   OASIS_In          = -1
   INTEGER                    ::   OASIS_Sent        = -1
   INTEGER                    ::   OASIS_SentOut     = -1
   INTEGER                    ::   OASIS_ToRest      = -1
   INTEGER                    ::   OASIS_ToRestOut   = -1
   INTEGER                    ::   OASIS_Recvd       = -1
   INTEGER                    ::   OASIS_RecvOut     = -1
   INTEGER                    ::   OASIS_FromRest    = -1
   INTEGER                    ::   OASIS_FromRestOut = -1
#endif

   INTEGER                    ::   nrcv         ! total number of fields received 
   INTEGER                    ::   nsnd         ! total number of fields sent 
   INTEGER                    ::   ncplmodel    ! Maximum number of models to/from which NEMO is potentialy sending/receiving data
   INTEGER, PUBLIC, PARAMETER ::   nmaxfld=60   ! Maximum number of coupling fields
   INTEGER, PUBLIC, PARAMETER ::   nmaxcat=5    ! Maximum number of coupling fields
   INTEGER, PUBLIC, PARAMETER ::   nmaxcpl=5    ! Maximum number of coupling fields
   LOGICAL, PARAMETER         ::   ltmp_wapatch = .TRUE.   ! patch to restore wraparound rows in cpl_send, cpl_rcv, cpl_define  
   INTEGER                    ::   nldi_save, nlei_save
   INTEGER                    ::   nldj_save, nlej_save
   
   TYPE, PUBLIC ::   FLD_CPL               !: Type for coupling field information
      LOGICAL               ::   laction   ! To be coupled or not
      CHARACTER(len = 8)    ::   clname    ! Name of the coupling field   
      CHARACTER(len = 1)    ::   clgrid    ! Grid type  
      REAL(wp)              ::   nsgn      ! Control of the sign change
      INTEGER, DIMENSION(nmaxcat,nmaxcpl) ::   nid   ! Id of the field (no more than 9 categories and 9 extrena models)
      INTEGER               ::   nct       ! Number of categories in field
      INTEGER               ::   ncplmodel ! Maximum number of models to/from which this variable may be sent/received
   END TYPE FLD_CPL

   TYPE(FLD_CPL), DIMENSION(nmaxfld), PUBLIC ::   srcv, ssnd   !: Coupling fields

   REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   exfld   ! Temporary buffer for receiving

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE cpl_init( cd_modname, kl_comm )
      !!-------------------------------------------------------------------
      !!             ***  ROUTINE cpl_init  ***
      !!
      !! ** Purpose :   Initialize coupled mode communication for ocean
      !!    exchange between AGCM, OGCM and COUPLER. (OASIS3 software)
      !!
      !! ** Method  :   OASIS3 MPI communication 
      !!--------------------------------------------------------------------
      CHARACTER(len = *), INTENT(in   ) ::   cd_modname   ! model name as set in namcouple file
      INTEGER           , INTENT(  out) ::   kl_comm      ! local communicator of the model
      !!--------------------------------------------------------------------

      ! WARNING: No write in numout in this routine
      !============================================

      !------------------------------------------------------------------
      ! 1st Initialize the OASIS system for the application
      !------------------------------------------------------------------
      CALL oasis_init_comp ( ncomp_id, TRIM(cd_modname), nerror )
      IF ( nerror /= OASIS_Ok ) &
         CALL oasis_abort (ncomp_id, 'cpl_init', 'Failure in oasis_init_comp')

      !------------------------------------------------------------------
      ! 3rd Get an MPI communicator for OPA local communication
      !------------------------------------------------------------------

      CALL oasis_get_localcomm ( kl_comm, nerror )
      IF ( nerror /= OASIS_Ok ) &
         CALL oasis_abort (ncomp_id, 'cpl_init','Failure in oasis_get_localcomm' )
      !
   END SUBROUTINE cpl_init


   SUBROUTINE cpl_define( krcv, ksnd, kcplmodel )
      !!-------------------------------------------------------------------
      !!             ***  ROUTINE cpl_define  ***
      !!
      !! ** Purpose :   Define grid and field information for ocean
      !!    exchange between AGCM, OGCM and COUPLER. (OASIS3 software)
      !!
      !! ** Method  :   OASIS3 MPI communication 
      !!--------------------------------------------------------------------
      INTEGER, INTENT(in) ::   krcv, ksnd     ! Number of received and sent coupling fields
      INTEGER, INTENT(in) ::   kcplmodel      ! Maximum number of models to/from which NEMO is potentialy sending/receiving data
      !
      INTEGER :: id_part
      INTEGER :: paral(5)       ! OASIS3 box partition
      INTEGER :: ishape(4)    ! shape of arrays passed to PSMILe
      INTEGER :: ji,jc,jm       ! local loop indicees
      CHARACTER(LEN=64) :: zclname
      CHARACTER(LEN=2) :: cli2
      !!--------------------------------------------------------------------

      ! patch to restore wraparound rows in cpl_send, cpl_rcv, cpl_define
      IF ( ltmp_wapatch ) THEN
         nldi_save = nldi   ;   nlei_save = nlei
         nldj_save = nldj   ;   nlej_save = nlej
         IF( nimpp           ==      1 ) nldi = 1
         IF( nimpp + jpi - 1 == jpiglo ) nlei = jpi
         IF( njmpp           ==      1 ) nldj = 1
         IF( njmpp + jpj - 1 == jpjglo ) nlej = jpj
      ENDIF 
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'cpl_define : initialization in coupled ocean/atmosphere case'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*)

      ncplmodel = kcplmodel
      IF( kcplmodel > nmaxcpl ) THEN
         CALL oasis_abort ( ncomp_id, 'cpl_define', 'ncplmodel is larger than nmaxcpl, increase nmaxcpl')   ;   RETURN
      ENDIF

      nrcv = krcv
      IF( nrcv > nmaxfld ) THEN
         CALL oasis_abort ( ncomp_id, 'cpl_define', 'nrcv is larger than nmaxfld, increase nmaxfld')   ;   RETURN
      ENDIF

      nsnd = ksnd
      IF( nsnd > nmaxfld ) THEN
         CALL oasis_abort ( ncomp_id, 'cpl_define', 'nsnd is larger than nmaxfld, increase nmaxfld')   ;   RETURN
      ENDIF
      !
      ! ... Define the shape for the area that excludes the halo
      !     For serial configuration (key_mpp_mpi not being active)
      !     nl* is set to the global values 1 and jp*glo.
      !
      ishape(1) = 1
      ishape(2) = nlei-nldi+1
      ishape(3) = 1
      ishape(4) = nlej-nldj+1
      !
      ! ... Allocate memory for data exchange
      !
      ALLOCATE(exfld(nlei-nldi+1, nlej-nldj+1), stat = nerror)
      IF( nerror > 0 ) THEN
         CALL oasis_abort ( ncomp_id, 'cpl_define', 'Failure in allocating exfld')   ;   RETURN
      ENDIF
      !
      ! -----------------------------------------------------------------
      ! ... Define the partition 
      ! -----------------------------------------------------------------
      
      paral(1) = 2                                              ! box partitioning
      paral(2) = jpiglo * (nldj-1+njmpp-1) + (nldi-1+nimpp-1)   ! NEMO lower left corner global offset    
      paral(3) = nlei-nldi+1                                    ! local extent in i 
      paral(4) = nlej-nldj+1                                    ! local extent in j
      paral(5) = jpiglo                                         ! global extent in x
      
      IF( ln_ctl ) THEN
         WRITE(numout,*) ' multiexchg: paral (1:5)', paral
         WRITE(numout,*) ' multiexchg: jpi, jpj =', jpi, jpj
         WRITE(numout,*) ' multiexchg: nldi, nlei, nimpp =', nldi, nlei, nimpp
         WRITE(numout,*) ' multiexchg: nldj, nlej, njmpp =', nldj, nlej, njmpp
      ENDIF
   
      CALL oasis_def_partition ( id_part, paral, nerror, jpiglo*jpjglo )
      !
      ! ... Announce send variables. 
      !
      ssnd(:)%ncplmodel = kcplmodel
      !
      DO ji = 1, ksnd
         IF ( ssnd(ji)%laction ) THEN

            IF( ssnd(ji)%nct > nmaxcat ) THEN
               CALL oasis_abort ( ncomp_id, 'cpl_define', 'Number of categories of '//   &
                  &              TRIM(ssnd(ji)%clname)//' is larger than nmaxcat, increase nmaxcat' )
               RETURN
            ENDIF
            
            DO jc = 1, ssnd(ji)%nct
               DO jm = 1, kcplmodel

                  IF ( ssnd(ji)%nct .GT. 1 ) THEN
                     WRITE(cli2,'(i2.2)') jc
                     zclname = TRIM(ssnd(ji)%clname)//'_cat'//cli2
                  ELSE
                     zclname = ssnd(ji)%clname
                  ENDIF
                  IF ( kcplmodel  > 1 ) THEN
                     WRITE(cli2,'(i2.2)') jm
                     zclname = 'model'//cli2//'_'//TRIM(zclname)
                  ENDIF
#if defined key_agrif
                  IF( agrif_fixed() /= 0 ) THEN 
                     zclname=TRIM(Agrif_CFixed())//'_'//TRIM(zclname)
                  END IF
#endif
                  IF( ln_ctl ) WRITE(numout,*) "Define", ji, jc, jm, " "//TRIM(zclname), " for ", OASIS_Out
                  CALL oasis_def_var (ssnd(ji)%nid(jc,jm), zclname, id_part   , (/ 2, 1 /),   &
                     &                OASIS_Out          , ishape , OASIS_REAL, nerror )
                  IF ( nerror /= OASIS_Ok ) THEN
                     WRITE(numout,*) 'Failed to define transient ', ji, jc, jm, " "//TRIM(zclname)
                     CALL oasis_abort ( ssnd(ji)%nid(jc,jm), 'cpl_define', 'Failure in oasis_def_var' )
                  ENDIF
                  IF( ln_ctl .AND. ssnd(ji)%nid(jc,jm) /= -1 ) WRITE(numout,*) "variable defined in the namcouple"
                  IF( ln_ctl .AND. ssnd(ji)%nid(jc,jm) == -1 ) WRITE(numout,*) "variable NOT defined in the namcouple"
               END DO
            END DO
         ENDIF
      END DO
      !
      ! ... Announce received variables. 
      !
      srcv(:)%ncplmodel = kcplmodel
      !
      DO ji = 1, krcv
         IF ( srcv(ji)%laction ) THEN 
            
            IF( srcv(ji)%nct > nmaxcat ) THEN
               CALL oasis_abort ( ncomp_id, 'cpl_define', 'Number of categories of '//   &
                  &              TRIM(srcv(ji)%clname)//' is larger than nmaxcat, increase nmaxcat' )
               RETURN
            ENDIF
            
            DO jc = 1, srcv(ji)%nct
               DO jm = 1, kcplmodel
                  
                  IF ( srcv(ji)%nct .GT. 1 ) THEN
                     WRITE(cli2,'(i2.2)') jc
                     zclname = TRIM(srcv(ji)%clname)//'_cat'//cli2
                  ELSE
                     zclname = srcv(ji)%clname
                  ENDIF
                  IF ( kcplmodel  > 1 ) THEN
                     WRITE(cli2,'(i2.2)') jm
                     zclname = 'model'//cli2//'_'//TRIM(zclname)
                  ENDIF
#if defined key_agrif
                  IF( agrif_fixed() /= 0 ) THEN 
                     zclname=TRIM(Agrif_CFixed())//'_'//TRIM(zclname)
                  END IF
#endif
                  IF( ln_ctl ) WRITE(numout,*) "Define", ji, jc, jm, " "//TRIM(zclname), " for ", OASIS_In
                  CALL oasis_def_var (srcv(ji)%nid(jc,jm), zclname, id_part   , (/ 2, 1 /),   &
                     &                OASIS_In           , ishape , OASIS_REAL, nerror )
                  IF ( nerror /= OASIS_Ok ) THEN
                     WRITE(numout,*) 'Failed to define transient ', ji, jc, jm, " "//TRIM(zclname)
                     CALL oasis_abort ( srcv(ji)%nid(jc,jm), 'cpl_define', 'Failure in oasis_def_var' )
                  ENDIF
                  IF( ln_ctl .AND. srcv(ji)%nid(jc,jm) /= -1 ) WRITE(numout,*) "variable defined in the namcouple"
                  IF( ln_ctl .AND. srcv(ji)%nid(jc,jm) == -1 ) WRITE(numout,*) "variable NOT defined in the namcouple"

               END DO
            END DO
         ENDIF
      END DO
      
      !------------------------------------------------------------------
      ! End of definition phase
      !------------------------------------------------------------------
      !     
#if defined key_agrif
      IF( agrif_fixed() == Agrif_Nb_Fine_Grids() ) THEN
#endif
      !  
#if defined key_iomput && defined key_xios25  
!needed for XIOS 2.5  
      CALL xios_oasis_enddef()  
#endif  
      ! 
      CALL oasis_enddef(nerror)
      IF( nerror /= OASIS_Ok )   CALL oasis_abort ( ncomp_id, 'cpl_define', 'Failure in oasis_enddef')
#if defined key_agrif
      ENDIF
#endif
      !
      IF ( ltmp_wapatch ) THEN
         nldi = nldi_save   ;   nlei = nlei_save
         nldj = nldj_save   ;   nlej = nlej_save
      ENDIF
   END SUBROUTINE cpl_define
   
   
   SUBROUTINE cpl_snd( kid, kstep, pdata, kinfo )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_snd  ***
      !!
      !! ** Purpose : - At each coupling time-step,this routine sends fields
      !!      like sst or ice cover to the coupler or remote application.
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   kid       ! variable index in the array
      INTEGER                   , INTENT(  out) ::   kinfo     ! OASIS3 info argument
      INTEGER                   , INTENT(in   ) ::   kstep     ! ocean time-step in seconds
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pdata
      !!
      INTEGER                                   ::   jc,jm     ! local loop index
      !!--------------------------------------------------------------------
      ! patch to restore wraparound rows in cpl_send, cpl_rcv, cpl_define
      IF ( ltmp_wapatch ) THEN
         nldi_save = nldi   ;   nlei_save = nlei
         nldj_save = nldj   ;   nlej_save = nlej
         IF( nimpp           ==      1 ) nldi = 1
         IF( nimpp + jpi - 1 == jpiglo ) nlei = jpi
         IF( njmpp           ==      1 ) nldj = 1
         IF( njmpp + jpj - 1 == jpjglo ) nlej = jpj
      ENDIF
      !
      ! snd data to OASIS3
      !
      DO jc = 1, ssnd(kid)%nct
         DO jm = 1, ssnd(kid)%ncplmodel
        
            IF( ssnd(kid)%nid(jc,jm) /= -1 ) THEN
               CALL oasis_put ( ssnd(kid)%nid(jc,jm), kstep, pdata(nldi:nlei, nldj:nlej,jc), kinfo )
               
               IF ( ln_ctl ) THEN        
                  IF ( kinfo == OASIS_Sent     .OR. kinfo == OASIS_ToRest .OR.   &
                     & kinfo == OASIS_SentOut  .OR. kinfo == OASIS_ToRestOut ) THEN
                     WRITE(numout,*) '****************'
                     WRITE(numout,*) 'oasis_put: Outgoing ', ssnd(kid)%clname
                     WRITE(numout,*) 'oasis_put: ivarid ', ssnd(kid)%nid(jc,jm)
                     WRITE(numout,*) 'oasis_put:  kstep ', kstep
                     WRITE(numout,*) 'oasis_put:   info ', kinfo
                     WRITE(numout,*) '     - Minimum value is ', MINVAL(pdata(nldi:nlei,nldj:nlej,jc))
                     WRITE(numout,*) '     - Maximum value is ', MAXVAL(pdata(nldi:nlei,nldj:nlej,jc))
                     WRITE(numout,*) '     -     Sum value is ',    SUM(pdata(nldi:nlei,nldj:nlej,jc))
                     WRITE(numout,*) '****************'
                  ENDIF
               ENDIF
               
            ENDIF
            
         ENDDO
      ENDDO
      IF ( ltmp_wapatch ) THEN
         nldi = nldi_save   ;   nlei = nlei_save
         nldj = nldj_save   ;   nlej = nlej_save
      ENDIF
      !
    END SUBROUTINE cpl_snd


   SUBROUTINE cpl_rcv( kid, kstep, pdata, pmask, kinfo )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_rcv  ***
      !!
      !! ** Purpose : - At each coupling time-step,this routine receives fields
      !!      like stresses and fluxes from the coupler or remote application.
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   kid       ! variable index in the array
      INTEGER                   , INTENT(in   ) ::   kstep     ! ocean time-step in seconds
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pdata     ! IN to keep the value if nothing is done
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pmask     ! coupling mask
      INTEGER                   , INTENT(  out) ::   kinfo     ! OASIS3 info argument
      !!
      INTEGER                                   ::   jc,jm     ! local loop index
      LOGICAL                                   ::   llaction, llfisrt
      !!--------------------------------------------------------------------
      ! patch to restore wraparound rows in cpl_send, cpl_rcv, cpl_define
      IF ( ltmp_wapatch ) THEN
         nldi_save = nldi   ;   nlei_save = nlei
         nldj_save = nldj   ;   nlej_save = nlej
      ENDIF
      !
      ! receive local data from OASIS3 on every process
      !
      kinfo = OASIS_idle
      !
      DO jc = 1, srcv(kid)%nct
         IF ( ltmp_wapatch ) THEN
            IF( nimpp           ==      1 ) nldi = 1
            IF( nimpp + jpi - 1 == jpiglo ) nlei = jpi
            IF( njmpp           ==      1 ) nldj = 1
            IF( njmpp + jpj - 1 == jpjglo ) nlej = jpj
         ENDIF
         llfisrt = .TRUE.

         DO jm = 1, srcv(kid)%ncplmodel

            IF( srcv(kid)%nid(jc,jm) /= -1 ) THEN

               CALL oasis_get ( srcv(kid)%nid(jc,jm), kstep, exfld, kinfo )         
               
               llaction =  kinfo == OASIS_Recvd   .OR. kinfo == OASIS_FromRest .OR.   &
                  &        kinfo == OASIS_RecvOut .OR. kinfo == OASIS_FromRestOut
               
               IF ( ln_ctl )   WRITE(numout,*) "llaction, kinfo, kstep, ivarid: " , llaction, kinfo, kstep, srcv(kid)%nid(jc,jm)
               
               IF ( llaction ) THEN
                  
                  kinfo = OASIS_Rcv
                  IF( llfisrt ) THEN 
                     pdata(nldi:nlei,nldj:nlej,jc) =                                 exfld(:,:) * pmask(nldi:nlei,nldj:nlej,jm)
                     llfisrt = .FALSE.
                  ELSE
                     pdata(nldi:nlei,nldj:nlej,jc) = pdata(nldi:nlei,nldj:nlej,jc) + exfld(:,:) * pmask(nldi:nlei,nldj:nlej,jm)
                  ENDIF
                  
                  IF ( ln_ctl ) THEN        
                     WRITE(numout,*) '****************'
                     WRITE(numout,*) 'oasis_get: Incoming ', srcv(kid)%clname
                     WRITE(numout,*) 'oasis_get: ivarid '  , srcv(kid)%nid(jc,jm)
                     WRITE(numout,*) 'oasis_get:   kstep', kstep
                     WRITE(numout,*) 'oasis_get:   info ', kinfo
                     WRITE(numout,*) '     - Minimum value is ', MINVAL(pdata(nldi:nlei,nldj:nlej,jc))
                     WRITE(numout,*) '     - Maximum value is ', MAXVAL(pdata(nldi:nlei,nldj:nlej,jc))
                     WRITE(numout,*) '     -     Sum value is ',    SUM(pdata(nldi:nlei,nldj:nlej,jc))
                     WRITE(numout,*) '****************'
                  ENDIF
                  
               ENDIF
               
            ENDIF
            
         ENDDO

         IF ( ltmp_wapatch ) THEN
            nldi = nldi_save   ;   nlei = nlei_save
            nldj = nldj_save   ;   nlej = nlej_save
         ENDIF
         !--- Fill the overlap areas and extra hallows (mpp)
         !--- check periodicity conditions (all cases)
         IF( .not. llfisrt ) THEN
            CALL lbc_lnk( 'cpl_oasis3', pdata(:,:,jc), srcv(kid)%clgrid, srcv(kid)%nsgn )   
         ENDIF
 
      ENDDO
      !
   END SUBROUTINE cpl_rcv


   INTEGER FUNCTION cpl_freq( cdfieldname )  
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_freq  ***
      !!
      !! ** Purpose : - send back the coupling frequency for a particular field
      !!----------------------------------------------------------------------
      CHARACTER(len = *), INTENT(in) ::   cdfieldname    ! field name as set in namcouple file
      !!
      INTEGER               :: id
      INTEGER               :: info
      INTEGER, DIMENSION(1) :: itmp
      INTEGER               :: ji,jm     ! local loop index
      INTEGER               :: mop
      !!----------------------------------------------------------------------
      cpl_freq = 0   ! defaut definition
      id = -1        ! defaut definition
      !
      DO ji = 1, nsnd
         IF (ssnd(ji)%laction ) THEN
            DO jm = 1, ncplmodel
               IF( ssnd(ji)%nid(1,jm) /= -1 ) THEN
                  IF( TRIM(cdfieldname) == TRIM(ssnd(ji)%clname) ) THEN
                     id = ssnd(ji)%nid(1,jm)
                     mop = OASIS_Out
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      DO ji = 1, nrcv
         IF (srcv(ji)%laction ) THEN
            DO jm = 1, ncplmodel
               IF( srcv(ji)%nid(1,jm) /= -1 ) THEN
                  IF( TRIM(cdfieldname) == TRIM(srcv(ji)%clname) ) THEN
                     id = srcv(ji)%nid(1,jm)
                     mop = OASIS_In
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      !
      IF( id /= -1 ) THEN
#if ! defined key_oa3mct_v1v2
         CALL oasis_get_freqs(id, mop, 1, itmp, info)
#else
         CALL oasis_get_freqs(id,      1, itmp, info)
#endif
         cpl_freq = itmp(1)
      ENDIF
      !
   END FUNCTION cpl_freq


   SUBROUTINE cpl_finalize
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_finalize  ***
      !!
      !! ** Purpose : - Finalizes the coupling. If MPI_init has not been
      !!      called explicitly before cpl_init it will also close
      !!      MPI communication.
      !!----------------------------------------------------------------------
      !
      DEALLOCATE( exfld )
      IF (nstop == 0) THEN
         CALL oasis_terminate( nerror )         
      ELSE
         CALL oasis_abort( ncomp_id, "cpl_finalize", "NEMO ABORT STOP" )
      ENDIF       
      !
   END SUBROUTINE cpl_finalize

#if ! defined key_oasis3

   !!----------------------------------------------------------------------
   !!   No OASIS Library          OASIS3 Dummy module...
   !!----------------------------------------------------------------------

   SUBROUTINE oasis_init_comp(k1,cd1,k2)
      CHARACTER(*), INTENT(in   ) ::  cd1
      INTEGER     , INTENT(  out) ::  k1,k2
      k1 = -1 ; k2 = -1
      WRITE(numout,*) 'oasis_init_comp: Error you sould not be there...', cd1
   END SUBROUTINE oasis_init_comp

   SUBROUTINE oasis_abort(k1,cd1,cd2)
      INTEGER     , INTENT(in   ) ::  k1
      CHARACTER(*), INTENT(in   ) ::  cd1,cd2
      WRITE(numout,*) 'oasis_abort: Error you sould not be there...', cd1, cd2
   END SUBROUTINE oasis_abort

   SUBROUTINE oasis_get_localcomm(k1,k2)
      INTEGER     , INTENT(  out) ::  k1,k2
      k1 = -1 ; k2 = -1
      WRITE(numout,*) 'oasis_get_localcomm: Error you sould not be there...'
   END SUBROUTINE oasis_get_localcomm

   SUBROUTINE oasis_def_partition(k1,k2,k3,k4)
      INTEGER     , INTENT(  out) ::  k1,k3
      INTEGER     , INTENT(in   ) ::  k2(5)
      INTEGER     , INTENT(in   ) ::  k4
      k1 = k2(1) ; k3 = k2(5)+k4
      WRITE(numout,*) 'oasis_def_partition: Error you sould not be there...'
   END SUBROUTINE oasis_def_partition

   SUBROUTINE oasis_def_var(k1,cd1,k2,k3,k4,k5,k6,k7)
      CHARACTER(*), INTENT(in   ) ::  cd1
      INTEGER     , INTENT(in   ) ::  k2,k3(2),k4,k5(2,2),k6
      INTEGER     , INTENT(  out) ::  k1,k7
      k1 = -1 ; k7 = -1
      WRITE(numout,*) 'oasis_def_var: Error you sould not be there...', cd1
   END SUBROUTINE oasis_def_var

   SUBROUTINE oasis_enddef(k1)
      INTEGER     , INTENT(  out) ::  k1
      k1 = -1
      WRITE(numout,*) 'oasis_enddef: Error you sould not be there...'
   END SUBROUTINE oasis_enddef
  
   SUBROUTINE oasis_put(k1,k2,p1,k3)
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::  p1
      INTEGER                 , INTENT(in   ) ::  k1,k2
      INTEGER                 , INTENT(  out) ::  k3
      k3 = -1
      WRITE(numout,*) 'oasis_put: Error you sould not be there...'
   END SUBROUTINE oasis_put

   SUBROUTINE oasis_get(k1,k2,p1,k3)
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::  p1
      INTEGER                 , INTENT(in   ) ::  k1,k2
      INTEGER                 , INTENT(  out) ::  k3
      p1(1,1) = -1. ; k3 = -1
      WRITE(numout,*) 'oasis_get: Error you sould not be there...'
   END SUBROUTINE oasis_get

   SUBROUTINE oasis_get_freqs(k1,k5,k2,k3,k4)
      INTEGER              , INTENT(in   ) ::  k1,k2
      INTEGER, DIMENSION(1), INTENT(  out) ::  k3
      INTEGER              , INTENT(  out) ::  k4,k5
      k3(1) = k1 ; k4 = k2 ; k5 = k2
      WRITE(numout,*) 'oasis_get_freqs: Error you sould not be there...'
   END SUBROUTINE oasis_get_freqs

   SUBROUTINE oasis_terminate(k1)
      INTEGER     , INTENT(  out) ::  k1
      k1 = -1
      WRITE(numout,*) 'oasis_terminate: Error you sould not be there...'
   END SUBROUTINE oasis_terminate
   
#endif

   !!=====================================================================
END MODULE cpl_oasis3
