MODULE prtctl
   !!======================================================================
   !!                       ***  MODULE prtctl   ***
   !! Ocean system : print all SUM trends for each processor domain
   !!======================================================================
   !! History :  9.0  !  05-07  (C. Talandier) original code
   !!            3.4  !  11-11  (C. Harris) decomposition changes for running with CICE
   !!----------------------------------------------------------------------
   USE dom_oce          ! ocean space and time domain variables
#if defined key_nemocice_decomp
   USE ice_domain_size, only: nx_global, ny_global
#endif
   USE in_out_manager   ! I/O manager
   USE lib_mpp          ! distributed memory computing

   IMPLICIT NONE
   PRIVATE

   INTEGER , DIMENSION(:), ALLOCATABLE, SAVE ::   numid
   INTEGER , DIMENSION(:), ALLOCATABLE, SAVE ::   nlditl , nldjtl    ! first, last indoor index for each i-domain
   INTEGER , DIMENSION(:), ALLOCATABLE, SAVE ::   nleitl , nlejtl    ! first, last indoor index for each j-domain
   INTEGER , DIMENSION(:), ALLOCATABLE, SAVE ::   nimpptl, njmpptl   ! i-, j-indexes for each processor
   INTEGER , DIMENSION(:), ALLOCATABLE, SAVE ::   nlcitl , nlcjtl    ! dimensions of every subdomain
   INTEGER , DIMENSION(:), ALLOCATABLE, SAVE ::   ibonitl, ibonjtl   !

   REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::   t_ctll , s_ctll    ! previous tracer trend values
   REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::   u_ctll , v_ctll    ! previous velocity trend values

   INTEGER ::   ktime   ! time step

   PUBLIC prt_ctl         ! called by all subroutines
   PUBLIC prt_ctl_info    ! called by all subroutines
   PUBLIC prt_ctl_init    ! called by opa.F90
   PUBLIC sub_dom         ! called by opa.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE prt_ctl (tab2d_1, tab3d_1, mask1, clinfo1, tab2d_2, tab3d_2,   &
      &                                  mask2, clinfo2, kdim, clinfo3 )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE prt_ctl  ***
      !!
      !! ** Purpose : - print sum control of 2D or 3D arrays over the same area 
      !!                in mono and mpp case. This way can be usefull when
      !!                debugging a new parametrization in mono or mpp. 
      !!
      !! ** Method  : 2 possibilities exist when setting the ln_ctl parameter to
      !!                .true. in the ocean namelist:
      !!              - to debug a MPI run .vs. a mono-processor one; 
      !!                the control print will be done over each sub-domain.
      !!                The nictl[se] and njctl[se] parameters in the namelist must 
      !!                be set to zero and [ij]splt to the corresponding splitted
      !!                domain in MPI along respectively i-, j- directions.
      !!              - to debug a mono-processor run over the whole domain/a specific area; 
      !!                in the first case the nictl[se] and njctl[se] parameters must be set
      !!                to zero else to the indices of the area to be controled. In both cases
      !!                isplt and jsplt must be set to 1.
      !!              - All arguments of the above calling sequence are optional so their
      !!                name must be explicitly typed if used. For instance if the 3D
      !!                array tn(:,:,:) must be passed through the prt_ctl subroutine, 
      !!                it must looks like: CALL prt_ctl(tab3d_1=tn).
      !!
      !!                    tab2d_1 : first 2D array
      !!                    tab3d_1 : first 3D array
      !!                    mask1   : mask (3D) to apply to the tab[23]d_1 array
      !!                    clinfo1 : information about the tab[23]d_1 array
      !!                    tab2d_2 : second 2D array
      !!                    tab3d_2 : second 3D array
      !!                    mask2   : mask (3D) to apply to the tab[23]d_2 array
      !!                    clinfo2 : information about the tab[23]d_2 array
      !!                    kdim    : k- direction for 3D arrays 
      !!                    clinfo3 : additional information 
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:)  , INTENT(in), OPTIONAL ::   tab2d_1
      REAL(wp), DIMENSION(:,:,:), INTENT(in), OPTIONAL ::   tab3d_1
      REAL(wp), DIMENSION(:,:,:), INTENT(in), OPTIONAL ::   mask1
      CHARACTER (len=*)         , INTENT(in), OPTIONAL ::   clinfo1
      REAL(wp), DIMENSION(:,:)  , INTENT(in), OPTIONAL ::   tab2d_2
      REAL(wp), DIMENSION(:,:,:), INTENT(in), OPTIONAL ::   tab3d_2
      REAL(wp), DIMENSION(:,:,:), INTENT(in), OPTIONAL ::   mask2
      CHARACTER (len=*)         , INTENT(in), OPTIONAL ::   clinfo2
      INTEGER                   , INTENT(in), OPTIONAL ::   kdim
      CHARACTER (len=*)         , INTENT(in), OPTIONAL ::   clinfo3
      !
      CHARACTER (len=15) :: cl2
      INTEGER ::  jn, sind, eind, kdir,j_id
      REAL(wp) :: zsum1, zsum2, zvctl1, zvctl2
      REAL(wp), DIMENSION(jpi,jpj)     :: ztab2d_1, ztab2d_2
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zmask1, zmask2, ztab3d_1, ztab3d_2
      !!----------------------------------------------------------------------

      ! Arrays, scalars initialization 
      kdir      = jpkm1
      cl2       = ''
      zsum1     = 0.e0
      zsum2     = 0.e0
      zvctl1    = 0.e0
      zvctl2    = 0.e0
      ztab2d_1(:,:)   = 0.e0
      ztab2d_2(:,:)   = 0.e0
      ztab3d_1(:,:,:) = 0.e0
      ztab3d_2(:,:,:) = 0.e0
      zmask1  (:,:,:) = 1.e0
      zmask2  (:,:,:) = 1.e0

      ! Control of optional arguments
      IF( PRESENT(clinfo2) )   cl2                  = clinfo2
      IF( PRESENT(kdim)    )   kdir                 = kdim
      IF( PRESENT(tab2d_1) )   ztab2d_1(:,:)        = tab2d_1(:,:)
      IF( PRESENT(tab2d_2) )   ztab2d_2(:,:)        = tab2d_2(:,:)
      IF( PRESENT(tab3d_1) )   ztab3d_1(:,:,1:kdir) = tab3d_1(:,:,1:kdir)
      IF( PRESENT(tab3d_2) )   ztab3d_2(:,:,1:kdir) = tab3d_2(:,:,1:kdir)
      IF( PRESENT(mask1)   )   zmask1  (:,:,:)      = mask1  (:,:,:)
      IF( PRESENT(mask2)   )   zmask2  (:,:,:)      = mask2  (:,:,:)

      IF( lk_mpp .AND. jpnij > 1 ) THEN       ! processor number
         sind = narea
         eind = narea
      ELSE                                    ! processors total number
         sind = 1
         eind = ijsplt
      ENDIF

      ! Loop over each sub-domain, i.e. the total number of processors ijsplt
      DO jn = sind, eind
         ! Set logical unit
         j_id = numid(jn - narea + 1)
         ! Set indices for the SUM control
         IF( .NOT. lsp_area ) THEN
            IF (lk_mpp .AND. jpnij > 1)   THEN
               nictls = MAX(  1, nlditl(jn) )
               nictle = MIN(jpi, nleitl(jn) )
               njctls = MAX(  1, nldjtl(jn) )
               njctle = MIN(jpj, nlejtl(jn) )
               ! Do not take into account the bound of the domain
               IF( ibonitl(jn) == -1 .OR. ibonitl(jn) == 2 ) nictls = MAX(2, nictls)
               IF( ibonjtl(jn) == -1 .OR. ibonjtl(jn) == 2 ) njctls = MAX(2, njctls)
               IF( ibonitl(jn) ==  1 .OR. ibonitl(jn) == 2 ) nictle = MIN(nictle, nleitl(jn) - 1)
               IF( ibonjtl(jn) ==  1 .OR. ibonjtl(jn) == 2 ) njctle = MIN(njctle, nlejtl(jn) - 1)
            ELSE
               nictls = MAX(  1, nimpptl(jn) - 1 + nlditl(jn) )
               nictle = MIN(jpi, nimpptl(jn) - 1 + nleitl(jn) )
               njctls = MAX(  1, njmpptl(jn) - 1 + nldjtl(jn) )
               njctle = MIN(jpj, njmpptl(jn) - 1 + nlejtl(jn) )
               ! Do not take into account the bound of the domain
               IF( ibonitl(jn) == -1 .OR. ibonitl(jn) == 2 ) nictls = MAX(2, nictls)
               IF( ibonjtl(jn) == -1 .OR. ibonjtl(jn) == 2 ) njctls = MAX(2, njctls)
               IF( ibonitl(jn) ==  1 .OR. ibonitl(jn) == 2 ) nictle = MIN(nictle, nimpptl(jn) + nleitl(jn) - 2)
               IF( ibonjtl(jn) ==  1 .OR. ibonjtl(jn) == 2 ) njctle = MIN(njctle, njmpptl(jn) + nlejtl(jn) - 2)
            ENDIF
         ENDIF

         IF( PRESENT(clinfo3)) THEN
            IF ( clinfo3 == 'tra' )  THEN
               zvctl1 = t_ctll(jn)
               zvctl2 = s_ctll(jn)
            ELSEIF ( clinfo3 == 'dyn' )   THEN
               zvctl1 = u_ctll(jn)
               zvctl2 = v_ctll(jn)
            ENDIF
         ENDIF

         ! Compute the sum control
         ! 2D arrays
         IF( PRESENT(tab2d_1) )   THEN
            zsum1 = SUM( ztab2d_1(nictls:nictle,njctls:njctle)*zmask1(nictls:nictle,njctls:njctle,1) )
            zsum2 = SUM( ztab2d_2(nictls:nictle,njctls:njctle)*zmask2(nictls:nictle,njctls:njctle,1) )
         ENDIF

         ! 3D arrays
         IF( PRESENT(tab3d_1) )   THEN
            zsum1 = SUM( ztab3d_1(nictls:nictle,njctls:njctle,1:kdir)*zmask1(nictls:nictle,njctls:njctle,1:kdir) )
            zsum2 = SUM( ztab3d_2(nictls:nictle,njctls:njctle,1:kdir)*zmask2(nictls:nictle,njctls:njctle,1:kdir) )
         ENDIF

         ! Print the result
         IF( PRESENT(clinfo3) )   THEN
            WRITE(j_id,FMT='(a,D23.16,3x,a,D23.16)')clinfo1, zsum1-zvctl1, cl2, zsum2-zvctl2
            SELECT CASE( clinfo3 )
            CASE ( 'tra-ta' ) 
               t_ctll(jn) = zsum1
            CASE ( 'tra' ) 
                t_ctll(jn) = zsum1
                s_ctll(jn) = zsum2
            CASE ( 'dyn' ) 
                u_ctll(jn) = zsum1
                v_ctll(jn) = zsum2 
            END SELECT
         ELSEIF ( PRESENT(clinfo2) .OR. PRESENT(tab2d_2) .OR. PRESENT(tab3d_2) )   THEN
            WRITE(j_id,FMT='(a,D23.16,3x,a,D23.16)')clinfo1, zsum1, cl2, zsum2
         ELSE
            WRITE(j_id,FMT='(a,D23.16)')clinfo1, zsum1
         ENDIF

      ENDDO
      !
   END SUBROUTINE prt_ctl


   SUBROUTINE prt_ctl_info (clinfo1, ivar1, clinfo2, ivar2, itime)
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE prt_ctl_info  ***
      !!
      !! ** Purpose : - print information without any computation
      !!
      !! ** Action  : - input arguments
      !!                    clinfo1 : information about the ivar1
      !!                    ivar1   : value to print
      !!                    clinfo2 : information about the ivar2
      !!                    ivar2   : value to print
      !!----------------------------------------------------------------------
      CHARACTER (len=*), INTENT(in)           ::   clinfo1
      INTEGER          , INTENT(in), OPTIONAL ::   ivar1
      CHARACTER (len=*), INTENT(in), OPTIONAL ::   clinfo2
      INTEGER          , INTENT(in), OPTIONAL ::   ivar2
      INTEGER          , INTENT(in), OPTIONAL ::   itime
      !
      INTEGER :: jn, sind, eind, iltime, j_id
      !!----------------------------------------------------------------------

      IF( lk_mpp .AND. jpnij > 1 ) THEN       ! processor number
         sind = narea
         eind = narea
      ELSE                                    ! total number of processors
         sind = 1
         eind = ijsplt
      ENDIF

      ! Set to zero arrays at each new time step
      IF( PRESENT(itime) )   THEN
         iltime = itime
         IF( iltime > ktime )   THEN
            t_ctll(:) = 0.e0   ;   s_ctll(:) = 0.e0
            u_ctll(:) = 0.e0   ;   v_ctll(:) = 0.e0
            ktime = iltime
         ENDIF
      ENDIF

      ! Loop over each sub-domain, i.e. number of processors ijsplt
      DO jn = sind, eind
         !
         j_id = numid(jn - narea + 1)         ! Set logical unit
         !
         IF( PRESENT(ivar1) .AND. PRESENT(clinfo2) .AND. PRESENT(ivar2) )   THEN
            WRITE(j_id,*)clinfo1, ivar1, clinfo2, ivar2
         ELSEIF ( PRESENT(ivar1) .AND. PRESENT(clinfo2) .AND. .NOT. PRESENT(ivar2) )   THEN
            WRITE(j_id,*)clinfo1, ivar1, clinfo2
         ELSEIF ( PRESENT(ivar1) .AND. .NOT. PRESENT(clinfo2) .AND. PRESENT(ivar2) )   THEN
            WRITE(j_id,*)clinfo1, ivar1, ivar2
         ELSEIF ( PRESENT(ivar1) .AND. .NOT. PRESENT(clinfo2) .AND. .NOT. PRESENT(ivar2) )   THEN
            WRITE(j_id,*)clinfo1, ivar1
         ELSE
            WRITE(j_id,*)clinfo1
         ENDIF
         !
      END DO
      !
   END SUBROUTINE prt_ctl_info


   SUBROUTINE prt_ctl_init
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE prt_ctl_init  ***
      !!
      !! ** Purpose :   open ASCII files & compute indices
      !!----------------------------------------------------------------------
      INTEGER ::   jn, sind, eind, j_id
      CHARACTER (len=28) :: clfile_out
      CHARACTER (len=23) :: clb_name
      CHARACTER (len=19) :: cl_run
      !!----------------------------------------------------------------------

      ! Allocate arrays
      ALLOCATE( nlditl(ijsplt) , nleitl(ijsplt) , nimpptl(ijsplt) , ibonitl(ijsplt) ,   &
         &      nldjtl(ijsplt) , nlejtl(ijsplt) , njmpptl(ijsplt) , ibonjtl(ijsplt) ,   &
         &      nlcitl(ijsplt) , t_ctll(ijsplt) , u_ctll (ijsplt) ,                     &
         &      nlcjtl(ijsplt) , s_ctll(ijsplt) , v_ctll (ijsplt)                       )

      ! Initialization 
      t_ctll(:) = 0.e0
      s_ctll(:) = 0.e0
      u_ctll(:) = 0.e0
      v_ctll(:) = 0.e0
      ktime = 1

      IF( lk_mpp .AND. jpnij > 1 ) THEN
         sind = narea
         eind = narea
         clb_name = "('mpp.output_',I4.4)"
         cl_run = 'MULTI processor run'
         ! use indices for each area computed by mpp_init subroutine
         nlditl(1:jpnij) = nldit(:) 
         nleitl(1:jpnij) = nleit(:) 
         nldjtl(1:jpnij) = nldjt(:) 
         nlejtl(1:jpnij) = nlejt(:) 
         !
         nimpptl(1:jpnij) = nimppt(:)
         njmpptl(1:jpnij) = njmppt(:)
         !
         nlcitl(1:jpnij) = nlcit(:)
         nlcjtl(1:jpnij) = nlcjt(:)
         !
         ibonitl(1:jpnij) = ibonit(:)
         ibonjtl(1:jpnij) = ibonjt(:)
      ELSE
         sind = 1
         eind = ijsplt
         clb_name = "('mono.output_',I4.4)"
         cl_run = 'MONO processor run '
         ! compute indices for each area as done in mpp_init subroutine
         CALL sub_dom
      ENDIF

      ALLOCATE( numid(eind-sind+1) )

      DO jn = sind, eind
         WRITE(clfile_out,FMT=clb_name) jn-1
         CALL ctl_opn( numid(jn -narea + 1), clfile_out, 'REPLACE', 'FORMATTED', 'SEQUENTIAL', 1, numout, .FALSE. )
         j_id = numid(jn -narea + 1)
         WRITE(j_id,*)
         WRITE(j_id,*) '                 L O D Y C - I P S L'
         WRITE(j_id,*) '                     O P A model'
         WRITE(j_id,*) '            Ocean General Circulation Model'
         WRITE(j_id,*) '               version OPA 9.0  (2005) '
         WRITE(j_id,*)
         WRITE(j_id,*) '                   PROC number: ', jn
         WRITE(j_id,*)
         WRITE(j_id,FMT="(19x,a20)")cl_run

         ! Print the SUM control indices
         IF( .NOT. lsp_area )   THEN
            nictls = nimpptl(jn) + nlditl(jn) - 1
            nictle = nimpptl(jn) + nleitl(jn) - 1
            njctls = njmpptl(jn) + nldjtl(jn) - 1
            njctle = njmpptl(jn) + nlejtl(jn) - 1
         ENDIF
         WRITE(j_id,*) 
         WRITE(j_id,*) 'prt_ctl :  Sum control indices'
         WRITE(j_id,*) '~~~~~~~'
         WRITE(j_id,*)
         WRITE(j_id,9000)'                                nlej   = ', nlejtl(jn), '              '
         WRITE(j_id,9000)'                  ------------- njctle = ', njctle, ' -------------'
         WRITE(j_id,9001)'                  |                                       |'
         WRITE(j_id,9001)'                  |                                       |'
         WRITE(j_id,9001)'                  |                                       |'
         WRITE(j_id,9002)'           nictls = ', nictls,  '                           nictle = ', nictle
         WRITE(j_id,9002)'           nldi   = ', nlditl(jn),  '                           nlei   = ', nleitl(jn)
         WRITE(j_id,9001)'                  |                                       |'
         WRITE(j_id,9001)'                  |                                       |'
         WRITE(j_id,9001)'                  |                                       |'
         WRITE(j_id,9004)'  njmpp  = ',njmpptl(jn),'   ------------- njctls = ', njctls, ' -------------'
         WRITE(j_id,9003)'           nimpp  = ', nimpptl(jn), '        nldj   = ', nldjtl(jn), '              '
         WRITE(j_id,*)
         WRITE(j_id,*)

9000     FORMAT(a41,i4.4,a14)
9001     FORMAT(a59)
9002     FORMAT(a20,i4.4,a36,i4.4)
9003     FORMAT(a20,i4.4,a17,i4.4)
9004     FORMAT(a11,i4.4,a26,i4.4,a14)
      END DO
      !
   END SUBROUTINE prt_ctl_init


   SUBROUTINE sub_dom
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sub_dom  ***
      !!                    
      !! ** Purpose :   Lay out the global domain over processors. 
      !!                CAUTION: 
      !!                This part has been extracted from the mpp_init
      !!                subroutine and names of variables/arrays have been 
      !!                slightly changed to avoid confusion but the computation
      !!                is exactly the same. Any modification about indices of
      !!                each sub-domain in the mppini.F90 module should be reported 
      !!                here.
      !!
      !! ** Method  :   Global domain is distributed in smaller local domains.
      !!                Periodic condition is a function of the local domain position
      !!                (global boundary or neighbouring domain) and of the global
      !!                periodic
      !!                Type :         jperio global periodic condition
      !!
      !! ** Action  : - set domain parameters
      !!                    nimpp     : longitudinal index 
      !!                    njmpp     : latitudinal  index
      !!                    narea     : number for local area
      !!                    nlcil      : first dimension
      !!                    nlcjl      : second dimension
      !!                    nbondil    : mark for "east-west local boundary"
      !!                    nbondjl    : mark for "north-south local boundary"
      !!
      !! History :
      !!        !  94-11  (M. Guyon)  Original code
      !!        !  95-04  (J. Escobar, M. Imbard)
      !!        !  98-02  (M. Guyon)  FETI method
      !!        !  98-05  (M. Imbard, J. Escobar, L. Colombet )  SHMEM and MPI versions
      !!   8.5  !  02-08  (G. Madec)  F90 : free form
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj, jn               ! dummy loop indices
      INTEGER ::   &
         ii, ij,                         &  ! temporary integers
         irestil, irestjl,               &  !    "          "
         ijpi  , ijpj, nlcil,            &  ! temporary logical unit
         nlcjl , nbondil, nbondjl,       &
         nrecil, nrecjl, nldil, nleil, nldjl, nlejl

      INTEGER, DIMENSION(jpi,jpj) ::   iimpptl, ijmpptl, ilcitl, ilcjtl   ! workspace
      REAL(wp) ::   zidom, zjdom            ! temporary scalars
      INTEGER ::   inum                     ! local logical unit
      !!----------------------------------------------------------------------

      !
      !
      !  1. Dimension arrays for subdomains
      ! -----------------------------------
      !  Computation of local domain sizes ilcitl() ilcjtl()
      !  These dimensions depend on global sizes isplt,jsplt and jpiglo,jpjglo
      !  The subdomains are squares leeser than or equal to the global
      !  dimensions divided by the number of processors minus the overlap
      !  array (cf. par_oce.F90).

#if defined key_nemocice_decomp
      ijpi = ( nx_global+2-2*nn_hls + (isplt-1) ) / isplt + 2*nn_hls
      ijpj = ( ny_global+2-2*nn_hls + (jsplt-1) ) / jsplt + 2*nn_hls 
#else
      ijpi = ( jpiglo-2*nn_hls + (isplt-1) ) / isplt + 2*nn_hls
      ijpj = ( jpjglo-2*nn_hls + (jsplt-1) ) / jsplt + 2*nn_hls
#endif


      nrecil  = 2 * nn_hls
      nrecjl  = 2 * nn_hls
      irestil = MOD( jpiglo - nrecil , isplt )
      irestjl = MOD( jpjglo - nrecjl , jsplt )

      IF(  irestil == 0 )   irestil = isplt
#if defined key_nemocice_decomp

      ! In order to match CICE the size of domains in NEMO has to be changed
      ! The last line of blocks (west) will have fewer points 
      DO jj = 1, jsplt 
         DO ji=1, isplt-1 
            ilcitl(ji,jj) = ijpi 
         END DO 
         ilcitl(isplt,jj) = jpiglo - (isplt - 1) * (ijpi - nrecil)
      END DO 

#else 

      DO jj = 1, jsplt
         DO ji = 1, irestil
            ilcitl(ji,jj) = ijpi
         END DO
         DO ji = irestil+1, isplt
            ilcitl(ji,jj) = ijpi -1
         END DO
      END DO

#endif
      
      IF( irestjl == 0 )   irestjl = jsplt
#if defined key_nemocice_decomp 

      ! Same change to domains in North-South direction as in East-West. 
      DO ji = 1, isplt 
         DO jj=1, jsplt-1 
            ilcjtl(ji,jj) = ijpj 
         END DO 
         ilcjtl(ji,jsplt) = jpjglo - (jsplt - 1) * (ijpj - nrecjl)
      END DO 

#else 

      DO ji = 1, isplt
         DO jj = 1, irestjl
            ilcjtl(ji,jj) = ijpj
         END DO
         DO jj = irestjl+1, jsplt
            ilcjtl(ji,jj) = ijpj -1
         END DO
      END DO

#endif
      zidom = nrecil
      DO ji = 1, isplt
         zidom = zidom + ilcitl(ji,1) - nrecil
      END DO
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*)' sum ilcitl(i,1) = ', zidom, ' jpiglo = ', jpiglo
      
      zjdom = nrecjl
      DO jj = 1, jsplt
         zjdom = zjdom + ilcjtl(1,jj) - nrecjl
      END DO
      IF(lwp) WRITE(numout,*)' sum ilcitl(1,j) = ', zjdom, ' jpjglo = ', jpjglo
      IF(lwp) WRITE(numout,*)
      

      !  2. Index arrays for subdomains
      ! -------------------------------

      iimpptl(:,:) = 1
      ijmpptl(:,:) = 1
      
      IF( isplt > 1 ) THEN
         DO jj = 1, jsplt
            DO ji = 2, isplt
               iimpptl(ji,jj) = iimpptl(ji-1,jj) + ilcitl(ji-1,jj) - nrecil
            END DO
         END DO
      ENDIF

      IF( jsplt > 1 ) THEN
         DO jj = 2, jsplt
            DO ji = 1, isplt
               ijmpptl(ji,jj) = ijmpptl(ji,jj-1)+ilcjtl(ji,jj-1)-nrecjl
            END DO
         END DO
      ENDIF
      
      ! 3. Subdomain description
      ! ------------------------

      DO jn = 1, ijsplt
         ii = 1 + MOD( jn-1, isplt )
         ij = 1 + (jn-1) / isplt
         nimpptl(jn) = iimpptl(ii,ij)
         njmpptl(jn) = ijmpptl(ii,ij)
         nlcitl (jn) = ilcitl (ii,ij)     
         nlcil       = nlcitl (jn)     
         nlcjtl (jn) = ilcjtl (ii,ij)     
         nlcjl       = nlcjtl (jn)
         nbondjl = -1                                    ! general case
         IF( jn   >  isplt          )   nbondjl = 0      ! first row of processor
         IF( jn   >  (jsplt-1)*isplt )  nbondjl = 1     ! last  row of processor
         IF( jsplt == 1             )   nbondjl = 2      ! one processor only in j-direction
         ibonjtl(jn) = nbondjl
         
         nbondil = 0                                     ! 
         IF( MOD( jn, isplt ) == 1 )   nbondil = -1      !
         IF( MOD( jn, isplt ) == 0 )   nbondil =  1      !
         IF( isplt            == 1 )   nbondil =  2      ! one processor only in i-direction
         ibonitl(jn) = nbondil
         
         nldil =  1   + nn_hls
         nleil = nlcil - nn_hls
         IF( nbondil == -1 .OR. nbondil == 2 )   nldil = 1
         IF( nbondil ==  1 .OR. nbondil == 2 )   nleil = nlcil
         nldjl =  1   + nn_hls
         nlejl = nlcjl - nn_hls
         IF( nbondjl == -1 .OR. nbondjl == 2 )   nldjl = 1
         IF( nbondjl ==  1 .OR. nbondjl == 2 )   nlejl = nlcjl
         nlditl(jn) = nldil
         nleitl(jn) = nleil
         nldjtl(jn) = nldjl
         nlejtl(jn) = nlejl
      END DO
      !
      ! Save processor layout in layout_prtctl.dat file 
      IF(lwp) THEN
         CALL ctl_opn( inum, 'layout_prtctl.dat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE., narea )
         WRITE(inum,'(a)') 'nproc nlcil nlcjl nldil nldjl nleil nlejl nimpptl njmpptl ibonitl ibonjtl'
         !
         DO jn = 1, ijsplt
            WRITE(inum,'(i5,6i6,4i8)') jn-1,nlcitl(jn),  nlcjtl(jn), &
               &                            nlditl(jn),  nldjtl(jn), &
               &                            nleitl(jn),  nlejtl(jn), &
               &                           nimpptl(jn), njmpptl(jn), &
               &                           ibonitl(jn), ibonjtl(jn)
         END DO
         CLOSE(inum)   
      END IF
      !
      !
   END SUBROUTINE sub_dom

   !!======================================================================
END MODULE prtctl
