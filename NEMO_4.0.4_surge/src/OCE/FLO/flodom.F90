MODULE flodom
   !!======================================================================
   !!                       ***  MODULE  flodom  ***
   !! Ocean floats :   domain
   !!======================================================================
   !! History :  OPA  ! 1998-07 (Y.Drillet, CLIPPER)  Original code
   !!  NEMO      3.3  ! 2011-09 (C.Bricaud,S.Law-Chune Mercator-Ocean): add ARIANE convention + comsecitc changes
   !!----------------------------------------------------------------------
   !!   flo_dom               : initialization of floats
   !!   add_new_floats        : add new floats (long/lat/depth)
   !!   add_new_ariane_floats : add new floats with araine convention (i/j/k)
   !!   findmesh              : compute index of position 
   !!   dstnce                : compute distance between face mesh and floats 
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE flo_oce         ! ocean drifting floats
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   flo_dom         ! routine called by floats.F90
   PUBLIC   flo_dom_alloc   ! Routine called in floats.F90

   CHARACTER (len=21) ::  clname1 = 'init_float'              ! floats initialisation filename
   CHARACTER (len=21) ::  clname2 = 'init_float_ariane'       ! ariane floats initialisation filename


   INTEGER , ALLOCATABLE, DIMENSION(:) ::   iimfl, ijmfl, ikmfl       ! index mesh of floats
   INTEGER , ALLOCATABLE, DIMENSION(:) ::   idomfl, ivtest, ihtest    !   -     
   REAL(wp), ALLOCATABLE, DIMENSION(:) ::   zgifl, zgjfl,  zgkfl      ! distances in indexes

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE flo_dom
      !! ---------------------------------------------------------------------
      !!                  ***  ROUTINE flo_dom  ***
      !!                 
      !!  ** Purpose :   Initialisation of floats
      !!
      !!  ** Method  :   We put the floats  in the domain with the latitude,
      !!               the longitude (degree) and the depth (m).
      !!----------------------------------------------------------------------      
      INTEGER            ::   jfl    ! dummy loop  
      INTEGER            ::   inum   ! logical unit for file read
      !!---------------------------------------------------------------------
      
      ! Initialisation with the geographical position or restart
      
      IF(lwp) WRITE(numout,*) 'flo_dom : compute initial position of floats'
      IF(lwp) WRITE(numout,*) '~~~~~~~~'
      IF(lwp) WRITE(numout,*) '           jpnfl = ',jpnfl
      
      !-------------------------!
      ! FLOAT RESTART FILE READ !
      !-------------------------!
      IF( ln_rstflo )THEN

         IF(lwp) WRITE(numout,*) '        float restart file read'
         
         ! open the restart file 
         !----------------------
         CALL ctl_opn( inum, 'restart_float', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp )

         ! read of the restart file
         READ(inum,*)   ( tpifl  (jfl), jfl=1, jpnrstflo),   & 
                        ( tpjfl  (jfl), jfl=1, jpnrstflo),   &
                        ( tpkfl  (jfl), jfl=1, jpnrstflo),   &
                        ( nisobfl(jfl), jfl=1, jpnrstflo),   &
                        ( ngrpfl (jfl), jfl=1, jpnrstflo)    
         CLOSE(inum)

         ! if we want a  surface drift  ( like PROVOR floats )
         IF( ln_argo ) nisobfl(1:jpnrstflo) = 0
         
         ! It is possible to add new floats.          
         !---------------------------------
         IF( jpnfl > jpnrstflo )THEN

            IF(lwp) WRITE(numout,*) '        add new floats'

            IF( ln_ariane )THEN  !Add new floats with ariane convention
                CALL flo_add_new_ariane_floats(jpnrstflo+1,jpnfl) 
            ELSE                 !Add new floats with long/lat convention
                CALL flo_add_new_floats(jpnrstflo+1,jpnfl)
            ENDIF
         ENDIF

      !--------------------------------------!
      ! FLOAT INITILISATION: NO RESTART FILE !
      !--------------------------------------!
      ELSE    !ln_rstflo

         IF( ln_ariane )THEN       !Add new floats with ariane convention
            CALL flo_add_new_ariane_floats(1,jpnfl)
         ELSE                      !Add new floats with long/lat convention
            CALL flo_add_new_floats(1,jpnfl)
         ENDIF

      ENDIF
            
   END SUBROUTINE flo_dom

   SUBROUTINE flo_add_new_floats(kfl_start, kfl_end)
      !! -------------------------------------------------------------
      !!                 ***  SUBROUTINE add_new_arianefloats  ***
      !!          
      !! ** Purpose :   
      !!
      !!       First initialisation of floats
      !!       the initials positions of floats are written in a file
      !!       with a variable to know if it is a isobar float a number 
      !!       to identified who want the trajectories of this float and 
      !!       an index for the number of the float         
      !!       open the init file 
      !!               
      !! ** Method  : 
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kfl_start, kfl_end
      !!
      INTEGER           :: inum ! file unit
      INTEGER           :: jfl,ji, jj, jk ! dummy loop indices
      INTEGER           :: itrash         ! trash var for reading
      INTEGER           :: ifl            ! number of floats to read
      REAL(wp)          :: zdxab, zdyad
      LOGICAL           :: llinmesh
      CHARACTER(len=80) :: cltmp
      !!---------------------------------------------------------------------
      ifl = kfl_end-kfl_start+1

      ! we get the init values 
      !-----------------------
      CALL ctl_opn( inum , clname1, 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp )
      DO jfl = kfl_start,kfl_end
         READ(inum,*) flxx(jfl),flyy(jfl),flzz(jfl), nisobfl(jfl),ngrpfl(jfl),itrash
         if(lwp)write(numout,*)'read:',jfl,flxx(jfl),flyy(jfl),flzz(jfl), nisobfl(jfl),ngrpfl(jfl),itrash ; call flush(numout)
      END DO
      CLOSE(inum)
            
      ! Test to find the grid point coordonate with the geographical position            
      !----------------------------------------------------------------------
      DO jfl = kfl_start,kfl_end
         ihtest(jfl) = 0
         ivtest(jfl) = 0
         ikmfl(jfl) = 0
# if   defined key_mpp_mpi
         DO ji = MAX(nldi,2), nlei
            DO jj = MAX(nldj,2), nlej   ! NO vector opt.
# else         
         DO ji = 2, jpi
            DO jj = 2, jpj   ! NO vector opt.
# endif                     
               ! For each float we find the indexes of the mesh                      
               CALL flo_findmesh(glamf(ji-1,jj-1),gphif(ji-1,jj-1),   &
                                 glamf(ji-1,jj  ),gphif(ji-1,jj  ),   &
                                 glamf(ji  ,jj  ),gphif(ji  ,jj  ),   &
                                 glamf(ji  ,jj-1),gphif(ji  ,jj-1),   &
                                 flxx(jfl)       ,flyy(jfl)       ,   &
                                 glamt(ji  ,jj  ),gphit(ji  ,jj  ), llinmesh)
               IF( llinmesh )THEN
                  iimfl(jfl) = ji
                  ijmfl(jfl) = jj
                  ihtest(jfl) = ihtest(jfl)+1
                  DO jk = 1, jpk-1
                     IF( (gdepw_n(ji,jj,jk) <= flzz(jfl)) .AND. (gdepw_n(ji,jj,jk+1) > flzz(jfl)) ) THEN
                        ikmfl(jfl) = jk
                        ivtest(jfl) = ivtest(jfl) + 1
                     ENDIF
                  END DO
               ENDIF
            END DO
         END DO

         ! If the float is in a mesh computed by an other processor we put iimfl=ijmfl=-1               
         IF( ihtest(jfl) ==  0 ) THEN
            iimfl(jfl) = -1
            ijmfl(jfl) = -1
         ENDIF
      END DO

      !Test if each float is in one and only one proc
      !----------------------------------------------
      IF( lk_mpp )   THEN 
         CALL mpp_sum('flodom', ihtest,jpnfl)
         CALL mpp_sum('flodom', ivtest,jpnfl)
      ENDIF
      DO jfl = kfl_start,kfl_end

         IF( (ihtest(jfl) > 1 ) .OR. ( ivtest(jfl) > 1) ) THEN
             WRITE(cltmp,'(A10,i4.4,A20)' )'THE FLOAT',jfl,' IS NOT IN ONLY ONE MESH'
             CALL ctl_stop('STOP',TRIM(cltmp) )
         ENDIF
         IF( (ihtest(jfl) == 0) ) THEN
             WRITE(cltmp,'(A10,i4.4,A20)' )'THE FLOAT',jfl,' IS IN NO MESH'
             CALL ctl_stop('STOP',TRIM(cltmp) )
         ENDIF
      END DO

      ! We compute the distance between the float and the face of the mesh            
      !-------------------------------------------------------------------
      DO jfl = kfl_start,kfl_end

         ! Made only if the float is in the domain of the processor              
         IF( (iimfl(jfl) >= 0) .AND. (ijmfl(jfl) >= 0) ) THEN

            ! TEST TO KNOW IF THE FLOAT IS NOT INITIALISED IN THE COAST
            idomfl(jfl) = 0
            IF( tmask(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)) == 0. ) idomfl(jfl) = 1

            ! Computation of the distance between the float and the faces of the mesh
            !            zdxab
            !             .
            !        B----.---------C
            !        |    .         |
            !        |<------>flo   |
            !        |        ^     |
            !        |        |.....|....zdyad
            !        |        |     |
            !        A--------|-----D
            !
            zdxab = flo_dstnce( flxx(jfl), flyy(jfl), glamf(iimfl(jfl)-1,ijmfl(jfl)-1), flyy(jfl) )
            zdyad = flo_dstnce( flxx(jfl), flyy(jfl), flxx(jfl), gphif(iimfl(jfl)-1,ijmfl(jfl)-1) )

            ! Translation of this distances (in meter) in indexes
            zgifl(jfl)= (iimfl(jfl)-0.5) + zdxab/e1u(iimfl(jfl)-1,ijmfl(jfl)) + (mig(1)-1)
            zgjfl(jfl)= (ijmfl(jfl)-0.5) + zdyad/e2v(iimfl(jfl),ijmfl(jfl)-1) + (mjg(1)-1)
            zgkfl(jfl) = (( gdepw_n(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)+1) - flzz(jfl) )* ikmfl(jfl))   &
               &                 / (  gdepw_n(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)+1)                              &
               &                    - gdepw_n(iimfl(jfl),ijmfl(jfl),ikmfl(jfl) ) )                             &
               &                 + (( flzz(jfl)-gdepw_n(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)) ) *(ikmfl(jfl)+1))   &
               &                 / (  gdepw_n(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)+1)                              &
               &                    - gdepw_n(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)) )
         ELSE
            zgifl(jfl) = 0.e0
            zgjfl(jfl) = 0.e0
            zgkfl(jfl) = 0.e0
         ENDIF

      END DO
                  
      ! The sum of all the arrays zgifl, zgjfl, zgkfl give 3 arrays with the positions of all the floats.
      IF( lk_mpp )   THEN 
         CALL mpp_sum( 'flodom', zgjfl, ifl )   ! sums over the global domain
         CALL mpp_sum( 'flodom', zgkfl, ifl )
      ENDIF
            
      DO jfl = kfl_start,kfl_end
         tpifl(jfl) = zgifl(jfl)
         tpjfl(jfl) = zgjfl(jfl)
         tpkfl(jfl) = zgkfl(jfl)
      END DO

      ! WARNING : initial position not in the sea         
      IF( .NOT. ln_rstflo ) THEN 
         DO jfl =  kfl_start,kfl_end
            IF( idomfl(jfl) == 1 ) THEN
               IF(lwp) WRITE(numout,*)'*****************************'
               IF(lwp) WRITE(numout,*)'!!!!!!!  WARNING   !!!!!!!!!!'
               IF(lwp) WRITE(numout,*)'*****************************'
               IF(lwp) WRITE(numout,*)'The float number',jfl,'is out of the sea.'
               IF(lwp) WRITE(numout,*)'geographical position',flxx(jfl),flyy(jfl),flzz(jfl)
               IF(lwp) WRITE(numout,*)'index position',tpifl(jfl),tpjfl(jfl),tpkfl(jfl)
            ENDIF
         END DO
      ENDIF

   END SUBROUTINE flo_add_new_floats

   SUBROUTINE flo_add_new_ariane_floats(kfl_start, kfl_end)
      !! -------------------------------------------------------------
      !!                 ***  SUBROUTINE add_new_arianefloats  ***
      !!          
      !! ** Purpose :   
      !!       First initialisation of floats with ariane convention
      !!       
      !!       The indexes are read directly from file (warning ariane
      !!       convention, are refered to 
      !!       U,V,W grids - and not T-) 
      !!       The isobar advection is managed with the sign of tpkfl ( >0 -> 3D
      !!       advection, <0 -> 2D) 
      !!       Some variables are not read, as - gl         : time index; 4th
      !!       column        
      !!                                       - transport  : transport ; 5th
      !!                                       column
      !!       and paste in the jtrash var
      !!       At the end, ones need to replace the indexes on T grid
      !!       RMQ : there is no float groups identification !
      !!
      !!               
      !! ** Method  : 
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kfl_start, kfl_end
      !!
      INTEGER  :: inum         ! file unit
      INTEGER  :: ierr, ifl
      INTEGER  :: jfl, jfl1    ! dummy loop indices
      INTEGER  :: itrash       ! trash var for reading  
      CHARACTER(len=80) :: cltmp

      !!----------------------------------------------------------------------
      nisobfl(kfl_start:kfl_end) = 1 ! we assume that by default we want 3D advection

      ifl = kfl_end - kfl_start + 1  ! number of floats to read  

      ! we check that the number of floats in the init_file are consistant with the namelist
      IF( lwp ) THEN

         jfl1=0
         ierr=0
         CALL ctl_opn( inum, clname2, 'OLD', 'FORMATTED', 'SEQUENTIAL',  1, numout, .TRUE., 1 )
         DO WHILE (ierr .EQ. 0)
            jfl1=jfl1+1
            READ(inum,*, iostat=ierr)
         END DO
         CLOSE(inum)
         IF( (jfl1-1) .NE. ifl )THEN 
            WRITE(cltmp,'(A25,A20,A3,i4.4,A10,i4.4)')"the number of floats in ",TRIM(clname2), &
                                                     " = ",jfl1," is not equal to jfl= ",ifl
            CALL ctl_stop('STOP',TRIM(cltmp) )
         ENDIF

      ENDIF
            
      ! we get the init values 
      CALL ctl_opn( inum, clname2, 'OLD', 'FORMATTED', 'SEQUENTIAL', 1, numout, .TRUE., 1 )
      DO jfl = kfl_start, kfl_end
          READ(inum,*) tpifl(jfl),tpjfl(jfl),tpkfl(jfl),itrash, itrash
              
          IF ( tpkfl(jfl) .LT. 0. ) nisobfl(jfl) = 0 !set the 2D advection according to init_float
          ngrpfl(jfl)=jfl
      END DO

      ! conversion from ariane index to T grid index
      tpkfl(kfl_start:kfl_end) = abs(tpkfl)-0.5 ! reversed vertical axis
      tpifl(kfl_start:kfl_end) = tpifl+0.5
      tpjfl(kfl_start:kfl_end) = tpjfl+0.5


   END SUBROUTINE flo_add_new_ariane_floats


   SUBROUTINE flo_findmesh( pax, pay, pbx, pby,   &
                            pcx, pcy, pdx, pdy,   &
                            px  ,py  ,ptx, pty, ldinmesh )
      !! -------------------------------------------------------------
      !!                ***  ROUTINE findmesh  ***
      !!      
      !! ** Purpose :   Find the index of mesh for the point spx spy
      !!
      !! ** Method  : 
      !!----------------------------------------------------------------------
      REAL(wp) ::   &
         pax, pay, pbx, pby,    &     ! ???
         pcx, pcy, pdx, pdy,    &     ! ???
         px, py,                &     ! longitude and latitude
         ptx, pty                     ! ???
      LOGICAL ::  ldinmesh            ! ???
      !!
      REAL(wp) ::   zabt, zbct, zcdt, zdat, zabpt, zbcpt, zcdpt, zdapt
      !!---------------------------------------------------------------------
      !! Statement function
      REAL(wp) ::   fsline
      REAL(wp) ::   psax, psay, psbx, psby, psx, psy
      fsline( psax, psay, psbx, psby, psx, psy ) = psy  * ( psbx - psax )   &
         &                                       - psx  * ( psby - psay )   &
         &                                       + psax *   psby - psay * psbx
      !!---------------------------------------------------------------------
      
      ! 4 semi plane defined by the 4 points and including the T point
      zabt = fsline(pax,pay,pbx,pby,ptx,pty)
      zbct = fsline(pbx,pby,pcx,pcy,ptx,pty)
      zcdt = fsline(pcx,pcy,pdx,pdy,ptx,pty)
      zdat = fsline(pdx,pdy,pax,pay,ptx,pty)
      
      ! 4 semi plane defined by the 4 points and including the extrememity
      zabpt = fsline(pax,pay,pbx,pby,px,py)
      zbcpt = fsline(pbx,pby,pcx,pcy,px,py)
      zcdpt = fsline(pcx,pcy,pdx,pdy,px,py)
      zdapt = fsline(pdx,pdy,pax,pay,px,py)
       
      ! We compare the semi plane T with the semi plane including the point
      ! to know if it is in this  mesh.
      ! For numerical reasons it is possible that for a point which is on
      ! the line we don't have exactly zero with fsline function. We want 
      ! that a point can't be in 2 mesh in the same time, so we put the 
      ! coefficient to zero if it is smaller than 1.E-12
      
      IF( ABS(zabpt) <= 1.E-12 ) zabpt = 0.
      IF( ABS(zbcpt) <= 1.E-12 ) zbcpt = 0.
      IF( ABS(zcdpt) <= 1.E-12 ) zcdpt = 0.
      IF( ABS(zdapt) <= 1.E-12 ) zdapt = 0.
      IF( (zabt*zabpt >  0.) .AND. (zbct*zbcpt >= 0. ) .AND. ( zcdt*zcdpt >= 0. ) .AND. ( zdat*zdapt > 0. )   &
         .AND. ( px <= MAX(pcx,pdx) ) .AND. ( px >= MIN(pax,pbx) )    &
         .AND. ( py <= MAX(pby,pcy) ) .AND. ( py >= MIN(pay,pdy) ) ) THEN
         ldinmesh=.TRUE.
      ELSE
         ldinmesh=.FALSE.
      ENDIF
      !
   END SUBROUTINE flo_findmesh


   FUNCTION flo_dstnce( pla1, phi1, pla2, phi2 )
      !! -------------------------------------------------------------
      !!                 ***  Function dstnce  ***
      !!          
      !! ** Purpose :   returns distance (in m) between two geographical
      !!                points
      !! ** Method  : 
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::   pla1, phi1, pla2, phi2   ! ???
      !!
      REAL(wp) :: dly1, dly2, dlx1, dlx2, dlx, dls, dld, dpi
      REAL(wp) :: flo_dstnce
      !!---------------------------------------------------------------------
      !
      dpi  = 2._wp * ASIN(1._wp)
      dls  = dpi / 180._wp
      dly1 = phi1 * dls
      dly2 = phi2 * dls
      dlx1 = pla1 * dls
      dlx2 = pla2 * dls
      !
      dlx = SIN(dly1) * SIN(dly2) + COS(dly1) * COS(dly2) * COS(dlx2-dlx1)
      !
      IF( ABS(dlx) > 1.0_wp ) dlx = 1.0_wp
      !
      dld = ATAN(SQRT( 1._wp * ( 1._wp-dlx )/( 1._wp+dlx ) )) * 222.24_wp / dls
      flo_dstnce = dld * 1000._wp
      !
   END FUNCTION flo_dstnce

   INTEGER FUNCTION flo_dom_alloc()
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION flo_dom_alloc  ***
      !!----------------------------------------------------------------------

      ALLOCATE( iimfl(jpnfl) , ijmfl(jpnfl) , ikmfl(jpnfl) ,                          &  
                idomfl(jpnfl), ivtest(jpnfl), ihtest(jpnfl),                 &
                zgifl(jpnfl) , zgjfl(jpnfl) , zgkfl(jpnfl)   , STAT=flo_dom_alloc )
      !
      CALL mpp_sum ( 'flodom', flo_dom_alloc )
      IF( flo_dom_alloc /= 0 )   CALL ctl_stop( 'STOP', 'flo_dom_alloc: failed to allocate arrays' )
   END FUNCTION flo_dom_alloc

   !!======================================================================
END MODULE flodom
