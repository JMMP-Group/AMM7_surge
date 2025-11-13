MODULE flowri
   !!======================================================================
   !!                       ***  MODULE  flowri  ***
   !!
   !! Ocean floats: write floats trajectory in ascii                    ln_flo_ascii = T
   !!                                    or in netcdf ( IOM or IOSPSL ) ln_flo_ascii = F           
   !!======================================================================
   !!  History :  OPA  !  1999-09  (Y. Drillet)    : Original code
   !!              -   !  2000-06  (J.-M. Molines) : Profiling floats for CLS 
   !!   NEMO      1.0  !  2002-10  (A. Bozec)  F90 : Free form and module
   !!             3.2  !  2010-08  (slaw, cbricaud): netcdf outputs and others 
   !!----------------------------------------------------------------------
   USE flo_oce         ! ocean drifting floats
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE lib_mpp         ! distribued memory computing library
   USE in_out_manager  ! I/O manager
   USE phycst          ! physic constants
   USE dianam          ! build name of file (routine)
   USE ioipsl
   USE iom             ! I/O library

   IMPLICIT NONE
   PRIVATE

   PUBLIC flo_wri         ! routine called by floats.F90
   PUBLIC flo_wri_alloc   ! routine called by floats.F90

   INTEGER :: jfl                            ! number of floats
   CHARACTER (len=80)  :: clname             ! netcdf output filename

   REAL(wp), ALLOCATABLE, DIMENSION(:) ::   zlon , zlat, zdep   ! 2D workspace
   REAL(wp), ALLOCATABLE, DIMENSION(:) ::   ztem , zsal, zrho   ! 2D workspace

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION flo_wri_alloc()
      !!-------------------------------------------------------------------
      !!                ***  FUNCTION flo_wri_alloc  ***
      !!-------------------------------------------------------------------
      ALLOCATE( ztem(jpnfl) , zsal(jpnfl) , zrho(jpnfl) , &
                zlon(jpnfl) , zlat(jpnfl) , zdep(jpnfl) , STAT=flo_wri_alloc)
      !  
      CALL mpp_sum ( 'flowri', flo_wri_alloc )
      IF( flo_wri_alloc /= 0 )   CALL ctl_stop( 'STOP', 'flo_wri_alloc: failed to allocate arrays.' )
   END FUNCTION flo_wri_alloc

   SUBROUTINE flo_wri( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE flo_wri ***
      !!             
      !! ** Purpose :   Write position of floats in "trajec_float.nc",according
      !!                to ARIANE TOOLS (http://stockage.univ-brest.fr/~grima/Ariane/ )  n
      !!                nomenclature
      !!    
      !!      
      !! ** Method  :   The frequency of  ??? is nwritefl
      !!      
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER  :: kt                               ! time step

      !! * Local declarations
      INTEGER  :: iafl , ibfl , icfl             ! temporary integer
      INTEGER  :: ia1fl, ib1fl, ic1fl            !   "
      INTEGER  :: iafloc,ibfloc,ia1floc,ib1floc  !   "
      INTEGER  :: irec, irecflo

      REAL(wp) :: zafl,zbfl,zcfl                 ! temporary real
      REAL(wp) :: ztime                          !   "

      INTEGER, DIMENSION(2)          :: icount
      INTEGER, DIMENSION(2)          :: istart
      INTEGER, DIMENSION(1)          :: ish
      INTEGER, DIMENSION(2)          :: ish2
      !!----------------------------------------------------------------------
      
      !-----------------------------------------------------
      ! I- Save positions, temperature, salinty and density 
      !-----------------------------------------------------
      zlon(:)=0.0 ; zlat(:)=0.0 ; zdep(:)=0.0 
      ztem(:)=0.0 ; zsal(:)=0.0 ; zrho(:)=0.0 

      DO jfl = 1, jpnfl

         iafl  = INT (tpifl(jfl))            ! I-index of the nearest point before
         ibfl  = INT (tpjfl(jfl))            ! J-index of the nearest point before
         icfl  = INT (tpkfl(jfl))            ! K-index of the nearest point before
         ia1fl = iafl + 1                    ! I-index of the nearest point after
         ib1fl = ibfl + 1                    ! J-index of the nearest point after
         ic1fl = icfl + 1                    ! K-index of the nearest point after
         zafl  = tpifl(jfl) - REAL(iafl,wp)  ! distance  ?????
         zbfl  = tpjfl(jfl) - REAL(ibfl,wp)  ! distance  ?????
         zcfl  = tpkfl(jfl) - REAL(icfl,wp)  ! distance  ?????

         IF( lk_mpp ) THEN
               
            iafloc = mi1( iafl )
            ibfloc = mj1( ibfl )
 
            IF( nldi <= iafloc .AND. iafloc <= nlei .AND. &
              & nldj <= ibfloc .AND. ibfloc <= nlej       ) THEN 

               !the float is inside of current proc's area
               ia1floc = iafloc + 1
               ib1floc = ibfloc + 1
     
               !save position of the float
               zlat(jfl) = (1.-zafl)*(1.-zbfl)*gphit(iafloc ,ibfloc ) + (1.-zafl) * zbfl * gphit(iafloc ,ib1floc)   &
                     +     zafl *(1.-zbfl)*gphit(ia1floc,ibfloc ) +     zafl  * zbfl * gphit(ia1floc,ib1floc)   
               zlon(jfl) = (1.-zafl)*(1.-zbfl)*glamt(iafloc ,ibfloc ) + (1.-zafl) * zbfl * glamt(iafloc ,ib1floc)   &
                     +     zafl *(1.-zbfl)*glamt(ia1floc,ibfloc ) +     zafl  * zbfl * glamt(ia1floc,ib1floc)
               zdep(jfl) = (1.-zcfl)*gdepw_n(iafloc,ibfloc,icfl ) + zcfl * gdepw_n(iafloc,ibfloc,ic1fl)     

               !save temperature, salinity and density at this position
               ztem(jfl) = tsn(iafloc,ibfloc,icfl,jp_tem)
               zsal (jfl) = tsn(iafloc,ibfloc,icfl,jp_sal)
               zrho (jfl) = (rhd(iafloc,ibfloc,icfl)+1)*rau0

            ENDIF

         ELSE  ! mono proc case  

            iafloc  = iafl
            ibfloc  = ibfl
            ia1floc = iafloc + 1
            ib1floc = ibfloc + 1

            !save position of the float               
            zlat(jfl) = (1.-zafl)*(1.-zbfl)*gphit(iafloc ,ibfloc ) + (1.-zafl) * zbfl * gphit(iafloc ,ib1floc)   &
                      +     zafl *(1.-zbfl)*gphit(ia1floc,ibfloc ) +     zafl  * zbfl * gphit(ia1floc,ib1floc)
            zlon(jfl) = (1.-zafl)*(1.-zbfl)*glamt(iafloc ,ibfloc ) + (1.-zafl) * zbfl * glamt(iafloc ,ib1floc)   &
                      +     zafl *(1.-zbfl)*glamt(ia1floc,ibfloc ) +     zafl  * zbfl * glamt(ia1floc,ib1floc)
            zdep(jfl) = (1.-zcfl)*gdepw_n(iafloc,ibfloc,icfl ) + zcfl * gdepw_n(iafloc,ibfloc,ic1fl)

            ztem(jfl) = tsn(iafloc,ibfloc,icfl,jp_tem)
            zsal(jfl) = tsn(iafloc,ibfloc,icfl,jp_sal)
            zrho(jfl) = (rhd(iafloc,ibfloc,icfl)+1)*rau0
          
         ENDIF

      END DO ! loop on float

      !Only proc 0 writes all positions : SUM of positions on all procs
      IF( lk_mpp ) THEN
         CALL mpp_sum( 'flowri', zlon, jpnfl )   ! sums over the global domain
         CALL mpp_sum( 'flowri', zlat, jpnfl )   ! sums over the global domain
         CALL mpp_sum( 'flowri', zdep, jpnfl )   ! sums over the global domain
         CALL mpp_sum( 'flowri', ztem, jpnfl )   ! sums over the global domain
         CALL mpp_sum( 'flowri', zsal, jpnfl )   ! sums over the global domain
         CALL mpp_sum( 'flowri', zrho, jpnfl )   ! sums over the global domain
      ENDIF


      !-------------------------------------!
      ! II- WRITE WRITE WRITE WRITE WRITE   !
      !-------------------------------------!

      !--------------------------!
      ! II-1 Write in ascii file !
      !--------------------------!

      IF( ln_flo_ascii )THEN

         IF( ( kt == nn_it000 .OR. MOD( kt,nn_writefl)== 0 ) .AND. lwp )THEN

            !II-1-a Open ascii file
            !----------------------
            IF( kt == nn_it000 ) THEN
               CALL ctl_opn( numflo, 'trajec_float', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE. )
               irecflo = NINT( (nitend-nn_it000) / FLOAT(nn_writefl) )
               WRITE(numflo,*) cexper, irecflo, jpnfl, nn_writefl
            ENDIF

            !II-1-b Write in ascii file
            !-----------------------------
            WRITE(numflo,*) zlon,zlat,zdep,nisobfl,ngrpfl,ztem,zsal, FLOAT(ndastp)


            !II-1-c Close netcdf file
            !-------------------------
            IF( kt == nitend )   CLOSE( numflo )

         ENDIF

      !-----------------------------------------------------
      ! II-2 Write in netcdf file
      !-----------------------------------------------------

      ELSE

      !II-2-a Write with IOM
      !----------------------

#if defined key_iomput
         CALL iom_put( "traj_lon"     , zlon )
         CALL iom_put( "traj_lat"     , zlat )
         CALL iom_put( "traj_dep"     , zdep )
         CALL iom_put( "traj_temp"    , ztem )
         CALL iom_put( "traj_salt"    , zsal  )
         CALL iom_put( "traj_dens"    , zrho )
         CALL iom_put( "traj_group"   , REAL(ngrpfl,wp) )
#else

      !II-2-b Write with IOIPSL
      !------------------------

         IF( ( kt == nn_it000 .OR. MOD( kt,nn_writefl)== 0 ) .AND. lwp )THEN


            !II-2-b-1 Open netcdf file
            !-------------------------
            IF( kt==nn_it000 )THEN   ! Create and open

               CALL dia_nam( clname, nn_writefl, 'trajec_float' )
               clname=TRIM(clname)//".nc"

               CALL fliocrfd( clname , (/'ntraj' , '    t' /), (/ jpnfl , -1/) , numflo )
   
               CALL fliodefv( numflo, 'traj_lon'    , (/1,2/), v_t=flio_r8, long_name="Longitude"           , units="degrees_east"  )
               CALL fliodefv( numflo, 'traj_lat'    , (/1,2/), v_t=flio_r8, long_name="Latitude"            , units="degrees_north" )
               CALL fliodefv( numflo, 'traj_depth'  , (/1,2/), v_t=flio_r8, long_name="Depth"               , units="meters" )
               CALL fliodefv( numflo, 'time_counter', (/2/)  , v_t=flio_r8, long_name="Time axis"           & 
                         & , units="seconds since start of the run " )
               CALL fliodefv( numflo, 'traj_temp'   , (/1,2/), v_t=flio_r8, long_name="Temperature"         , units="C" )
               CALL fliodefv( numflo, 'traj_salt'   , (/1,2/), v_t=flio_r8, long_name="Salinity"            , units="PSU" )
               CALL fliodefv( numflo, 'traj_dens'   , (/1,2/), v_t=flio_r8, long_name="Density"             , units="kg/m3" )
               CALL fliodefv( numflo, 'traj_group'  , (/1/)  , v_t=flio_r8, long_name="number of the group" , units="no unit" )

               CALL flioputv( numflo , 'traj_group' , REAL(ngrpfl,wp) )
  
            ELSE  ! Re-open
       
               CALL flioopfd( TRIM(clname), numflo , "WRITE" )

            ENDIF

            !II-2-b-2 Write in  netcdf file
            !-------------------------------
            irec =  INT( (kt-nn_it000+1)/nn_writefl ) +1
            ztime = ( kt-nn_it000 + 1 ) * rdt

            CALL flioputv( numflo , 'time_counter', ztime , start=(/irec/) )

            DO jfl = 1, jpnfl

               istart = (/jfl,irec/)

               CALL flioputv( numflo , 'traj_lon'    , zlon(jfl), start=istart )
               CALL flioputv( numflo , 'traj_lat'    , zlat(jfl), start=istart )  
               CALL flioputv( numflo , 'traj_depth'  , zdep(jfl), start=istart )  
               CALL flioputv( numflo , 'traj_temp'   , ztem(jfl), start=istart )  
               CALL flioputv( numflo , 'traj_salt'   , zsal(jfl), start=istart )  
               CALL flioputv( numflo , 'traj_dens'   , zrho(jfl), start=istart )  

            ENDDO

            !II-2-b-3 Close netcdf file
            !---------------------------
            CALL flioclo( numflo )

         ENDIF

#endif
      ENDIF ! netcdf writing
   
   END SUBROUTINE flo_wri

   !!=======================================================================
END MODULE flowri
