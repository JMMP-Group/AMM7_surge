MODULE iscplhsb
   !!======================================================================
   !!                       ***  MODULE  iscplhsb  ***
   !! Ocean forcing: ice sheet/ocean coupling (conservation)
   !!=====================================================================
   !! History :  NEMO  ! 2015-01 P. Mathiot: original 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   iscpl_alloc    : variable allocation
   !!   iscpl_hsb      : compute and store the input of heat/salt/volume 
   !!                    into the system due to the coupling process
   !!   iscpl_div      : correction of divergence to keep volume conservation
   !!----------------------------------------------------------------------
   USE oce             ! global tra/dyn variable
   USE dom_oce         ! ocean space and time domain
   USE domwri          ! ocean space and time domain
   USE domngb          ! 
   USE phycst          ! physical constants
   USE sbc_oce         ! surface boundary condition variables
   USE iscplini        ! 
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE lib_fortran     ! MPP library
   USE lbclnk          !

   IMPLICIT NONE
   PRIVATE
   
   PUBLIC   iscpl_div   
   PUBLIC   iscpl_cons       
   !! * Substitutions  
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE iscpl_cons(ptmask_b, psmask_b, pe3t_b, pts_flx, pvol_flx, prdt_iscpl)
      !!---------------------------------------------------------------------- 
      !!                   ***  ROUTINE iscpl_cons  ***
      !! 
      !! ** Purpose :   compute input into the system during the coupling step
      !!                compute the correction term
      !!                compute where the correction have to be applied
      !! 
      !! ** Method  :   compute tsn*e3t-tsb*e3tb and e3t-e3t_b
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:  ), INTENT(in ) :: ptmask_b    !! mask before
      REAL(wp), DIMENSION(:,:,:  ), INTENT(in ) :: pe3t_b      !! scale factor before
      REAL(wp), DIMENSION(:,:    ), INTENT(in ) :: psmask_b    !! mask before
      REAL(wp), DIMENSION(:,:,:,:), INTENT(out) :: pts_flx     !! corrective flux to have tracer conservation
      REAL(wp), DIMENSION(:,:,:  ), INTENT(out) :: pvol_flx    !! corrective flux to have volume conservation
      REAL(wp),                     INTENT(in ) :: prdt_iscpl  !! coupling period 
      !
      INTEGER  ::   ji  , jj  , jk           ! loop index
      INTEGER  ::   jip1, jim1, jjp1, jjm1
      REAL(wp) ::   summsk, zsum , zsumn, zjip1_ratio  , zjim1_ratio, zdtem, zde3t, z1_rdtiscpl
      REAL(wp) ::   zarea , zsum1, zsumb, zjjp1_ratio  , zjjm1_ratio, zdsal
      REAL(wp), DIMENSION(jpi,jpj)        ::   zdssh   ! workspace
      REAL(wp), DIMENSION(:), ALLOCATABLE ::   zlon, zlat
      REAL(wp), DIMENSION(:), ALLOCATABLE ::   zcorr_vol, zcorr_tem, zcorr_sal
      INTEGER , DIMENSION(:), ALLOCATABLE ::   ixpts, iypts, izpts, inpts
      INTEGER :: jpts, npts
      !!----------------------------------------------------------------------

      ! get imbalance (volume heat and salt)
      ! initialisation difference
      zde3t = 0._wp   ;   zdsal = 0._wp   ;   zdtem = 0._wp

      ! initialisation correction term
      pvol_flx(:,:,:  ) = 0._wp
      pts_flx (:,:,:,:) = 0._wp
      
      z1_rdtiscpl = 1._wp / prdt_iscpl 

      ! mask tsn and tsb 
      tsb(:,:,:,jp_tem) = tsb(:,:,:,jp_tem) * ptmask_b(:,:,:)
      tsn(:,:,:,jp_tem) = tsn(:,:,:,jp_tem) *  tmask  (:,:,:)
      tsb(:,:,:,jp_sal) = tsb(:,:,:,jp_sal) * ptmask_b(:,:,:)
      tsn(:,:,:,jp_sal) = tsn(:,:,:,jp_sal) *  tmask  (:,:,:)

      !==============================================================================
      ! diagnose the heat, salt and volume input and compute the correction variable
      !==============================================================================

      ! 
      zdssh(:,:) = sshn(:,:) * ssmask(:,:) - sshb(:,:) * psmask_b(:,:)
      IF (.NOT. ln_linssh ) zdssh = 0.0_wp ! already included in the levels by definition
      
      DO jk = 1,jpk-1
         DO jj = 2,jpj-1
            DO ji = fs_2,fs_jpim1
               IF (tmask_h(ji,jj) == 1._wp) THEN

                  ! volume differences
                  zde3t = e3t_n(ji,jj,jk) * tmask(ji,jj,jk) - pe3t_b(ji,jj,jk) * ptmask_b(ji,jj,jk)

                  ! heat diff
                  zdtem = tsn(ji,jj,jk,jp_tem) * e3t_n(ji,jj,jk) *  tmask  (ji,jj,jk)   &
                        - tsb(ji,jj,jk,jp_tem) * pe3t_b (ji,jj,jk) * ptmask_b(ji,jj,jk)
                  ! salt diff
                  zdsal = tsn(ji,jj,jk,jp_sal) * e3t_n(ji,jj,jk) *  tmask  (ji,jj,jk)   &
                        - tsb(ji,jj,jk,jp_sal) * pe3t_b (ji,jj,jk) * ptmask_b(ji,jj,jk)
               
                  ! shh changes
                  IF ( ptmask_b(ji,jj,jk) == 1._wp .OR. tmask(ji,jj,jk) == 1._wp ) THEN 
                     zde3t = zde3t + zdssh(ji,jj) ! zdssh = 0 if vvl
                     zdssh(ji,jj) = 0._wp
                  END IF

                  ! volume, heat and salt differences in each cell 
                  pvol_flx(ji,jj,jk)       =   pvol_flx(ji,jj,jk)        + zde3t * z1_rdtiscpl
                  pts_flx (ji,jj,jk,jp_sal)=   pts_flx (ji,jj,jk,jp_sal) + zdsal * z1_rdtiscpl 
                  pts_flx (ji,jj,jk,jp_tem)=   pts_flx (ji,jj,jk,jp_tem) + zdtem * z1_rdtiscpl

                  ! case where we close a cell: check if the neighbour cells are wet 
                  IF ( tmask(ji,jj,jk) == 0._wp .AND. ptmask_b(ji,jj,jk) == 1._wp ) THEN

                     jip1=ji+1 ; jim1=ji-1 ; jjp1=jj+1 ; jjm1=jj-1 ;

                     zsum =   e1e2t(ji  ,jjp1) * tmask(ji  ,jjp1,jk) + e1e2t(ji  ,jjm1) * tmask(ji  ,jjm1,jk) &
                       &    + e1e2t(jim1,jj  ) * tmask(jim1,jj  ,jk) + e1e2t(jip1,jj  ) * tmask(jip1,jj  ,jk)

                     IF ( zsum /= 0._wp ) THEN
                        zjip1_ratio   = e1e2t(jip1,jj  ) * tmask(jip1,jj  ,jk) / zsum
                        zjim1_ratio   = e1e2t(jim1,jj  ) * tmask(jim1,jj  ,jk) / zsum
                        zjjp1_ratio   = e1e2t(ji  ,jjp1) * tmask(ji  ,jjp1,jk) / zsum
                        zjjm1_ratio   = e1e2t(ji  ,jjm1) * tmask(ji  ,jjm1,jk) / zsum

                        pvol_flx(ji  ,jjp1,jk       ) = pvol_flx(ji  ,jjp1,jk       ) + pvol_flx(ji,jj,jk       ) * zjjp1_ratio
                        pvol_flx(ji  ,jjm1,jk       ) = pvol_flx(ji  ,jjm1,jk       ) + pvol_flx(ji,jj,jk       ) * zjjm1_ratio
                        pvol_flx(jip1,jj  ,jk       ) = pvol_flx(jip1,jj  ,jk       ) + pvol_flx(ji,jj,jk       ) * zjip1_ratio
                        pvol_flx(jim1,jj  ,jk       ) = pvol_flx(jim1,jj  ,jk       ) + pvol_flx(ji,jj,jk       ) * zjim1_ratio
                        pts_flx (ji  ,jjp1,jk,jp_sal) = pts_flx (ji  ,jjp1,jk,jp_sal) + pts_flx (ji,jj,jk,jp_sal) * zjjp1_ratio
                        pts_flx (ji  ,jjm1,jk,jp_sal) = pts_flx (ji  ,jjm1,jk,jp_sal) + pts_flx (ji,jj,jk,jp_sal) * zjjm1_ratio
                        pts_flx (jip1,jj  ,jk,jp_sal) = pts_flx (jip1,jj  ,jk,jp_sal) + pts_flx (ji,jj,jk,jp_sal) * zjip1_ratio
                        pts_flx (jim1,jj  ,jk,jp_sal) = pts_flx (jim1,jj  ,jk,jp_sal) + pts_flx (ji,jj,jk,jp_sal) * zjim1_ratio
                        pts_flx (ji  ,jjp1,jk,jp_tem) = pts_flx (ji  ,jjp1,jk,jp_tem) + pts_flx (ji,jj,jk,jp_tem) * zjjp1_ratio
                        pts_flx (ji  ,jjm1,jk,jp_tem) = pts_flx (ji  ,jjm1,jk,jp_tem) + pts_flx (ji,jj,jk,jp_tem) * zjjm1_ratio
                        pts_flx (jip1,jj  ,jk,jp_tem) = pts_flx (jip1,jj  ,jk,jp_tem) + pts_flx (ji,jj,jk,jp_tem) * zjip1_ratio
                        pts_flx (jim1,jj  ,jk,jp_tem) = pts_flx (jim1,jj  ,jk,jp_tem) + pts_flx (ji,jj,jk,jp_tem) * zjim1_ratio

                        ! set to 0 the cell we distributed over neigbourg cells
                        pvol_flx(ji,jj,jk       ) = 0._wp
                        pts_flx (ji,jj,jk,jp_sal) = 0._wp
                        pts_flx (ji,jj,jk,jp_tem) = 0._wp

                     ELSE IF (zsum == 0._wp ) THEN
                        ! case where we close a cell and no adjacent cell open
                        ! check if the cell beneath is wet
                        IF ( tmask(ji,jj,jk+1) == 1._wp ) THEN
                           pvol_flx(ji,jj,jk+1)       =  pvol_flx(ji,jj,jk+1)        + pvol_flx(ji,jj,jk)
                           pts_flx (ji,jj,jk+1,jp_sal)=  pts_flx (ji,jj,jk+1,jp_sal) + pts_flx (ji,jj,jk,jp_sal)
                           pts_flx (ji,jj,jk+1,jp_tem)=  pts_flx (ji,jj,jk+1,jp_tem) + pts_flx (ji,jj,jk,jp_tem)

                           ! set to 0 the cell we distributed over neigbourg cells
                           pvol_flx(ji,jj,jk       ) = 0._wp
                           pts_flx (ji,jj,jk,jp_sal) = 0._wp
                           pts_flx (ji,jj,jk,jp_tem) = 0._wp
                        ELSE
                        ! case no adjacent cell on the horizontal and on the vertical
                           IF ( lwp ) THEN   ! JMM : cAution this warning may occur on any mpp subdomain but numout is only
                                             ! open for narea== 1 (lwp=T)
                           WRITE(numout,*) 'W A R N I N G iscpl: no adjacent cell on the vertical and horizontal'
                           WRITE(numout,*) '                     ',mig(ji),' ',mjg(jj),' ',jk
                           WRITE(numout,*) '                     ',ji,' ',jj,' ',jk,' ',narea
                           WRITE(numout,*) ' we are now looking for the closest wet cell on the horizontal '
                           ENDIF
                        ! We deal with these points later.
                        END IF
                     END IF
                  END IF
               END IF
            END DO
         END DO
      END DO

!!gm  ERROR !!!!
!!    juste use tmask_i  or in case of ISF smask_i (to be created to compute the sum without halos)
!
!      CALL lbc_sum(pvol_flx(:,:,:       ),'T',1.)
!      CALL lbc_sum(pts_flx (:,:,:,jp_sal),'T',1.)
!      CALL lbc_sum(pts_flx (:,:,:,jp_tem),'T',1.)
      CALL ctl_stop( 'STOP', ' iscpl_cons:   please modify this MODULE !' )
!!gm end
      ! if no neighbour wet cell in case of 2close a cell", need to find the nearest wet point 
      ! allocation and initialisation of the list of problematic point
      ALLOCATE( inpts(jpnij) )
      inpts(:) = 0

      ! fill narea location with the number of problematic point
      DO jk = 1,jpk-1
         DO jj = 2,jpj-1
            DO ji = fs_2,fs_jpim1
               IF (     ptmask_b(ji,jj,jk) == 1._wp .AND. tmask(ji,jj,jk+1)  == 0._wp .AND. tmask_h(ji,jj) == 1._wp  &
                  .AND. SUM(tmask(ji-1:ji+1,jj,jk)) + SUM(tmask(ji,jj-1:jj+1,jk)) == 0._wp) THEN
                  inpts(narea) = inpts(narea) + 1 
               END IF
            END DO
         END DO
      END DO

      ! build array of total problematic point on each cpu (share to each cpu)
      CALL mpp_max('iscplhsb', inpts,jpnij) 

      ! size of the new variable
      npts  = SUM(inpts)    
      
      ! allocation of the coordinates, correction, index vector for the problematic points 
      ALLOCATE(ixpts(npts), iypts(npts), izpts(npts), zcorr_vol(npts), zcorr_sal(npts), zcorr_tem(npts), zlon(npts), zlat(npts))
      ixpts(:) = -9999 ; iypts(:) = -9999 ; izpts(:) = -9999 ; zlon(:) = -1.0e20_wp ; zlat(:) = -1.0e20_wp
      zcorr_vol(:) = -1.0e20_wp
      zcorr_sal(:) = -1.0e20_wp
      zcorr_tem(:) = -1.0e20_wp

      ! fill new variable
      jpts = SUM(inpts(1:narea-1))
      DO jk = 1,jpk-1
         DO jj = 2,jpj-1
            DO ji = fs_2,fs_jpim1
               IF (     ptmask_b(ji,jj,jk) == 1._wp .AND. tmask(ji,jj,jk+1)  == 0._wp .AND. tmask_h(ji,jj) == 1._wp  &
                  .AND. SUM(tmask(ji-1:ji+1,jj,jk)) + SUM(tmask(ji,jj-1:jj+1,jk)) == 0._wp) THEN
                  jpts = jpts + 1  ! positioning in the inpts vector for the area narea
                  ixpts(jpts) = ji           ; iypts(jpts) = jj ; izpts(jpts) = jk
                  zlon (jpts) = glamt(ji,jj) ; zlat (jpts) = gphit(ji,jj)
                  zcorr_vol(jpts) = pvol_flx(ji,jj,jk)
                  zcorr_sal(jpts) = pts_flx (ji,jj,jk,jp_sal)
                  zcorr_tem(jpts) = pts_flx (ji,jj,jk,jp_tem)

                  ! set flx to 0 (safer)
                  pvol_flx(ji,jj,jk       ) = 0.0_wp
                  pts_flx (ji,jj,jk,jp_sal) = 0.0_wp
                  pts_flx (ji,jj,jk,jp_tem) = 0.0_wp
               END IF
            END DO
         END DO
      END DO

      ! build array of total problematic point on each cpu (share to each cpu)
      ! point coordinates
      CALL mpp_max('iscplhsb', zlat ,npts)
      CALL mpp_max('iscplhsb', zlon ,npts)
      CALL mpp_max('iscplhsb', izpts,npts)

      ! correction values 
      CALL mpp_max('iscplhsb', zcorr_vol,npts)
      CALL mpp_max('iscplhsb', zcorr_sal,npts)
      CALL mpp_max('iscplhsb', zcorr_tem,npts)

      ! put correction term in the closest cell          
      DO jpts = 1,npts
         CALL dom_ngb(zlon(jpts), zlat(jpts), ixpts(jpts), iypts(jpts),'T', izpts(jpts))
         DO jj = mj0(iypts(jpts)),mj1(iypts(jpts))
            DO ji = mi0(ixpts(jpts)),mi1(ixpts(jpts))
               jk = izpts(jpts)

               IF (tmask_h(ji,jj) == 1._wp) THEN
                  ! correct the vol_flx in the closest cell
                  pvol_flx(ji,jj,jk)        =  pvol_flx(ji,jj,jk       ) + zcorr_vol(jpts)
                  pts_flx (ji,jj,jk,jp_sal) =  pts_flx (ji,jj,jk,jp_sal) + zcorr_sal(jpts)
                  pts_flx (ji,jj,jk,jp_tem) =  pts_flx (ji,jj,jk,jp_tem) + zcorr_tem(jpts)

                  ! set correction to 0
                  zcorr_vol(jpts) = 0.0_wp
                  zcorr_sal(jpts) = 0.0_wp
                  zcorr_tem(jpts) = 0.0_wp
               END IF
            END DO
         END DO
      END DO

      ! deallocate variables 
      DEALLOCATE(inpts)
      DEALLOCATE(ixpts, iypts, izpts, zcorr_vol, zcorr_sal, zcorr_tem, zlon, zlat)
    
      ! add contribution store on the hallo (lbclnk remove one of the contribution)
      pvol_flx(:,:,:       ) = pvol_flx(:,:,:       ) * tmask(:,:,:)
      pts_flx (:,:,:,jp_sal) = pts_flx (:,:,:,jp_sal) * tmask(:,:,:)
      pts_flx (:,:,:,jp_tem) = pts_flx (:,:,:,jp_tem) * tmask(:,:,:)

!!gm  ERROR !!!!
!!    juste use tmask_i  or in case of ISF smask_i (to be created to compute the sum without halos)
!
!      ! compute sum over the halo and set it to 0.
!      CALL lbc_sum(pvol_flx(:,:,:       ),'T',1._wp)
!      CALL lbc_sum(pts_flx (:,:,:,jp_sal),'T',1._wp)
!      CALL lbc_sum(pts_flx (:,:,:,jp_tem),'T',1._wp)
!!gm end

      !
   END SUBROUTINE iscpl_cons


   SUBROUTINE iscpl_div( phdivn )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE iscpl_div  ***
      !!
      !! ** Purpose :   update the horizontal divergenc
      !!
      !! ** Method  :
      !!                CAUTION : iscpl is positive (inflow) and expressed in m/s
      !!
      !! ** Action  :   phdivn   increase by the iscpl correction term
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   phdivn   ! horizontal divergence
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               phdivn(ji,jj,jk) = phdivn(ji,jj,jk) + hdiv_iscpl(ji,jj,jk) / e3t_n(ji,jj,jk)
            END DO
         END DO
      END DO
      !
   END SUBROUTINE iscpl_div

END MODULE iscplhsb
