MODULE obs_sort
   !!=====================================================================
   !!                       ***  MODULE obs_sort  ***
   !! Observation diagnostics: Various tools for sorting etc.
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   sort_dp_indx : Get indicies for ascending order for a double prec. array
   !!   index_sort   : Get indicies for ascending order for a double prec. array
   !!---------------------------------------------------------------------
   !! * Modules used
   USE par_kind, ONLY : & ! Precision variables
      & dp
  
   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE index_sort    ! Get indicies for ascending order for a double prec. array
   
   PUBLIC sort_dp_indx   ! Get indicies for ascending order for a double prec. array
  
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE sort_dp_indx( kvals, pvals, kindx )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE sort_dp_indx  ***
      !!          
      !! ** Purpose : Get indicies for ascending order for a double precision array
      !!
      !! ** Method  : Call index_sort routine
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  06-05  (K. Mogensen)  Original code
      !!        !  06-10  (A. Weaver) Cleaning
      !!----------------------------------------------------------------------

      !! * Arguments
      INTEGER, INTENT(IN) :: kvals     ! Number of elements to be sorted
      REAL(KIND=dp), DIMENSION(kvals), INTENT(IN) :: &
         & pvals            ! Array to be sorted
      INTEGER, DIMENSION(kvals), INTENT(OUT) ::  &
         & kindx            ! Indices for ordering of array

      !! * Local declarations

      !-----------------------------------------------------------------------
      ! Call qsort routine
      !-----------------------------------------------------------------------
      IF (kvals>=1) THEN

         CALL index_sort( pvals, kindx, kvals )

      ENDIF

   END SUBROUTINE sort_dp_indx

   SUBROUTINE index_sort( pval, kindx, kvals )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE index_sort  ***
      !!          
      !! ** Purpose : Get indicies for ascending order for a double precision array
      !!
      !! ** Method  : Heapsort
      !!
      !! ** Action  : 
      !!
      !! References : http://en.wikipedia.org/wiki/Heapsort
      !!
      !! History :
      !!        !  06-05  (K. Mogensen)  Original code
      !!        !  06-10  (A. Weaver) Cleaning
      !!----------------------------------------------------------------------

      !! * Arguments
      INTEGER, INTENT(IN) :: kvals         ! Number of values
      REAL(KIND=dp), DIMENSION(kvals), INTENT(IN) :: &
         & pval                            ! Array to be sorted
      INTEGER, DIMENSION(kvals), INTENT(INOUT) :: &
         & kindx                           ! Indicies for ordering

      !! * Local declarations
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jt
      INTEGER :: jn
      INTEGER :: jparent
      INTEGER :: jchild

      DO ji = 1, kvals
         kindx(ji) = ji
      END DO
      
      ji = kvals/2 + 1
      jn = kvals

      main_loop : DO

         IF ( ji > 1 ) THEN
            ji = ji-1
            jt = kindx(ji)
         ELSE
            jt = kindx(jn)
            kindx(jn) = kindx(1)
            jn = jn-1
            IF ( jn <= 1 ) THEN
               kindx(1) = jt
               EXIT main_loop
            ENDIF
         ENDIF

         jparent = ji
         jchild =  2 * ji

         inner_loop : DO

            IF ( jchild > jn ) EXIT inner_loop
            IF ( jchild < jn ) THEN
               IF ( pval(kindx(jchild)) < pval(kindx(jchild+1)) ) THEN
                 jchild = jchild+1
               ENDIF
            ENDIF
            IF  ( pval(jt) < pval(kindx(jchild))) THEN
               kindx(jparent) = kindx(jchild)
               jparent = jchild
               jchild  = jchild*2
            ELSE 
               jchild = jn + 1 
            ENDIF

         END DO inner_loop

         kindx(jparent) = jt

      END DO main_loop
      
   END SUBROUTINE index_sort

END MODULE obs_sort
 
