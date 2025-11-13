PROGRAM nemo
   !!======================================================================
   !!                     ***  PROGRAM nemo  ***
   !!
   !! ** Purpose :   encapsulate nemo_gcm so that it can also be called
   !!              together with the linear tangent and adjoint models
   !!======================================================================
   !! History :   OPA  ! 2001-02  (M. Imbard, A. Weaver)  Original code
   !!   NEMO      1.0  ! 2003-10  (G. Madec) F90
   !!----------------------------------------------------------------------
   USE nemogcm   ! NEMO system   (nemo_gcm routine)
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
   !
   CALL nemo_gcm           ! NEMO direct code
   ! 
   !!======================================================================
END PROGRAM nemo
