MODULE errioipsl
!-
!$Id: errioipsl.f90 2281 2010-10-15 14:21:13Z smasson $
!-
! This software is governed by the CeCILL license
! See IOIPSL/IOIPSL_License_CeCILL.txt
!---------------------------------------------------------------------
IMPLICIT NONE
!-
PRIVATE
!-
PUBLIC :: ipslnlf, ipslerr, ipslerr_act, ipslerr_inq, histerr, ipsldbg
!-
  INTEGER :: n_l=6, ilv_cur=0, ilv_max=0
  LOGICAL :: ioipsl_debug=.FALSE., lact_mode=.TRUE.
!-
!===
CONTAINS
!===
SUBROUTINE ipslnlf (new_number,old_number)
!!--------------------------------------------------------------------
!! The "ipslnlf" routine allows to know and modify
!! the current logical number for the messages.
!!
!! SUBROUTINE ipslnlf (new_number,old_number)
!!
!! Optional INPUT argument
!!
!! (I) new_number : new logical number of the file
!!
!! Optional OUTPUT argument
!!
!! (I) old_number : current logical number of the file
!!--------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,OPTIONAL,INTENT(IN)  :: new_number
  INTEGER,OPTIONAL,INTENT(OUT) :: old_number
!---------------------------------------------------------------------
  IF (PRESENT(old_number)) THEN
    old_number = n_l
  ENDIF
  IF (PRESENT(new_number)) THEN
    n_l = new_number
  ENDIF
!---------------------
END SUBROUTINE ipslnlf
!===
SUBROUTINE ipslerr (plev,pcname,pstr1,pstr2,pstr3)
!---------------------------------------------------------------------
!! The "ipslerr" routine
!! allows to handle the messages to the user.
!!
!! INPUT
!!
!! plev   : Category of message to be reported to the user
!!          1 = Note to the user
!!          2 = Warning to the user
!!          3 = Fatal error
!! pcname : Name of subroutine which has called ipslerr
!! pstr1   
!! pstr2  : Strings containing the explanations to the user
!! pstr3
!---------------------------------------------------------------------
   IMPLICIT NONE
!-
   INTEGER :: plev
   CHARACTER(LEN=*) :: pcname,pstr1,pstr2,pstr3
!-
   CHARACTER(LEN=30),DIMENSION(3) :: pemsg = &
  &  (/ "NOTE TO THE USER FROM ROUTINE ", &
  &     "WARNING FROM ROUTINE          ", &
  &     "FATAL ERROR FROM ROUTINE      " /)
!---------------------------------------------------------------------
   IF ( (plev >= 1).AND.(plev <= 3) ) THEN
     ilv_cur = plev
     ilv_max = MAX(ilv_max,plev)
     WRITE(n_l,'(/,A," ",A)') TRIM(pemsg(plev)),TRIM(pcname)
     WRITE(n_l,'(3(" --> ",A,/))') TRIM(pstr1),TRIM(pstr2),TRIM(pstr3)
   ENDIF
   IF ( (plev == 3).AND.lact_mode) THEN
     WRITE(n_l,'("Fatal error from IOIPSL. STOP in ipslerr with code")')
     STOP 1
   ENDIF
!---------------------
END SUBROUTINE ipslerr
!===
SUBROUTINE ipslerr_act (new_mode,old_mode)
!!--------------------------------------------------------------------
!! The "ipslerr_act" routine allows to know and modify
!! the current "action mode" for the error messages,
!! and reinitialize the error level values.
!!
!! SUBROUTINE ipslerr_act (new_mode,old_mode)
!!
!! Optional INPUT argument
!!
!! (I) new_mode : new error action mode
!!                .TRUE.  -> STOP     in case of fatal error
!!                .FALSE. -> CONTINUE in case of fatal error
!!
!! Optional OUTPUT argument
!!
!! (I) old_mode : current error action mode
!!--------------------------------------------------------------------
  IMPLICIT NONE
!-
  LOGICAL,OPTIONAL,INTENT(IN)  :: new_mode
  LOGICAL,OPTIONAL,INTENT(OUT) :: old_mode
!---------------------------------------------------------------------
  IF (PRESENT(old_mode)) THEN
    old_mode = lact_mode
  ENDIF
  IF (PRESENT(new_mode)) THEN
    lact_mode = new_mode
  ENDIF
  ilv_cur = 0
  ilv_max = 0
!-------------------------
END SUBROUTINE ipslerr_act
!===
SUBROUTINE ipslerr_inq (current_level,maximum_level)
!!--------------------------------------------------------------------
!! The "ipslerr_inq" routine allows to know
!! the current level of the error messages
!! and the maximum level encountered since the
!! last call to "ipslerr_act".
!!
!! SUBROUTINE ipslerr_inq (current_level,maximum_level)
!!
!! Optional OUTPUT argument
!!
!! (I) current_level : current error level
!! (I) maximum_level : maximum error level
!!--------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,OPTIONAL,INTENT(OUT) :: current_level,maximum_level
!---------------------------------------------------------------------
  IF (PRESENT(current_level)) THEN
    current_level = ilv_cur
  ENDIF
  IF (PRESENT(maximum_level)) THEN
    maximum_level = ilv_max
  ENDIF
!-------------------------
END SUBROUTINE ipslerr_inq
!===
SUBROUTINE histerr (plev,pcname,pstr1,pstr2,pstr3)
!---------------------------------------------------------------------
!- INPUT
!- plev   : Category of message to be reported to the user
!-          1 = Note to the user
!-          2 = Warning to the user
!-          3 = Fatal error
!- pcname : Name of subroutine which has called histerr
!- pstr1   
!- pstr2  : String containing the explanations to the user
!- pstr3
!---------------------------------------------------------------------
   IMPLICIT NONE
!-
   INTEGER :: plev
   CHARACTER(LEN=*) :: pcname,pstr1,pstr2,pstr3
!-
   CHARACTER(LEN=30),DIMENSION(3) :: pemsg = &
  &  (/ "NOTE TO THE USER FROM ROUTINE ", &
  &     "WARNING FROM ROUTINE          ", &
  &     "FATAL ERROR FROM ROUTINE      " /)
!---------------------------------------------------------------------
   IF ( (plev >= 1).AND.(plev <= 3) ) THEN
     WRITE(*,'("     ")')
     WRITE(*,'(A," ",A)') TRIM(pemsg(plev)),TRIM(pcname)
     WRITE(*,'(" --> ",A)') pstr1
     WRITE(*,'(" --> ",A)') pstr2
     WRITE(*,'(" --> ",A)') pstr3
   ENDIF
   IF (plev == 3) THEN
     STOP 'Fatal error from IOIPSL. See stdout for more details'
   ENDIF
!---------------------
END SUBROUTINE histerr
!===
SUBROUTINE ipsldbg (new_status,old_status)
!!--------------------------------------------------------------------
!! The "ipsldbg" routine
!! allows to activate or deactivate the debug,
!! and to know the current status of the debug.
!!
!! SUBROUTINE ipsldbg (new_status,old_status)
!!
!! Optional INPUT argument
!!
!! (L) new_status : new status of the debug
!!
!! Optional OUTPUT argument
!!
!! (L) old_status : current status of the debug
!!--------------------------------------------------------------------
  IMPLICIT NONE
!-
  LOGICAL,OPTIONAL,INTENT(IN)  :: new_status
  LOGICAL,OPTIONAL,INTENT(OUT) :: old_status
!---------------------------------------------------------------------
  IF (PRESENT(old_status)) THEN
    old_status = ioipsl_debug
  ENDIF
  IF (PRESENT(new_status)) THEN
    ioipsl_debug = new_status
  ENDIF
!---------------------
END SUBROUTINE ipsldbg
!===
!-------------------
END MODULE errioipsl
