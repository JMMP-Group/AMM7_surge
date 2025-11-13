MODULE stringop
!-
!$Id: stringop.f90 2281 2010-10-15 14:21:13Z smasson $
!-
! This software is governed by the CeCILL license
! See IOIPSL/IOIPSL_License_CeCILL.txt
!---------------------------------------------------------------------
CONTAINS
!=
SUBROUTINE cmpblank (str)
!---------------------------------------------------------------------
!- Compact blanks
!---------------------------------------------------------------------
  CHARACTER(LEN=*),INTENT(inout) :: str
!-
  INTEGER :: lcc,ipb
!---------------------------------------------------------------------
  lcc = LEN_TRIM(str)
  ipb = 1
  DO
    IF (ipb >= lcc)   EXIT
    IF (str(ipb:ipb+1) == '  ') THEN
      str(ipb+1:) = str(ipb+2:lcc)
      lcc = lcc-1
    ELSE
      ipb = ipb+1
    ENDIF
  ENDDO
!----------------------
END SUBROUTINE cmpblank
!===
INTEGER FUNCTION cntpos (c_c,l_c,c_r,l_r)
!---------------------------------------------------------------------
!- Finds number of occurences of c_r in c_c
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*),INTENT(in) :: c_c
  INTEGER,INTENT(IN) :: l_c
  CHARACTER(LEN=*),INTENT(in) :: c_r
  INTEGER,INTENT(IN) :: l_r
!-
  INTEGER :: ipos,indx
!---------------------------------------------------------------------
  cntpos = 0
  ipos   = 1
  DO
    indx = INDEX(c_c(ipos:l_c),c_r(1:l_r))
    IF (indx > 0) THEN
      cntpos = cntpos+1
      ipos   = ipos+indx+l_r-1
    ELSE
      EXIT
    ENDIF
  ENDDO
!------------------
END FUNCTION cntpos
!===
INTEGER FUNCTION findpos (c_c,l_c,c_r,l_r)
!---------------------------------------------------------------------
!- Finds position of c_r in c_c
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*),INTENT(in) :: c_c
  INTEGER,INTENT(IN) :: l_c
  CHARACTER(LEN=*),INTENT(in) :: c_r
  INTEGER,INTENT(IN) :: l_r
!---------------------------------------------------------------------
  findpos = INDEX(c_c(1:l_c),c_r(1:l_r))
  IF (findpos == 0)  findpos=-1
!-------------------
END FUNCTION findpos
!===
SUBROUTINE find_str (str_tab,str,pos)
!---------------------------------------------------------------------
!- This subroutine looks for a string in a table
!---------------------------------------------------------------------
!- INPUT
!-   str_tab  : Table  of strings
!-   str      : Target we are looking for
!- OUTPUT
!-   pos      : -1 if str not found, else value in the table
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*),DIMENSION(:),INTENT(in) :: str_tab
  CHARACTER(LEN=*),INTENT(in) :: str
  INTEGER,INTENT(out) :: pos
!-
  INTEGER :: nb_str,i
!---------------------------------------------------------------------
  pos = -1
  nb_str=SIZE(str_tab)
  IF ( nb_str > 0 ) THEN
    DO i=1,nb_str
      IF ( TRIM(str_tab(i)) == TRIM(str) ) THEN
        pos = i
        EXIT
      ENDIF
    ENDDO
  ENDIF
!----------------------
END SUBROUTINE find_str
!===
SUBROUTINE nocomma (str)
!---------------------------------------------------------------------
!- Replace commas with blanks
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: str
!-
  INTEGER :: i
!---------------------------------------------------------------------
  DO i=1,LEN_TRIM(str)
    IF (str(i:i) == ',')   str(i:i) = ' '
  ENDDO
!---------------------
END SUBROUTINE nocomma
!===
SUBROUTINE strlowercase (str)
!---------------------------------------------------------------------
!- Converts a string into lowercase
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: str
!-
  INTEGER :: i,ic
!---------------------------------------------------------------------
  DO i=1,LEN_TRIM(str)
    ic = IACHAR(str(i:i))
    IF ( (ic >= 65).AND.(ic <= 90) )  str(i:i) = ACHAR(ic+32)
  ENDDO
!--------------------------
END SUBROUTINE strlowercase
!===
SUBROUTINE struppercase (str)
!---------------------------------------------------------------------
!- Converts a string into uppercase
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: str
!-
  INTEGER :: i,ic
!---------------------------------------------------------------------
  DO i=1,LEN_TRIM(str)
    ic = IACHAR(str(i:i))
    IF ( (ic >= 97).AND.(ic <= 122) )  str(i:i) = ACHAR(ic-32)
  ENDDO
!--------------------------
END SUBROUTINE struppercase
!===
SUBROUTINE str_xfw (c_string,c_word,l_ok)
!---------------------------------------------------------------------
!- Given a character string "c_string", of arbitrary length,
!- returns a logical flag "l_ok" if a word is found in it,
!- the first word "c_word" if found and the new string "c_string"
!- without the first word "c_word"
!---------------------------------------------------------------------
  CHARACTER(LEN=*),INTENT(INOUT) :: c_string
  CHARACTER(LEN=*),INTENT(OUT) :: c_word
  LOGICAL,INTENT(OUT) :: l_ok
!-
  INTEGER :: i_b,i_e
!---------------------------------------------------------------------
  l_ok = (LEN_TRIM(c_string) > 0)
  IF (l_ok) THEN
    i_b = VERIFY(c_string,' ')
    i_e = INDEX(c_string(i_b:),' ')
    IF (i_e == 0) THEN
      c_word = c_string(i_b:)
      c_string = ""
    ELSE
      c_word = c_string(i_b:i_b+i_e-2)
      c_string = ADJUSTL(c_string(i_b+i_e-1:))
    ENDIF
  ENDIF
!---------------------
END SUBROUTINE str_xfw
!===
!------------------
END MODULE stringop
