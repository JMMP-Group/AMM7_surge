MODULE declare
! -*- Mode: f90 -*-
!$Id: ncunderflow.f90 2281 2010-10-15 14:21:13Z smasson $
!-
! This software is governed by the CeCILL license
! See IOIPSL/IOIPSL_License_CeCILL.txt
!-
! f90 -L/usr/local/lib -lnetcdf -align dcommons -g
!     -ladebug -check format -check bounds
!     -check output_conversion -fpe1
!     -I/usr/local/include -free -arch host -tune host
!     -warn declarations -warn argument_checking
!     ncunderflow.f  -o ncunderflow
!
! ifc -FR -cl,ncunderflow.pcl -o ncunderflow ncunderflow.f
!     -L/usr/local/install/netcdf/lib/libnetcdf.a -lPEPCF90
!
  IMPLICIT NONE
  INTEGER, PARAMETER :: r4 = 4, r8 = 8, i4 = 4, i8 = 8
  INTEGER, PARAMETER :: il = KIND(1)
  LOGICAL :: ldebug = .FALSE.
  INTEGER (kind = il) :: nout = 0, nerr = 0         ! Standard output, standard error
  CHARACTER (LEN=4), PARAMETER :: cerror = 'VOID'
END MODULE declare
!!
MODULE mod_nfdiag
CONTAINS
  SUBROUTINE nfdiag ( kios, clmess, lcd)
    !!
    !! Imprime un message d'erreur NetCDF
    !!
    USE declare
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
    !!
    INTEGER (kind=i4), INTENT (in) :: kios
    CHARACTER (len = *), INTENT (in) :: clmess
    LOGICAL, INTENT (in), OPTIONAL :: lcd 
    CHARACTER (len = 80) :: clt
    LOGICAL :: ld
    !!
    IF ( PRESENT ( lcd)) THEN
        ld = lcd
    ELSE
        ld = ldebug
    ENDIF
    !!
    clt = TRIM ( NF_STRERROR ( kios) )
    !!
    IF ( ld ) THEN 
        IF ( kios == NF_NOERR ) THEN 
            WRITE ( unit = nout, fmt = * ) "OK : ", TRIM (clmess)
        ELSE
            WRITE ( unit = nout, fmt = * ) "ERROR : ", TRIM (clmess), " : ", TRIM ( clt), " : ", kios
            IF ( .NOT. ld ) STOP
        END IF
    ELSE
        IF ( kios /= NF_NOERR ) THEN 
            WRITE ( unit = nout, fmt = * ) "ERROR : ", TRIM (clmess), " : ", TRIM ( clt), " : ", kios
            STOP
        END IF
    ENDIF
    !!
    RETURN
    !!
  END SUBROUTINE nfdiag
  !!
END MODULE mod_nfdiag

MODULE mod_lec
CONTAINS
  !!
  SUBROUTINE lec (chaine, cval, c_c)
    !!
    USE declare
    IMPLICIT NONE
    !!
    CHARACTER (len = *), INTENT ( inout) :: chaine
    CHARACTER (len = *), INTENT ( inout) :: cval
    CHARACTER (len=*), OPTIONAL :: c_c
    INTEGER (kind = il) :: ji, ji1, ji2, ji3, jl, jb
    INTEGER (kind = i4) :: index
    !!
    !! Read character string up to ':' or ',', or in c_c if present
    !! Returns the real before the character (xerror if not available)
    !! Reduce the string
    !!
    jl = LEN (chaine) ; jb = LEN_TRIM (chaine)
    IF ( ldebug) WRITE ( nout, *) 'Lec : jl, jb ', jl, jb
    IF ( jb == 0 ) THEN 
        cval = cerror
    ELSE
        ji1 = INDEX (chaine, ':') ; ji2 = INDEX (chaine, ',')
        IF ( PRESENT (c_c)) THEN 
            ji3 = INDEX (chaine, c_c) ; ji = MAX (ji1, ji2, ji3)
        ELSE
            ji = MAX (ji1, ji2)
        ENDIF
        IF ( ji == 0 ) THEN 
            READ ( chaine (1:jb) , fmt = * ) cval
            chaine (1:jl-jb) = chaine (jb+1:jl)
        ELSE IF ( ji == 1 ) THEN 
            cval = cerror
            chaine (1:jl-1) = chaine (2:jl)
        ELSE
            cval = chaine (1:ji-1)
            chaine (1:jl-ji) = chaine (ji+1:jl )
        END IF
    END IF
    !!
  END SUBROUTINE lec
END MODULE mod_lec

PROGRAM ncunderflow

  ! Ce programme ouvre un fichier de donnees au format netcdf
  ! et met a zero toutes les valeurs trop petites pour etre
  ! representees par un reel sur 4 octets au format IEEE
  !
  ! Revision 2.0  2004/04/05 14:47:50  adm
  ! JB+MAF+AC: switch to IOIPSL 2.0 (1)
  !
  ! Revision 1.1  2003/04/09 15:21:56  adm
  ! add ncunderflow in IOIPSL
  ! and modify AA_make to take it into account
  ! SD + MAF
  !
  ! Revision 1.1  2001/02/07  14:36:07  jypeter
  ! J-Y Peterschmitt / LMCE / 07/02/2001
  ! Initial revision
  !
  USE declare
  USE mod_nfdiag
  USE mod_lec
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  INTEGER (kind=il), EXTERNAL :: iargc

  ! Nombre maximal de dimensions : 6

  INTEGER (kind=il), PARAMETER :: jpmaxdim = 6, jpmaxvar = 1024

  CHARACTER (len = 128) :: clnomprog, clnomfic
  CHARACTER (len = 1024) :: clistvar, clecline
  CHARACTER (len = 128), DIMENSION(jpmaxdim) :: clnomdim
  CHARACTER (len = 128), DIMENSION(jpmaxvar) :: clvarcmd, clvarfic, clvar ! Nom des variables dans le fichier est sur la ligne de commande.
  LOGICAL :: lrever = .FALSE. ! Si .true., on traite toutes les variables sauf celle de la ligne de commande
  LOGICAL :: lnocoord = .FALSE. ! Si .truee., on exclu les variables coordonnées
  LOGICAL :: lverbose = .TRUE.
  
  INTEGER (kind=il) :: incid, ircode, ivarid, ivartype, inbdim, inbatt
  INTEGER (kind=il) :: nvarcmd, nvarfic, nvar, nfile, jvarcmd, jvarfic, jvar, jfile, ierr
  INTEGER (kind=il) :: ji, jdim3, jdim4, jdim5, jdim6, j1, j2, j3, jarg, ncumul
  INTEGER (kind=il), DIMENSION(jpmaxdim) :: idimid, idimsize, istart, icount
  REAL (kind=r4), DIMENSION(:,:), ALLOCATABLE :: zdatacorr
  REAL (kind=r8), DIMENSION(:,:), ALLOCATABLE :: zdata
  REAL (kind=r4) :: reps = TINY (1.0_r4) * 10.0_r4
  LOGICAL :: lok

  ! Verification du nombre de parametres
  IF(iargc() .LT. 2) THEN
      CALL usage
      STOP
  ENDIF

  ! Aide
  jarg = 1
  Lab1: DO WHILE ( jarg <= 3 )
    IF (ldebug) WRITE(nout,*) 'lecture ligne commande ', jarg
    CALL getarg (jarg,clecline)
    IF ( clecline(1:1) /= '-' ) EXIT Lab1
    IF ( clecline(1:2) == '-h' .OR. clecline(1:2) == '-?'  ) THEN
        CALL usage
        STOP 
    ELSE IF ( clecline(1:2) == '-x' ) THEN
        lrever = .TRUE. 
    ELSE IF ( clecline(1:2) == '-d' ) THEN
        ldebug = .TRUE. 
    ELSE IF ( clecline(1:2) == '-V' ) THEN
        lverbose = .FALSE. 
    ELSE IF ( clecline(1:2) == '-v' ) THEN
        jarg = jarg + 1
        ! Recuperation des noms de variables
        IF (ldebug) WRITE(nout,*) 'lecture liste vriables ', jarg
        CALL getarg (jarg,clistvar)
        clistvar = TRIM(ADJUSTL(clistvar))
        jvarcmd = 0 ; nvarcmd = 0
        SeekVar: DO WHILE ( .TRUE. )
          CALL lec ( clistvar, clvarcmd(jvarcmd+1)(:) )
          IF ( TRIM(clvarcmd(jvarcmd+1)(:)) == cerror ) EXIT SeekVar
          jvarcmd = jvarcmd + 1
          nvarcmd = jvarcmd
          IF (ldebug) WRITE(nout,*) 'affecte variable ', jvarcmd, TRIM(clvarcmd(jvarcmd))
        END DO SeekVar
    ENDIF
    jarg = jarg + 1
  END DO Lab1

  ! Boucle sur les fichiers
  FileLoop: DO jfile = jarg, iargc()
    
    ! Recuperation du nom du fichier a traiter
    CALL getarg ( jfile, clnomfic)
    
    ! Ouverture du fichier
    CALL nfdiag ( NF_OPEN ( TRIM(clnomfic), NF_WRITE, incid ), "Opening " // TRIM(clnomfic) )
    WRITE (nout,*) TRIM(clnomfic)

    ! Recuparation de la liste des variables du fichier
    nvarfic = 0
    DO jvarfic = 1, jpmaxvar
      j3 = NF_INQ_VAR ( incid, jvarfic, clvarfic(jvarfic)(:), ivartype, inbdim, idimid, inbatt)
      IF ( j3 /= NF_NOERR ) EXIT
      nvarfic = jvarfic
    END DO

    ! Liste des variables a traiter
    IF ( lrever ) THEN
        IF ( nvarcmd == 0) THEN 
            clvar = clvarfic
            nvar  = nvarfic
        ELSE 
            jvar = 0
            DO jvarfic = 1, nvarfic
              lok = .TRUE.
              DO jvarcmd = 1, nvarcmd
                IF ( TRIM(clvarfic(jvarfic)(:)) == TRIM(clvarcmd(jvarcmd)(:)) ) THEN
                    lok = .FALSE.
                END IF
              END DO
              IF ( lok) THEN 
                  jvar = jvar + 1
                  clvar(jvar) = clvarfic(jvarfic)
              END IF
            END DO
            nvar = jvar
        END IF
    ELSE
        clvar = clvarcmd
        nvar  = nvarcmd
    END IF

    ncumul = 0
    VarLoop: DO jvar = 1, nvar 
      
      IF (lverbose) &
         & WRITE(nout, FMT='("Correction de ", A, " dans ", A, " : ", $)') TRIM(clvar(jvar)(:)), TRIM(clnomfic)

      ! Passage de netcdf en mode 'erreurs non fatales'
      !  CALL ncpopt(NCVERBOS)
      ! En fait, on reste dans le mode par defaut, dans lequel une erreur
      ! netcdf cause un arret du programme. Du coup, il n'est pas
      ! necessaire de tester la valeur de la variable ircode
      ! ATTENTION! Si jamais on veut arreter le programme a cause d'une
      ! erreur ne provenant pas de netcdf, il faut penser a fermer
      ! manuellement le fichier avec un appel a ncclos
      
      ! Recuperation de l'identificateur de la variable
      CALL nfdiag ( NF_INQ_VARID ( incid, TRIM(clvar(jvar)(:)), ivarid), "Get var id " // TRIM(clvar(jvar)(:)))

      ivartype = 0 ; idimid = 0 ; inbdim = 0 ; inbatt = 0
      ! Recuperation du nombre de dimensions de la variable
      CALL nfdiag ( NF_INQ_VAR ( incid, ivarid, clvar(jvar)(:), ivartype, inbdim, idimid, inbatt),  &
         & "Get var info " // TRIM(clvar(jvar)(:)))
      
      IF(inbdim .GT. jpmaxdim) THEN
          WRITE(nout,*)
          WRITE(nout, *) 'La variable ', TRIM(clvar(jvar)(:)), ' a trop de dimensions'
          CALL nfdiag ( NF_CLOSE (incid), "Closing file")
          STOP
      ENDIF
      
      ! Recuperation des dimensions effectives
      idimsize(3:jpmaxdim) = 1 ! Au cas ou la variable n'ait que
      ! 2 ou 3 dims, on initialise ces valeurs
      ! qui serviront dans le controle des boucles
      ! et qui auraient une valeur indefinie sinon
      DO ji = 1, inbdim
        CALL nfdiag ( NF_INQ_DIM ( incid, idimid(ji), clnomdim(ji), idimsize(ji)), "NF_INQ_DIM")
        IF (lverbose) WRITE(nout, '(A,A,A,I3,$)') '    ', TRIM(clnomdim(ji)), ' = ', idimsize(ji)
        IF ( idimsize(ji) == 0 ) THEN
            WRITE(nout, '(A,A,A,A,I3)') TRIM(clvar(jvar)(:)), ', ', TRIM(clnomdim(ji)), ' = ', idimsize(ji)
            CYCLE VarLoop
        END IF
      ENDDO
      IF (lverbose) WRITE(nout,*)
      idimsize = MAX ( idimsize, 1)
      ncumul = ncumul + 1
      
      ! Determination du type de la variable, en fonction du nom de
      ! la premiere dimension
!$$$  IF(INDEX(TRIM(clnomdim(1)),'ongitude') .NE. 0) THEN
!$$$      ! var de type map ou 3d
!$$$      write(nout, *) '  --> MAP/3D'
!$$$  ELSE IF(INDEX(TRIM(clnomdim(1)),'atitude') .NE. 0) THEN
!$$$      ! var de type xsec 
!$$$      write(nout, *) '  --> XSEC'
!$$$  ELSE
!$$$      WRITE(nout, *)
!$$$      WRITE(nout, *) 'Bizarre, la premiere dimension n''est ni "longitude" ni "latitude"'
!$$$      CALL ncclos(incid, ircode)
!$$$      STOP
!$$$  ENDIF

      ! Reservation de memoire pour charger et traiter
      ! une grille idimsize(1)*idimsize(2) de la variable
      ALLOCATE(zdata(idimsize(1), idimsize(2)), stat=ierr)
      IF(ierr .NE. 0) THEN
          WRITE(nout, *) 'Erreur d''allocation memoire pour zdata'
          CALL nfdiag ( NF_CLOSE (incid), "NF_CLOSE")
          STOP
      ENDIF
      ALLOCATE(zdatacorr(idimsize(1), idimsize(2)), stat=ierr)
      IF(ierr .NE. 0) THEN
          WRITE(nout, *) 'Erreur d''allocation memoire pour zdatacorr'
          CALL nfdiag ( NF_CLOSE (incid), "NF_CLOSE")
          STOP
      ENDIF

      ! Parametrisation de la partie de la variable a charger en memoire
      ! (une 'grille' que l'on lira autant de fois qu'il y a de niveaux et
      !  de pas de temps)
      ! Rappel : seuls les elements 1..inbdim des tableaux sont
      ! significatifs et utiles

      icount = 0
      
      DO jdim6 = 1, idimsize(6)
        DO jdim5 = 1, idimsize(5)
          DO jdim4 = 1, idimsize(4)
            DO jdim3 = 1, idimsize(3)
              istart = (/   1        ,      1      , jdim3, jdim4, jdim5, jdim6 /)
              icount = (/ idimsize(1),  idimsize(2),   1  ,  1   ,  1   ,  1    /)

              ! Chargement d'une 'grille' de donnees, en real*8
              CALL nfdiag ( NF_GET_VARA_DOUBLE(incid, ivarid, istart(1:inbdim), icount(1:inbdim), zdata), &
                 & "NF_GET_VARA_DOUBLE")
              ! Mise a zero de toutes les valeurs trop petites pour etre
              ! representees par un reel sur 4 octets au format IEEE.
              ! Le truc est de faire une operation nulle (addition de 0)
              ! sur des donnees qui posent problemes, EN AYANT COMPILE LE PROG
              ! AVEC l'OPTION "-fpe1". Dans ce cas, les valeurs trop petites
              ! sont remplacees par zero (0.0) et le programme continue,
              ! au lieu de planter.
              ! Il est possible de faire afficher le nb de valeurs qui ont pose
              ! un pb en utilisant en plus l'option "-check underflow"
              zdata = zdata + 0.0_r8
              zdatacorr = REAL(zdata, KIND=r4)
              WHERE ( ABS (zdatacorr) < reps) zdatacorr = 0.0_r4
              
              ! Sauvegarde de la grille corrigee dans le fichier
              ! (a la place de la grille initiale), en real*4
              CALL nfdiag ( NF_PUT_VARA_REAL(incid, ivarid, istart, icount, zdatacorr), "NF_PUT_VARA_REAL" )
              
            END DO
          END DO
        END DO
      END DO
      
      DEALLOCATE ( zdata)
      DEALLOCATE ( zdatacorr)
      
    END DO VarLoop
    
    WRITE (nout,*) 'ncunderflow, nombre de variables corrigees : ', ncumul

    ! Fermeture du fichier
    CALL nfdiag ( NF_CLOSE (incid), "Closing" )
    
  END DO FileLoop

CONTAINS
  SUBROUTINE usage
    IMPLICIT NONE
    CALL getarg (0, clnomprog)

    WRITE(nout, FMT='("Command : ", A)') TRIM(clnomprog)
    WRITE(nout, FMT='("Removes underflows in NetCDF files") ')
    WRITE(nout, FMT='("Usage : ", A, " [-x] [-V] [-d] -v nomvar[,nomvar] nomfic [nomfic]")' ) TRIM(clnomprog)
    WRITE(nout, FMT='("Options : ")' )
    WRITE(nout, FMT='("  -V : mode verbose off. Default is verbose on.")' )
    WRITE(nout, FMT='("  -d : debug mode on. Default is debug off.")' )
    WRITE(nout, FMT='("  -v : gives list of variables to be corrected, separated by a coma.")' )
    WRITE(nout, FMT='("  -x : reverses meaning of -v : given variable are not corrected")' )
    WRITE(nout, FMT='("       if -x is given, and not -v, all variables are corrected.")' )
    
    
    STOP
  END SUBROUTINE usage
  
END PROGRAM ncunderflow
