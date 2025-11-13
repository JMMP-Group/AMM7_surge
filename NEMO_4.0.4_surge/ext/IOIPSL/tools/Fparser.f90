PROGRAM fparser
!-
!$Id: Fparser.f90 2281 2010-10-15 14:21:13Z smasson $
!-
! This software is governed by the CeCILL license
! See IOIPSL/IOIPSL_License_CeCILL.txt

  USE stringop

  IMPLICIT NONE
  !
  !
  !     Parses the code to create the Config.in Config.default and Config.help
  !     which are used by the tk shell.
  !
  !
  INTEGER nbkeymax, nbhelpmax, nbcasemax, nbsourmax, nbelmax
  PARAMETER (nbkeymax=100, nbhelpmax=50, nbcasemax=50, nbsourmax=20,nbelmax=nbhelpmax+10)
  INTEGER nbfilesmax
  PARAMETER (nbfilesmax=150)

  !
  CHARACTER*120 :: configs(nbkeymax,nbelmax)
  CHARACTER*120 :: tmp_help, tmp_key, tmp_desc, tmp_def
  INTEGER :: keylen(nbkeymax), nbkeys
  INTEGER :: key_pos(nbkeymax), help_pos(nbkeymax,2), def_pos(nbkeymax,2)
  INTEGER :: des_pos(nbkeymax), IF_pos(nbkeymax)
  CHARACTER*6 TYPE_op(nbkeymax)
  !
  CHARACTER*120 :: def_out(nbkeymax, nbhelpmax)
  INTEGER :: nbdef_out(nbkeymax)
  !
  CHARACTER*120 :: tke
  !
  CHARACTER*2 :: nbstr
  !
  CHARACTER*80 :: files(nbfilesmax), source(nbsourmax), filetmp
  CHARACTER*80 :: tmp, main_name
  CHARACTER*120 :: keycase(nbcasemax), tmp_CASE
  INTEGER :: nbcase, ii, find, nbsource
  LOGICAL :: next_source, next_name, last_or

  LOGICAL :: is_main, cont

  CHARACTER*1 :: backslash, simplequote, doublequote

  INTEGER :: ia, iread, iret, IFF, ih, nb_line, iv, id
  INTEGER :: ind_space, ind_comma, ind_USE
  INTEGER :: nbfiles, nb_key, nb_key_file
  !
  INTEGER, EXTERNAL ::  iargc, getarg 
  !
  !
  next_source = .FALSE.
  next_name = .FALSE.
  is_main = .FALSE.
  nbsource = 0
  nbfiles = 0
  main_name = 'IPSL'
  !
  backslash = ACHAR(92)
  simplequote = ACHAR(39)
  doublequote = ACHAR(34)
  !
  !
  !
  !     Analyse command line
  !
  !
  !     Get the number of arguments, that is the options and the
  !     files to be parsed.
  !
  !

  iread = iargc()
  !
  DO ia=1,iread
     !
     iret = getarg(ia,tmp)
     !
     IF (next_source) THEN

        nbsource = nbsource + 1
        IF ( nbsource .GT. nbsourmax) THEN
           WRITE(*,*) 'Too many files to source in the arguments.' 
           WRITE(*,*) 'Increase nbsourmax'
           STOP
        ELSE
           source(nbsource) = tmp(1:LEN_TRIM(tmp))
        ENDIF
        next_source = .FALSE.

     ELSE IF (next_name) THEN
        main_name = tmp(1:LEN_TRIM(tmp))
        next_name = .FALSE.

     ELSE
        !     
        IF ( INDEX(tmp,'-m') .GT. 0) THEN
           is_main = .TRUE.
        ELSE IF ( INDEX(tmp,'-n') .GT. 0) THEN
           next_name = .TRUE.
        ELSE IF ( INDEX(tmp,'-s') .GT. 0) THEN
           next_source = .TRUE.
        ELSE IF ( INDEX(tmp,'-h') .GT. 0) THEN
           WRITE(*,*) 'USAGE : Fparse [-name NAME] '
           WRITE(*,*) '               [-source file_to_source]'
           WRITE(*,*) '               [-main] FORTAN_files'
        ELSE
           nbfiles = nbfiles + 1
           IF ( nbfiles .GT. nbfilesmax) THEN
              WRITE(*,*) 'Too many files to include in &
                   &                 the arguments.' 
              WRITE(*,*) 'Increase nbfilesmax'
              STOP
           ELSE
              files(nbfiles) = tmp(1:LEN_TRIM(tmp))
           ENDIF
        ENDIF

     ENDIF

  ENDDO
  !
  IF ( nbfiles .LT. 1 ) THEN
     WRITE(*,*) 'No files provided'
     STOP
  ENDIF
  !
  !
  !     1.0 Read files and extract the lines which we need
  !
  !
  nb_key = 0
  !
  DO IFF=1,nbfiles
     !
     filetmp = files(IFF)
     CALL READ_from_file(filetmp, nbkeymax, nbelmax, configs, nb_key, keylen)
     !
  ENDDO
  !
  !     2.0 Scan the information we have extracted from the file for the elements we need
  !
  !
  CALL analyse_configs(nbkeymax, nb_key, nbelmax, keylen, configs, key_pos, help_pos, def_pos, des_pos, IF_pos, TYPE_op)
  !
  !
  !  3.0 Prepare the default values to put them in an array
  !
  !
  DO ia = 1,nb_key
     !
     !  3.1 Go to blank delimited lines
     !
     nbdef_out(ia) = 0
     !
     DO ii=def_pos(ia,1), def_pos(ia,2)
        !
        tmp_help = configs(ia,ii)
        ind_comma = INDEX(tmp_help(1:len_TRIM(tmp_help)),',')
        DO WHILE (ind_comma .GT. 0) 
           tmp_help(ind_comma:ind_comma) = ' '
           ind_comma = INDEX(tmp_help,',')
        ENDDO
        CALL cmpblank(tmp_help)
        configs(ia,ii) = tmp_help
        !
        !  3.2 extract the values
        !
        tmp_help = TRIM(ADJUSTL(configs(ia,ii)))
        ind_space= INDEX(tmp_help(1:LEN_TRIM(tmp_help)),' ')
        ! Get the first one (there is no space in between)
        IF ( ind_space .EQ. 0) THEN
           nbdef_out(ia) = nbdef_out(ia) + 1
           def_out(ia, nbdef_out(ia)) = tmp_help(1:LEN_TRIM(tmp_help))
        ELSE
           ! Get all those which are before spaces
           DO WHILE (ind_space .GT. 0) 
              nbdef_out(ia) = nbdef_out(ia) + 1
              def_out(ia, nbdef_out(ia)) = tmp_help(1:ind_space)
              tmp_help = ADJUSTL(tmp_help(ind_space+1:LEN_TRIM(tmp_help)))
              ind_space= INDEX(tmp_help(1:LEN_TRIM(tmp_help)),' ')
           ENDDO
           !    Get the last one which does not have a space behind
           IF ( LEN_TRIM(tmp_help) .GT. 0) THEN
              nbdef_out(ia) = nbdef_out(ia) + 1
              def_out(ia, nbdef_out(ia)) = tmp_help(1:LEN_TRIM(tmp_help))
           ENDIF
           !
        ENDIF
     ENDDO
     !
  ENDDO
  !
  !
  !
  !      4.0 OPEN Config.in Defaults and Help files
  !
  !
  OPEN (16, FILE='Config.in')
  OPEN (17, FILE='Config.help')
  OPEN (18, FILE='Config.defaults')
  !
  !     Some explantation
  !
  DO IFF=16,18
     WRITE(IFF,'(1a)') '# '
     WRITE(IFF,'(1a)') '# File created by Fparser, DO NOT EDIT'
     WRITE(IFF,'(2a)') '# ', main_name(1:LEN_TRIM(main_name))
     WRITE(IFF,'(1a)') '# '
     WRITE(IFF,'(1a)') '# '
  ENDDO
  !
  WRITE(17,'(2a)') '# Format of this file: description<nl>', &
       &     ' variable<nl>helptext<nl><nl>.'
  WRITE(17,'(2a)') '# If the question being documented is of', &
       &     ' type "choice", we list'
  WRITE(17,'(2a)') '# only the first occurring config variable.', &
       &     ' The help texts'
  WRITE(17,'(2a)') '# must not contain empty lines. No variable', &
       &     ' should occur twice; if it'
  WRITE(17,'(2a)') '# does, only the first occurrence will be', &
       &     ' used by Configure. The lines'
  WRITE(17,'(2a)') '# in a help text should be indented two', &
       &     ' positions. Lines starting with'
  WRITE(17,'(2a)') '# "#" are ignored. To be nice to menuconfig,', &
       &     ' limit your lines to 70'
  WRITE(17,'(2a)') '# characters. Use emacs" kfill.el to edit', &
       &     ' this file or you lose.'
  WRITE(17,'(2a)') '#'
  !
  IF ( is_main ) THEN
     WRITE(16,'(3a)') 'mainmenu_name "Configuration of model ', &
          &        main_name(1:LEN_TRIM(main_name)), '"'
     WRITE(16,'(1a)') '# '
  ENDIF
  !
  WRITE(16,'(1a)') 'mainmenu_option next_comment'
  WRITE(16,'(3a)') 'comment "', main_name(1:LEN_TRIM(main_name)), '"'
  WRITE(16,'(1a)') '# '
  !
  !   5.0 Loop through the KEYWORDS to prepare the output
  !
  DO IFF =1,nb_key
     !
     !     Config.in file
     !

     !
     !     Is it a conditional option ?
     !
     IF ( IF_pos(IFF) .GE. 0)  THEN
        tmp_help = configs(IFF,IF_pos(IFF))
        !
        IF ( (index(tmp_help,'||') .LE. 0) .AND. (index(tmp_help,'&&') .LE. 0) ) THEN
           IF ( tmp_help(1:1) .EQ. '!') THEN
              WRITE(16,'(3a)') 'if [ "$', tmp_help(2:LEN_TRIM(tmp_help)),  '" = "n" ]; then'
           ELSE
              WRITE(16,'(3a)') 'if [ "$', tmp_help(1:LEN_TRIM(tmp_help)),  '" = "y" ]; then'
           ENDIF
        ELSE
           !
           last_or = .TRUE.
           nbcase = 0
           !
           DO WHILE( INDEX(tmp_help,'||') .GT. 0)
              ii = INDEX(tmp_help,'||')
              nbcase = nbcase + 1
              if ( nbcase .EQ. 1 ) THEN
                 IF ( tmp_help(1:1) .EQ. '!') THEN
                    WRITE(16,'(3a)') 'if [ "$', tmp_help(2:ii-1),  '" = "n" \\'
                 ELSE
                    WRITE(16,'(3a)') 'if [ "$', tmp_help(1:ii-1),  '" = "y" \\'
                 ENDIF
              ELSE
                 IF ( tmp_help(1:1) .EQ. '!') THEN
                    WRITE(16,'(3a)') '-o "$', tmp_help(2:ii-1),  '" = "n" \\'
                 ELSE
                    WRITE(16,'(3a)') '-o "$', tmp_help(1:ii-1),  '" = "y" \\'
                 ENDIF
              ENDIF
              tmp_help = TRIM(ADJUSTL(tmp_help(ii+2:LEN_TRIM(tmp_help))))
           ENDDO
           !
           DO WHILE( INDEX(tmp_help,'&&') .GT. 0)
              ii = INDEX(tmp_help,'&&')
              nbcase = nbcase + 1
              if ( nbcase .EQ. 1 ) THEN
                 IF ( tmp_help(1:1) .EQ. '!') THEN
                    WRITE(16,'(3a)') 'if [ "$', tmp_help(2:ii-1),  '" = "n" \\'
                 ELSE
                    WRITE(16,'(3a)') 'if [ "$', tmp_help(1:ii-1),  '" = "y" \\'
                 ENDIF
              ELSE
                 IF ( tmp_help(1:1) .EQ. '!') THEN
                    WRITE(16,'(3a)') '-a "$', tmp_help(2:ii-1),  '" = "n" \\'
                 ELSE
                    WRITE(16,'(3a)') '-a "$', tmp_help(1:ii-1),  '" = "y" \\'
                 ENDIF
              ENDIF
              tmp_help = TRIM(ADJUSTL(tmp_help(ii+2:LEN_TRIM(tmp_help))))
              last_or = .FALSE.
           ENDDO
           !
           IF ( last_or ) THEN
              IF ( tmp_help(1:1) .EQ. '!') THEN
                 WRITE(16,'(3a)') '-o "$', tmp_help(2:LEN_TRIM(tmp_help)),  '" = "n" ]; then'
              ELSE
                 WRITE(16,'(3a)') '-o "$', tmp_help(1:LEN_TRIM(tmp_help)),  '" = "y" ]; then'
              ENDIF
           ELSE
              IF ( tmp_help(1:1) .EQ. '!') THEN
                 WRITE(16,'(3a)') '-a "$', tmp_help(2:LEN_TRIM(tmp_help)),  '" = "n" ]; then'
              ELSE
                 WRITE(16,'(3a)') '-a "$', tmp_help(1:LEN_TRIM(tmp_help)),  '" = "y" ]; then'
              ENDIF
           ENDIF
        ENDIF
        WRITE(16,'(1a)') '       '
     ENDIF
     !
     !      Extract the information from configs
     !
     DO iv = 1,nbdef_out(IFF)

        IF (nbdef_out(IFF) .EQ. 1) THEN
           tmp_key = configs(IFF,key_pos(IFF))
           tmp_desc = configs(IFF,des_pos(IFF))
           tmp_def = def_out(IFF,iv)
        ELSE
           tmp_key = configs(IFF,key_pos(IFF))
           WRITE(nbstr,'(I2.2)') iv
           tmp_key = tmp_key(1:LEN_TRIM(tmp_key))//'__'//nbstr
           tmp_desc = configs(IFF,des_pos(IFF))
           IF ( iv .EQ. 1) THEN
              tmp_desc = tmp_desc(1:LEN_TRIM(tmp_desc))//' (Vector)'
           ELSE
              tmp_desc = 'Cont...    '//tmp_key(1:LEN_TRIM(tmp_key))
           ENDIF
           tmp_def = def_out(IFF,iv)
        ENDIF
        !
        !
        !
        IF (INDEX(TYPE_op(IFF),'bool') .GT. 0) THEN
           !
           WRITE(16,'(4a)') 'bool "', tmp_desc(1:LEN_TRIM(tmp_desc)), &
                &              '" ',tmp_key(1:LEN_TRIM(tmp_key))
           !
        ELSE IF (INDEX(TYPE_op(IFF),'hex') .GT. 0) THEN
           !
           WRITE(16,'(6a)') 'hex "', tmp_desc(1:LEN_TRIM(tmp_desc)) &
                &              ,'" ',tmp_key(1:LEN_TRIM(tmp_key)) &
                &              ,' ',tmp_def(1:LEN_TRIM(tmp_def))
           !
        ELSE IF (INDEX(TYPE_op(IFF),'choice') .GT. 0) THEN
           !
           !   Get number of options
           !
           nbcase = 0
           DO WHILE( INDEX(tmp_key,'||') .GT. 0)
              ii = INDEX(tmp_key,'||')
              nbcase = nbcase + 1
              keycase(nbcase) = tmp_key(1:ii-1)
              tmp_key=tmp_key(ii+2:LEN_TRIM(tmp_key))
           ENDDO
           nbcase = nbcase + 1
           keycase(nbcase) = tmp_key(1:LEN_TRIM(tmp_key))

           WRITE(16,'(4a)') "choice '", tmp_desc(1:LEN_TRIM(tmp_desc))," '",backslash
           !
           !   List options
           !
           tmp_CASE = keycase(1)
           WRITE(16,'(5a)') '        "', tmp_CASE(1:LEN_TRIM(tmp_CASE)), "          "&
                &,tmp_CASE(1:LEN_TRIM(tmp_CASE)), backslash
           !
           DO ii=2,nbcase-1
              tmp_CASE = keycase(ii)
              WRITE(16,'(5a)') '         ',   tmp_CASE(1:LEN_TRIM(tmp_CASE)),  '          ',&
                   & tmp_CASE(1:LEN_TRIM(tmp_CASE)),  backslash
           ENDDO
           !
           tmp_CASE = keycase(nbcase)
           WRITE(16,'(6a)') '         ', &
                &              tmp_CASE(1:LEN_TRIM(tmp_CASE)), &
                &              '          ', tmp_CASE(1:LEN_TRIM(tmp_CASE)), &
                &              '"  ',tmp_def(1:LEN_TRIM(tmp_def)) 
           !
        ELSE
           WRITE(*,'(2a)') 'Uniplemented operation : ', TYPE_op(IFF)
           STOP
        ENDIF
        !
        !     Config.help file
        !
        tmp_key = configs(IFF,key_pos(IFF))
        IF (INDEX(TYPE_op(IFF),'choice') .GT. 0) THEN
           ii = INDEX(tmp_key,'||')-1
        ELSE
           ii = LEN_TRIM(tmp_key)
        ENDIF

        IF ( nbdef_out(IFF) .GT. 1) THEN
           WRITE(17,'(1a)') tmp_desc(1:LEN_TRIM(tmp_desc))
           WRITE(nbstr,'(I2.2)') iv
           tke = tmp_key(1:ii)//'__'//nbstr
           WRITE(17,'(1a)') tke(1:LEN_TRIM(tke))
           WRITE(17,'(1a)') '  (Vector)'
        ELSE
           WRITE(17,'(1a)') tmp_desc(1:LEN_TRIM(tmp_desc))
           WRITE(17,'(1a)') tmp_key(1:ii)
        ENDIF
        !
        DO ih=help_pos(IFF,1),help_pos(IFF,2)
           tmp_help = configs(IFF,ih)
           WRITE(17,'("  ",1a)') tmp_help(1:LEN_TRIM(tmp_help))
        ENDDO
        !
        !     Config.default file
        !
        IF (INDEX(TYPE_op(IFF),'choice') .GT. 0) THEN

           WRITE(18,'(2a)') tmp_def(1:LEN_TRIM(tmp_def)),'=y'

        ELSE

           WRITE(18,'(3a)') tmp_key(1:LEN_TRIM(tmp_key)),'=', &
                &              tmp_def(1:LEN_TRIM(tmp_def))

        ENDIF
        !     
        !     Add some empty line to all files
        !     
        WRITE(16,'(1a)') '       '
        WRITE(17,'(1a)') '       '
        WRITE(17,'(1a)') '       '
     ENDDO
     !
     !
     !     Close the IF if needed
     !

     IF ( IF_pos(IFF) .GT. 0) THEN
        WRITE(16,'(1a)') 'fi'
        WRITE(16,'(1a)') '       '
     ENDIF

     !
  ENDDO
  !
  WRITE(16,'(1a)') 'endmenu'
  WRITE(16,'(1a)') '       '
  IF ( nbsource .GT. 0) THEN
     DO ih=1,nbsource
        tmp = source(ih)
        WRITE(16,'(1a)') '              '
        WRITE(16,'(3a)') 'source ',tmp(1:LEN_TRIM(tmp)), &
             &           '/Config.in'
     ENDDO
  ENDIF
  !
  !
  CLOSE(16)
  CLOSE(17)
  CLOSE(18)
  !
  !
  !
  STOP

END PROGRAM fparser
!
!
!==========================================================
!
!
SUBROUTINE READ_from_file(file, nbkeymax, nbelmax, configs, nbitems, itemlen)
  !
  USE stringop
  !
  IMPLICIT NONE
  !
  !
  !     This routine reads the file and adds the config info it finds to the configs array.
  !     Thus the nbitems is an imput variable as it can be increased as we go through the files.
  !
  !
  CHARACTER*(*) :: file
  INTEGER :: nbkeymax, nbelmax
  CHARACTER*120 :: configs(nbkeymax, nbelmax)
  INTEGER ::  nbitems, itemlen(nbkeymax)
  !
  INTEGER :: conf_pos, ip
  CHARACTER*250 line
  LOGICAL :: cont, conf_END
  !
  cont = .TRUE.
  conf_END = .TRUE.
  !
  OPEN (12, file=file)
  !
  !   1.0 Loop over all the lines of a given file to extract all the configuration line
  !
  DO WHILE (cont)
    READ(12,'(a)',END=9999) line
    !
    !   1.0  A configuration line is detected by the line below.
    !
    IF ( INDEX(line,'Config') .EQ. 1 .OR. INDEX(line,'!'//'Config') .GE. 1 ) THEN
        !
        IF ( conf_END ) THEN
            nbitems = nbitems + 1
            IF ( nbitems .GT. nbkeymax) THEN
                WRITE(*,*) 'read_from_file : The number of keys in the input array is too small for this file'
                STOP
            ENDIF
            itemlen(nbitems) = 0
            conf_END = .FALSE.
        ENDIF
        !
        itemlen(nbitems) = itemlen(nbitems) + 1
        IF ( itemlen(nbitems) .GT. nbelmax ) THEN
            WRITE(*,*) 'read_from_file : The number of elements per key in the input array is too small'
            STOP
        ENDIF
        !
        !  The detected line is shaved !
        !
        IF ( INDEX(line,'Config') .EQ. 1) THEN
            conf_pos = 7
        ELSE
            conf_pos = INDEX(line,'!'//'Config') +7
        ENDIF
        line = line(conf_pos:LEN_TRIM(line))
        line = TRIM(ADJUSTL(line))
        CALL cmpblank(line)
        !
        configs(nbitems,itemlen(nbitems)) = line
        !
    ELSE
        !
        !   Look for the end of a configuration structure.
        !  It is determined by a call to the getin subroutine
        !
        CALL strlowercase(line)
        CALL cmpblank(line)
        ip = INDEX(line,' (')
        DO WHILE (ip .GT. 0)
          line = line(1:ip-1)//line(ip+1:LEN_TRIM(line))
          ip = INDEX(line,' (')
        ENDDO
        IF ( INDEX(line, 'call getin(') .GT. 0 .OR. INDEX(line, 'call setvar(') .GT. 0) THEN
            conf_END = .TRUE.
        ENDIF
        !
    ENDIF
    !
    cont = .TRUE.
    GOTO 8888
9999 cont = .FALSE.
8888 CONTINUE
    
    ENDDO
  !
  CLOSE(12)
  !
  END SUBROUTINE READ_from_file
  !
  !==========================================================
  !
  !
  SUBROUTINE analyse_configs(nbkmax, nb_key, nbelmax, keylen, configs, key_pos, help_pos, def_pos, des_pos, IF_pos, TYPE_op)
    !
  USE stringop
  !
    IMPLICIT NONE
    !
    !
    !    This subroutine will localize the KEYWORDS in the configs array
    !    and extract all their arguments. For the moment 5 arguments are recognized :
    !    KEY  : The keyword by which the all is identified
    !    HELP : This identifies the help text
    !    DEF : The default value of for this KEYWORD
    !    DESC : A short description, not more than one line
    !    IF : Specifies the other Keyword it depend on. This is a nice features for the menus as it can hide 
    !            things we do not need
    !
    !    The DEF and HELP keywords can be multi line
    !
    INTEGER :: nbkmax, nb_key, nbelmax
    INTEGER :: keylen(nbkmax)
    INTEGER :: key_pos(nbkmax), help_pos(nbkmax,2), def_pos(nbkmax,2), des_pos(nbkmax), IF_pos(nbkmax)
    CHARACTER*120 :: configs(nbkmax,nbelmax)
    CHARACTER*6 :: TYPE_op(nbkmax)
    !
    !   This is the number of arguments we need to find an end for and the total number of arguments  we can have.
    !   Thus these parameters needs to be updated when the list of arguments to the routine is changed
    !
    INTEGER, PARAMETER :: toendlen=2, indexlen=5
    !
    INTEGER :: toend(toendlen), foundend(toendlen), kindex(indexlen)
    INTEGER :: ik, il, ieq
    CHARACTER*120 :: tmp_str, tmp_str2
    !
    !
    key_pos(1:nb_key)=-1
    help_pos(1:nb_key,1:2)=-1
    def_pos(1:nb_key,1:2)=-1
    des_pos(1:nb_key)=-1
    IF_pos(1:nb_key)=-1
    TYPE_op(1:nb_key)='hex'
    !
    DO ik=1,nb_key
      !
      !
      DO il=1,keylen(ik)
        !
        ieq = INDEX(configs(ik,il),'=')
        tmp_str = configs(ik,il)
        tmp_str = tmp_str(1:ieq)
        CALL struppercase(tmp_str)
        !
        !      Decide if this is a reserved name and where it fits
        !
        !      At the same time we clean up the configs array
        !
        IF ( INDEX(tmp_str,'KEY') .GT. 0) THEN
            IF ( key_pos(ik) .GT. 0) THEN
                WRITE(*,*) 'analyse_config : Already have a KEYWORD, check that you have a call to getin'
                WRITE(*,*) 'analyse_config : ', configs(ik,il)
                STOP
            ENDIF
            key_pos(ik) = il
            tmp_str2 = configs(ik,il)
            tmp_str2 = tmp_str2(ieq+1:LEN_TRIM(tmp_str2))
            configs(ik,il) = TRIM(ADJUSTL(tmp_str2))
            !
            !   Here we have to check that we are not in an 'choice' case
            !
            IF ( INDEX(tmp_str2,'||') .GT. 0) THEN
                TYPE_op(ik) = 'choice'
            ENDIF
            !
        ENDIF
        !
        IF ( INDEX(tmp_str,'DEF') .GT. 0) THEN
            IF ( def_pos(ik,1) .GT. 0) THEN
                WRITE(*,*) 'analyse_config : Already have a DEF, check that you have a call to getin'
                WRITE(*,*) 'analyse_config : ', configs(ik,il)
                STOP
            ENDIF
            def_pos(ik,1) = il
            tmp_str2 = configs(ik,il)
            tmp_str2 = tmp_str2(ieq+1:LEN_TRIM(tmp_str2))
            tmp_str2 = TRIM(ADJUSTL(tmp_str2))
            configs(ik,il) = tmp_str2 
            !
            !  Here we can check if we have a boolean operation
            !  We also wish to standardise the value of booleans
            !
            CALL struppercase(tmp_str2)
            IF (INDEX(tmp_str2,'Y') .EQ. 1 .AND. LEN_TRIM(tmp_str2) .EQ. 1 .OR.&
               & INDEX(tmp_str2,'T') .EQ. 1 .AND. LEN_TRIM(tmp_str2) .EQ. 1 .OR.&
               & INDEX(tmp_str2,'YES') .EQ. 1 .AND. LEN_TRIM(tmp_str2) .EQ. 3 .OR.&
               & INDEX(tmp_str2,'TRUE') .EQ. 1 .AND. LEN_TRIM(tmp_str2) .EQ. 4 .OR.&
               & INDEX(tmp_str2,'.TRUE.') .EQ. 1) THEN
                configs(ik,il) = 'y'
                TYPE_op(ik) = 'bool'
            ENDIF
            !
            IF (INDEX(tmp_str2,'N') .EQ. 1 .AND. LEN_TRIM(tmp_str2) .EQ. 1 .OR.&
               & INDEX(tmp_str2,'F') .EQ. 1 .AND. LEN_TRIM(tmp_str2) .EQ. 1 .OR.&
               & INDEX(tmp_str2,'NO') .EQ. 1 .AND. LEN_TRIM(tmp_str2) .EQ. 2 .OR.&
               & INDEX(tmp_str2,'FALSE') .EQ. 1  .AND. LEN_TRIM(tmp_str2) .EQ. 5 .OR.&
               & INDEX(tmp_str2,'.FALSE.') .EQ. 1) THEN
                configs(ik,il) = 'n'
                TYPE_op(ik) = 'bool'
            ENDIF
            !
            ! Here we check if we have a default behavior and put a standard name
            !
            IF (INDEX(tmp_str2,'DEF') .EQ. 1 .OR. INDEX(tmp_str2,'NONE') .EQ. 1) THEN
                configs(ik,il) = 'default'
            ENDIF
            !
        ENDIF
        !
        IF ( INDEX(tmp_str,'DESC') .GT. 0) THEN
            IF ( des_pos(ik) .GT. 0) THEN
                WRITE(*,*) 'analyse_config : Already have a DESC, check that you have a call to getin'
                WRITE(*,*) 'analyse_config : ', configs(ik,il)
                STOP
            ENDIF
            des_pos(ik) = il
            tmp_str2 = configs(ik,il)
            tmp_str2 = tmp_str2(ieq+1:LEN_TRIM(tmp_str2))
            configs(ik,il) = TRIM(ADJUSTL(tmp_str2))
        ENDIF
        !
        IF ( INDEX(tmp_str,'IF') .GT. 0) THEN
            IF ( IF_pos(ik) .GT. 0) THEN
                WRITE(*,*) 'analyse_config : Already have a IF, check that you have a call to getin'
                WRITE(*,*) 'analyse_config : ', configs(ik,il)
                STOP
            ENDIF
            IF_pos(ik) = il
            tmp_str2 = configs(ik,il)
            tmp_str2 = tmp_str2(ieq+1:LEN_TRIM(tmp_str2))
            configs(ik,il) = TRIM(ADJUSTL(tmp_str2))
        ENDIF 
        !
        IF ( INDEX(tmp_str,'HELP') .GT. 0) THEN
            help_pos(ik,1) = il
            tmp_str2 = configs(ik,il)
            tmp_str2 = tmp_str2(ieq+1:LEN_TRIM(tmp_str2))
            configs(ik,il) = TRIM(ADJUSTL(tmp_str2))
        ENDIF
        !
      ENDDO
      !
      !     Check if we not missing some important informations as for instance
      !
      !     THE KEYWORD 
      !
      IF ( key_pos(ik) .LT. 1) THEN
          WRITE(*,*) 'analyse_configs : Could not find a keyword in the following entry :'
          DO il=1,keylen(ik)
            WRITE(*,'(a70)') configs(ik,il)
          ENDDO
          STOP
      ENDIF
      !
      !      THE DEFAULT VALUE
      !
      IF ( def_pos(ik,1) .LT. 1) THEN
          WRITE(*,*) 'analyse_configs : Could not find a default value in the following entry :'
          DO il=1,keylen(ik)
            WRITE(*,'(a70)') configs(ik,il)
          ENDDO
          STOP
      ENDIF
      !
      !   Get the end of all the multi line arguments
      !
      toend(1) = MAX(def_pos(ik,1),1)
      toend(2) = MAX(help_pos(ik,1),1)
      foundend(:) = keylen(ik)
      kindex(1) = MAX(key_pos(ik),1)
      kindex(2) = MAX(des_pos(ik),1)
      kindex(3) = MAX(def_pos(ik,1),1)
      kindex(4) = MAX(IF_pos(ik),1)
      kindex(5) = MAX(help_pos(ik,1),1)
      CALL find_ends(toendlen, toend, indexlen, kindex, foundend)
      def_pos(ik,2) = foundend(1)
      help_pos(ik,2) = foundend(2)
      !
    ENDDO
    !
  END SUBROUTINE analyse_configs
  !
  SUBROUTINE find_ends(toendlen, toend, indexlen, kindex, foundend)
    !
    IMPLICIT NONE
    !
    !
    !  We find the end of the text for all the elements in the key which are multi line
    !   This subroutine aims at providing a flexible way to determine this so that other
    !   elements in the Keyword can be multi line. For the moment it is only the Help and Ded 
    !  which are allowed to be multi line.
    !
    !  Foundend need to be initialized to the maximum value of the elements
    !
    !
    INTEGER :: toendlen, toend(toendlen), indexlen, kindex(indexlen), foundend(toendlen)
    !
    INTEGER :: whmin(1), ie, ii
    !
    DO ie=1,toendlen
      !
      whmin = MINLOC(toend(1:toendlen))
      !
      DO ii=1,indexlen
        IF ( kindex(ii) .GT. toend(whmin(1)) .AND. foundend(whmin(1)) .GE. kindex(ii)) THEN
            foundend(whmin(1)) = kindex(ii)-1
            toend(whmin(1)) = 100000
        ENDIF
      ENDDO
      !
    ENDDO
    !
  END SUBROUTINE find_ends
