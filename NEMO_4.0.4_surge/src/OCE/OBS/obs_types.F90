MODULE obs_types
   !!=====================================================================
   !!                       ***  MODULE  obs_types  ***
   !! Observation diagnostics: Observation type integer to character 
   !!                          translation
   !!=====================================================================

   !!---------------------------------------------------------------------
   !! 
   !!   The NetCDF variable CWMO_INST_TYP_COR is used to identify the 
   !!   different instrument types for coriolis data.
   !!
   !!    WMO  NEMOVAR  TYPE                   DESCRIPTION
   !!    ---  -------  ----   --------------------------------------------
   !!    800     0     MBT    (1941-) mechanical bathythermograph data
   !!    401     1     XBT    (1967-) expendable bathythermograph data
   !!    830     2     CTD    (1967-) high resolution CTD data
   !!    820     3     MRB    (1990-) moored buoy data
   !!    831     4     PFL    (1994-) profiling float data
   !!    995     5     DRB    (1998-) drifting buoy data
   !!    997     6     APB    (1997-) autonomous pinniped bathythermograph
   !!    996     7     UOR    (1992-) undulating oceanographic recorder
   !!    741     8     OSD    (1800-) low resolution (bottle) CTD data
  !!
   !! History :
   !!        !  06-03  (K. Mogensen) Original code
   !!        !  06-10  (A. Weaver) Cleanup
   !!---------------------------------------------------------------------

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE 

   !! * Shared Module variables

   INTEGER, PUBLIC, PARAMETER :: ntyp1770 = 1023
!RBbug useless ?   CHARACTER(LEN=4), PUBLIC, DIMENSION(0:ntyp1770) :: cwmotyp1770
   CHARACTER(LEN=80), PUBLIC, DIMENSION(0:ntyp1770) :: cwmonam1770
   CHARACTER(LEN=3), PUBLIC, DIMENSION(0:ntyp1770) :: ctypshort

   INTEGER, PUBLIC, PARAMETER :: ntypalt = 8
   CHARACTER(LEN=40), PUBLIC, DIMENSION(0:ntypalt) :: calttyp

   PUBLIC obs_typ_init
   PUBLIC obs_wmo_init
   PUBLIC obs_alt_typ_init

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS 
   
   SUBROUTINE obs_typ_init
      !!---------------------------------------------------------------------
      !!
      !!                   *** ROUTINE obs_wmo_init ***
      !!
      !! ** Purpose : Initialize code tables
      !!
      !! ** Method  : 
      !!
      !! ** Action  : 
      !!
      !! References : 
      !!
      !! History :  
      !!      ! :  2007-06 (K. Mogensen) Original code
      !!----------------------------------------------------------------------
      !! * Modules used
      
      !! * Arguments

      !! * Local declarations

      CALL obs_wmo_init

      CALL obs_alt_typ_init

   END SUBROUTINE obs_typ_init

   SUBROUTINE obs_wmo_init
      !!---------------------------------------------------------------------
      !!
      !!                   *** ROUTINE obs_wmo_init ***
      !!
      !! ** Purpose : Initialize WMO code 1770 code tables
      !!
      !! ** Method  : Direct initialisation of variables
      !!
      !! ** Action  : 
      !!
      !! References : WORLD  METEOROLOGICAL  ORGANIZATION
      !!              Manual on Codes
      !!              International Codes
      !!              VOLUME I.1 (Annex II to WMO Technical Regulations)
      !!              Part A -- Alphanumeric Codes
      !!              1995 edition
      !!              WMO-No. 306 
      !!              Secretariat of the World Meteorological Organization 
      !!              Geneva, Switzerland
      !!
      !! History :  
      !!      ! :  2007-04 (K. Mogensen) Original code
      !!----------------------------------------------------------------------
      !! * Modules used
      
      !! * Arguments

      !! * Local declarations
      INTEGER :: ji

      DO ji = 0, ntyp1770

         cwmonam1770(ji) = 'Not defined'
         ctypshort(ji) = '---'

!         IF ( ji < 1000 ) THEN
!            WRITE(cwmotyp1770(ji),'(1X,I3.3)') ji
!         ELSE
!            WRITE(cwmotyp1770(ji),'(I4.4)') ji
!         ENDIF

      END DO

      cwmonam1770(   1) = 'Sippican T-4'
      cwmonam1770(   2) = 'Sippican T-4'
      cwmonam1770(  11) = 'Sippican T-5'
      cwmonam1770(  21) = 'Sippican Fast Deep'
      cwmonam1770(  31) = 'Sippican T-6'
      cwmonam1770(  32) = 'Sippican T-6'
      cwmonam1770(  41) = 'Sippican T-7'
      cwmonam1770(  42) = 'Sippican T-7'
      cwmonam1770(  51) = 'Sippican Deep Blue'
      cwmonam1770(  52) = 'Sippican Deep Blue'
      cwmonam1770(  61) = 'Sippican T-10'
      cwmonam1770(  71) = 'Sippican T-11'
      cwmonam1770( 201) = 'TSK T-4'
      cwmonam1770( 202) = 'TSK T-4'
      cwmonam1770( 211) = 'TSK T-6'
      cwmonam1770( 212) = 'TSK T-6'
      cwmonam1770( 221) = 'TSK T-7'
      cwmonam1770( 222) = 'TSK T-7'
      cwmonam1770( 231) = 'TSK T-5'
      cwmonam1770( 241) = 'TSK T-10'
      cwmonam1770( 251) = 'TSK Deep Blue'
      cwmonam1770( 252) = 'TSK Deep Blue'
      cwmonam1770( 261) = 'TSK AXBT '
      cwmonam1770( 401) = 'Sparton XBT-1'
      cwmonam1770( 411) = 'Sparton XBT-3'
      cwmonam1770( 421) = 'Sparton XBT-4'
      cwmonam1770( 431) = 'Sparton XBT-5'
      cwmonam1770( 441) = 'Sparton XBT-5DB'
      cwmonam1770( 451) = 'Sparton XBT-6'
      cwmonam1770( 461) = 'Sparton XBT-7'
      cwmonam1770( 462) = 'Sparton XBT-7'
      cwmonam1770( 471) = 'Sparton XBT-7DB'
      cwmonam1770( 481) = 'Sparton XBT-10'
      cwmonam1770( 491) = 'Sparton XBT-20'
      cwmonam1770( 501) = 'Sparton XBT-20DB'
      cwmonam1770( 510) = 'Sparton 536 AXBT'
      cwmonam1770( 700) = 'Sippican XCTD standard'
      cwmonam1770( 710) = 'Sippican XCTD deep'
      cwmonam1770( 720) = 'Sippican AXCTD'
      cwmonam1770( 730) = 'Sippican SXCTD'
      cwmonam1770( 741) = 'TSK XCTD'
      cwmonam1770( 742) = 'TSK XCTD-2 '
      cwmonam1770( 743) = 'TSK XCTD-2F '
      cwmonam1770( 751) = 'TSK AXCTD '
      cwmonam1770( 800) = 'Mechanical BT'
      cwmonam1770( 810) = 'Hydrocast'
      cwmonam1770( 820) = 'Thermistor Chain'
      cwmonam1770( 825) = 'Temperature (sonic) and pressure probes'
      cwmonam1770( 830) = 'CTD'
      cwmonam1770( 831) = 'CTD-P-ALACE float'
      cwmonam1770( 840) = 'PROVOR, No conductivity sensor '
      cwmonam1770( 841) = 'PROVOR, Seabird conductivity sensor '
      cwmonam1770( 842) = 'PROVOR, FSI conductivity sensor '
      cwmonam1770( 845) = 'Web Research, No conductivity sensor '
      cwmonam1770( 846) = 'Web Research, Seabird conductivity sensor '
      cwmonam1770( 847) = 'Web Research. FSI conductivity sensor'
      cwmonam1770( 850) = 'SOLO, No conductivity sensor '
      cwmonam1770( 851) = 'SOLO, Seabird conductivity sensor '
      cwmonam1770( 852) = 'SOLO, FSI conductivity sensor'
      cwmonam1770( 855) = 'Profiling float, NINJA, no conductivity sensor'
      cwmonam1770( 856) = 'Profiling float, NINJA, SBE conductivity sensor'
      cwmonam1770( 857) = 'Profiling float, NINJA, FSI conductivity sensor'
      cwmonam1770( 858) = 'Profiling float, NINJA, TSK conductivity sensor'
      cwmonam1770( 900) = 'Sippican T-12 XBT'
      cwmonam1770(1023) = 'Missing value'

      DO ji = 853, 854
         cwmonam1770(ji) = 'Reserved'
      END DO

      DO ji = 859, 899
         cwmonam1770(ji) = 'Reserved'
      END DO

      DO ji = 901, 999
         cwmonam1770(ji) = 'Reserved'
      END DO

      DO ji = 1000, 1022
         cwmonam1770(ji) = 'Reserved'
      END DO

      ctypshort(800) = 'MBT'
      ctypshort(401) = 'XBT'
      ctypshort(830) = 'CTD'
      ctypshort(820) = 'MRB'
      ctypshort(831) = 'PFL'
      ctypshort(995) = 'DRB'
      ctypshort(997) = 'APB'
      ctypshort(996) = 'UOR'
      ctypshort(700:799) = 'OSD'

   END SUBROUTINE obs_wmo_init

   SUBROUTINE obs_alt_typ_init
      !!---------------------------------------------------------------------
      !!
      !!                   *** ROUTINE obs_alt_typ_init ***
      !!
      !! ** Purpose : Initialize CLS altimeter code tables
      !!
      !! ** Method  : Direct initialisation of variables
      !!
      !! ** Action  : 
      !!
      !! References : CLS
      !1              SSALTO/DUACS User Handbook
      !!              (M)SLA and (M)ADT Near-Real Time and
      !!              Delayed time products
      !!              CLS-DOS-NT-06-034
      !!              2006
      !!              CLS
      !!              8-10 Rue Hermes
      !!              Parc Technologique du Canal
      !!              31526 Ramonville St-Agne
      !!              France
      !!
      !! History :  
      !!      ! :  2007-06 (K. Mogensen) Original code
      !!----------------------------------------------------------------------
      !! * Modules used
      
      !! * Arguments

      !! * Local declarations

      calttyp(0) = 'Unknown'
      calttyp(1) = 'ERS-1'
      calttyp(2) = 'ERS-2'
      calttyp(3) = 'Topex/Poseidon'
      calttyp(4) = 'Topex/Poseidon on its new orbit'
      calttyp(5) = 'GFO'
      calttyp(6) = 'Jason-1'
      calttyp(7) = 'Envisat'
      calttyp(8) = 'Jason-2'

   END SUBROUTINE obs_alt_typ_init

END MODULE obs_types
