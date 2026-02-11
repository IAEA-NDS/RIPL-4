C 29/10/06                                                      ECIS06  LDIS-000
      SUBROUTINE LDIS(WV,IPP,PIP,NVAL,VAL,LO)                           LDIS-001
C INPUT OF DISPERSION PARAMETERS.                                       LDIS-002
C INPUT:     WV(J,I): CENTRE OF MASS ENERGY OF THE LEVEL I FOR J=3.     LDIS-003
C                     LABORATORY ENERGY OF THE LEVEL I FOR J=13.        LDIS-004
C            IPP(I,J):FIRST LEVEL USING POTENTIAL J FOR I=1 (TEMPORARY).LDIS-005
C                     -1 TO READ DISPERSION PARAMETERS FOR I=2.         LDIS-006
C            NVAL,VAL:OPTICAL POTENTIALS.                               LDIS-007
C            LO(I):   LOGICAL CONTROLS:                                 LDIS-008
C               LO(7)  =.TRUE. MATRIX ELEMENT AND FORM FACTORS READ.    LDIS-009
C               LO(109)=.TRUE. FOR DIRAC POTENTIALS.                    LDIS-010
C OUTPUT:    IPP,PIP: EQUIVALENT BY CALL, PARAMETERS OF THE DISPERSION  LDIS-011
C                     RELATIONS (SEE INPUT DESCRIPTION).                LDIS-012
C                                                                       LDIS-013
C FOR THE COMMON  /COUPL/ SEE CALC.                                     LDIS-014
C                                                                       LDIS-015
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /COUPL/:                     LDIS-016
C  NPP:       NUMBER OF OPTICAL POTENTIALS.                             LDIS-017
C   USED:     NPP.                                                      LDIS-018
C                                                                       LDIS-019
C***********************************************************************LDIS-020
      IMPLICIT REAL*8 (A-H,O-Z)                                         LDIS-021
      LOGICAL LO(150)                                                   LDIS-022
      CHARACTER*4 CLEG(6)                                               LDIS-023
      DIMENSION WV(22,*),IPP(34,*),PIP(17,*),NVAL(2,*),VAL(42,*),VV(6)  LDIS-024
      COMMON /COUPL/ IQM,IQMAX,NBT1,NFA,NPP,NSPIN,NVA                   LDIS-025
      COMMON /INOUT/ MR,MW,MS                                           LDIS-026
      DATA CLEG /' VOL','UME ',' SCA','LAR ',' VEC','TOR '/             LDIS-027
      DO 13 IJ=1,NPP                                                    LDIS-028
C COLLECT SIX DEPTHS OF POTENTIALS.                                     LDIS-029
      IP=IPP(1,IJ)                                                      LDIS-030
      IF (LO(7)) GO TO 2                                                LDIS-031
      DO 1 I=1,6                                                        LDIS-032
    1 VV(I)=VAL(4*I-3,IP)                                               LDIS-033
      GO TO 4                                                           LDIS-034
    2 DO 3 I=1,6                                                        LDIS-035
    3 VV(I)=VAL(NVAL(1,8*IP+I-6),1)                                     LDIS-036
C INPUT OF DISPERSION PARAMETERS.                                       LDIS-037
    4 READ (MR,1000) IPPIP,(IPP(J,IJ),J=2,6)                            LDIS-038
      IF (IPP(6,IJ).NE.1) GO TO 5                                       LDIS-039
      WRITE (MW,1001) IJ                                                LDIS-040
      GO TO 13                                                          LDIS-041
    5 IF (IPPIP.NE.1) IPP(1,IJ)=-IPP(1,IJ)                              LDIS-042
      IF (MOD(IPP(2,IJ),2).NE.0) GO TO 14                               LDIS-043
      READ (MR,1002) (PIP(J,IJ),J=4,17)                                 LDIS-044
      DO 6 I=1,3                                                        LDIS-045
      IF (MOD(IPP(I+2,IJ),2).NE.0) GO TO 14                             LDIS-046
      IF ((IPP(I+2,IJ).NE.0).AND.(VV(2*I).EQ.0.D0)) GO TO 15            LDIS-047
    6 CONTINUE                                                          LDIS-048
      IF ((PIP(8,IJ).NE.0.D0).AND.(VV(5).EQ.0.D0)) GO TO 16             LDIS-049
      IF ((PIP(9,IJ).NE.0.D0).AND.(VV(6).EQ.0.D0)) GO TO 17             LDIS-050
      IF ((PIP(17,IJ).NE.0.D0).AND.(VV(1).EQ.0.D0)) GO TO 18            LDIS-051
      IF (IPP(1,IJ).GT.0) GO TO 7                                       LDIS-052
      WRITE (MW,1003)                                                   LDIS-053
      IZ=3                                                              LDIS-054
      GO TO 8                                                           LDIS-055
    7 WRITE (MW,1004)                                                   LDIS-056
      IZ=13                                                             LDIS-057
    8 IF (PIP(4,IJ).EQ.0.D0) PIP(4,IJ)=WV(IZ,IP)                        LDIS-058
      IF (PIP(5,IJ).EQ.0.D0) PIP(5,IJ)=-6.8D0                           LDIS-059
      IF (PIP(6,IJ).EQ.0.D0) PIP(6,IJ)=PIP(5,IJ)                        LDIS-060
      WRITE (MW,1005) (PIP(I,IP),I=4,6)                                 LDIS-061
      IF (PIP(8,IJ).NE.0.D0) WRITE (MW,1006) VV(5),PIP(8,IJ),PIP(4,IJ)  LDIS-062
      IF (PIP(9,IJ).NE.0.D0) WRITE (MW,1007) VV(6),PIP(9,IJ),PIP(4,IJ)  LDIS-063
      IF (PIP(17,IJ).NE.0.D0) WRITE (MW,1008) VV(1),PIP(17,IJ),PIP(4,IJ)LDIS-064
      I1=1                                                              LDIS-065
      IF (LO(109)) I1=3                                                 LDIS-066
      IF (IPP(3,IJ).LT.0) GO TO 9                                       LDIS-067
      WRITE (MW,1009) CLEG(I1),CLEG(I1+1),IPP(3,IJ),PIP(10,IJ)          LDIS-068
      IF ((IPP(2,IJ).NE.0).OR.(PIP(11,IJ).NE.0.D0)) WRITE (MW,1010) PIP(LDIS-069
     15,IP),PIP(7,IJ),IPP(2,IJ),PIP(11,IJ),PIP(12,IJ)                   LDIS-070
      GO TO 10                                                          LDIS-071
    9 WRITE (MW,1011) CLEG(I1),CLEG(I1+1),IPP(3,IJ),PIP(10,IJ),PIP(11,IPLDIS-072
     1),PIP(12,IJ)                                                      LDIS-073
   10 IF (LO(109)) GO TO 12                                             LDIS-074
      IF (IPP(4,IJ).LT.0) GO TO 11                                      LDIS-075
      WRITE (MW,1012) IPP(4,IJ),PIP(13,IJ),PIP(14,IJ),PIP(15,IJ)        LDIS-076
      GO TO 13                                                          LDIS-077
   11 IF (PIP(14,IJ).EQ.0.D0) PIP(14,IJ)=2.D0*PIP(13,IJ)                LDIS-078
      WRITE (MW,1013) IPP(4,IJ),PIP(13,IJ),PIP(14,IJ)                   LDIS-079
      GO TO 13                                                          LDIS-080
   12 IF (IPP(4,IJ).GE.0) WRITE (MW,1009) CLEG(5),CLEG(6),IPP(4,IJ),PIP(LDIS-081
     113,IJ)                                                            LDIS-082
      IF (IPP(4,IJ).LT.0) WRITE (MW,1011) CLEG(5),CLEG(6),IPP(4,IJ),PIP(LDIS-083
     113,IJ),PIP(14,IJ),PIP(15,IJ)                                      LDIS-084
   13 CONTINUE                                                          LDIS-085
      RETURN                                                            LDIS-086
   14 WRITE (MW,1014) IJ,(IPP(I,IJ),I=2,5)                              LDIS-087
      GO TO 20                                                          LDIS-088
   15 WRITE (MW,1015) IJ,(IPP(I+2,IJ),VV(2*I),I=1,3)                    LDIS-089
      GO TO 19                                                          LDIS-090
   16 WRITE (MW,1016) IJ,IP,PIP(7,IJ),VV(5)                             LDIS-091
      GO TO 19                                                          LDIS-092
   17 WRITE (MW,1017) IJ,IP,PIP(8,IJ),VV(6)                             LDIS-093
      GO TO 19                                                          LDIS-094
   18 WRITE (MW,1018) IJ,IP,PIP(15,IJ),VV(1)                            LDIS-095
   19 WRITE (MW,1019)                                                   LDIS-096
   20 WRITE (MW,1020)                                                   LDIS-097
      STOP                                                              LDIS-098
 1000 FORMAT (14I5)                                                     LDIS-099
 1001 FORMAT (' DISPERSION COEFFICIENTS READ FOR EACH LEVEL: DATA=',I5) LDIS-100
 1002 FORMAT (7F10.5)                                                   LDIS-101
 1003 FORMAT (' USE OF CENTRE OF MASS ENERGIES.')                       LDIS-102
 1004 FORMAT (' USE OF LABORATORY ENERGIES.')                           LDIS-103
 1005 FORMAT (' IMAGINARY DEPTHS READ FOR',F12.6,' MEV   FERMI ENERGY:',LDIS-104
     1F12.6,' MEV   TRESHOLD ENERGY:',F12.6)                            LDIS-105
 1006 FORMAT (' EXPONENTIAL DECREASE OF REAL SPIN-ORBIT POTENTIAL:',3X,FLDIS-106
     112.6,'*EXP(-',F10.6,'(E-',F10.6,')) MEV')                         LDIS-107
 1007 FORMAT (' LINEAR DECREASE OF IMAGINARY SPIN-ORBIT POTENTIAL:',3X,FLDIS-108
     112.6,'-',F10.6,'*(E-',F10.6,') MEV')                              LDIS-109
 1008 FORMAT (' EXPONENTIAL DECREASE OF HARTREE-FOCK POTENTIAL:',6X,F12.LDIS-110
     16,'*EXP(-',F10.6,'(E-',F10.6,')) MEV')                            LDIS-111
 1009 FORMAT (2A4,':   POWER =',I3,5X,'B =',F10.4)                      LDIS-112
 1010 FORMAT (' FOR ENERGIES ''E'' SUCH THAT |',F12.6,'-''E''| >',F12.6,LDIS-113
     1' MEV'/' CORRECTIONS WITH POWER',I3,' FOR NEGATIVES ENERGIES AND CLDIS-114
     2OEFFICIENT',F12.6,' FOR POSITIVE ENERGIES'/15X,'DAMPING FACTOR',F1LDIS-115
     32.6)                                                              LDIS-116
 1011 FORMAT (2A4,':   POWER =',I3,5X,'SUM OF B1 =',F10.4,5X,'AND B2 =',LDIS-117
     1F10.4,5X,'CONTRIBUTION OF B1 =',F10.5)                            LDIS-118
 1012 FORMAT (' SURFACE:   POWER =',I3,5X,'B =',F10.4,5X,' DECREASING RALDIS-119
     1TE =',F10.6,5X,' NON LOCALITY PARAMETER =',F10.4)                 LDIS-120
 1013 FORMAT (' SURFACE:   POWER =',I3,5X,'DIFFERENCE OF B1 =',F10.4,5X,LDIS-121
     1'AND B2 =',F10.4)                                                 LDIS-122
 1014 FORMAT (' FOR DISPERSION RELATIONS OF POTENTIAL',I3,' ALL THE FOLLLDIS-123
     1OWING INTEGERS MUST BE EVEN'/10X,'N2 =',I4,10X,'NV =',I4,10X,'NS =LDIS-124
     2',I4,10X,'NL =',I4)                                               LDIS-125
 1015 FORMAT (' THE POTENTIAL',I4,' CANNOT BE USED FOR DISPERSION RELATILDIS-126
     1ONS BECAUSE AN IMAGINARY STRENGTH IS 0:'/10X,'NV =',I3,'  WV =',1PLDIS-127
     2,D13.6,10X,'NS =',I3,'  WS =',D13.6,10X,'NL =',I3,'  LS =',D13.6) LDIS-128
 1016 FORMAT (' THE POTENTIAL',I4,' CANNOT BE USED WITH VARIATION OF THELDIS-129
     1 REAL SPIN-ORBIT OF WHICH THE STRENGTH IS 0':/10X,'PIP(7,',I2,') =LDIS-130
     2',1P,D13.6,'  VLS =',D13.6)                                       LDIS-131
 1017 FORMAT (' THE POTENTIAL',I4,' CANNOT BE USED WITH VARIATION OF THELDIS-132
     1 IMAGINARY SPIN-ORBIT OF WHICH THE STRENGTH IS 0':/10X,'PIP(8,',I2LDIS-133
     2,') =',1P,D13.6,'  WLS =',D13.6)                                  LDIS-134
 1018 FORMAT (' THE POTENTIAL',I4,' CANNOT BE USED WITH VARIATION OF THELDIS-135
     1 HARTREE-FOCK POTENTIAL OF WHICH THE STRENGTH IS 0':/10X,'PIP(15,'LDIS-136
     2,I2,') =',1P,D13.6,'  V =',D13.6)                                 LDIS-137
 1019 FORMAT (' GIVE VALUES FOR ANOTHER ENERGY'/)                       LDIS-138
 1020 FORMAT (//' IN LDIS  ...  STOP  ...')                             LDIS-139
      END                                                               LDIS-140
