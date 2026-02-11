C 28/02/07                                                      ECIS06  LECD-000
      SUBROUTINE LECD(WV,NP,KFIT,NESSAI,YY,NJY,NMAX,LO,MF,MFM,FM,DONN,DWLECD-001
     1,NW,NMX)                                                          LECD-002
C INPUT OF EXPERIMENTAL DATA AND SEARCH INFORMATIONS.                   LECD-003
C INPUT:     WV(J,*): ENERGY IN THE CENTRE OF MASS IN MEV FOR J=3.      LECD-004
C            NP:      INDICATIONS FOR PLOTS OF CROSS-SECTION            LECD-005
C            KFIT:    NUMBER OF FITS STORED IN A SEARCH.                LECD-006
C            NESSAI:  MAXIMUM NUMBER OF EVALUATIONS IN THE SEARCH.      LECD-007
C            YY(3):   PARAMETERS FOR THE SEARCH.                        LECD-008
C            NJY:     MAXIMUM INDEX OF NON STANDARD OBSERVABLE.         LECD-009
C            NMAX:    MAXIMUM NUMBER OF DATA FOR THE WORKING SPACE.     LECD-010
C            LO(I):   LOGICAL CONTROLS:                                 LECD-011
C               LO(32) =.TRUE. AUTOMATIC SEARCH ON SOME PARAMETERS.     LECD-012
C               LO(81) =.TRUE. HAUSER-FESHBACH CORRECTIONS.             LECD-013
C               LO(72) =.TRUE. NO OUTPUT OF EXPERIMENTAL DATA WHEN THEY LECD-014
C                              ARE READ.                                LECD-015
C               LO(85) =.TRUE. FISSION TRANSMISSION COEFFICIENTS.       LECD-016
C               LO(86) =.TRUE. GAMMA EMISSION IN COMPOUND NUCLEUS.      LECD-017
C OUTPUT:    NMX:     NUMBER OF INDEXES.                                LECD-018
C            MF:      CONTINUATION OF THE SECOND PART OF MF AS DESCRIBEDLECD-019
C                     IN SUBROUTINE DEPH.                               LECD-020
C            MFM,FM:  IN EQUIVALENCE BY CALL:                           LECD-021
C                     MFM(1,*) CHANNEL,                                 LECD-022
C                     MFM(2,*), MFM(3,*) BEGINNING AND END OF DATA,     LECD-023
C                     MFM(4,*) INDICATION CENTRE OF MASS OR LABORATORY  LECD-024
C                              SYSTEM BY 0 OR 1,                        LECD-025
C                     FM(3,*) WEIGHT,                                   LECD-026
C                     FM(4,*) AND FM(5,*) NORM AND ITS ERROR,           LECD-027
C                     FM(6,*) PLACE FOR CALCULATED NORMALISATION,       LECD-028
C                     FM(7,*) PLACE FOR CALCULATED CHI2.                LECD-029
C            DONN:    EXPERIMENTAL DATA: ANGLE, VALUE, ERROR,           LECD-030
C                     ANGULAR WIDTH AND PLACE FOR CORRECTED ERROR.      LECD-031
C            DW:      ACCURACY OF PARAMETERS IN SEARCH.                 LECD-032
C            NW:      INDEXES OF PARAMETERS IN SEARCH.                  LECD-033
C            NJY:     MAXIMUM INDEX OF NON STANDARD OBSERVABLE.         LECD-034
C                                                                       LECD-035
C      IF THE NUMBER OF DATA FOR AN OBSERVABLE IS 0, IT IS SUMMED       LECD-036
C      WITH THE NEXT OBSERVABLE WHICH MUST BE OF THE SAME KIND.         LECD-037
C                                                                       LECD-038
C FOR THE COMMON  /INTEG/ SEE CALC.                                     LECD-039
C                                                                       LECD-040
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /INTEG/:                     LECD-041
C  NCOLS:     NUMBER OF CHANNELS WITH ANGULAR DISTRIBUTIONS.            LECD-042
C  NCOLT:     NUMBER OF CHANNELS INCLUDING UNCOUPLED STATES.            LECD-043
C  JTH:       MAXIMUM NUMBER OF ANGLES FOR A PLOT.                      LECD-044
C  NCOLR:     NUMBER OF EXPERIMENTAL ANGULAR DISTRIBUTIONS.             LECD-045
C  NREC:      NUMBER OF VARIABLES IN SEARCH.                            LECD-046
C  NTOT:      NUMBER OF EXPERIMENTAL DATA.                              LECD-047
C   DEFINED:  NTOT.                                                     LECD-048
C   USED:                                                               LECD-049
C                                                                       LECD-050
C***********************************************************************LECD-051
      IMPLICIT REAL*8 (A-H,O-Z)                                         LECD-052
      LOGICAL LO(150),LT,LX                                             LECD-053
      DIMENSION WV(22,*),MF(10,*),MFM(14,*),FM(7,*),DONN(6,*),NP(2),DW(*LECD-054
     1),NW(*),YY(3)                                                     LECD-055
      COMMON /INOUT/ MR,MW,MS                                           LECD-056
      COMMON /INTEG/ IDMT,NPLACE,NCOLL,NJMAX,ITERM,JDM,JIT,KMIN,KMAX,NCOLECD-057
     1LS,NCOLT,NBET,LMX,LMAX1,NLT,ISM,NJC,JTX,KCC,MS1,MS2,KBA,KAB,KBC,JTLECD-058
     2H,NCOLR,NREC,NTOT,LMAX2,KE,ITEMM,KXT,LMAX3,NRZ,NTZ,IPM,IPK,MCM(2),LECD-059
     3NCT(6)                                                            LECD-060
      LT=.TRUE.                                                         LECD-061
      NTOT=0                                                            LECD-062
      NCLR=0                                                            LECD-063
      NMX=NREC+1                                                        LECD-064
      DO 5 IV=1,NCOLR                                                   LECD-065
C INPUT OF THE DEFINITION OF THE ANGULAR DISTRIBUTION.                  LECD-066
      LX=.FALSE.                                                        LECD-067
      READ (MR,1000) LX,MFM(4,IV),NT,MFM(1,IV),MF(2,IV),(FM(J,IV),J=3,5)LECD-068
      IF (MFM(1,IV).GT.NCOLS) GO TO 8                                   LECD-069
      IF (WV(3,MFM(1,IV)).LE.0.D0) GO TO 9                              LECD-070
      IF (MF(2,IV).GT.19) GO TO 10                                      LECD-071
      IF (FM(3,IV).EQ.0.D0) FM(3,IV)=1.D0                               LECD-072
      IF (FM(4,IV).EQ.0.D0) FM(4,IV)=1.D0                               LECD-073
      FM(6,IV)=1.D0                                                     LECD-074
      FM(7,IV)=0.D0                                                     LECD-075
      MF(1,IV)=MFM(1,IV)                                                LECD-076
      JTH=MAX0(JTH,NT)                                                  LECD-077
      NJY=MAX0(NJY,-MF(2,IV))                                           LECD-078
      MF(5,IV)=1                                                        LECD-079
      MF2=MIN0(2,MFM(1,IV))                                             LECD-080
      IF (MF(2,IV).EQ.0.OR.MF(2,IV).EQ.1) MF(5,IV)=NP(MF2)              LECD-081
      IF (.NOT.LO(72).AND.NT.EQ.0) WRITE (MW,1001) IV,NT,MFM(1,IV),MF(2,LECD-082
     1IV)                                                               LECD-083
      IF (IV.EQ.1) GO TO 1                                              LECD-084
      IF ((.NOT.LT).AND.(MF(2,IV).NE.MF(2,IV-1))) GO TO 11              LECD-085
    1 LT=NT.NE.0                                                        LECD-086
      NS=NTOT+1                                                         LECD-087
      NTOT=NTOT+NT                                                      LECD-088
      MFM(2,IV)=NS                                                      LECD-089
      MFM(3,IV)=NTOT                                                    LECD-090
      IF (.NOT.LT) GO TO 5                                              LECD-091
      NCLR=NCLR+1                                                       LECD-092
      IF (.NOT.LO(72)) WRITE (MW,1002) IV,NT,MFM(1,IV),MF(2,IV),(FM(I,IVLECD-093
     1),I=3,5)                                                          LECD-094
      IF (.NOT.LO(72).AND.(MF(2,IV).EQ.19)) WRITE (MW,1003)             LECD-095
      LX=LX.AND.(MF(2,IV).EQ.0.OR.MF(2,IV).EQ.1.OR.MF(2,IV).EQ.19)      LECD-096
      IF (.NOT.LO(72).AND.LX) WRITE (MW,1004)                           LECD-097
      IF (.NOT.LO(72).AND.MFM(4,IV).EQ.1) WRITE (MW,1005)               LECD-098
      IF (3*NTOT+3.GE.NMAX) CALL MEMO('LECD',NMAX,3*NTOT+3)             LECD-099
C INPUT OF THE ANGULAR DISTRIBUTION DATA.                               LECD-100
      DO 2 I=NS,NTOT                                                    LECD-101
      READ (MR,1006) (DONN(J,I),J=1,5)                                  LECD-102
      IF (LX) DONN(3,I)=DONN(2,I)*DONN(3,I)*.01D0                       LECD-103
      DONN(6,I)=DONN(3,I)                                               LECD-104
      IF (.NOT.LO(72)) WRITE (MW,1007) I,(DONN(J,I),J=1,5)              LECD-105
    2 CONTINUE                                                          LECD-106
      IF (MF(2,IV).LT.19) GO TO 4                                       LECD-107
C CHECK OF TOTAL CROSS-SECTION DATA.                                    LECD-108
      DO 3 I=NS,NTOT                                                    LECD-109
      J=IDINT(DONN(1,I)+.1D0)                                           LECD-110
      IF (J.LT.0.OR.J.GT.NCOLT+2) GO TO 12                              LECD-111
      IF (J.GT.0.AND.J.LE.NCOLT.AND.WV(3,J).LE.0.D0) GO TO 13           LECD-112
      IF (LO(81).AND.J.GT.NCOLS) GO TO 12                               LECD-113
      IF ((.NOT.LO(85)).AND.J.EQ.NCOLT+1) GO TO 12                      LECD-114
      IF ((.NOT.LO(86)).AND.J.EQ.NCOLT+2) GO TO 12                      LECD-115
    3 CONTINUE                                                          LECD-116
    4 IF (FM(5,IV).EQ.0.D0) GO TO 5                                     LECD-117
      NTOT=NTOT+1                                                       LECD-118
C NORMALISATION AS A DATA FOR THE SEARCH.                               LECD-119
      DONN(2,NTOT)=FM(4,IV)                                             LECD-120
      DONN(3,NTOT)=FM(5,IV)                                             LECD-121
      IF (.NOT.LO(72)) WRITE (MW,1008) NTOT,DONN(2,NTOT),DONN(3,NTOT)   LECD-122
    5 CONTINUE                                                          LECD-123
      IF ((.NOT.LO(72)).AND.(NCOLR.NE.NCLR)) WRITE (MW,1009) NCLR,NCOLR LECD-124
      IF ((.NOT.LO(32)).OR.(NREC.EQ.0)) GO TO 7                         LECD-125
      WRITE (MW,1010) NREC,KFIT,NESSAI,YY(1),YY(2)                      LECD-126
      NISE=12*NTOT+2*NREC                                               LECD-127
      IF (NISE+NMX.GE.2*NMAX) CALL MEMO('LECD',NMAX,(NISE+NMX)/2)       LECD-128
      READ (MR,1006) (DW(6*NTOT+I),I=1,NREC)                            LECD-129
      READ (MR,1011) (NW(NISE+I),I=1,NREC)                              LECD-130
      WRITE (MW,1012) (I,NW(NISE+I),DW(6*NTOT+I),I=1,NREC)              LECD-131
      DO 6 I=1,NREC                                                     LECD-132
C A NEGATIVE VALUE -K INSTEAD OF INDEXES OF PARAMETERS MEANS THAT K     LECD-133
C PARAMETERS WILL BE KEPT PROPORTIONAL - INPUT OF THEIR INDEXES.        LECD-134
      K=-NW(NISE+I)                                                     LECD-135
      IF (K.LT.0) GO TO 6                                               LECD-136
      NW(NISE+NMX)=K                                                    LECD-137
      IF (NISE+NMX+K.GE.2*NMAX) CALL MEMO('LECD',NMAX,(NISE+NMX+K)/2)   LECD-138
      READ (MR,1011) (NW(NISE+NMX+J),J=1,K)                             LECD-139
      WRITE (MW,1013) I,(NW(NISE+NMX+J),J=1,K)                          LECD-140
      NW(NISE+I)=-NMX                                                   LECD-141
      NMX=NMX+K+1                                                       LECD-142
    6 CONTINUE                                                          LECD-143
    7 IF (LT) RETURN                                                    LECD-144
      WRITE (MW,1014)                                                   LECD-145
      GO TO 14                                                          LECD-146
    8 WRITE (MW,1015) MFM(1,IV)                                         LECD-147
      GO TO 14                                                          LECD-148
    9 WRITE (MW,1016) MFM(1,IV)                                         LECD-149
      GO TO 14                                                          LECD-150
   10 WRITE (MW,1017)                                                   LECD-151
      GO TO 14                                                          LECD-152
   11 WRITE (MW,1018) MF(2,IV),MF(2,IV-1)                               LECD-153
      GO TO 14                                                          LECD-154
   12 WRITE (MW,1019) DONN(1,I),I                                       LECD-155
      GO TO 14                                                          LECD-156
   13 WRITE (MW,1020) DONN(1,I),I                                       LECD-157
   14 WRITE (MW,1021)                                                   LECD-158
      STOP                                                              LECD-159
 1000 FORMAT (L1,I1,I3,2I5,5X,3F10.5)                                   LECD-160
 1001 FORMAT (/'  ANG. DISTR.',I3,I9,' DATA   LEVEL =',I3,5X,'KIND =',I3LECD-161
     1,' UNRESOLVED WITH THE FOLLOWING ONE.'/)                          LECD-162
 1002 FORMAT (/'  ANG. DISTR.',I3,I9,' DATA   LEVEL =',I3,5X,'KIND =',I3LECD-163
     1//5X,'WEIGHT',1P,D12.4,5X,'NORM',D12.4,'  WITH ERROR',D12.4//23X,'LECD-164
     2ANGLE',14X,'VALUE',15X,'ERROR',12X,'ANG. WIDTH',10X,'ANG. ERROR'/)LECD-165
 1003 FORMAT (' THESE DATA ARE TOTAL CROSS-SECTIONS.')                  LECD-166
 1004 FORMAT (' INPUT OF PERCENTAGE ERRORS.')                           LECD-167
 1005 FORMAT (' ANGLES IN THE LABORATORY SYSTEM.')                      LECD-168
 1006 FORMAT (7F10.5)                                                   LECD-169
 1007 FORMAT (5X,I5,1P,5D20.6)                                          LECD-170
 1008 FORMAT (5X,I5,' DATA WHICH IS A NORMALISATION',1P,D20.6,' WITH ERRLECD-171
     1OR',D20.6)                                                        LECD-172
 1009 FORMAT (/5X,I5,' DIFFERENT ANGULAR DISTRIBUTIONS INSTEAD OF',I5/) LECD-173
 1010 FORMAT (//5X,I5,' PARAMETERS IN SEARCH'/5X,I5,' RESULTS STORED'/5XLECD-174
     1,I5,' RUNS    STARTING SCALE',F10.2/12X,'MULTIPLICATION FACTOR',F1LECD-175
     20.2/)                                                             LECD-176
 1011 FORMAT (14I5)                                                     LECD-177
 1012 FORMAT (5X,I5,5X,I5,5X,F15.8)                                     LECD-178
 1013 FORMAT (5X,I5,' VARIABLE DEFINED AS',20I5/(18X,20I5))             LECD-179
 1014 FORMAT (//'  THE LAST ANGULAR DISTRIBUTION INCLUDE NO DATA.')     LECD-180
 1015 FORMAT (//' LEVEL NUMBER',I4,' TOO LARGE.')                       LECD-181
 1016 FORMAT (//' LEVEL',I4,' CLOSED.')                                 LECD-182
 1017 FORMAT (//' THERE ARE ONLY 20 KIND OF PROGRAMMED OBSERVABLES'//'  LECD-183
     1FOR OTHER KINDS, INTRODUCE A NEGATIVE NUMBER AND, LATER, THEIR DESLECD-184
     2CRIPTION.')                                                       LECD-185
 1018 FORMAT (//' THE KIND =',I3,' OF THIS OBSERVABLE IS NOT THE SAME ASLECD-186
     1 THE KIND =',I3,' OF THE PREVIOUS ONE WHICH WAS EMPTY.')          LECD-187
 1019 FORMAT (' ANGULAR DATA',F10.5,' FOR',I4,' DATA NOT ALLOWED FOR TOTLECD-188
     1AL CROSS-SECTION.')                                               LECD-189
 1020 FORMAT (' ANGULAR DATA',F10.5,' FOR',I4,' DATA NOT ALLOWED FOR CLOLECD-190
     1SED CHANNEL.')                                                    LECD-191
 1021 FORMAT (//' IN LECD  ...  STOP  ...')                             LECD-192
      END                                                               LECD-193
