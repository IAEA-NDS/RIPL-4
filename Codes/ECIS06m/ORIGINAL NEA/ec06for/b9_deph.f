C 31/08/06                                                      ECIS06  DEPH-000
      SUBROUTINE DEPH(NCOL,WV,IPI,NGR,NPR,NMAX,LO,MF,KTGR,NIV,NJX,NJY)  DEPH-001
C HELICITY QUANTUM NUMBERS AND CHOICE OF OBSERVABLES FOR THE OUTPUT.    DEPH-002
C INPUT:     NCOL:    NUMBER OF COUPLED NUCLEAR LEVELS FOR DIMENSION OF DEPH-003
C                     ARRAY NIV, EQUAL TO NCOLL.                        DEPH-004
C            WV(I,*): MASS OF THE INCIDENT PARTICLE FOR J=1,            DEPH-005
C                     ENERGY IN THE CENTRE OF MASS SYSTEM FOR J=3.      DEPH-006
C            IPI(I,*):PARITY (0 FOR + AND 1 FOR -) FOR J=1,             DEPH-007
C                     MULTIPLICITY OF INCIDENT PARTICLE FOR J=2,        DEPH-008
C                     MULTIPLICITY OF THE TARGET FOR J=3,               DEPH-009
C                     PRODUCT OF CHARGES FOR J=4.                       DEPH-010
C            NGR:     INDICATION FOR PLOTS OF CROSS-SECTION: NUMBER     DEPH-011
C                     OF POWERS OF 10 BY 100 POINTS.                    DEPH-012
C            NPR:     INDICATIONS FOR PLOTS OF POLARISATIONS:           DEPH-013
C                     1 FIRST ONE,2 SECOND,3 FIRST AND SECOND ..ETC     DEPH-014
C                     ON A BINARY BASIS. THE FIRST VALUE OF NGR AND NPR DEPH-015
C                     IS FOR THE GROUND STATE AND THE SECOND ONE FOR THEDEPH-016
C                     EXCITED STATES                                    DEPH-017
C            NMAX:    MAXIMUM NUMBER OF AMPLITUDES.                     DEPH-018
C            LO(I):   LOGICAL CONTROLS:                                 DEPH-019
C               LO(44) =.TRUE. COULOMB CORRECTIONS.                     DEPH-020
C               LO(46) =.TRUE. RESTRICTED COULOMB CORRECTIONS.          DEPH-021
C               LO(66) =.TRUE. NO CALCULATION AT EQUIDISTANT ANGLES.    DEPH-022
C               LO(81) =.TRUE. HAUSER-FESHBACH CORRECTIONS.             DEPH-023
C               LO(92) =.TRUE. NON STANDARD OBSERVABLES AT EQUIDISTANT  DEPH-024
C                              ANGLES.                                  DEPH-025
C               LO(100)=.TRUE. DIRAC EQUATION.                          DEPH-026
C OUTPUT:    IPI(I,*):BEGINNING AND END FOR EACH LEVEL IN THE TWO PARTS DEPH-027
C                     OF THE TABLE MF FOR J=5 TO 8.                     DEPH-028
C            MF:      TABLE OF QUANTUM NUMBERS AND OBSERVABLES.         DEPH-029
C            KTGR:    LENGTH OF THE TABLE MF.                           DEPH-030
C            NIV:     ADDRESS OF COULOMB INTEGRALS IN NIV(*,*,3).       DEPH-031
C            NJX:     NUMBER OF SETS OF COULOMB FUNCTIONS AND INTEGRALS.DEPH-032
C            NJY:     MAXIMUM LABEL OF NON STANDARD OBSERVABLE.         DEPH-033
C                                                                       DEPH-034
C FOR THE COMMON  /INTEG/ SEE CALC.                                     DEPH-035
C                                                                       DEPH-036
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /INTEG/:                     DEPH-037
C  NCOLS:     NUMBER OF CHANNELS WITH ANGULAR DISTRIBUTIONS.            DEPH-038
C  NJC:       MAXIMUM NUMBER OF OBSERVABLES AT EQUIDISTANT ANGLES.      DEPH-039
C  JTX:       MAXIMUM NUMBER OF CALCULATED VALUES FOR A PLOT.           DEPH-040
C  KCC:       NUMBER OF INDEPENDENT AMPLITUDES WITH UNCOUPLED STATES.   DEPH-041
C  MS1:       LARGEST PARTICLE MULTIPLICITY.                            DEPH-042
C  MS2:       LARGEST TARGET MULTIPLICITY.                              DEPH-043
C  KBA:       NUMBER OF INDEPENDENT AMPLITUDES WITHOUT UNCOUPLED STATES.DEPH-044
C  KAB:       MAXIMUM NUMBER OF EQUATIONS.                              DEPH-045
C  KBC:       MAXIMUM NUMBER OF SOLUTIONS.                              DEPH-046
C  NCT(1):    NUMBER OF EQUATIONS FOR POSITIVE PARITY.                  DEPH-047
C  NCT(2):    NUMBER OF EQUATIONS FOR NEGATIVE PARITY.                  DEPH-048
C  NCT(3):    NUMBER OF SOLUTIONS FOR POSITIVE PARITY.                  DEPH-049
C  NCT(4):    NUMBER OF SOLUTIONS FOR NEGATIVE PARITY.                  DEPH-050
C  NCT(5):    NUMBER OF COMPOUND NUCLEUS COEFFICIENTS, POSITIVE PARITY. DEPH-051
C  NCT(6):    NUMBER OF COMPOUND NUCLEUS COEFFICIENTS, NEGATIVE PARITY. DEPH-052
C   DEFINED:  NJC,JTX,KCC,MS1,MS2,KBA,KAB,KBC,NCT.                      DEPH-053
C   USED:     NCOLS.                                                    DEPH-054
C                                                                       DEPH-055
C   **** TABLE MF **** FIRST PART                                       DEPH-056
C FOR EACH INDEPENDENT AMPLITUDE (WITH RESPECT TO PARITY ONLY):         DEPH-057
C  MF(1,*)  HELICITY OF THE OUTGOING PARTICLE,                          DEPH-058
C  MF(2,*)  HELICITY OF THE RESIDUAL TARGET,                            DEPH-059
C  MF(3,*)  HELICITY OF THE INCOMING PARTICLE,                          DEPH-060
C  MF(4,*)  HELICITY OF THE TARGET,                                     DEPH-061
C    (THESE HELICITIES ARE NUMBERED FROM THE LOWEST VALUE)              DEPH-062
C  MF(5,*)  TWICE THE MAGNETIC QUANTUM NUMBER OF THE ROTATION MATRIX    DEPH-063
C    ELEMENTS RELATED TO THE INITIAL STATE,                             DEPH-064
C  MF(6,*)  TWICE THE MAGNETIC QUANTUM NUMBER OF THE ROTATION MATRIX    DEPH-065
C    ELEMENTS FOR THE FINAL STATE,                                      DEPH-066
C    (HOWEVER, WITH ABSOLUTE VALUE 99999 IT INDICATES THAT THE LAST     DEPH-067
C    COMPUTED MATRIX ELEMENTS CAN BE USED WITH THE SIGN OF MF(10,*).    DEPH-068
C    A RE-ORDERING OF THE HELICITIES HAS BEEN DONE TO BE ABLE TO DO SO) DEPH-069
C  MF(7,*)  DIRECT ADDRESS OF THE AMPLITUDE IN AN ONE-ROW MATRIX,       DEPH-070
C  MF(8,*)  DIRECT ADDRESS OF PARITY RELATED AMPLITUDE OR 0 IF THERE IS DEPH-071
C    NONE,                                                              DEPH-072
C  MF(9,*)  RELATIVE SIGN BETWEEN THE AMPLITUDES MF(7,*) AND MF(8,*),   DEPH-073
C  MF(10,*) RELATIVE SIGN FOR ROTATION MATRIX ELEMENTS IF MF(6,*)=99999.DEPH-074
C***********************************************************************DEPH-075
      IMPLICIT REAL*8 (A-H,O-Z)                                         DEPH-076
      LOGICAL LO(150)                                                   DEPH-077
      DIMENSION WV(22,*),IPI(11,*),MF(10,*),NGR(2),NPR(2),NIV(NCOL,NCOL,DEPH-078
     13)                                                                DEPH-079
      COMMON /INOUT/ MR,MW,MS                                           DEPH-080
      COMMON /INTEG/ IDMT,NPLACE,NCOLL,NJMAX,ITERM,JDM,JIT,KMIN,KMAX,NCODEPH-081
     1LS,NCOLT,NBET,LMX,LMAX1,NLT,ISM,NJC,JTX,KCC,MS1,MS2,KBA,KAB,KBC,JTDEPH-082
     2H,NCOLR,NREC,NTOT,LMAX2,KE,ITEMM,KXT,LMAX3,NRZ,NTZ,IPM,IPK,MCM(2),DEPH-083
     3NCT(6)                                                            DEPH-084
      NJC=2                                                             DEPH-085
      NJY=0                                                             DEPH-086
      JTX=0                                                             DEPH-087
      NCT(5)=0                                                          DEPH-088
      NCT(6)=0                                                          DEPH-089
      KCC=0                                                             DEPH-090
      MS1=0                                                             DEPH-091
      MS2=0                                                             DEPH-092
      NI=IPI(2,1)                                                       DEPH-093
      MI=IPI(3,1)                                                       DEPH-094
C LOOP ON THE NUCLEAR LEVELS.                                           DEPH-095
      DO 13 I=1,NCOLS                                                   DEPH-096
      IP=MOD(IPI(1,I)+IPI(1,1),2)                                       DEPH-097
      KI=KCC+1                                                          DEPH-098
      IPI(6,I)=KI                                                       DEPH-099
      NJ=IPI(2,I)                                                       DEPH-100
      MJ=IPI(3,I)                                                       DEPH-101
      MS1=MAX0(MS1,NJ)                                                  DEPH-102
      MS2=MAX0(MS2,MJ)                                                  DEPH-103
      DO 10 I1=1,NJ                                                     DEPH-104
      J1=NJ+1-I1                                                        DEPH-105
      DO 9 I2=1,MJ                                                      DEPH-106
      J2=MJ+1-I2                                                        DEPH-107
      DO 8 I3=1,NI                                                      DEPH-108
      J3=NI+1-I3                                                        DEPH-109
      DO 7 I4=1,MI                                                      DEPH-110
      J4=MI+1-I4                                                        DEPH-111
      IF (KCC.LT.KI) GO TO 3                                            DEPH-112
C SEARCH FOR PARITY CONJUGATE AMONG THE AMPLITUDES ALREADY OBTAINED.    DEPH-113
      DO 1 J=KI,KCC                                                     DEPH-114
      IF ((MF(1,J).EQ.J1).AND.(MF(2,J).EQ.J2).AND.(MF(3,J).EQ.J3).AND.(MDEPH-115
     1F(4,J).EQ.J4)) GO TO 2                                            DEPH-116
    1 CONTINUE                                                          DEPH-117
      GO TO 3                                                           DEPH-118
    2 MF(8,J)=I1+NJ*(I2-1+MJ*(I4-1+MI*(I3-1)))                          DEPH-119
      MF(9,J)=1-2*MOD(IP+J1+I2+J3+I4,2)                                 DEPH-120
      GO TO 7                                                           DEPH-121
C NEW AMPLITUDE.                                                        DEPH-122
    3 KCC=KCC+1                                                         DEPH-123
      IF (5*KCC.GE.NMAX) CALL MEMO('DEPH',NMAX,5*KCC)                   DEPH-124
      MF(1,KCC)=I1                                                      DEPH-125
      MF(2,KCC)=I2                                                      DEPH-126
      MF(3,KCC)=I3                                                      DEPH-127
      MF(4,KCC)=I4                                                      DEPH-128
      MF(5,KCC)=2*(MF(3,KCC)-MF(4,KCC))-IPI(2,1)+IPI(3,1)               DEPH-129
      MF(6,KCC)=2*(MF(1,KCC)-MF(2,KCC))-IPI(2,I)+IPI(3,I)               DEPH-130
      MF(7,KCC)=I1+NJ*(I2-1+MJ*(I4-1+MI*(I3-1)))                        DEPH-131
      MF(8,KCC)=0                                                       DEPH-132
      MF(9,KCC)=0                                                       DEPH-133
      MF(10,KCC)=0                                                      DEPH-134
      IF (KCC.LE.KI) GO TO 7                                            DEPH-135
C SEARCH FOR RELATED ROTATION MATRIX ELEMENTS.                          DEPH-136
      KK=KCC-1                                                          DEPH-137
      M1=IABS(MF(5,KCC)+MF(6,KCC))/2                                    DEPH-138
      M2=IABS(MF(5,KCC)-MF(6,KCC))/2                                    DEPH-139
      DO 6 J=KI,KK                                                      DEPH-140
      IF (MF(6,J).EQ.99999) GO TO 6                                     DEPH-141
      M3=IABS(MF(5,J)+MF(6,J))/2                                        DEPH-142
      M4=IABS(MF(5,J)-MF(6,J))/2                                        DEPH-143
      IF (M1.NE.M3.OR.M2.NE.M4) GO TO 6                                 DEPH-144
      MF(10,KCC)=1-MOD(IABS(MF(5,J)-MF(5,KCC)+MF(6,KCC)-MF(6,J))/2,4)   DEPH-145
      MF(6,KCC)=99999                                                   DEPH-146
      JA=J+1                                                            DEPH-147
      IF (JA.EQ.KCC) GO TO 7                                            DEPH-148
C PERMUTATION OF AMPLITUDES.                                            DEPH-149
      DO 5 M1=1,10                                                      DEPH-150
      M2=MF(M1,KCC)                                                     DEPH-151
      DO 4 M3=JA,KCC                                                    DEPH-152
      M4=KCC+JA-M3                                                      DEPH-153
    4 MF(M1,M4)=MF(M1,M4-1)                                             DEPH-154
    5 MF(M1,JA)=M2                                                      DEPH-155
      GO TO 7                                                           DEPH-156
    6 CONTINUE                                                          DEPH-157
    7 CONTINUE                                                          DEPH-158
    8 CONTINUE                                                          DEPH-159
    9 CONTINUE                                                          DEPH-160
   10 CONTINUE                                                          DEPH-161
      IPI(7,I)=KCC                                                      DEPH-162
C COMPUTATION OF THE NUMBER OF COUPLED EQUATIONS.                       DEPH-163
      L=NJ*MJ                                                           DEPH-164
      KK1=L/2                                                           DEPH-165
      KK2=L-KK1                                                         DEPH-166
      IF (KK1.EQ.KK2) GO TO 11                                          DEPH-167
      N=MJ/2+NJ/2+IPI(1,I)                                              DEPH-168
      IF (2*(N/2).EQ.N) GO TO 11                                        DEPH-169
      KK1=KK1+1                                                         DEPH-170
      KK2=KK1-1                                                         DEPH-171
   11 NCT(5)=NCT(5)+KK2                                                 DEPH-172
      NCT(6)=NCT(6)+KK1                                                 DEPH-173
      IF (I.NE.1) GO TO 12                                              DEPH-174
      NCT(3)=NCT(5)                                                     DEPH-175
      NCT(4)=NCT(6)                                                     DEPH-176
   12 IF (I.NE.NCOLL) GO TO 13                                          DEPH-177
      NCT(1)=NCT(5)                                                     DEPH-178
      NCT(2)=NCT(6)                                                     DEPH-179
      KAB=MAX0(NCT(1),NCT(2))                                           DEPH-180
      KBC=MAX0(NCT(3),NCT(4))                                           DEPH-181
      IF (NCT(4).EQ.0) KAB=NCT(1)                                       DEPH-182
      IF (NCT(3).EQ.0) KAB=NCT(2)                                       DEPH-183
      KBA=KCC                                                           DEPH-184
   13 CONTINUE                                                          DEPH-185
      IF (KBA.EQ.KCC) WRITE (MW,1000) KBC,KAB,KBA                       DEPH-186
      IF (KBA.NE.KCC) WRITE (MW,1001) KBC,KAB,KBA,KCC                   DEPH-187
      KTGR=KCC                                                          DEPH-188
      IPI(9,NCOLS)=IPI(7,NCOLS)                                         DEPH-189
      IF (LO(66)) GO TO 25                                              DEPH-190
C***********************************************************************DEPH-191
C   **** TABLE MF-FM **** SECOND PART                                   DEPH-192
C  THIS PART OF THE TABLE IS PROLONGATED IN LECD FOR EACH ANGULAR       DEPH-193
C  DISTRIBUTION AND WILL BE UPDATED BY OBSE. MF(1,*), MF(2,*) AND       DEPH-194
C  MF(5,*) ONLY ARE DEFINED HERE.                                       DEPH-195
C                                                                       DEPH-196
C  - FOR EACH OBSERVABLE:                                               DEPH-197
C  MF(1,*)             LEVEL.                                           DEPH-198
C  MF(2,*)             KIND: SEE BELOW **** STANDARD DESCRIPTIONS ****  DEPH-199
C  MF(3,I), MF(4,I)    BEGINNING AND END OF THE DESCRIPTION IN THE ARRAYDEPH-200
C                      AM DEFINED IN SUBROUTINE OBSE,                   DEPH-201
C  MF(5,*)             INDICATION FOR PLOTS,                            DEPH-202
C  MF(6,*) TO M(10,*)  CONTAINS THE LEGEND (ADDED IN SUBROUTINE OBSE).  DEPH-203
C                                                                       DEPH-204
C  ****STANDARD OPTIONS****                                             DEPH-205
C FOR ALL STATES : CROSS SECTION.                                       DEPH-206
C FOR GROUND STATE WITH CHARGED PARTICLES: CROSS SECTION DIVIDED BY     DEPH-207
C       RUTHERFORD'S CROSS-SECTION.                                     DEPH-208
C FOR SPIN 1/2 - GROUND STATE  POLARISATION.                            DEPH-209
C FOR SPIN 1/2 EXCITED STATES  VECTOR ANALYSING POWER, POLARISATION     DEPH-210
C       AND SPIN-FLIP.                                                  DEPH-211
C FOR SPIN LARGER THAN 1/2  GROUND STATE IT11, T20, T21 AND T22.        DEPH-212
C FOR EXCITED STATES, IT11, VECTOR POLARISATION, T20, T21 AND T22.      DEPH-213
C                                                                       DEPH-214
C  *** NON-STANDARD OPTIONS ***                                         DEPH-215
C  THE FIRST ONE MUST BE THE CROSS-SECTION.                             DEPH-216
C ANY OBSERVABLE NOT DESCRIBED BELOW IS IDENTIFIED BY A NEGATIVE NUMBER DEPH-217
C  OF WHICH THE POSITIVE VALUE MUST BE FOUND BY THE SUBROUTINE OBSE     DEPH-218
C  FOLLOWED BY THE DESCRIPTION OF THE OBSERVABLE.                       DEPH-219
C                                                                       DEPH-220
C **** STANDARD DESCRIPTIONS ****                                       DEPH-221
C  0  CROSS SECTION.                                                    DEPH-222
C  1  CROSS SECTION DIVIDED BY RUTHERFORD'S CROSS SECTION.              DEPH-223
C  2  VECTOR ANALYSING POWER.                                           DEPH-224
C  3  VECTOR POLARISATION.                                              DEPH-225
C  4  T20.                                                              DEPH-226
C  5  T21.                                                              DEPH-227
C  6  T22.                                                              DEPH-228
C  7  KYY OR D  DEFINED AS -A(1100 1100)-A(1100 1-100).                 DEPH-229
C  8  KXX OR R  DEFINED AS  A(1100 1100)-A(1100 1-100).                 DEPH-230
C  9  KZZ OR A' DEFINED AS  A(1000 1000).                               DEPH-231
C 10  KXZ OR R' DEFINED AS -SQRT(2.) A(1100 1000).                      DEPH-232
C 11  KZX OR A  DEFINED AS -SQRT(2.) A(1000 1100).                      DEPH-233
C 12  SPIN-FLIP :  (A(0000,0000)+A(1100,1100)+A(1100,1-100))/2.         DEPH-234
C 13  VECTOR ANALYSING POWER OF THE TARGET.                             DEPH-235
C   (NOTE A RATIO SQRT(2.) WITH IT11 FOR SPIN 1/2 FOR 2, 3 AND 13).     DEPH-236
C 14  AYY  DEFINED AS -A(1111 0000)-A(111-1 0000).                      DEPH-237
C 15  AXX  DEFINED AS  A(1111 0000)-A(111-1 0000).                      DEPH-238
C 16  AZZ  DEFINED AS  A(1010 0000).                                    DEPH-239
C 17  AXZ  DEFINED AS -SQRT(2.) A(1110 0000).                           DEPH-240
C 18  AZX  DEFINED AS -SQRT(2.) A(1011 0000).                           DEPH-241
C 19  RESERVED FOR A SET OF EXPERIMENTAL DATA WHICH ARE REACTION        DEPH-242
C          CROSS-SECTIONS.                                              DEPH-243
C***********************************************************************DEPH-244
      KX=KCC                                                            DEPH-245
      IY=1                                                              DEPH-246
      IF (.NOT.LO(92)) GO TO 14                                         DEPH-247
      READ (MR,1002) (IPI(9,I),I=1,NCOLL)                               DEPH-248
      WRITE (MW,1003) (I,IPI(9,I),I=1,NCOLL)                            DEPH-249
   14 DO 24 I=1,NCOLL                                                   DEPH-250
      IF (LO(92)) GO TO 19                                              DEPH-251
C THERE MUST BE AT LEAST PLACE TO STORE SIX OBSERVABLES                 DEPH-252
      IF (5*KX+30.GE.NMAX) CALL MEMO('DEPH',NMAX,5*KX+30)               DEPH-253
      KX=KX+1                                                           DEPH-254
      IPI(9,I)=1                                                        DEPH-255
      MF(2,KX)=0                                                        DEPH-256
      IF ((I.NE.1).OR.(IPI(4,1).EQ.0)) GO TO 15                         DEPH-257
      IPI(9,1)=2                                                        DEPH-258
      KX=KX+1                                                           DEPH-259
      MF(2,KX)=1                                                        DEPH-260
   15 MF(5,KX)=NGR(IY)                                                  DEPH-261
      IF (NI.EQ.1) GO TO 19                                             DEPH-262
      IPI(9,I)=IPI(9,I)+1                                               DEPH-263
      KX=KX+1                                                           DEPH-264
      MF(2,KX)=2                                                        DEPH-265
      IX=NPR(IY)/2                                                      DEPH-266
      MF(5,KX)=NPR(1)-2*IX                                              DEPH-267
      IF (NI.GT.2) GO TO 17                                             DEPH-268
      IF (I.EQ.1) GO TO 19                                              DEPH-269
      IF (IPI(2,I).GT.3) GO TO 16                                       DEPH-270
      IF (IPI(2,I).EQ.1) GO TO 19                                       DEPH-271
      IPI(9,I)=IPI(9,I)+1                                               DEPH-272
      KX=KX+1                                                           DEPH-273
      MF(2,KX)=3                                                        DEPH-274
      MF(5,KX)=MOD(IX,2)                                                DEPH-275
      IX=IX/2                                                           DEPH-276
   16 KX=KX+1                                                           DEPH-277
      MF(2,KX)=12                                                       DEPH-278
      MF(5,KX)=MOD(IX,2)                                                DEPH-279
      IPI(9,I)=IPI(9,I)+1                                               DEPH-280
      GO TO 19                                                          DEPH-281
   17 IPI(9,I)=IPI(9,I)+3                                               DEPH-282
      DO 18 L=4,6                                                       DEPH-283
      KX=KX+1                                                           DEPH-284
      MF(2,KX)=L                                                        DEPH-285
      MF(5,KX)=MOD(IX,2)                                                DEPH-286
   18 IX=IX/2                                                           DEPH-287
   19 K1=KTGR+1                                                         DEPH-288
      IPI(8,I)=K1                                                       DEPH-289
      KTGR=KTGR+IPI(9,I)                                                DEPH-290
      IPI(9,I)=KTGR                                                     DEPH-291
      IF (5*KTGR.GE.NMAX) CALL MEMO('DEPH',NMAX,5*KTGR)                 DEPH-292
      DO 20 J=K1,KTGR                                                   DEPH-293
   20 MF(1,J)=I                                                         DEPH-294
      IF (.NOT.LO(92)) GO TO 22                                         DEPH-295
      READ (MR,1002) (MF(2,J),J=K1,KTGR)                                DEPH-296
      IF (MF(2,K1).NE.0) WRITE (MW,1004)                                DEPH-297
      READ (MR,1002) (MF(5,J),J=K1,KTGR)                                DEPH-298
      MF(2,K1)=0                                                        DEPH-299
      WRITE (MW,1005) I,(MF(2,J),J=K1,KTGR)                             DEPH-300
      WRITE (MW,1006) (MF(5,J),J=K1,KTGR)                               DEPH-301
      DO 21 J=K1,KTGR                                                   DEPH-302
      IF (MF(2,J).GT.18) GO TO 32                                       DEPH-303
   21 CONTINUE                                                          DEPH-304
   22 JT=0                                                              DEPH-305
      DO 23 J=K1,KTGR                                                   DEPH-306
      NJY=MAX0(NJY,-MF(2,J))                                            DEPH-307
      IF (MF(2,J).EQ.0.OR.MF(2,J).EQ.1) GO TO 23                        DEPH-308
      JT=JT+1                                                           DEPH-309
   23 CONTINUE                                                          DEPH-310
      JTX=MAX0(JTX,JT)                                                  DEPH-311
      NJC=MAX0(NJC,KTGR-K1+2)                                           DEPH-312
   24 IY=2                                                              DEPH-313
      IF (LO(81)) NJC=MAX0(NJC+2,6)                                     DEPH-314
   25 IY=0                                                              DEPH-315
      IF (LO(44)) IY=1                                                  DEPH-316
      IF (LO(46).AND.LO(44)) WRITE (MW,1007)                            DEPH-317
      DO 28 I=1,NCOLL                                                   DEPH-318
      IF (LO(46).AND.LO(44)) GO TO 27                                   DEPH-319
      DO 26 J=1,NCOLL                                                   DEPH-320
   26 NIV(I,J,3)=IY                                                     DEPH-321
      IF (LO(108)) NIV(I,I,3)=1                                         DEPH-322
      GO TO 28                                                          DEPH-323
   27 READ (MR,1002) (NIV(I,J,3),J=1,NCOLL)                             DEPH-324
      WRITE (MW,1008) I,(J,NIV(I,J,3),J=1,NCOLL)                        DEPH-325
   28 CONTINUE                                                          DEPH-326
      NJ=NJX                                                            DEPH-327
      DO 30 I=1,NCOLL                                                   DEPH-328
      DO 29 J=1,I                                                       DEPH-329
      L=NIV(I,J,3)                                                      DEPH-330
      NIV(I,J,3)=0                                                      DEPH-331
C NO COUL. CORR. WITH DIFF. PARTICLE MASS, CHARGES OR A CLOSED CHANNEL  DEPH-332
      IF (L.EQ.0.OR.IPI(4,I).NE.IPI(4,J).OR.(DABS(WV(1,I)-WV(1,J)).GT.1DDEPH-333
     10).OR.WV(3,I).LT.0.D0.OR.WV(3,J).LT.0.D0) GO TO 29                DEPH-334
      NJX=NJX+1                                                         DEPH-335
      NIV(I,J,3)=NJX                                                    DEPH-336
   29 NIV(J,I,3)=NIV(I,J,3)                                             DEPH-337
   30 CONTINUE                                                          DEPH-338
      IF (NJ.EQ.NJX) RETURN                                             DEPH-339
      WRITE (MW,1009)                                                   DEPH-340
      DO 31 I=1,NCOLL                                                   DEPH-341
   31 WRITE (MW,1010) I,(NIV(I,J,3),J,J=1,NCOLL)                        DEPH-342
      RETURN                                                            DEPH-343
   32 WRITE (MW,1011)                                                   DEPH-344
      STOP                                                              DEPH-345
 1000 FORMAT (/5X,I3,' SOLUTIONS',I10,' COUPLED EQUATIONS',I10,' INDEPENDEPH-346
     1DENT AMPLITUDES.'/)                                               DEPH-347
 1001 FORMAT (/5X,I3,' SOLUTIONS',I10,' COUPLED EQUATIONS',I10,' INDEPENDEPH-348
     1DENT AMPLITUDES',I10,' SETS OF COMPOUND COEFFICIENTS.'/)          DEPH-349
 1002 FORMAT (14I5)                                                     DEPH-350
 1003 FORMAT (/' EQUIDISTANT ANGLES OUTPUT:'//(10X,'CHANNEL',I3,5X,I3,' DEPH-351
     1OBSERVABLES.'))                                                   DEPH-352
 1004 FORMAT (' THE FIRST OBSERVABLE MUST BE THE CROSS SECTION: LABEL 0.DEPH-353
     1')                                                                DEPH-354
 1005 FORMAT (' OBSERVABLES FOR CHANNEL',I3,'  : ',18I5/(24X,18I5))     DEPH-355
 1006 FORMAT (12X,'GRAPH INFORMATION: ',18I5/(24X,18I5))                DEPH-356
 1007 FORMAT (/' COULOMB CORRECTIONS ( 1=YES, 0=NO):')                  DEPH-357
 1008 FORMAT (' I =',I2,4X,10(' J =',I2,':',I1,2X)/(10X,10(' J =',I2,':'DEPH-358
     1,I1,2X)))                                                         DEPH-359
 1009 FORMAT (/' STORAGE OF COULOMB CORRECTIONS:')                      DEPH-360
 1010 FORMAT (' I =',I2,2X,6(I6,' FOR J =',I2)/(10X,6(I6,' FOR J =',I2))DEPH-361
     1)                                                                 DEPH-362
 1011 FORMAT (//' THE STANDARD OBSERVABLES ARE ONLY 19 FOR ANY OTHER GIVDEPH-363
     1E A NEGATIVE INTEGER'///' IN DEPH  ...  STOP  ...')               DEPH-364
      END                                                               DEPH-365
