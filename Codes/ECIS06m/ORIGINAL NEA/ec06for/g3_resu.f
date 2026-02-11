C 23/03/07                                                      ECIS06  RESU-000
      SUBROUTINE RESU(IPI,SR,TX,SRX,JMAX,KMAX,NCOLL,NCOLS,MF,CQ,MFM,FM,IRESU-001
     1PJ,IPK,DONN,NCOLR,NCO,COE,WV,WVM,FCN,NOI,XD,JMIN,NRZ,NJC,RES,AM,EXRESU-002
     2,THE,DXX,SPG,NZZ,LO)                                              RESU-003
C COMPUTES CROSS-SECTIONS AND POLARISATIONS - COMPARE TO EXPERIMENTAL   RESU-004
C RESULTS -   OBTAINS EXPERIMENTAL NORMALISATIONS AND PARTIAL CHI2.     RESU-005
C INPUT:     IPI(J,I):PARITY OF THE NUCLEAR STATES (+/-) FOR J=1,       RESU-006
C                     MULTIPLICITY OF PARTICLE AND TARGET J=2,3,        RESU-007
C                     FIRST/LAST AMPL. AND OBSERVABLE FOR EACH          RESU-008
C                     LEVEL FOR J=6 TO 9. SEE CALX.                     RESU-009
C            SR:      HELICITY SCATTERING COEFFICIENTS.                 RESU-010
C            TX:      INELASTIC CROSS-SECTIONS IN MILLIBARNS            RESU-011
C                     FOLLOWED BY HAUSER-FESHBACH COEFFICIENTS.         RESU-012
C            SRX:     EQUIVALENT WITH SR BY CALL TO SAVE SR AND TX FOR  RESU-013
C                     MINIMUM CHI2 IN A SEARCH.                         RESU-014
C            JMAX:    MAXIMUM NUMBER OF CHANNEL SPINS.                  RESU-015
C            KMAX:    MAXIMUM NUMBER OF L FOR COMPOUND NUCLEUS.         RESU-016
C            NCOLL:   NUMBER OF COUPLED CHANNELS.                       RESU-017
C            NCOLS:   NUMBER OF CHANNELS WITH ANGULAR DISTRIBUTION.     RESU-018
C            MF,CQ:   TABLES OF HELICITY, DESCRIPTION OF OBSERVABLES    RESU-019
C                     ... ETC.  SEE DEPH AND OBSE.                      RESU-020
C            MFM,FM:  DESCRIPTION OF EXPERIMENTAL DATA.  SEE LECD.      RESU-021
C            IPJ:     NUMBER OF CHANNEL SPINS USED.                     RESU-022
C            IPK:     NUMBER OF L FOR COMPOUND NUCLEUS.                 RESU-023
C            DONN:    EXPERIMENTAL DATA:ANGLE, VALUE, EXPERIMENTAL      RESU-024
C                     ERROR ,ANGULAR WIDTH, ANGULAR ERROR AND           RESU-025
C                     CALCULATED ERROR.                                 RESU-026
C            NCOLR:   NUMBER OF EXPERIMENTAL ANGULAR DISTRIBUTIONS.     RESU-027
C            NCO,COE: INDICATIONS FOR OBSERVABLES.  SEE OBSE.           RESU-028
C            WV(I,J): WAVE NUMBER AND COULOMB PARAMETER.                RESU-029
C            WVM:     SAME AS WV FOR THE CONTINUUM.                     RESU-030
C            FCN:     COMPOUND NUCLEUS COEFFICIENTS.                    RESU-031
C            NOI:     STARTING AND FINAL ADDRESSES FOR CONTINUA.        RESU-032
C            XD:      ENERGY STEP FOR THE CONTINUA.                     RESU-033
C            JMIN:    TWICE MINIMUM CHANNEL SPIN.                       RESU-034
C            NRZ:     LENGTH OF SR+TX TO BE SAVED FOR MINIMUM CHI2.     RESU-035
C            NJC:     FIRST DIMENSION OF WORKING AREA EX.               RESU-036
C            NZZ:     TOTAL LENGTH OF WORKING FIELD,RETURNS SPACE USED. RESU-037
C            LO:      LOGICAL CONTROLS:                                 RESU-038
C               LO(8)  =.TRUE. RELATIVISTIC KINEMATICS.                 RESU-039
C               LO(31) =.TRUE. INPUT OF EXPERIMENTAL DATA AND CHI2      RESU-040
C                              CALCULATION.                             RESU-041
C               LO(32) =.TRUE. AUTOMATIC SEARCH ON SOME PARAMETERS.     RESU-042
C               LO(33) =.TRUE. SYMMETRISED CHI2 FOR CROSS-SECTIONS.     RESU-043
C               LO(59) =.TRUE. PRINT RESULTS ON FILES 58 AND 59.        RESU-044
C               LO(64) =.TRUE. PRINT RESULTS ON FILES 64 AND 66.        RESU-045
C               LO(66) =.TRUE. NO CALCULATION AT EQUIDISTANT ANGLES.    RESU-046
C               LO(67) =.TRUE. NO PLOT OF EXPERIMENTAL DATA.            RESU-047
C               LO(68) =.TRUE. NO PLOT OF EQUIDISTANT CROSS-SECTIONS.   RESU-048
C               LO(69) =.TRUE. NO PLOT OF EQUIDISTANT POLARISATIONS.    RESU-049
C               LO(74) =.TRUE. OUTPUT OF TIME IN DIFFERENT STEPS.       RESU-050
C               LO(81) =.TRUE. HAUSER-FESHBACH CORRECTIONS.             RESU-051
C               LO(84) =.TRUE. UNCOUPLED LEVELS FOR COMPOUND NUCLEUS.   RESU-052
C               LO(85) =.TRUE. FISSION TRANSMISSION COEFFICIENTS.       RESU-053
C               LO(86) =.TRUE. GAMMA EMISSION IN COMPOUND NUCLEUS.      RESU-054
C               LO(91) =.TRUE. ANGLES IN THE LABORATORY SYSTEM.         RESU-055
C               LO(116)=.TRUE. NO OUTPUT.                               RESU-056
C               LO(118)=.TRUE. FOR LAST RESULTS.                        RESU-057
C               LO(126)=.TRUE. SOME OBSERVABLES IN LABORATORY SYSTEM.   RESU-058
C OUTPUT:    RES:     DIFFERENCE BETWEEN EXPERIMENTAL AND CALCULATED    RESU-059
C                     VALUE DIVIDED BY ERROR.                           RESU-060
C WORKING AREAS:                                                        RESU-061
C            AM:      FOR GENERAL PURPOSES.                             RESU-062
C            EX:      FOR OBSERVABLES.                                  RESU-063
C            THE:     ANGLES FOR PLOT.                                  RESU-064
C            DXX:     CROSS-SECTIONS FOR PLOTS.                         RESU-065
C            SPG:     POLARISATIONS FOR PLOTS                           RESU-066
C                                                                       RESU-067
C FOR THE COMMON  /ANGUL/ SEE LECT.                                     RESU-068
C FOR THE COMMONS /DCONS/ AND /TITRE/ SEE CALC.                         RESU-069
C FOR THE COMMONS /DCHI2/ AND /NCOMP/  SEE CALX.                        RESU-070
C FOR THE COMMON  /RESCT/ SEE SCAT.                                     RESU-071
C                                                                       RESU-072
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /ANGUL/:                     RESU-073
C  THETA1:    FIRST ANGLE FOR COMPUTATION AT EQUIDISTANT ANGLES.        RESU-074
C  DTHETA:    STEP FOR COMPUTATION AT EQUIDISTANT ANGLES.               RESU-075
C  THETA2:    LAST ANGLE FOR COMPUTATION AT EQUIDISTANT ANGLES.         RESU-076
C  DTHE:      AVERAGING ANGLE.                                          RESU-077
C  NCJ:       NUMBER OF FACTORISATIONS OF 1/(1-COS(THETA)) IN AMPLITUDE.RESU-078
C  JMM(1):    NUMBER OF CHANNEL SPINS USED FOR MINIMUM CHI2.            RESU-079
C  JMM(2):    NUMBER OF L FOR COMPOUND NUCLEUS FOR MINIMUM CHI2.        RESU-080
C   DEFINED:  JMM.                                                      RESU-081
C   USED:     THETA1,DTHETA,THETA2,DTHE,NCJ,JMM.                        RESU-082
C                                                                       RESU-083
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCHI2/:                     RESU-084
C  CHI2:      CHI-SQUARE COMPUTED IN SUBROUTINE RESU.                   RESU-085
C  CHI2M:     MINIMUM CHI-SQUARE IN THE SEARCH.                         RESU-086
C   DEFINED:  CHI2,CHI2M.                                               RESU-087
C   USED:     CHI2,CHI2M.                                               RESU-088
C                                                                       RESU-089
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     RESU-090
C  CMB:       ATOMIC MASS CM DIVIDED BY H-BAR*C.                        RESU-091
C   USED:     CMB.                                                      RESU-092
C                                                                       RESU-093
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /NCOMP/:                     RESU-094
C  NSP(1):    NUMBER OF UNCOUPLED LEVELS FOR COMPOUND NUCLEUS           RESU-095
C             INCLUDING DISCRETISATION OF CONTINUUM.                    RESU-096
C  NSP(3):    NUMBER OF THESE LEVELS WITHOUT ANGULAR DISTRIBUTION.      RESU-097
C  NCONT:     NUMBER OF CONTINUUM FOR COMPOUND NUCLEUS.                 RESU-098
C  NIE:       NUMBER OF UNCOUPLED STATES ADDED FOR DISCRETISATION.      RESU-099
C  NCOLX:     TOTAL NUMBER OF LEVELS WITHOUT DISCRETISATION.            RESU-100
C  NDP:       ADDRESS OF FISSION AND GAMMA FINAL RESULTS.               RESU-101
C   DEFINED:  NDP.                                                      RESU-102
C   USED:     NSP,NCONT,NIE,NCOLX,NDP.                                  RESU-103
C                                                                       RESU-104
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /RESCT/:                     RESU-105
C  EXCN:      COMPOUND-NUCLEUS CROSS-SECTION.                           RESU-106
C  XA,XB:     CONSTANTS FOR CHANGE TO LABORATORY SYSTEM.                RESU-107
C  JN1,JN2:   FIRST AND LAST AMPLITUDE IN THE TABLE.                    RESU-108
C  NTT:       TOTAL NUMBER OF AMPLITUDES FOR THIS STATE AND A GIVEN J.  RESU-109
C  MTT:       SIZE OF THE SCATTERING MATRIX FOR THIS STATE AND ANGLE.   RESU-110
C  MT1:       MULTIPLICITY OF THE OUTGOING PARTICLE.                    RESU-111
C  MT2:       MULTIPLICITY OF THE RESIDUAL TARGET.                      RESU-112
C  MT3:       MULTIPLICITY OF THE INCIDENT PARTICLE.                    RESU-113
C  MT4:       MULTIPLICITY OF THE INITIAL TARGET.                       RESU-114
C  NOUT:      NUCLEAR STATE CONSIDERED.                                 RESU-115
C  NIX,NFX:   FIRST AND LAST OBSERVABLE IN THE TABLE.                   RESU-116
C  KLT:       1 FOR EXPERIMENTAL DATA.                                  RESU-117
C             2 FOR EQUIDISTANT ANGLES.                                 RESU-118
C             3 FOR PURE COMPOUND NUCLEUS.                              RESU-119
C  LKT:       LOGICAL=.TRUE. FOR CENTRE OF MASS SYSTEM.                 RESU-120
C   DEFINED:  XA,XB,JN1,JN2,NTT,MTT,MT1,MT2;MT3,MT4,NOUT,NIX,NFX,KLT,   RESU-121
C             LKT.                                                      RESU-122
C   USED:     EXCN,XA,XB,JN1,JN2,NTT,MTT,MT1,MT2;MT3,MT4,NOUT,NIX,NFX,  RESU-123
C             KLT,LKT.                                                  RESU-124
C                                                                       RESU-125
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /TITRE/:                     RESU-126
C  TITLE(18): TITLE OF THE RUN PRINTED AS HEADING OF RESULTS.           RESU-127
C   USED:     TITLE.                                                    RESU-128
C                                                                       RESU-129
C INTERNAL LOGICAL LT(3):                                               RESU-130
C            XD:      ENERGY STEP FOR THE CONTINUA.                     RESU-131
C            LT(1)=.TRUE. FOR CROSS-SECTION IN MILLIBARNS.              RESU-132
C            LT(2)=.TRUE. FOR SYMMETRISED CHI2, CROSS-SECTION AND       RESU-133
C                         NORMALISATION ERROR.                          RESU-134
C            LT(3)=.TRUE. FOR GRAPH OF POLARISATION.                    RESU-135
C            LKT  =.FALSE. FOR ANGLE IN THE LABORATORY SYSTEM.          RESU-136
C***********************************************************************RESU-137
      IMPLICIT REAL*8 (A-H,O-Z)                                         RESU-138
      LOGICAL LO(150),LT(3),LKT                                         RESU-139
      DIMENSION IPI(11,*),SR(2,JMAX,*),TX(*),SRX(*),MF(10,*),MFM(14,*),FRESU-140
     1M(7,*),DONN(6,*),NCO(*),COE(*),WV(22,*),WVM(22,*),FCN(KMAX,*),NOI(RESU-141
     22,*),XD(3,*),RES(*),AM(*),EX(NJC,4),THE(*),DXX(*),SPG(*),ZX(3),ZY(RESU-142
     33),SP(6)                                                          RESU-143
      CHARACTER*1 SIGM(2)                                               RESU-144
      CHARACTER*4 LG(10),TITLE,CQ(10,1)                                 RESU-145
      COMMON /ANGUL/ THETA1,THETA2,DTHETA,DTHE,NCJ,NL(3),JMM(2)         RESU-146
      COMMON /DCHI2/ CHI2,CHI2M,YY(3)                                   RESU-147
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XZ                            RESU-148
      COMMON /INOUT/ MR,MW,MS                                           RESU-149
      COMMON /NCOMP/ NSP(3),NRD(2),NCONT,NCOJ,NCONS,NIE,NCOLX,NDP,NDQ,ACRESU-150
     1N(20)                                                             RESU-151
      COMMON /RESCT/ EXCN,XA,XB,JN1,JN2,NTT,MTT,MT1,MT2,MT3,MT4,NOUT,NIXRESU-152
     1,NFX,KLT,LKT                                                      RESU-153
      COMMON /TITRE/ TITLE(18)                                          RESU-154
      DATA SIGM,LG,I3,I4,I5,I7,LT /'+','-',' COM','POUN','D NU','CLEU','RESU-155
     1S ',2*'    ','   D','IREC','T ',4*0,3*.FALSE./                    RESU-156
      NESP=0                                                            RESU-157
      MT3=IPI(2,1)                                                      RESU-158
      MT4=IPI(3,1)                                                      RESU-159
      IF (.NOT.LO(118)) GO TO 2                                         RESU-160
C COPY THE SCATTERING COEFFICIENTS OBTAINED FOR MINIMUM CHI2.           RESU-161
      IPJ=JMM(1)                                                        RESU-162
      IPK=JMM(2)                                                        RESU-163
      DO 1 I=1,NRZ                                                      RESU-164
    1 SRX(I)=SRX(I+NRZ)                                                 RESU-165
    2 IF (.NOT.LO(31)) GO TO 29                                         RESU-166
C CALCULATION AT EXPERIMENTAL ANGLES.                                   RESU-167
      KLT=1                                                             RESU-168
      CHI2=0.D0                                                         RESU-169
      NOUT=0                                                            RESU-170
      IF (LO(64)) WRITE (64,1000) WV(1,1),WV(13,1),WV(2,1),IPI(4,1),NCOLRESU-171
     1R                                                                 RESU-172
      KI=IPI(9,NCOLL)                                                   RESU-173
      JI=1                                                              RESU-174
      KA=0                                                              RESU-175
C PSEUDO LOOP ON THE ANGULAR DISTRIBUTIONS.                             RESU-176
    3 AA=0.D0                                                           RESU-177
      BB=0.D0                                                           RESU-178
      CC=0.D0                                                           RESU-179
      NIX=KI+JI                                                         RESU-180
      NFX=NIX                                                           RESU-181
      LT(1)=MF(2,NIX).EQ.0                                              RESU-182
      LT(2)=LO(33).AND.((IABS(MF(2,NIX)).LE.1.OR.MF(2,NIX).EQ.19).AND.FMRESU-183
     1(5,JI).EQ.0.D0)                                                   RESU-184
      JIM=JI                                                            RESU-185
      JIN=JI                                                            RESU-186
    4 J1=MFM(2,JI)                                                      RESU-187
      J2=MFM(3,JI)                                                      RESU-188
      LKT=MFM(4,JI).NE.1                                                RESU-189
      IF (J2.GE.J1) GO TO 5                                             RESU-190
      JI=JI+1                                                           RESU-191
      GO TO 4                                                           RESU-192
    5 JIF=JI                                                            RESU-193
      J=J1                                                              RESU-194
C PSEUDO LOOP ON THE EXPERIMENTAL DATA.                                 RESU-195
    6 DONN(6,J)=DONN(3,J)                                               RESU-196
      KC=1                                                              RESU-197
      THETA=DONN(1,J)                                                   RESU-198
    7 JI=JIM                                                            RESU-199
      ZX(KC)=0.D0                                                       RESU-200
      ZY(KC)=0.D0                                                       RESU-201
      IF (MF(2,KI+JI).NE.19) GO TO 63                                   RESU-202
C DATA ARE TOTAL CROSS SECTIONS.                                        RESU-203
      I=IDINT(DONN(1,J)*1.000001D0)                                     RESU-204
      IF (I.GT.NCOLX-NCONT) GO TO 8                                     RESU-205
      IF (I.EQ.-1) ZX(1)=TX(1)                                          RESU-206
      IF (I.EQ.0) ZX(1)=TX(1)-TX(2)                                     RESU-207
      IF (I.EQ.0.AND.NCOLL.NE.NCOLX) ZX(1)=ZX(1)-TX(NCOLL+2)            RESU-208
      IF (I.LE.0) GO TO 18                                              RESU-209
      K=MAX0(I,1)+1+NCOLL                                               RESU-210
      ZX(1)=0.D0                                                        RESU-211
      IF (NCOLL.NE.NCOLX) ZX(1)=TX(NCOLL+I+1)                           RESU-212
      IF (I.LE.NCOLL) ZX(1)=ZX(1)+TX(I+1)                               RESU-213
      GO TO 18                                                          RESU-214
    8 IF (I.GT.NCOLX) GO TO 10                                          RESU-215
      IK=NCOLL+NCOLX-NCONT+1                                            RESU-216
      II=I+NCONT-NCOLX                                                  RESU-217
      IJ=NOI(1,II)                                                      RESU-218
      JI=NOI(2,II)                                                      RESU-219
      ZX(1)=0.D0                                                        RESU-220
      IF (IJ.GT.JI) GO TO 18                                            RESU-221
      DO 9 II=IJ,JI                                                     RESU-222
    9 ZX(1)=ZX(1)+TX(IK+II)                                             RESU-223
      GO TO 18                                                          RESU-224
   10 ZX(1)=TX(NDP+I-NCOLX)                                             RESU-225
      GO TO 18                                                          RESU-226
   11 IF (LT(1)) EX(2,KC)=1.D0                                          RESU-227
      ZX(KC)=(ZX(KC)*ZY(KC)+EX(2,KC)*EX(1,KC))/(ZY(KC)+EX(1,KC))        RESU-228
      ZY(KC)=ZY(KC)+EX(1,KC)                                            RESU-229
      IF (JI.EQ.JIF) GO TO 12                                           RESU-230
      JI=JI+1                                                           RESU-231
      GO TO 63                                                          RESU-232
   12 IF (LT(1)) ZX(KC)=ZX(KC)*ZY(KC)                                   RESU-233
      GO TO ( 13 , 14 , 15 ),KC                                         RESU-234
   13 IF (DONN(4,J).EQ.0.D0) GO TO 18                                   RESU-235
      KC=2                                                              RESU-236
      THETA=DONN(1,J)-DONN(4,J)                                         RESU-237
      GO TO 7                                                           RESU-238
   14 KC=3                                                              RESU-239
      THETA=DONN(1,J)+DONN(4,J)                                         RESU-240
      GO TO 7                                                           RESU-241
   15 IF (LT(1)) GO TO 16                                               RESU-242
      ZX(1)=(ZX(1)*ZY(1)+ZX(2)*ZY(2)+ZX(3)*ZY(3))/(ZY(1)+ZY(2)+ZY(3))   RESU-243
      GO TO 17                                                          RESU-244
   16 ZX(1)=(ZX(1)+ZX(2)+ZX(3))/3.D0                                    RESU-245
   17 IF (DONN(5,J).NE.0.D0) DONN(6,J)=DSQRT(DONN(3,J)*DONN(3,J)+((ZX(2)RESU-246
     1-ZX(3))*DONN(2,J)*DONN(5,J)/(2.D0*DONN(4,J)*ZX(1)))**2)           RESU-247
   18 RES(J)=ZX(1)                                                      RESU-248
      IF (LT(2)) DONN(6,J)=DONN(6,J)*DSQRT(FM(4,JI)*RES(J)/DONN(2,J))   RESU-249
      XC=DONN(6,J)**2/FM(3,JI)                                          RESU-250
      AA=AA+RES(J)**2/XC                                                RESU-251
      BB=BB+RES(J)*DONN(2,J)/XC                                         RESU-252
      CC=CC+DONN(2,J)**2/XC                                             RESU-253
      J=J+1                                                             RESU-254
C END OF THE PSEUDO LOOP ON EXPERIMENTAL DATA.                          RESU-255
      IF (J.LE.J2) GO TO 6                                              RESU-256
C COMPUTATION OF NORMALISATION AND CHI2.                                RESU-257
      IF (FM(5,JI).EQ.0.D0) GO TO 21                                    RESU-258
      XA=FM(3,JI)/FM(5,JI)**2                                           RESU-259
      AA=AA+XA                                                          RESU-260
      XB=FM(4,JI)                                                       RESU-261
      BB=BB+XA*XB                                                       RESU-262
      CC=CC+XA*XB*XB                                                    RESU-263
      IF (JI.EQ.NCOLR) GO TO 19                                         RESU-264
      IF (FM(4,JI).NE.FM(4,JI+1).OR.FM(5,JI).NE.FM(5,JI+1)) GO TO 19    RESU-265
      IF ((2*IABS(MF(2,NIX))-3)*(2*IABS(MF(2,NIX+1))-3).LT.0) GO TO 19  RESU-266
      JI=JI+1                                                           RESU-267
      JIM=JI                                                            RESU-268
      GO TO 4                                                           RESU-269
   19 IF (LT(2)) GO TO 20                                               RESU-270
      FM(6,JI)=BB/AA                                                    RESU-271
      GO TO 22                                                          RESU-272
   20 FM(6,JI)=DSQRT(CC/AA)                                             RESU-273
      GO TO 22                                                          RESU-274
   21 FM(6,JI)=FM(4,JI)                                                 RESU-275
   22 DO 27 JJ=JIN,JI                                                   RESU-276
      NIX=KI+JJ                                                         RESU-277
      J1=MFM(2,JJ)                                                      RESU-278
      J2=MFM(3,JJ)                                                      RESU-279
      J3=J2-J1+1                                                        RESU-280
      IF (LO(64)) WRITE (64,1001) (MF(J,NIX),J=1,2),J3,(CQ(J,NIX),J=6,10RESU-281
     1)                                                                 RESU-282
      IF (J3.LE.0) GO TO 27                                             RESU-283
      IF (JJ.NE.JI) FM(6,JJ)=FM(6,JI)                                   RESU-284
      IF (LO(116)) GO TO 23                                             RESU-285
C OUTPUT OF THE CALCULATED AND THE EXPERIMENTAL VALUES.                 RESU-286
      M=0                                                               RESU-287
      WRITE (MW,1002) TITLE                                             RESU-288
      IF (MF(2,NIX).NE.19) WRITE (MW,1003) MFM(1,JJ)                    RESU-289
      WRITE (MW,1004) (CQ(J,NIX),J=6,10)                                RESU-290
   23 A1=FM(6,JJ)                                                       RESU-291
      A6=DSQRT(FM(3,JJ))                                                RESU-292
      AA=0.D0                                                           RESU-293
      IF (FM(5,JJ).EQ.0.D0) GO TO 24                                    RESU-294
      A2=(A1-FM(4,JJ))/FM(5,JJ)                                         RESU-295
      AA=A2**2                                                          RESU-296
      RES(J2+1)=A2*A6                                                   RESU-297
      CHI2=CHI2+RES(J2+1)**2                                            RESU-298
   24 DO 26 J=J1,J2                                                     RESU-299
      A2=DONN(2,J)/A1                                                   RESU-300
      A3=DONN(6,J)/A1                                                   RESU-301
      A4=((RES(J)-A2)/A3)**2                                            RESU-302
      IF (LO(116)) GO TO 25                                             RESU-303
      M=M+1                                                             RESU-304
      A5=DONN(3,J)/A1                                                   RESU-305
      THE(M)=DONN(1,J)                                                  RESU-306
      DXX(M)=A2                                                         RESU-307
      SPG(M)=RES(J)                                                     RESU-308
      WRITE (MW,1005) DONN(1,J),RES(J),A2,A5,A3,A4                      RESU-309
      IF (LO(64)) WRITE (64,1005) DONN(1,J),RES(J),A2,A5,A3,A4          RESU-310
   25 AA=AA+A4                                                          RESU-311
      RES(J)=(RES(J)-A2)*A6/A3                                          RESU-312
   26 CHI2=CHI2+RES(J)**2                                               RESU-313
      FM(7,JJ)=AA                                                       RESU-314
      IF (LO(116)) GO TO 27                                             RESU-315
      WRITE (MW,1006) (FM(J,JJ),J=3,7)                                  RESU-316
      IF (MF(1,NIX).EQ.19) GO TO 27                                     RESU-317
      IF (LO(74)) CALL HORA                                             RESU-318
      IF (LO(67)) GO TO 27                                              RESU-319
      LT(3)=.NOT.(LT(1).OR.MF(2,NIX).EQ.1.OR.MF(2,NIX).EQ.19)           RESU-320
      CALL GRAL(THE,SPG,DXX,M,MF(1,NIX),CQ(1,NIX),1,LT(3),.FALSE.)      RESU-321
   27 CONTINUE                                                          RESU-322
      JI=JI+1                                                           RESU-323
C END OF THE PSEUDO LOOP ON ANGULAR DISTRIBUTIONS.                      RESU-324
      IF (JI.LE.NCOLR) GO TO 3                                          RESU-325
      IF (.NOT.LO(116)) WRITE (MW,1007) CHI2                            RESU-326
      IF (.NOT.((CHI2.LE.CHI2M).AND.LO(32))) GO TO 29                   RESU-327
      JMM(1)=IPJ                                                        RESU-328
      JMM(2)=IPK                                                        RESU-329
      CHI2M=CHI2                                                        RESU-330
C IF THE CHI2 DECREASED, SAVE THE SCATTERING COEFFICIENTS.              RESU-331
      DO 28 I=1,NRZ                                                     RESU-332
   28 SRX(I+NRZ)=SRX(I)                                                 RESU-333
   29 IF (LO(116)) GO TO 62                                             RESU-334
      LKT=.NOT.LO(91)                                                   RESU-335
      IF (LO(66).AND.(.NOT.(LO(81).OR.LO(59)))) GO TO 62                RESU-336
C COMPUTATION AT EQUIDISTANT ANGLES.                                    RESU-337
      WRITE (MW,1002) TITLE                                             RESU-338
      JG=IDINT((THETA2-THETA1)/DTHETA+1.5D0)                            RESU-339
      IF (LO(64)) WRITE (66,1008) WV(1,1),WV(13,1),WV(2,1),IPI(4,1),NCOLRESU-340
     1S                                                                 RESU-341
      IF (.NOT.LO(59)) GO TO 30                                         RESU-342
      ND=1                                                              RESU-343
      IF (WV(5,1).EQ.0.D0) ND=3                                         RESU-344
      IF (LO(81)) ND=ND+1                                               RESU-345
      WRITE (58,1009) WV(1,1),WV(13,1),WV(2,1),IPI(4,1),ND              RESU-346
      IF (NCOLS.NE.1) WRITE (59,1010) WV(1,1),WV(13,1),WV(2,1),IPI(4,1),RESU-347
     1NCOLS-1                                                           RESU-348
   30 IF (WV(5,1).NE.0.D0) GO TO 31                                     RESU-349
      WRITE (MW,1011) TX(1)                                             RESU-350
      IF (LO(59)) WRITE (58,1012) TX(1)                                 RESU-351
   31 RX=TX(1)-TX(2)                                                    RESU-352
      IF (LO(59)) WRITE (58,1012) RX                                    RESU-353
      IF (LO(81)) GO TO 32                                              RESU-354
      WRITE (MW,1013) RX                                                RESU-355
      GO TO 42                                                          RESU-356
C COMPOUND NUCLEUS RESULTS.                                             RESU-357
   32 WRITE (MW,1014) RX                                                RESU-358
      NDP=2*NCOLL+NSP(1)+1                                              RESU-359
      RX=RX-TX(NCOLL+2)                                                 RESU-360
      WRITE (MW,1015) RX                                                RESU-361
      IF (LO(85)) WRITE (MW,1016) TX(NDP+1)                             RESU-362
      IF (LO(86)) WRITE (MW,1017) TX(NDP+2)                             RESU-363
      RX=TX(NDP+1)+TX(NDP+2)                                            RESU-364
      WRITE (MW,1018)                                                   RESU-365
      RY=0.D0                                                           RESU-366
      DO 33 I=1,NCOLL                                                   RESU-367
      II=IPI(1,I)+1                                                     RESU-368
      SP2=0.5D0*DFLOAT(IPI(3,I)-1)                                      RESU-369
      RY=RY+TX(NCOLL+I+1)                                               RESU-370
   33 WRITE (MW,1019) I,SP2,SIGM(II),WV(3,I),TX(NCOLL+I+1)              RESU-371
      WRITE (MW,1020) RY                                                RESU-372
      RX=RX+RY                                                          RESU-373
      IF (.NOT.LO(84)) GO TO 41                                         RESU-374
      IF (NCOLL.EQ.NCOLS) GO TO 35                                      RESU-375
      WRITE (MW,1021)                                                   RESU-376
      RY=0.D0                                                           RESU-377
      N1=NCOLL+1                                                        RESU-378
      DO 34 I=N1,NCOLS                                                  RESU-379
      II=IPI(1,I)+1                                                     RESU-380
      SP2=0.5D0*DFLOAT(IPI(3,I)-1)                                      RESU-381
      RY=RY+TX(NCOLL+I+1)                                               RESU-382
   34 WRITE (MW,1019) I,SP2,SIGM(II),WV(3,I),TX(NCOLL+I+1)              RESU-383
      WRITE (MW,1020) RY                                                RESU-384
      RX=RX+RY                                                          RESU-385
   35 NSP1=NSP(3)                                                       RESU-386
      IF (NCONT.NE.0) NSP1=NSP1-NIE                                     RESU-387
      IF (NSP1.LT.1) GO TO 37                                           RESU-388
      WRITE (MW,1022)                                                   RESU-389
      RY=0.D0                                                           RESU-390
      DO 36 I=1,NSP1                                                    RESU-391
      J=I+NCOLS                                                         RESU-392
      II=IPI(1,J)+1                                                     RESU-393
      SP2=0.5D0*DFLOAT(IPI(3,J)-1)                                      RESU-394
      RY=RY+TX(NCOLL+J+1)                                               RESU-395
   36 WRITE (MW,1019) J,SP2,SIGM(II),WV(3,J),TX(NCOLL+J+1)              RESU-396
      WRITE (MW,1020) RY                                                RESU-397
      RX=RX+RY                                                          RESU-398
   37 IF (NCONT.EQ.0) GO TO 41                                          RESU-399
      IK=NCOLL+NCOLX-NCONT+1                                            RESU-400
      DO 40 I=1,NCONT                                                   RESU-401
      IJ=NOI(1,I)                                                       RESU-402
      JI=NOI(2,I)                                                       RESU-403
      RY=0.D0                                                           RESU-404
      IF (IJ.GT.JI) GO TO 40                                            RESU-405
      DO 38 II=IJ,JI                                                    RESU-406
   38 RY=RY+TX(IK+II)                                                   RESU-407
      WRITE (MW,1023) I,RY                                              RESU-408
      RX=RX+RY                                                          RESU-409
      DO 39 II=IJ,JI                                                    RESU-410
      RY=TX(IK+II)/XD(2,II)                                             RESU-411
   39 WRITE (MW,1024) II,WVM(3,II),XD(2,II),TX(IK+II),RY                RESU-412
   40 CONTINUE                                                          RESU-413
   41 WRITE (MW,1025) RX                                                RESU-414
      WRITE (MW,1002) TITLE                                             RESU-415
      IF (LO(59)) WRITE (58,1012) RX                                    RESU-416
   42 INIV=1                                                            RESU-417
      SP2=0.5D0*DFLOAT(IPI(3,INIV)-1)                                   RESU-418
      WRITE (MW,1026) SP2,SIGM(IPI(1,INIV)+1)                           RESU-419
      NOUT=0                                                            RESU-420
      KLT=2                                                             RESU-421
C PSEUDO DO LOOP ON LEVELS.                                             RESU-422
   43 IF (WV(3,INIV).GT.0.D0) GO TO 44                                  RESU-423
      WRITE (MW,1027) INIV                                              RESU-424
      IF (LO(64)) WRITE (66,1028) INIV,SP2,SIGM(IPI(1,INIV)+1)          RESU-425
      IF (LO(59)) WRITE (59,1029) INIV-1                                RESU-426
      GO TO 57                                                          RESU-427
   44 RX=TX(INIV+1)                                                     RESU-428
      IF (LO(81)) RX=RX+TX(NCOLL+INIV+1)                                RESU-429
      IF (INIV.EQ.1) GO TO 45                                           RESU-430
      WRITE (MW,1002) TITLE                                             RESU-431
      WRITE (MW,1030) INIV,SP2,SIGM(IPI(1,INIV)+1)                      RESU-432
      WRITE (MW,1031) RX                                                RESU-433
      GO TO 46                                                          RESU-434
   45 IF (WV(5,1).NE.0.D0) GO TO 47                                     RESU-435
      WRITE (MW,1032) RX                                                RESU-436
   46 IF (LO(81)) WRITE (MW,1033) TX(INIV+1)                            RESU-437
      IF (.NOT.LO(59)) GO TO 47                                         RESU-438
      IF (INIV.EQ.1) WRITE (58,1012) RX                                 RESU-439
      IF (INIV.NE.1) WRITE (59,1034) RX,INIV-1                          RESU-440
   47 IF (LO(81)) WRITE (MW,1035) TX(NCOLL+INIV+1)                      RESU-441
      IF (LO(66)) GO TO 62                                              RESU-442
      IF (JG.LE.0) GO TO 57                                             RESU-443
      THETB=THETA1                                                      RESU-444
      IJ=0                                                              RESU-445
      II=1                                                              RESU-446
C PSEUDO DO LOOP ON ANGLES.                                             RESU-447
   48 THE(II)=THETB                                                     RESU-448
      THETA=THETB                                                       RESU-449
      KC=1                                                              RESU-450
      GO TO 64                                                          RESU-451
   49 DXX(II)=EX(I5,1)                                                  RESU-452
      IF (LO(81)) GO TO 50                                              RESU-453
      I1=MIN0(I3,7)                                                     RESU-454
      WRITE (MW,1036) THETB,(EX(K,1),K=2,I1)                            RESU-455
      GO TO 51                                                          RESU-456
   50 I1=MIN0(I3,5)                                                     RESU-457
      IF (I1.LT.3) GO TO 52                                             RESU-458
      WRITE (MW,1036) THETB,EX(2,1),EX(I3+1,1),EX(I4,1),(EX(K,1),K=3,I1)RESU-459
   51 IF (I1.GE.I3) GO TO 53                                            RESU-460
      I2=I1+1                                                           RESU-461
      I1=MIN0(I1+5,I3)                                                  RESU-462
      WRITE (MW,1037) (EX(K,1),K=I2,I1)                                 RESU-463
      GO TO 51                                                          RESU-464
   52 WRITE (MW,1036) THETB,EX(2,1),EX(I3+1,1),EX(I4,1)                 RESU-465
   53 IF (.NOT.LO(64)) GO TO 54                                         RESU-466
      WRITE (66,1038) (MF(2,NIX+K-2),THETA,EX(K,1),(CQ(L,NIX+K-2),L=6,10RESU-467
     1),K=2,I3)                                                         RESU-468
      IF (LO(81)) WRITE (66,1039) THETA,EX(I3+1,1),(LG(K),K=1,5),THETA,ERESU-469
     1X(I4,1),(LG(K),K=6,10)                                            RESU-470
   54 IF (I7.LE.0) GO TO 56                                             RESU-471
      DO 55 K=1,I7                                                      RESU-472
      IJ=IJ+1                                                           RESU-473
   55 SPG(IJ)=EX(K+I5,1)                                                RESU-474
   56 THETB=THETB+DTHETA                                                RESU-475
      II=II+1                                                           RESU-476
C END OF THE PSEUDO DO LOOP ON ANGLES.                                  RESU-477
      IF (II.LE.JG) GO TO 48                                            RESU-478
      IF (.NOT.LO(68)) CALL GRAL(THE,DXX,DXX,JG,MF(1,NFX-I7),CQ(1,NFX-I7RESU-479
     1),1,.FALSE.,.TRUE.)                                               RESU-480
      IF (.NOT.LO(69).AND.I7.GT.0) CALL GRAL(THE,SPG,DXX,JG,MF(1,NIX+I5-RESU-481
     11),CQ(1,NIX+I5-1),I7,.TRUE.,.TRUE.)                               RESU-482
      IF (LO(74)) CALL HORA                                             RESU-483
   57 INIV=INIV+1                                                       RESU-484
      SP2=0.5D0*DFLOAT(IPI(3,INIV)-1)                                   RESU-485
C END OF THE PSEUDO DO LOOP ON LEVELS.                                  RESU-486
      IF (INIV.LE.NCOLL) GO TO 43                                       RESU-487
      IF (.NOT.LO(81)) GO TO 62                                         RESU-488
      JN2=NCOLL                                                         RESU-489
      KLT=3                                                             RESU-490
   58 JN1=JN2+1                                                         RESU-491
      JN2=MIN0(JN1+5,NCOLS)                                             RESU-492
      IS=JN1-1                                                          RESU-493
      IF (JN1.GT.JN2) GO TO 62                                          RESU-494
      DO 59 I=JN1,JN2                                                   RESU-495
   59 SP(I-IS)=0.5D0*DFLOAT(IPI(3,I)-1)                                 RESU-496
      WRITE (MW,1002) TITLE                                             RESU-497
      WRITE (MW,1040) JN1,JN2,(I,SP(I-IS),SIGM(IPI(1,I)+1),I=JN1,JN2)   RESU-498
      THETB=THETA1                                                      RESU-499
      GO TO 61                                                          RESU-500
   60 WRITE (MW,1041) THETB,(EX(I-IS,1),I=JN1,JN2)                      RESU-501
      THETB=THETB+DTHETA                                                RESU-502
   61 THETA=THETB                                                       RESU-503
      KC=1                                                              RESU-504
      IF (THETB.LE.THETA2) GO TO 71                                     RESU-505
      GO TO 58                                                          RESU-506
   62 NZZ=2*NESP                                                        RESU-507
      RETURN                                                            RESU-508
C FOR EXPERIMENTAL DATA.                                                RESU-509
   63 IF (JI.EQ.KA) GO TO 71                                            RESU-510
      INIV=MFM(1,JI)                                                    RESU-511
      KA=JI                                                             RESU-512
      GO TO 65                                                          RESU-513
C FOR EQUIDISTANT ANGLES.                                               RESU-514
   64 IF (INIV.EQ.NOUT) GO TO 71                                        RESU-515
   65 NOUT=INIV                                                         RESU-516
C CHANGE OF LEVEL.                                                      RESU-517
      JN1=IPI(6,NOUT)                                                   RESU-518
      JN2=IPI(7,NOUT)                                                   RESU-519
      MT1=IPI(2,NOUT)                                                   RESU-520
      MT2=IPI(3,NOUT)                                                   RESU-521
      NTT=JN2+1-JN1                                                     RESU-522
      MTT=MT1*MT2*MT3*MT4                                               RESU-523
      M2=1+2*MTT                                                        RESU-524
      M3=M2+2*NTT                                                       RESU-525
      M4=M3+IPJ+NCJ                                                     RESU-526
      IF (LO(126)) M4=MAX0(M4,6*MTT)                                    RESU-527
      NESP=MAX0(NESP,M4)                                                RESU-528
      IF (2*M4.GT.NZZ) CALL MEMO('RESU',NZZ,2*M4)                       RESU-529
      IF (LKT) GO TO 67                                                 RESU-530
      IF (LO(8)) GO TO 66                                               RESU-531
      XA=1.D0                                                           RESU-532
      XB=WV(4,1)*WV(1,NOUT)/(WV(4,NOUT)*WV(2,1))                        RESU-533
      GO TO 67                                                          RESU-534
   66 XA=DSQRT(1.D0+(WV(4,1)/(CMB*WV(2,1)))**2)                         RESU-535
      XB=DSQRT(WV(1,NOUT)**2+(WV(4,NOUT)/CMB)**2)*WV(4,1)/(WV(4,NOUT)*WVRESU-536
     1(2,1))                                                            RESU-537
   67 IF (KLT.NE.2) GO TO 71                                            RESU-538
C FOR EQUIDISTANT ANGLES.                                               RESU-539
      NIX=IPI(8,INIV)                                                   RESU-540
      NFX=IPI(9,INIV)                                                   RESU-541
      I4=2+NFX-NIX                                                      RESU-542
      I3=I4                                                             RESU-543
      I5=2                                                              RESU-544
      I8=I4-1                                                           RESU-545
      IF (I4.GT.2.AND.MF(2,NIX+1).EQ.1) I5=3                            RESU-546
      I7=I4-I5                                                          RESU-547
      IF (LO(81)) GO TO 68                                              RESU-548
      NGX=MIN0(NFX,NIX+5)                                               RESU-549
      WRITE (MW,1042) ((CQ(L,K),L=6,10),K=NIX,NGX)                      RESU-550
      GO TO 70                                                          RESU-551
   68 I4=I4+2                                                           RESU-552
      I8=I8+2                                                           RESU-553
      NGX=MIN0(NFX,NIX+3)                                               RESU-554
      NJX=NIX+1                                                         RESU-555
      IF (NGX.GE.NJX) GO TO 69                                          RESU-556
      WRITE (MW,1042) (CQ(L,NIX),L=6,10),LG                             RESU-557
      GO TO 70                                                          RESU-558
   69 WRITE (MW,1042) (CQ(L,NIX),L=6,10),LG,((CQ(L,K),L=6,10),K=NJX,NGX)RESU-559
   70 IF (LO(64)) WRITE (66,1043) INIV,SP2,SIGM(IPI(1,INIV)+1),I8,JG    RESU-560
      IF (NGX.EQ.NFX) GO TO 71                                          RESU-561
      NJX=NGX+1                                                         RESU-562
      NGX=MIN0(NFX,NJX+4)                                               RESU-563
      WRITE (MW,1044) ((CQ(L,K),L=6,10),K=NJX,NGX)                      RESU-564
      IF (NGX.LT.NFX) GO TO 70                                          RESU-565
   71 CALL SCAT(SR,MF,JMAX,KMAX,THETA,FCN,IPJ,IPK,NCO,COE,IPI,AM,AM(M2),RESU-566
     1AM(M3),EX(1,KC),WV,NCJ,JMIN,LO)                                   RESU-567
      GO TO ( 11 , 72 , 73 ) , KLT                                      RESU-568
   72 IF (.NOT.LO(81)) GO TO 73                                         RESU-569
      EX(I4-1,KC)=EXCN                                                  RESU-570
      EX(I4,KC)=EX(2,KC)-EXCN                                           RESU-571
   73 IF (DTHE.EQ.0.D0) GO TO ( 49 , 49 , 60 ) , KLT                    RESU-572
      GO TO ( 74 , 75 , 76 ),KC                                         RESU-573
   74 KC=2                                                              RESU-574
      THETA=THETB-DTHE                                                  RESU-575
      GO TO 71                                                          RESU-576
   75 KC=3                                                              RESU-577
      THETA=THETB+DTHE                                                  RESU-578
      GO TO 71                                                          RESU-579
   76 IF (KLT.EQ.3) GO TO 79                                            RESU-580
      AA=EX(1,1)+EX(1,2)+EX(1,3)                                        RESU-581
      DO 78 K=2,I4                                                      RESU-582
      IF ((K.GT.I5).AND.(K.LE.I3)) GO TO 77                             RESU-583
      EX(K,1)=(EX(K,1)+EX(K,2)+EX(K,3))/3.D0                            RESU-584
      GO TO 78                                                          RESU-585
   77 EX(K,1)=(EX(K,1)*EX(1,1)+EX(K,2)*EX(1,2)+EX(K,3)*EX(1,3))/AA      RESU-586
   78 CONTINUE                                                          RESU-587
      GO TO 49                                                          RESU-588
   79 DO 80 K=JN1,JN2                                                   RESU-589
   80 EX(K-IS,1)=(EX(K-IS,1)+EX(K-IS,2)+EX(K-IS,3))/3.D0                RESU-590
      GO TO 60                                                          RESU-591
 1000 FORMAT ('<EXP.DAT.>',F10.2,1P,D20.8,0P,F10.2,2I5)                 RESU-592
 1001 FORMAT (3I5,5A4)                                                  RESU-593
 1002 FORMAT ('1',5X,18A4//)                                            RESU-594
 1003 FORMAT (//30X,'**********   STATE',I5,'    **********'//)         RESU-595
 1004 FORMAT (//42X,5A4//10X,'ANGLE',10X,'CALC. VAL.',11X,'EXP. VAL.',10RESU-596
     1X,'EXP. ERROR',10X,'COR. ERROR',13X,'CHI2'/)                      RESU-597
 1005 FORMAT (6X,F10.3,5D20.5)                                          RESU-598
 1006 FORMAT (//' WEIGHT IN THE CHI2',8X,D15.5/' EXPERIMENTAL NORMALISATRESU-599
     1ION',D15.5/' NORMALISATION ERROR',7X,D15.5/' CALCULATED NORMALISATRESU-600
     2ION',2X,D15.5//' ***** CHI2 =',D15.6,'   *****'/)                 RESU-601
 1007 FORMAT (/' ************ CHI2 **********',D20.10//)                RESU-602
 1008 FORMAT ('<ANG.DIS.>',F10.2,1P,D20.8,0P,F10.2,2I5)                 RESU-603
 1009 FORMAT ('<CROSS-S.>',F10.2,1P,D20.8,0P,F10.2,2I5)                 RESU-604
 1010 FORMAT ('<INE.C.S.>',F10.2,1P,D20.8,0P,F10.2,2I5)                 RESU-605
 1011 FORMAT (6X,'==> TOTAL CROSS SECTION =',F14.6,' MILLIBARNS.')      RESU-606
 1012 FORMAT (1P,D12.5)                                                 RESU-607
 1013 FORMAT (' TOTAL REACTION CROSS SECTION =',F14.6,' MILLIBARNS.')   RESU-608
 1014 FORMAT (' TOTAL REACTION CROSS SECTION =',F14.6,' MILLIBARNS',11X,RESU-609
     1'( INCLUDING COMPOUND ELASTIC ).')                                RESU-610
 1015 FORMAT (' TOTAL DIRECT REACTION CROSS SECTION =',F14.6,' MILLIBARNRESU-611
     1S',4X,'( WITHOUT COMPOUND ELASTIC ).')                            RESU-612
 1016 FORMAT (8X,'FISSION CROSS SECTION =',F14.6,' MILLIBARNS.')        RESU-613
 1017 FORMAT (6X,'GAMMA RAY CROSS SECTION =',F14.6,' MILLIBARNS.')      RESU-614
 1018 FORMAT (/'  COMPOUND CROSS SECTION FOR SCATTERING TO COUPLED LEVELRESU-615
     1S.'/'  ======================================================='//'RESU-616
     2   LEVEL     SPIN       ENERGY       CROSS SECTION'/)             RESU-617
 1019 FORMAT (I5,F9.1,A1,'   AT',F10.5,' MEV',F14.5)                    RESU-618
 1020 FORMAT (/' SUM OF COMPOUND CONTRIBUTIONS',F14.6,' MILLIBARNS.')   RESU-619
 1021 FORMAT (/'  COMPOUND CROSS SECTION FOR SCATTERING TO UNCOUPLED LEVRESU-620
     1ELS WITH ANGULAR DISTRIBUTION.'/'  ===============================RESU-621
     2===================================================='//'   LEVEL  RESU-622
     3   SPIN       ENERGY       CROSS SECTION'/)                       RESU-623
 1022 FORMAT (/'  COMPOUND CROSS SECTION FOR SCATTERING TO LEVELS WITHOURESU-624
     1T ANGULAR DISTRIBUTION.'/'  ======================================RESU-625
     2======================================'//'   LEVEL     SPIN       RESU-626
     3ENERGY       CROSS SECTION'/)                                     RESU-627
 1023 FORMAT (/' CONTINUUM',I4/' TOTAL COMPOUND REACTION CROSS SECTION =RESU-628
     1',F14.6,' MILLIBARNS.'//8X,' ENERGY     STEP     CONTRIBUTION  VALRESU-629
     2UE (MILLIBARNS/MEV)'/)                                            RESU-630
 1024 FORMAT (2X,I3,2F10.5,2F15.5)                                      RESU-631
 1025 FORMAT (/' TOTAL SUM OF COMPOUND CONTRIBUTIONS',F14.6,' MILLIBARNSRESU-632
     1.')                                                               RESU-633
 1026 FORMAT (/' ELASTIC SCATTERING ON THE TARGET STATE OF SPIN =',F5.1,RESU-634
     1A1)                                                               RESU-635
 1027 FORMAT (//' CLOSED CHANNEL FOR THE TARGET STATE',I3)              RESU-636
 1028 FORMAT (I5,F5.1,A1,3X,'0',5X,'CLOSED CHANNEL.')                   RESU-637
 1029 FORMAT ('0.D0',8X,I3)                                             RESU-638
 1030 FORMAT (/' INELASTIC SCATTERING TO THE TARGET STATE',I3,'  SPIN ='RESU-639
     1,F5.1,A1)                                                         RESU-640
 1031 FORMAT (/6X,'INELASTIC CROSS SECTION =',F14.6,' MILLIBARNS.')     RESU-641
 1032 FORMAT (/'  TOTAL ELASTIC CROSS SECTION =',F14.6,' MILLIBARNS.')  RESU-642
 1033 FORMAT (9X,'DIRECT CROSS SECTION =',F14.6,' MILLIBARNS.')         RESU-643
 1034 FORMAT (1P,D12.5,I3)                                              RESU-644
 1035 FORMAT (7X,'COMPOUND CROSS SECTION =',F14.6,' MILLIBARNS.')       RESU-645
 1036 FORMAT (1X,F10.3,D16.5,2X,5F18.7)                                 RESU-646
 1037 FORMAT (29X,5F18.7)                                               RESU-647
 1038 FORMAT (I3,1P,2D12.5,5X,4A4,A2)                                   RESU-648
 1039 FORMAT (' -4',1P,2D12.5,5X,4A4,A2/' -5',2D12.5,5X,4A4,A2)         RESU-649
 1040 FORMAT (/' ANGULAR DISTRIBUTION OF COMPOUND SCATTERING ON LEVELS',RESU-650
     1I3,' TO',I3//6X,'ANGLE',6(I6,'-LEVEL',F5.1,A1))                   RESU-651
 1041 FORMAT (1X,F10.3,6F18.7)                                          RESU-652
 1042 FORMAT (/6X,'ANGLE ',6(4A4,A2))                                   RESU-653
 1043 FORMAT (I5,F5.1,A1,I4,I5)                                         RESU-654
 1044 FORMAT (30X,5(4A4,A2))                                            RESU-655
      END                                                               RESU-656
