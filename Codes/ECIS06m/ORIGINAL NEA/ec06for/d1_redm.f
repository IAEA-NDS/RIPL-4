C 27/02/07                                                      ECIS06  REDM-000
      SUBROUTINE REDM(IPI,NCOLL,NBETA,BETA,IPH,NVAR,VAR,FAC,IDT,LO,NIV,IREDM-001
     1Q,JQ,T,IT,IM)                                                     REDM-002
C REDUCED NUCLEAR MATRIX ELEMENTS.                                      REDM-003
C INPUT:     IPI(J,*):PARITY OF NUCLEAR STATES FOR J=1,                 REDM-004
C                     MULTIPLICITY OF THE PARTICLE FOR J=2,             REDM-005
C                     MULTIPLICITY OF THE TARGET FOR J=3.               REDM-006
C            NCOLL:   NUMBER OF COUPLED STATES                          REDM-007
C            NBETA:   QUANTUM NUMBERS OF DEFORMATIONS.                  REDM-008
C            BETA:    DEFORMATIONS, EQUIVALENT BY CALL WITH NBETA.      REDM-009
C            IPH:     DESCRIPTION OF VIBRATIONAL MODEL (SEE VIBM).      REDM-010
C            NVAR,VAR:EQUIVALENT BY CALL, PARAMETERS FOR SOME MODELS.   REDM-011
C            NVA:     NUMBER OF THESE PARAMETERS.                       REDM-012
C            FAC:     TABLE AND NUMBER OF LOGARITHMS OF FACTORIALS.     REDM-013
C            IDT:     LENGTH AVAILABLE IN THIS SUBROUTINE.              REDM-014
C            LO(I):   LOGICAL CONTROLS:                                 REDM-015
C               LO(1)  =.TRUE. ROTATIONAL MODEL-(.F.:VIBRATIONAL MODEL).REDM-016
C               LO(3)  =.TRUE. ANHARMONIC VIBRATIONAL OR ASYMMETRIC     REDM-017
C                              ROTATIONAL MODEL.                        REDM-018
C               LO(7)  =.TRUE. MATRIX ELEMENT AND FORM FACTORS READ.    REDM-019
C               LO(11) =.TRUE. DEFORMED COULOMB POTENTIAL.              REDM-020
C               LO(12) =.TRUE. DEFORMED IMAGINARY POTENTIAL.            REDM-021
C               LO(13) =.TRUE. DEFORMED REAL SPIN-ORBIT OR TENSOR.      REDM-022
C               LO(14) =.TRUE. DEFORMED IMAGINARY SPIN-ORBIT OR TENSOR. REDM-023
C               LO(15) =.TRUE. REDUCED NUCLEAR MATRIX ELEMENTS READ.    REDM-024
C               LO(19) =.TRUE. DEFORMED COULOMB SPIN-ORBIT POTENTIAL.   REDM-025
C               LO(20) =.TRUE. DISPERSION RELATIONS FOR TRANSITION      REDM-026
C                              FORM-FACTORS.                            REDM-027
C               LO(52) =.TRUE. REDUCED NUCLEAR MATRIX ELEMENTS OUTPUT.  REDM-028
C               LO(61) =.TRUE. REDUCED NUCLEAR MATRIX ELEMENTS WRITTEN  REDM-029
C                              ON FILE 61.                              REDM-030
C               LO(100)=.TRUE. DIRAC EQUATION.                          REDM-031
C               LO(101)=.TRUE. THERE IS A REAL SPIN-ORBIT POTENTIAL.    REDM-032
C               LO(103)=.TRUE. THERE IS A COULOMB SPIN-ORBIT POTENTIAL. REDM-033
C               LO(117)=.TRUE. FOR ALL CALCULATIONS EXCEPT THE FIRST.   REDM-034
C               LO(120)=.TRUE. OUTPUT AND LAST CALCULATION BEST ONE.    REDM-035
C               LO(121)=.TRUE. OPTICAL MODEL WITHOUT COUPLING.          REDM-036
C               LO(133)=.TRUE. STORE SCALAR AND COULOMB POTENTIAL       REDM-037
C                              INDEPENDENTLY.                           REDM-038
C OUTPUT:    NIV:     NUMBER AND REFERENCES OF INTERACTION FORM FACTORS REDM-039
C                     BETWEEN EACH CHANNELS. NIV(I1,I2,K): FIRST I OF   REDM-040
C                     T(3,I) FOR THE PAIR OF NUCLEAR STATES I1,I2 FOR   REDM-041
C                     K=1 AND LAST ONE FOR K=2.                         REDM-042
C            IQ,JQ,T: EQUIVALENT BY CALL OF IQ(6,I), JQ(I) AND T(3,I).  REDM-043
C            IT:      NUMBER OF NUCLEAR MATRIX ELEMENTS IN IQ,T. FOR A  REDM-044
C                     GIVEN VALUE OF I:                                 REDM-045
C                     IQ(1,I): REFERENCE TO TABLE OF FORM FACTORS,      REDM-046
C                     IQ(2,I): REFERENCE TO TABLE OF ANGULAR MOMENTA,   REDM-047
C                     IQ(3,I): ADDRESS OF THE ASSOCIATED SPIN-ORBIT FORMREDM-048
C                              FACTOR OR 0,                             REDM-049
C                     IQ(4,I): UNUSED,                                  REDM-050
C                     T(3,I):  REDUCED NUCLEAR MATRIX ELEMENT MULTIPLIEDREDM-051
C                              BY (-)**(L/2) WHERE L IS THE TRANSFERRED REDM-052
C                              ANGULAR MOMENTUM                         REDM-053
C            IM:      NUMBER OF ANGULAR MOMENTA.                        REDM-054
C***********************************************************************REDM-055
C BEYOND IQ(6,IT) ARE STORED IN JQ(*) THE TABLES:                       REDM-056
C            IVQ:     3*"IM" DATA, DESCRIPTION OF ANGULAR MOMENTA:      REDM-057
C                     1) L  TRANSFERRED ANGULAR MOMENTUM,               REDM-058
C                     2) 2*S WHERE S IS THE TRANSFER OF SPIN,           REDM-059
C                     3) 2*J WHERE J IS THE TRANSFER OF TOTAL           REDM-060
C                            SPIN OR 0 IF S=0.                          REDM-061
C            IVY:     7*"INTC", TABLE OF FORM FACTOR INDEPENDENT OF     REDM-062
C                     DISPERSION RELATIONS:                             REDM-063
C                     1) FORM FACTOR CONTROL NUMBER,                    REDM-064
C                     2) REFERENCE TO THE TABLE OF MULTIPOLES,          REDM-065
C                     3) 0 OR ADDRESS OF FIRST SPIN-ORBIT FORM FACTOR,  REDM-066
C                     4) 0 OR ADDRESS OF TEMPORARY COULOMB FORM FACTOR, REDM-067
C                     5) 0 OR ADDRESS OF TEMPORARY COULOMB SPIN-ORBIT,  REDM-068
C                     6) 0 OR ADDRESS OF CORRECTION TERM (POSITIVE FOR  REDM-069
C                        CORRECTED TERM, NEGATIVE FOR CORRECTION TERM), REDM-070
C                     7) ANGULAR MOMENTUM.                              REDM-071
C            IVZ:     4*"INTV" TABLE FORM FACTOR DEPENDENT OF DISPERSIONREDM-072
C                     RELATIONS:                                        REDM-073
C                     1) ADDRESS OF COMPUTATION (WITHOUT DISPERSION),   REDM-074
C                     2) FORM FACTOR CONTROL NUMBER FOR DISPERSION,     REDM-075
C                     3) 0 OR ADDRESS OF THE FIRST SPIN-ORBIT FORM      REDM-076
C                        FACTOR,                                        REDM-077
C                     4) ANGULAR MOMENTUM.                              REDM-078
C "INTC" AND "INTC" ARE IN COMMON /POTE1/. THE ADDRESSES OF THE TABLES  REDM-079
C "IVQ", "IVY" AND "IVZ" ARE DEFINED AFTER RETURN FROM THIS SUBROUTINE. REDM-080
C IF LO(15)=.TRUE. THE REDUCED MATRIX ELEMENTS ARE READ AT THE FIRST    REDM-081
C CALL AND THE PARAMETERS "VAR" ARE NOT USED - IN ANY CASE WHERE "VAR"  REDM-082
C ARE NOT USED,THE SUBROUTINE IS SKIPPED.                               REDM-083
C IF LO(61)=.TRUE. THE MATRIX ELEMENTS ARE PUNCHED AT THE LAST CALL     REDM-084
C OF A SEARCH IN THIS SUBROUTINE.                                       REDM-085
C***********************************************************************REDM-086
C                                                                       REDM-087
C THE COMMON /POTE2/ IS USED IN CALC, LECL, REDM, CAL1, POTE, QUAN,     REDM-088
C                               MTCH, INTI, INSH, INSI, INCH AND REST.  REDM-089
C THE COMMON /POTE1/ IS USED IN CALC, REDM, EXTP, POTE, ROTP, STDP,     REDM-090
C                               FOLD AND REST.                          REDM-091
C FOR THE COMMON  /COUPL/ SEE CALC.                                     REDM-092
C                                                                       REDM-093
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /COUPL/:                     REDM-094
C  NBT1:      NUMBER OF PHONONS (VIBRATIONS).                           REDM-095
C  NPP:       NUMBER OF OPTICAL POTENTIALS.                             REDM-096
C  NVA:       NUMBER OF NUCLEAR PARAMETERS.                             REDM-097
C   USED:     NBT1,NPP,NVA.                                             REDM-098
C                                                                       REDM-099
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /POTE1/:                     REDM-100
C  ITX(16):   STARTING ADDRESS OF DIFFERENT FORM FACTORS.               REDM-101
C             FOR SCHROEDINGER EQUATION, ITX(I)+1 IS THE STARTING       REDM-102
C             ADDRESS OF THE FORM FACTOR READ IN EXTP WITH ITYP=I       REDM-103
C             (POTENTIAL FOR I=1 TO 8, TRANSITION FOR I=9 TO 16).       REDM-104
C             FOR DIRAC EQUATIONS, ITX(1)=0,                            REDM-105
C             ITX(2)+1=ADDRESS OFF FIRST TRANSITION FORM FACTOR,        REDM-106
C             ITX(7)=ADDRESS OF LAST TRANSITION FORM FACTOR,            REDM-107
C             ITX(3)=ADDRESS OF LAST TEMPORARY CENTRAL POTENTIAL,       REDM-108
C             ITX(4)=ITX(7)-24,ITX(5)=ITX(3)-11,ITX(6)=ITX(2)-4.        REDM-109
C             ALL ARE USED FOR SCHROEDINGER, THE FIRST 8 FOR DIRAC.     REDM-110
C  IMAX:      MAXIMUM ANGULAR MOMENTUM.                                 REDM-111
C  INTC:      NUMBER OF FORM FACTORS WITHOUT DEFORMED SPIN-ORBIT        REDM-112
C             INCLUDING CORRECTION TERMS.                               REDM-113
C  INLS:      NUMBER OF SPIN-ORBIT FORM FACTORS NOT TAKING INTO ACCOUNT REDM-114
C             MULTIPLICATION BY 2.                                      REDM-115
C  INVC:      NUMBER OF COULOMB TRANSITION FORM FACTORS.                REDM-116
C  INVD:      IDEM FOR COULOMB SPIN-ORBIT.                              REDM-117
C  ITXM:      TOTAL NUMBER OF FORM FACTORS.                             REDM-118
C   DEFINED:  ITX,IMAX,INTC,INLS,INVC,INVD,ITXM.                        REDM-119
C                                                                       REDM-120
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /POTE2/:                     REDM-121
C  ITY(1):    STARTING ADDRESS OF REAL CENTRAL POTENTIAL (IT IS 0).     REDM-122
C  ITY(2):    STARTING ADDRESS OF IMAGINARY CENTRAL POTENTIAL.          REDM-123
C  ITY(3):    STARTING ADDRESS OF REAL SPIN-ORBIT POTENTIAL.            REDM-124
C  ITY(4):    STARTING ADDRESS OF IMAGINARY SPIN-ORBIT POTENTIAL.       REDM-125
C  ITY(5):    STARTING ADDRESS OF REAL CENTRAL TRANSITION.              REDM-126
C  ITY(6):    STARTING ADDRESS OF IMAGINARY CENTRAL TRANSITION.         REDM-127
C  ITY(7):    STARTING ADDRESS OF REAL SPIN-ORBIT TRANSITION.           REDM-128
C  ITY(8):    STARTING ADDRESS OF IMAGINARY SPIN-ORBIT TRANSITION.      REDM-129
C  ITY(9):    STARTING ADDRESS OF COULOMB CENTRAL POTENTIAL.            REDM-130
C  ITY(10):   STARTING ADDRESS OF COULOMB SPIN-ORBIT POTENTIAL.         REDM-131
C  ITY(11):   STARTING ADDRESS OF COULOMB CENTRAL TRANSITION.           REDM-132
C  ITY(12):   STARTING ADDRESS OF COULOMB SPIN-ORBIT TRANSITION.        REDM-133
C        ITY(2)=14*NCOLL AND ITY(5)=0 ONLY ARE USED FOR DIRAC EQUATIONS.REDM-134
C  INVT:      NUMBER OF TRANSITION FORM FACTORS WITHOUT SPIN-ORBIT.     REDM-135
C  INTV:      SAME AS INVT, TAKING INTO ACCOUNT DISPERSION.             REDM-136
C  INSL:      NUMBER OF SPIN-ORBIT FORM TRANSITION FACTORS NOT TAKING   REDM-137
C                  INTO ACCOUNT MULTIPLICATION BY 2.                    REDM-138
C  NPX:       NUMBER OF POTENTIALS TAKING INTO ACCOUNT DISPERSION.      REDM-139
C   DEFINED:  ITY,INVT,INTV,INSL.                                       REDM-140
C   USED:     ITY,INVT,INTV,INSL,NPX.                                   REDM-141
C                                                                       REDM-142
C***********************************************************************REDM-143
      IMPLICIT REAL*8 (A-H,O-Z)                                         REDM-144
      LOGICAL LO(150),LT                                                REDM-145
      DIMENSION NIV(NCOLL,NCOLL,*),IQ(6,*),JQ(*),T(3,*),IPI(11,*),NBETA(REDM-146
     118,*),BETA(9,*),IPH(2,*),NVAR(2,*),VAR(*),FAC(*)                  REDM-147
      COMMON /COUPL/ IQM,IQMAX,NBT1,NFA,NPP,NSPIN,NVA                   REDM-148
      COMMON /INOUT/ MR,MW,MS                                           REDM-149
      COMMON /POTE1/ ITX(16),IMAX,INTC,INLS,INVC,INVD,ITXM              REDM-150
      COMMON /POTE2/ ITY(12),INVT,INTV,INSL,NPX                         REDM-151
      IF (LO(120)) GO TO 42                                             REDM-152
      IF (LO(7).OR.LO(15)) GO TO 3                                      REDM-153
C STANDARD COMPUTATION OF THE REDUCED MATRIX ELEMENTS IN T(3,I)         REDM-154
C WITH FORM FACTOR IDENTIFICATION IN IQ(1,I), MULTIPOLE IN IQ(2,I)      REDM-155
C AND 0 OR 1 IN IQ(3,I) IF SPIN-ORBIT IS NOT OR IS DEFORMED.            REDM-156
C THEY ARE CALLED AT THE FIRST RUN OR IF "VAR" HAS BEEN CHANGED.        REDM-157
      IF (LO(1)) GO TO 1                                                REDM-158
      CALL VIBM(NIV,IQ,T,IPI,NCOLL,IT,IPH,NVAR,VAR,NBETA,FAC,IDT,LO)    REDM-159
      GO TO 17                                                          REDM-160
    1 IF (LO(3)) GO TO 2                                                REDM-161
      CALL ROTM(NIV,IQ,T,IPI,NCOLL,IPH,NBETA,NVAR,VAR,IT,FAC,IDT,LO)    REDM-162
      GO TO 17                                                          REDM-163
    2 CALL ROAM(NIV,IQ,T,IPI,NCOLL,IT,BETA,IPH,VAR,VAR(NVA+1),FAC,IDT,LOREDM-164
     1)                                                                 REDM-165
      GO TO 17                                                          REDM-166
C INPUT OF NUCLEAR MATRIX ELEMENTS FROM CARDS.                          REDM-167
    3 IF (LO(117)) GO TO 42                                             REDM-168
      IT=0                                                              REDM-169
      I=0                                                               REDM-170
      DO 16 I1=1,NCOLL                                                  REDM-171
      DO 15 I2=I1,NCOLL                                                 REDM-172
      NIV(I1,I2,1)=IT+1                                                 REDM-173
      READ (MR,1000) J1,J2,K                                            REDM-174
      WRITE (MW,1001) J1,J2,K                                           REDM-175
      IF ((J1.NE.I1).OR.(J2.NE.I2)) GO TO 46                            REDM-176
      IF (K.EQ.0) GO TO 15                                              REDM-177
      IJL=IABS(IPI(3,I1)-IPI(3,I2))                                     REDM-178
      JIL=(IPI(3,I1)+IPI(3,I2))                                         REDM-179
      IJS=IABS(IPI(2,I1)-IPI(2,I2))                                     REDM-180
      JIS=IPI(2,I1)+IPI(2,I2)                                           REDM-181
      IF (2*(IT+K).GT.IDT) CALL MEMO('REDM',IDT,2*(IT+K))               REDM-182
      DO 14 K1=1,K                                                      REDM-183
      IT=IT+1                                                           REDM-184
      READ (MR,1002) J,IQ(2,IT),NN,M,T(3,IT)                            REDM-185
      IF (NN.EQ.-1.AND.(.NOT.LO(11))) GO TO 47                          REDM-186
      N=NN                                                              REDM-187
      IF (NN.EQ.-1) N=0                                                 REDM-188
      IF (MOD(IQ(2,IT)+IPI(1,I1)+IPI(1,I2),2).NE.0) GO TO 48            REDM-189
      IF ((N.GT.JIS).OR.(N.LT.IJS).OR.(MOD(JIS+N,2).NE.0)) GO TO 49     REDM-190
      IF (LO(7)) GO TO 7                                                REDM-191
      IF (LO(1)) GO TO 5                                                REDM-192
      IF (LO(3)) GO TO 4                                                REDM-193
      IQP=J/(NBT1+1)                                                    REDM-194
      IPQ=MOD(J,NBT1+1)                                                 REDM-195
      IF ((IPQ.LE.0).OR.(MAX0(IPQ,IQP).GT.NBT1).OR.(IQP.LT.0)) GO TO 50 REDM-196
      GO TO 11                                                          REDM-197
    4 IF (J.LT.0.OR.J.GT.3) GO TO 51                                    REDM-198
      GO TO 11                                                          REDM-199
    5 IQP=J/1000                                                        REDM-200
      IF (LO(3)) GO TO 6                                                REDM-201
      IF (IQP.NE.1+IQ(2,IT)) GO TO 52                                   REDM-202
      GO TO 11                                                          REDM-203
    6 IQY=IQ(2,IT)/2+1                                                  REDM-204
      IQZ=(IQY*(IQY-1))/2                                               REDM-205
      IF (IQP.LE.IQZ.OR.IQP.GT.IQZ+IQY) GO TO 53                        REDM-206
      GO TO 11                                                          REDM-207
    7 IF (J.EQ.0) GO TO 10                                              REDM-208
      IF (IABS(J)-I-1) 8 , 10 , 54                                      REDM-209
    8 IF (IT.EQ.1) GO TO 11                                             REDM-210
      DO 9 II=2,IT                                                      REDM-211
      IF (J.EQ.IQ(1,II-1)) GO TO 11                                     REDM-212
    9 CONTINUE                                                          REDM-213
      GO TO 55                                                          REDM-214
   10 I=I+1                                                             REDM-215
      IF (J.EQ.0) J=I                                                   REDM-216
   11 IQ(1,IT)=J                                                        REDM-217
      IF (N.EQ.0) GO TO 12                                              REDM-218
      IF (2*IQ(2,IT).LT.IABS(N-M).OR.2*IQ(2,IT).GT.N+M.OR.MOD(N+M,2).NE.REDM-219
     10) GO TO 56                                                       REDM-220
      NM=M                                                              REDM-221
      IQ(3,IT)=10*M+10000*N                                             REDM-222
      GO TO 13                                                          REDM-223
   12 IF (M.LT.0) GO TO 57                                              REDM-224
      IQ(3,IT)=10*M                                                     REDM-225
      NM=2*IQ(2,IT)                                                     REDM-226
   13 IF (NM.LT.IJL.OR.NM.GT.JIL) GO TO 58                              REDM-227
      IF (NN.EQ.-1) IQ(3,IT)=IQ(3,IT)+1                                 REDM-228
   14 WRITE (MW,1003) IQ(1,IT),IQ(2,IT),NN,M,T(3,IT)                    REDM-229
   15 NIV(I1,I2,2)=IT                                                   REDM-230
   16 CONTINUE                                                          REDM-231
   17 K=NIV(NCOLL,NCOLL,2)                                              REDM-232
      DO 19 I1=1,NCOLL                                                  REDM-233
      DO 18 I2=I1,NCOLL                                                 REDM-234
      NIV(I2,I1,1)=NIV(I1,I2,1)                                         REDM-235
      NIV(I2,I1,2)=NIV(I1,I2,2)                                         REDM-236
   18 CONTINUE                                                          REDM-237
   19 CONTINUE                                                          REDM-238
C CHECK THAT THERE IS ENOUGH PLACE.                                     REDM-239
      IF (7*IT.GT.IDT) CALL MEMO('REDM',IDT,7*IT)                       REDM-240
      IM=0                                                              REDM-241
      INVT=0                                                            REDM-242
      INTV=0                                                            REDM-243
      INLS=0                                                            REDM-244
      INSL=0                                                            REDM-245
      INVC=0                                                            REDM-246
      INVD=0                                                            REDM-247
      INTC=0                                                            REDM-248
      IMAX=0                                                            REDM-249
C TABLE OF MULTIPOLES.                                                  REDM-250
      IF (IT.EQ.0) GO TO 40                                             REDM-251
      ITN=6*IT                                                          REDM-252
      DO 23 J=1,IT                                                      REDM-253
      IS=IQ(3,J)/10000                                                  REDM-254
      IJ=MOD(IQ(3,J),1000)/10                                           REDM-255
      IL=MOD(IQ(3,J),10)                                                REDM-256
      IF (IS.EQ.0) IJ=2*IQ(2,J)                                         REDM-257
      IF (IJ.EQ.0) IS=2*IQ(2,J)                                         REDM-258
      IF (IL.EQ.1) IS=-2                                                REDM-259
      IF (IS.NE.0) IQ(3,J)=0                                            REDM-260
      IF (IM.EQ.0) GO TO 21                                             REDM-261
      DO 20 K=1,IM,3                                                    REDM-262
      IF ((JQ(ITN+K).EQ.IQ(2,J)).AND.(JQ(ITN+K+1).EQ.IS).AND.(JQ(ITN+K+2REDM-263
     1).EQ.IJ)) GO TO 22                                                REDM-264
   20 CONTINUE                                                          REDM-265
   21 JQ(ITN+IM+1)=IQ(2,J)                                              REDM-266
      JQ(ITN+IM+2)=IS                                                   REDM-267
      JQ(ITN+IM+3)=IJ                                                   REDM-268
      IMAX=MAX0(IMAX,IQ(2,J))                                           REDM-269
      IM=IM+3                                                           REDM-270
      IQ(2,J)=IM/3                                                      REDM-271
      GO TO 23                                                          REDM-272
   22 IQ(2,J)=1+K/3                                                     REDM-273
   23 CONTINUE                                                          REDM-274
      ITM=ITN+IM+MOD(IM,2)                                              REDM-275
C TABLE OF FORM FACTORS.                                                REDM-276
C LT IS TRUE IF THE MULTIPOLE ORDER MATTERS FOR FORM FACTORS.           REDM-277
      LT=(LO(1).OR.LO(11).OR.LO(17).OR.LO(19)).AND.(.NOT.LO(7))         REDM-278
      DO 26 J=1,IT                                                      REDM-279
      IF (INVT.EQ.0) GO TO 25                                           REDM-280
      DO 24 L=1,INVT,7                                                  REDM-281
      IF (JQ(ITM+L).NE.IQ(1,J)) GO TO 24                                REDM-282
      IF (LT.AND.(JQ(ITM+L+1).NE.IQ(2,J))) GO TO 24                     REDM-283
      IQ(1,J)=1+L/7                                                     REDM-284
      JQ(ITM+L+2)=MAX0(IQ(3,J),JQ(ITM+L+2))                             REDM-285
      GO TO 26                                                          REDM-286
   24 CONTINUE                                                          REDM-287
   25 K=IQ(2,J)                                                         REDM-288
      K1=JQ(ITN+3*K-2)                                                  REDM-289
      JQ(ITM+INVT+1)=IQ(1,J)                                            REDM-290
      JQ(ITM+INVT+2)=K                                                  REDM-291
      JQ(ITM+INVT+3)=IQ(3,J)                                            REDM-292
      JQ(ITM+INVT+4)=0                                                  REDM-293
      JQ(ITM+INVT+5)=0                                                  REDM-294
      JQ(ITM+INVT+6)=0                                                  REDM-295
      JQ(ITM+INVT+7)=K1                                                 REDM-296
      IF ((IQ(1,J).LT.0).OR.((IQ(1,J).GT.1000).AND.(K1.LE.1).AND.(MOD(IQREDM-297
     1(1,J),1000).NE.0))) JQ(ITM+INVT+6)=1                              REDM-298
      IF (LO(11).AND.(MOD(JQ(ITN+3*K-1),2).EQ.0)) JQ(ITM+INVT+4)=1      REDM-299
      IF (LO(19)) JQ(ITM+INVT+5)=JQ(ITM+INVT+3)                         REDM-300
      INVT=INVT+7                                                       REDM-301
      IQ(1,J)=INVT/7                                                    REDM-302
   26 CONTINUE                                                          REDM-303
C SEARCH ON CORRECTION TERMS.                                           REDM-304
      INTC=INVT                                                         REDM-305
      DO 28 J=1,INVT,7                                                  REDM-306
      IF (JQ(ITM+J+5).EQ.0) GO TO 28                                    REDM-307
      K1=ITM+J-1                                                        REDM-308
      DO 27 L=1,7                                                       REDM-309
   27 JQ(ITM+INTC+L)=JQ(K1+L)                                           REDM-310
      JQ(ITM+INTC+6)=-(1+J/7)                                           REDM-311
      INTC=INTC+7                                                       REDM-312
      JQ(K1+6)=INTC/7                                                   REDM-313
   28 CONTINUE                                                          REDM-314
      DO 30 J=1,INTC,7                                                  REDM-315
      IF (JQ(ITM+J+3).EQ.0) GO TO 29                                    REDM-316
      INVC=INVC+1                                                       REDM-317
      JQ(ITM+J+3)=INVC                                                  REDM-318
   29 IF ((JQ(ITM+J+2).EQ.0).AND.(JQ(ITM+J+4).EQ.0)) GO TO 30           REDM-319
      INLS=INLS+1                                                       REDM-320
      JQ(ITM+J+2)=INLS                                                  REDM-321
      IF (JQ(ITM+J+4).EQ.0) GO TO 30                                    REDM-322
      INVD=INVD+1                                                       REDM-323
      JQ(ITM+J+4)=INVD                                                  REDM-324
   30 CONTINUE                                                          REDM-325
      ITP=ITM+INTC+MOD(INTC,2)                                          REDM-326
      NTC=INVT/7                                                        REDM-327
      DO 31 J=1,NTC                                                     REDM-328
      JQ(ITP+4*J-3)=J                                                   REDM-329
      JQ(ITP+4*J-2)=0                                                   REDM-330
      JQ(ITP+4*J-1)=JQ(ITM+7*J-4)                                       REDM-331
      JQ(ITP+4*J)=JQ(ITM+7*J)                                           REDM-332
   31 CONTINUE                                                          REDM-333
      INTV=4*NTC                                                        REDM-334
      IF (.NOT.LO(20)) GO TO 37                                         REDM-335
      DO 36 I1=1,NCOLL                                                  REDM-336
      DO 35 I2=I1,NCOLL                                                 REDM-337
      K1=NIV(I1,I2,1)                                                   REDM-338
      K2=NIV(I1,I2,2)                                                   REDM-339
      IF (K1.GT.K2) GO TO 35                                            REDM-340
      KK=I1*(NCOLL+1)+I2                                                REDM-341
      DO 34 K=K1,K2                                                     REDM-342
      J=IQ(1,K)                                                         REDM-343
      IF (JQ(ITP+4*J-2).NE.0) GO TO 32                                  REDM-344
      JQ(ITP+4*J-2)=KK                                                  REDM-345
      GO TO 34                                                          REDM-346
   32 DO 33 L=1,4                                                       REDM-347
   33 JQ(ITP+INTV+L)=JQ(ITP+4*J+L-4)                                    REDM-348
      JQ(ITP+INTV+2)=KK                                                 REDM-349
      INTV=INTV+4                                                       REDM-350
      IQ(1,K)=INTV/4                                                    REDM-351
   34 CONTINUE                                                          REDM-352
   35 CONTINUE                                                          REDM-353
   36 CONTINUE                                                          REDM-354
   37 IM=IM/3                                                           REDM-355
      INTV=INTV/4                                                       REDM-356
      INVT=INVT/7                                                       REDM-357
      INTC=INTC/7                                                       REDM-358
      INSL=0                                                            REDM-359
      DO 39 J=1,INTV                                                    REDM-360
      L=JQ(ITP+4*J-1)                                                   REDM-361
      INSL=MAX0(INSL,L)                                                 REDM-362
      DO 38 K=1,IT                                                      REDM-363
      IF (IQ(1,K).EQ.J) IQ(3,K)=L                                       REDM-364
   38 CONTINUE                                                          REDM-365
   39 CONTINUE                                                          REDM-366
C NUMBER OF REAL AND IMAGINARY FORM FACTORS.                            REDM-367
   40 NV=INTC                                                           REDM-368
      IF (LO(13).OR.LO(19)) NV=NV+2*INLS                                REDM-369
      ITY(1)=0                                                          REDM-370
      LO(121)=IM.EQ.0                                                   REDM-371
      IF (LO(100)) GO TO 41                                             REDM-372
      ITY(3)=NPX                                                        REDM-373
      ITY(5)=ITY(3)+NPX                                                 REDM-374
      IF (.NOT.(LO(101).OR.LO(103))) ITY(5)=ITY(3)                      REDM-375
      ITY(7)=ITY(5)+INTV                                                REDM-376
      ITY(2)=ITY(7)+2*INSL                                              REDM-377
      ITY(4)=ITY(2)+NPX                                                 REDM-378
      ITY(6)=ITY(4)+NPX                                                 REDM-379
      IF (.NOT.(LO(101).OR.LO(103))) ITY(6)=ITY(4)                      REDM-380
      ITY(8)=ITY(6)+INTV                                                REDM-381
      IF (.NOT.LO(12)) ITY(8)=ITY(6)                                    REDM-382
      ITX1=ITY(8)+2*INSL                                                REDM-383
      IF (.NOT.LO(14)) ITX1=ITY(8)                                      REDM-384
      INTXC=ITX1                                                        REDM-385
      IF (.NOT.LO(133)) INTXC=0                                         REDM-386
      ITY(9)=ITY(1)+INTXC                                               REDM-387
      ITY(10)=ITY(3)+INTXC                                              REDM-388
      ITY(11)=ITY(5)+INTXC                                              REDM-389
      ITY(12)=ITY(7)+INTXC                                              REDM-390
      IF (.NOT.(LO(103).OR.LO(11).OR.LO(19))) ITY(11)=ITY(10)           REDM-391
      IF (.NOT.(LO(11).OR.LO(19))) ITY(12)=ITY(11)                      REDM-392
      ITX2=ITY(12)+2*INSL                                               REDM-393
      IF (.NOT.LO(19)) ITX2=ITY(12)                                     REDM-394
      ITX(1)=MAX0(ITX1,ITX2)                                            REDM-395
      ITX(5)=ITX(1)+NPP                                                 REDM-396
      ITX(9)=ITX(5)+NPP                                                 REDM-397
      IF (.NOT.(LO(101).OR.LO(103))) ITX(9)=ITX(5)                      REDM-398
      ITX(13)=ITX(9)+INTC                                               REDM-399
      ITX(2)=ITX(13)+2*INLS                                             REDM-400
      ITX(6)=ITX(2)+NPP                                                 REDM-401
      ITX(10)=ITX(6)+NPP                                                REDM-402
      IF (.NOT.(LO(101).OR.LO(103))) ITX(10)=ITX(6)                     REDM-403
      ITX(14)=ITX(10)+INTC                                              REDM-404
      ITX(7)=ITX(14)+2*INLS                                             REDM-405
      IF (.NOT.LO(14)) ITX(7)=ITX(14)                                   REDM-406
      IF (.NOT.LO(12)) ITX(7)=ITX(10)                                   REDM-407
      ITX(8)=ITX(7)+NPP                                                 REDM-408
      ITX(3)=ITX(8)+NPP                                                 REDM-409
      ITX(4)=ITX(3)+NPP                                                 REDM-410
      ITX(11)=ITX(4)+NPP                                                REDM-411
      ITX(12)=ITX(11)+INTC                                              REDM-412
      ITX(15)=ITX(12)+INTC                                              REDM-413
      IF (.NOT.LO(12)) ITX(15)=ITX(12)                                  REDM-414
      ITX(16)=ITX(15)+INVC                                              REDM-415
      ITXM=ITX(16)+2*INVD+INTC                                          REDM-416
      GO TO 42                                                          REDM-417
   41 ITY(5)=0                                                          REDM-418
      ITX(1)=0                                                          REDM-419
      ITX(2)=14*NCOLL                                                   REDM-420
      ITY(2)=ITX(2)                                                     REDM-421
      ITX(7)=ITX(2)+4*(INTV+INSL)                                       REDM-422
      ITX(3)=ITX(7)+24*NPP                                              REDM-423
      ITX(4)=ITX(7)-24                                                  REDM-424
      ITX(5)=ITX(3)-11                                                  REDM-425
      ITX(6)=ITX(2)-4                                                   REDM-426
      ITXM=ITX(3)+11*INTC                                               REDM-427
   42 IF (.NOT.(LO(52).OR.LO(61))) RETURN                               REDM-428
C REDUCED MATRIX ELEMENTS PUNCHED ON CARDS ON REQUEST.                  REDM-429
C OUTPUT OF REDUCED MATRIX ELEMENTS ON REQUEST.                         REDM-430
      IT=NIV(NCOLL,NCOLL,2)                                             REDM-431
      ITM=6*IT-3                                                        REDM-432
      ITN=ITM+3*IM-4+MOD(IM,2)                                          REDM-433
      IF (LO(52)) WRITE (MW,1004)                                       REDM-434
      DO 45 I1=1,NCOLL                                                  REDM-435
      DO 44 I2=I1,NCOLL                                                 REDM-436
      K1=NIV(I1,I2,1)                                                   REDM-437
      K2=NIV(I1,I2,2)                                                   REDM-438
      IF (LO(52)) WRITE (MW,1005) I1,I2,K1,K2                           REDM-439
      K=K2-K1+1                                                         REDM-440
      IF (LO(61)) WRITE (61,1000) I1,I2,K                               REDM-441
      IF (K.EQ.0) GO TO 44                                              REDM-442
      DO 43 K=K1,K2                                                     REDM-443
      IF (LO(52)) WRITE (MW,1006) K,(IQ(J,K),J=1,3),T(3,K)              REDM-444
      IF (.NOT.LO(61)) GO TO 43                                         REDM-445
      M=IQ(2,K)                                                         REDM-446
      NM=JQ(ITN+7*M+1)                                                  REDM-447
      NS=JQ(ITM+3*M+2)                                                  REDM-448
      NJ=JQ(ITM+3*M+3)                                                  REDM-449
      IF (NS.LE.0) NJ=0                                                 REDM-450
      IF (NS.EQ.-2) NS=-1                                               REDM-451
      IF (NJ.EQ.0.AND.IQ(3,K).NE.0) NJ=1                                REDM-452
      WRITE (61,1007) NM,JQ(ITM+3*M+1),NS,NJ,T(3,K)                     REDM-453
   43 CONTINUE                                                          REDM-454
   44 CONTINUE                                                          REDM-455
   45 CONTINUE                                                          REDM-456
      IF (IT.EQ.0.OR.(.NOT.LO(52))) RETURN                              REDM-457
      WRITE (MW,1008) (I,(JQ(ITM+3*I+J),J=1,3),I=1,IM)                  REDM-458
      WRITE (MW,1009) (I,(JQ(ITN+7*I+J),J=1,7),I=1,INTC)                REDM-459
      ITP=ITN+7*INTC+3+MOD(INTC,2)                                      REDM-460
      WRITE (MW,1010) (I,(JQ(ITP+4*I+J),J=1,4),I=1,INTV)                REDM-461
      RETURN                                                            REDM-462
   46 WRITE (MW,1011) J1,J2,I1,I2                                       REDM-463
      GO TO 59                                                          REDM-464
   47 WRITE (MW,1012) NN                                                REDM-465
      GO TO 59                                                          REDM-466
   48 WRITE (MW,1013) IQ(2,IT),IPI(1,I1),IPI(1,I2)                      REDM-467
      GO TO 59                                                          REDM-468
   49 WRITE (MW,1014)                                                   REDM-469
      GO TO 59                                                          REDM-470
   50 WRITE (MW,1015)                                                   REDM-471
      GO TO 59                                                          REDM-472
   51 WRITE (MW,1016)                                                   REDM-473
      GO TO 59                                                          REDM-474
   52 WRITE (MW,1017)                                                   REDM-475
      GO TO 59                                                          REDM-476
   53 WRITE (MW,1018) N,IJS,JIS                                         REDM-477
      GO TO 59                                                          REDM-478
   54 WRITE (MW,1019) I,J                                               REDM-479
      GO TO 59                                                          REDM-480
   55 WRITE (MW,1020) J                                                 REDM-481
      GO TO 59                                                          REDM-482
   56 WRITE (MW,1021) M,IQ(2,IT),N                                      REDM-483
      GO TO 59                                                          REDM-484
   57 WRITE (MW,1022) M                                                 REDM-485
      GO TO 59                                                          REDM-486
   58 WRITE (MW,1023) NM,IJL,JIL                                        REDM-487
      GO TO 59                                                          REDM-488
   59 WRITE (MW,1024)                                                   REDM-489
      STOP                                                              REDM-490
 1000 FORMAT (3I5)                                                      REDM-491
 1001 FORMAT (' FOR I =',I4,'  AND IP =',I4,I8,' REDUCED MATRIX ELEMENTSREDM-492
     1')                                                                REDM-493
 1002 FORMAT (4I5,F20.12)                                               REDM-494
 1003 FORMAT (10X,'FORM FACTOR =',I4,' L =',I4,4X,'2*S =',I4,4X,'2*J =',REDM-495
     1I4,4X,' REDUCED MATRIX ELEMENT',D20.6)                            REDM-496
 1004 FORMAT (//' NUCLEAR REDUCED MATRIX ELEMENTS:'/)                   REDM-497
 1005 FORMAT (/' STATE',I3,'  WITH STATE',I3,10X,' FROM',I4,' TO',I4/)  REDM-498
 1006 FORMAT (10X,'N =',I3,3X,'FORM FACTOR =',I3,3X,'MULT. =',I3,3X,'SP.REDM-499
     1-O. =',I3,6X,'MATRIX ELEMENT',1P,D15.6)                           REDM-500
 1007 FORMAT (4I5,F20.12)                                               REDM-501
 1008 FORMAT (//5X,' CORRESPONDENCE TO MULTIPOLES:'//(20X,'N =',I3,5X,'LREDM-502
     1 =',I3,5X,'2*S =',I3,5X,'2*J =',I3))                              REDM-503
 1009 FORMAT (//4X,' CORRESPONDENCE TO FORM FACTORS IN THEIR COMPUTATIONREDM-504
     1:'//(10X,'N =',I3,4X,'F.F. =',I6,4X,'MULT. =',I4,4X,'SP.-O. =',I4,REDM-505
     24X,'COUL =',I4,4X,'SPDO =',I4,4X,'COR =',I4,4X,'L =',I4))         REDM-506
 1010 FORMAT (//4X,' CORRESPONDENCE TO FORM FACTORS IN THEIR USE:'//(10XREDM-507
     1,'N =',I3,4X,'ABOVE =',I4,4X,'DISP. =',I4,4X,'SP.-O. =',I4,4X,'L =REDM-508
     2',I4))                                                            REDM-509
 1011 FORMAT (//' INCORRECT ORDER OF INPUT FOR REDUCED MATRIX ELEMENTS :REDM-510
     1',2I6,' INSTEAD OF',2I6)                                          REDM-511
 1012 FORMAT (//' A MAGNETIC COULOMB INTERACTION CANNOT BE USED IF THE CREDM-512
     1OULOMB INTERACTION IS NOT DEFORMED:',I3)                          REDM-513
 1013 FORMAT (/' MULTIPOLE ORDER',I4,' NOT OF THE SAME PARITY AS',2I4)  REDM-514
 1014 FORMAT (' ERROR FOR HARMONIC VIBRATIONAL MODEL.')                 REDM-515
 1015 FORMAT (' ERROR FOR ANHARMONIC VIBRATIONAL MODEL.')               REDM-516
 1016 FORMAT (' ERROR FOR SYMMETRIC ROTATIONAL MODEL.')                 REDM-517
 1017 FORMAT (' ERROR FOR ASYMMETRIC ROTATIONAL MODEL.')                REDM-518
 1018 FORMAT (/' 2*S =',I4,' INCORRECT BETWEEN CHANNELS WITH 2*S =',2I4)REDM-519
 1019 FORMAT (/' LAST FORM FACTOR IDENTIFICATION',I3,'  NEW FORM FACTOR REDM-520
     1IDENTIFICATION',I3,'  TOO LARGE.')                                REDM-521
 1020 FORMAT (/' LAST FORM FACTOR IDENTIFICATION',I3,'  WAS NOT ALREADY REDM-522
     1READ WITH THIS SIGN.')                                            REDM-523
 1021 FORMAT (/' 2*J =',I4,' INCORRECT WITH L =',I4,' AND 2*S =',I4)    REDM-524
 1022 FORMAT (/' J-VALUE =',I4,'  USED AS DEFORMED SPIN-ORBIT CONTROL INREDM-525
     1CORRECT. LIMIT:',I2)                                              REDM-526
 1023 FORMAT (/' 2*TRANSFER OF ANG. MOMENTUM',I4,' NOT BETWEEN',I4,' ANDREDM-527
     1',I4)                                                             REDM-528
 1024 FORMAT (//' IN REDM  ...  STOP  ...')                             REDM-529
      END                                                               REDM-530
