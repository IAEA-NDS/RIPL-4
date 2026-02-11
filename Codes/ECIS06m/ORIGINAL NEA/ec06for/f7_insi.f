C 08/03/07                                                      ECIS06  INSI-000
      SUBROUTINE INSI(W,P,PD,PH,FAM,X,PAD1,PAD2,KITER,KAB,ISM,IPD,I,V,WWINSI-001
     1W,NAT,AT,LMD,NVI,NC,FAR,FAI,MC,LO,LT,Z)                           INSI-002
C  E. C. I. S. METHOD: INTEGRATION OF A SINGLE INHOMOGENEOUS EQUATION   INSI-003
C  BY THE NUMEROV METHOD       - SCHROEDINGER EQUATION -                INSI-004
C INPUT:     P:       COUPLED SOLUTION.                                 INSI-005
C            PD:      DERIVATIVE OF THE COUPLED SOLUTION.               INSI-006
C            PH:      HOMOGENEOUS SOLUTIONS.                            INSI-007
C            FAM(I):  WAVE NUMBER FOR I=1, COEFFICIENT OF CENTRAL       INSI-008
C                     POTENTIAL FOR I=3.                                INSI-009
C            KITER:   CURRENT ITERATION NUMBER.                         INSI-010
C            KAB:     MAXIMUM NUMBER OF EQUATIONS.                      INSI-011
C            ISM:     NUMBER OF RADIAL POINTS.                          INSI-012
C            IPD(IC): FUNCTION IC NEGLIGIBLE FOR THE FIRST IPD POINTS.  INSI-013
C            I:       CHANNEL NUMBER OF THE EQUATION.                   INSI-014
C            V:       REAL, IMAGINARY POTENTIALS AND COUPLINGS.         INSI-015
C            WWW:     COUPLING BETWEEN EQUATIONS COMPUTED IN INTI.      INSI-016
C            NAT,AT:  TABLE OF COUPLING COEFFICIENTS.                   INSI-017
C            LMD:     FIRST DIMENSION OF TABLES NAT AND AT.             INSI-018
C            NVI:     ADDRESSES OF COUPLING COEFFICIENTS.               INSI-019
C            NC:      NUMBER OF COUPLED CHANNELS.                       INSI-020
C            FAR,FAI: PHASE-SHIFTS TO UPDATE.                           INSI-021
C            MC:      NUCLEAR STATE NUMBER, ANGULAR MOMENTA....         INSI-022
C            LT:      .TRUE. TO COMPUTE ONLY THE DERIVATIVE.            INSI-023
C            Z:       COULOMB INTEGRAL FOR CORRECTIONS.                 INSI-024
C            LO(I):   LOGICAL CONTROLS:                                 INSI-025
C               LO(22) =.TRUE. NO USE OF PADE APPROXIMANTS.             INSI-026
C               LO(29) =.TRUE. NO DIAGONAL TERMS IN SECOND MEMBER.      INSI-027
C               LO(30) =.TRUE. PURE DWBA CALCULATION.                   INSI-028
C               LO(44) =.TRUE. COULOMB CORRECTIONS.                     INSI-029
C               LO(57) =.TRUE. OUTPUT PHASE-SHIFTS AT EACH ITERATION.   INSI-030
C               LO(104)=.TRUE. CONVERGENCE IS OBTAINED IN THE ITERATION.INSI-031
C               LO(105)=.TRUE. CONVERGENCE OBTAINED FOR THIS EQUATION.  INSI-032
C               LO(107)=.TRUE. ALL THE COUPLINGS CALCULATED BEFOREHAND. INSI-033
C               LO(110)=.TRUE. DERIVATIVES ARE NEEDED.                  INSI-034
C OUTPUT:    P:       SOLUTION FOR EQUATION I.                          INSI-035
C            PD:      DERIVATIVE OF THE SOLUTION FOR EQUATION I.        INSI-036
C            FAR,FAI: SCATTERING COEFFICIENT.                           INSI-037
C WORKING AREAS:                                                        INSI-038
C            W(2,ISM):REAL/IMAGINARY SECOND MEMBER.                     INSI-039
C            PAD1:    ITERATION RESULTS FOR THE CHANNEL I.              INSI-040
C            PAD2:    WORKING FIELD OF PADE,TWICE LONGER THAN PAD1.     INSI-041
C            X:       INTEGRAL OF REGULAR SOLUTION WITH SECOND MEMBER.  INSI-042
C                                                                       INSI-043
C FOR THE COMMON  /CONVE/ SEE CALX.                                     INSI-044
C FOR THE COMMON  /POTE2/ SEE REDM.                                     INSI-045
C                                                                       INSI-046
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /CONVE/:                     INSI-047
C  H:         STEP SIZE FOR INTEGRATION.                                INSI-048
C  BJM:       CONVERGENCE COEFFICIENT OF IMAGINARY POTENTIAL.           INSI-049
C  EITER:     CONVERGENCE CRITERION FOR S-MATRIX.                       INSI-050
C  ACONV:     CONVERGENCE CRITERION FOR POTENTIAL AND FUNCTION.         INSI-051
C  HCONV:     ACONV*H*H CONVERGENCE CRITERION FOR SECOND MEMBERS.       INSI-052
C   USED:     H,BJM,EITER,ACONV,HCONV.                                  INSI-053
C                                                                       INSI-054
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /POTE2/:                     INSI-055
C  ITY(2):    STARTING ADDRESS OF IMAGINARY CENTRAL POTENTIAL.          INSI-056
C  ITY(9):    STARTING ADDRESS OF COULOMB CENTRAL POTENTIAL.            INSI-057
C   USED:     ITY(2),ITY(9).                                            INSI-058
C                                                                       INSI-059
C***********************************************************************INSI-060
      IMPLICIT REAL*8 (A-H,O-Z)                                         INSI-061
      LOGICAL LO(150),LT                                                INSI-062
      DIMENSION W(2,*),P(2,ISM,*),PD(2,ISM,*),PH(ISM,4,*),FAM(KAB,10),X(INSI-063
     12,*),PAD1(2,*),PAD2(2,*),IPD(*),V(ISM,*),WWW(ISM,*),NAT(2*LMD,*),AINSI-064
     2T(LMD,*),NVI(KAB,KAB,4),FAR(*),FAI(*),MC(KAB,6),Z(4)              INSI-065
      COMMON /CONVE/ H,BJM,EITER,ACONV,CONJ,HCONV                       INSI-066
      COMMON /INOUT/ MR,MW,MS                                           INSI-067
      COMMON /POTE2/ ITY(12),INVT,INTV,INSL,NPX                         INSI-068
      IF (LT) GO TO 39                                                  INSI-069
      ISM1=ISM+1                                                        INSI-070
C PUT ZERO IN THE INHOMOGENEOUS TERMS.                                  INSI-071
      DO 1 IS=1,ISM                                                     INSI-072
      W(1,IS)=0.D0                                                      INSI-073
    1 W(2,IS)=0.D0                                                      INSI-074
      IDP=ISM1                                                          INSI-075
C NON DERIVATIVE TERM OF COUPLING POTENTIALS.                           INSI-076
      DO 27 IC=1,NC                                                     INSI-077
      IF (IPD(IC).GT.ISM) GO TO 27                                      INSI-078
      IF (LO(107)) GO TO 9                                              INSI-079
C NO PREVIOUS CALCULATION OF COUPLING POTENTIALS.                       INSI-080
      IF (IC.NE.I) GO TO 3                                              INSI-081
      IF (BJM.EQ.0.D0) GO TO 3                                          INSI-082
      KT=MC(IC,4)+ITY(2)                                                INSI-083
      DO 2 IS=1,ISM                                                     INSI-084
      W(1,IS)=W(1,IS)+FAM(I,7)*V(IS,KT)*P(2,IS,IC)*BJM                  INSI-085
    2 W(2,IS)=W(2,IS)-FAM(I,7)*V(IS,KT)*P(1,IS,IC)*BJM                  INSI-086
    3 IF (LO(29).AND.(IC.EQ.I)) GO TO 27                                INSI-087
C NO PREVIOUS CALCULATION OF COUPLING POTENTIALS.                       INSI-088
      K1=NVI(IC,I,1)                                                    INSI-089
      K2=NVI(IC,I,2)                                                    INSI-090
      IF (K1.GT.K2) GO TO 15                                            INSI-091
      K3=IPD(IC)                                                        INSI-092
      IDP=MIN0(IDP,K3)                                                  INSI-093
      DO 8 K=K1,K2                                                      INSI-094
      KT=NAT(1,K)                                                       INSI-095
      DO 4 IS=K3,ISM                                                    INSI-096
      W(1,IS)=W(1,IS)+AT(2,K)*V(IS,KT)*P(1,IS,IC)                       INSI-097
    4 W(2,IS)=W(2,IS)+AT(2,K)*V(IS,KT)*P(2,IS,IC)                       INSI-098
      IF (NAT(2,K).EQ.0) GO TO 6                                        INSI-099
      KY=NAT(2,K)+ITY(2)                                                INSI-100
      DO 5 IS=K3,ISM                                                    INSI-101
      W(1,IS)=W(1,IS)-AT(2,K)*V(IS,KY)*P(2,IS,IC)                       INSI-102
    5 W(2,IS)=W(2,IS)+AT(2,K)*V(IS,KY)*P(1,IS,IC)                       INSI-103
    6 IF (.NOT.LO(133)) GO TO 8                                         INSI-104
      IF (AT(3,K).EQ.0.D0) GO TO 8                                      INSI-105
      KZ=NAT(1,K)+ITY(9)                                                INSI-106
      DO 7 IS=K3,ISM                                                    INSI-107
      W(1,IS)=W(1,IS)+AT(3,K)*V(IS,KZ)*P(1,IS,IC)                       INSI-108
    7 W(2,IS)=W(2,IS)+AT(3,K)*V(IS,KZ)*P(2,IS,IC)                       INSI-109
    8 CONTINUE                                                          INSI-110
      GO TO 15                                                          INSI-111
C COUPLING POTENTIALS ALREADY CALCULATED.                               INSI-112
    9 L1=NVI(IC,I,1)                                                    INSI-113
      L2=NVI(IC,I,2)                                                    INSI-114
      K3=IPD(IC)                                                        INSI-115
      IDP=MIN0(IDP,K3)                                                  INSI-116
      IF (L1.EQ.0) GO TO 13                                             INSI-117
      IF (L2.GT.0) GO TO 11                                             INSI-118
      DO 10 IS=K3,ISM                                                   INSI-119
      W(1,IS)=W(1,IS)+WWW(IS,L1)*P(1,IS,IC)                             INSI-120
   10 W(2,IS)=W(2,IS)+WWW(IS,L1)*P(2,IS,IC)                             INSI-121
      GO TO 15                                                          INSI-122
   11 DO 12 IS=K3,ISM                                                   INSI-123
      W(1,IS)=W(1,IS)+WWW(IS,L1)*P(1,IS,IC)-WWW(IS,L2)*P(2,IS,IC)       INSI-124
   12 W(2,IS)=W(2,IS)+WWW(IS,L1)*P(2,IS,IC)+WWW(IS,L2)*P(1,IS,IC)       INSI-125
      GO TO 15                                                          INSI-126
   13 IF (L2.EQ.0) GO TO 15                                             INSI-127
      DO 14 IS=K3,ISM                                                   INSI-128
      W(1,IS)=W(1,IS)-WWW(IS,L2)*P(2,IS,IC)                             INSI-129
   14 W(2,IS)=W(2,IS)+WWW(IS,L2)*P(1,IS,IC)                             INSI-130
   15 IF (.NOT.LO(110)) GO TO 27                                        INSI-131
C DERIVATIVE TERM OF COUPLING POTENTIALS.                               INSI-132
      IF (LO(107)) GO TO 21                                             INSI-133
C NO PREVIOUS CALCULATION OF COUPLING POTENTIALS.                       INSI-134
      K1=K2+1                                                           INSI-135
      K2=NVI(IC,I,3)                                                    INSI-136
      IF (K1.GT.K2) GO TO 27                                            INSI-137
      DO 20 K=K1,K2                                                     INSI-138
      KT=NAT(1,K)                                                       INSI-139
      DO 16 IS=K3,ISM                                                   INSI-140
      W(1,IS)=W(1,IS)+AT(2,K)*V(IS,KT)*PD(1,IS,IC)                      INSI-141
   16 W(2,IS)=W(2,IS)+AT(2,K)*V(IS,KT)*PD(2,IS,IC)                      INSI-142
      IF (NAT(2,K).EQ.0) GO TO 18                                       INSI-143
      KT=NAT(2,K)+ITY(2)                                                INSI-144
      DO 17 IS=K3,ISM                                                   INSI-145
      W(1,IS)=W(1,IS)-AT(2,K)*V(IS,KT)*PD(2,IS,IC)                      INSI-146
   17 W(2,IS)=W(2,IS)+AT(2,K)*V(IS,KT)*PD(1,IS,IC)                      INSI-147
   18 IF (.NOT.LO(133)) GO TO 20                                        INSI-148
      IF (AT(3,K).EQ.0.D0) GO TO 20                                     INSI-149
      KT=NAT(2,K)+ITY(9)                                                INSI-150
      DO 19 IS=K3,ISM                                                   INSI-151
      W(1,IS)=W(1,IS)+AT(3,K)*V(IS,KT)*P(1,IS,IC)                       INSI-152
   19 W(2,IS)=W(2,IS)+AT(3,K)*V(IS,KT)*P(2,IS,IC)                       INSI-153
   20 CONTINUE                                                          INSI-154
      GO TO 27                                                          INSI-155
C COUPLING POTENTIALS ALREADY CALCULATED.                               INSI-156
   21 L1=NVI(IC,I,3)                                                    INSI-157
      L2=NVI(IC,I,4)                                                    INSI-158
      IF (L1.LE.0) GO TO 25                                             INSI-159
      IF (L2.GT.0) GO TO 23                                             INSI-160
      DO 22 IS=K3,ISM                                                   INSI-161
      W(1,IS)=W(1,IS)+WWW(IS,L1)*PD(1,IS,IC)                            INSI-162
   22 W(2,IS)=W(2,IS)+WWW(IS,L1)*PD(2,IS,IC)                            INSI-163
      GO TO 27                                                          INSI-164
   23 DO 24 IS=K3,ISM                                                   INSI-165
      W(1,IS)=W(1,IS)+WWW(IS,L1)*PD(1,IS,IC)-WWW(IS,L2)*PD(2,IS,IC)     INSI-166
   24 W(2,IS)=W(2,IS)+WWW(IS,L1)*PD(2,IS,IC)+WWW(IS,L2)*PD(1,IS,IC)     INSI-167
      GO TO 27                                                          INSI-168
   25 IF (L2.EQ.0) GO TO 27                                             INSI-169
      DO 26 IS=K3,ISM                                                   INSI-170
      W(1,IS)=W(1,IS)-WWW(IS,L2)*PD(2,IS,IC)                            INSI-171
   26 W(2,IS)=W(2,IS)+WWW(IS,L2)*PD(1,IS,IC)                            INSI-172
   27 CONTINUE                                                          INSI-173
      IPD(I)=IDP                                                        INSI-174
C SEARCH FOR THE FIRST NON NEGLIGIBLE VALUE.                            INSI-175
      IF (IDP.GT.ISM) GO TO 42                                          INSI-176
      DO 28 IS=IDP,ISM                                                  INSI-177
      IF ((DABS(W(1,IS))+DABS(W(2,IS))).GT.HCONV) GO TO 29              INSI-178
   28 CONTINUE                                                          INSI-179
   29 IPD(I)=IS                                                         INSI-180
      KT=MAX0(IPD(I),1)                                                 INSI-181
      IF (IPD(I).GE.ISM) GO TO 42                                       INSI-182
      IF (.NOT.LO(44)) GO TO 30                                         INSI-183
      W(1,ISM)=0.D0                                                     INSI-184
      W(2,ISM)=0.D0                                                     INSI-185
C INTEGRAL OF THE REGULAR FUNCTION WITH THE SECOND MEMBER.              INSI-186
   30 X(1,KT)=-PH(KT,2,I)*W(2,KT)+PH(KT,1,I)*W(1,KT)                    INSI-187
      X(2,KT)=PH(KT,1,I)*W(2,KT)+PH(KT,2,I)*W(1,KT)                     INSI-188
      KKT=KT+1                                                          INSI-189
      IF (KKT.GT.ISM) GO TO 32                                          INSI-190
      DO 31 IS=KKT,ISM                                                  INSI-191
      X(1,IS)=X(1,IS-1)+PH(IS,1,I)*W(1,IS)-PH(IS,2,I)*W(2,IS)           INSI-192
   31 X(2,IS)=X(2,IS-1)+PH(IS,1,I)*W(2,IS)+PH(IS,2,I)*W(1,IS)           INSI-193
C SCATTERING COEFFICIENTS.                                              INSI-194
   32 BRE=X(1,ISM)/(H*FAM(I,5))-Z(1)                                    INSI-195
      BIM=X(2,ISM)/(H*FAM(I,5))-Z(2)                                    INSI-196
      IF (LO(30)) GO TO 38                                              INSI-197
      HX=Z(3)*H                                                         INSI-198
      HY=Z(4)*H                                                         INSI-199
      IST=ISM+KT                                                        INSI-200
C INTEGRAL OF THE IRREGULAR FUNCTION WITH THE SECOND MEMBER IN HX/HY    INSI-201
C AND COMPUTATION OF THE SOLUTION WITH THE CORRECTION TERM W/12.        INSI-202
      DO 33 IS=KT,ISM                                                   INSI-203
      JS=IST-IS                                                         INSI-204
      P(1,JS,I)=(X(1,JS)*PH(JS,3,I)-X(2,JS)*PH(JS,4,I)-HX*PH(JS,1,I)+HY*INSI-205
     1PH(JS,2,I))/H-W(1,JS)/12.D0                                       INSI-206
      P(2,JS,I)=(X(1,JS)*PH(JS,4,I)+X(2,JS)*PH(JS,3,I)-HY*PH(JS,1,I)-HX*INSI-207
     1PH(JS,2,I))/H-W(2,JS)/12.D0                                       INSI-208
      HX=HX-PH(JS,3,I)*W(1,JS)+PH(JS,4,I)*W(2,JS)                       INSI-209
   33 HY=HY-PH(JS,3,I)*W(2,JS)-PH(JS,4,I)*W(1,JS)                       INSI-210
      IF (KT.EQ.1) GO TO 35                                             INSI-211
      K=KT-1                                                            INSI-212
      DO 34 IS=1,K                                                      INSI-213
      P(1,IS,I)=-(HX*PH(IS,1,I)-HY*PH(IS,2,I))/H                        INSI-214
   34 P(2,IS,I)=-(HY*PH(IS,1,I)+HX*PH(IS,2,I))/H                        INSI-215
   35 DO 36 KT=1,ISM                                                    INSI-216
      IF ((DABS(P(1,KT,I))+DABS(P(2,KT,I))).GT.ACONV) GO TO 37          INSI-217
   36 CONTINUE                                                          INSI-218
   37 IPD(I)=KT                                                         INSI-219
      IF (LO(22)) GO TO 38                                              INSI-220
      PAD1(1,KITER)=BRE                                                 INSI-221
      PAD1(2,KITER)=BIM                                                 INSI-222
C  TEST OF CONVERGENCE.                                                 INSI-223
   38 LO(105)=(DABS(BIM-FAI(1)).LE.EITER.AND.(DABS(BRE-FAR(1)).LE.EITER)INSI-224
     1)                                                                 INSI-225
      IF ((.NOT.LO(22)).AND.LO(104).AND.KITER.GT.3.AND.(.NOT.LO(105))) CINSI-226
     1ALL PADE(PAD1,PAD2,KITER,BRE,BIM,EITER,1.D0,0.D0,LO)              INSI-227
      LO(104)=LO(104).AND.LO(105)                                       INSI-228
      FAR(1)=BRE                                                        INSI-229
      FAI(1)=BIM                                                        INSI-230
      IF (LO(57)) WRITE (MW,1000) I,FAR(1),FAI(1),KITER,KT              INSI-231
      IF (.NOT.LO(110).OR.LO(92)) RETURN                                INSI-232
C COMPUTATION OF R*(D/DR) OF THE SOLUTIONS.                             INSI-233
   39 IST=ISM-3                                                         INSI-234
      DO 41 J=1,2                                                       INSI-235
      PD(J,1,I)=2.5D0*P(J,2,I)-.25D0*P(J,5,I)+(5.D0*P(J,4,I)-7.7D0*P(J,1INSI-236
     1,I)-10.D0*P(J,3,I)+.2D0*P(J,6,I))/6.D0                            INSI-237
      PD(J,2,I)=-.8D0*P(J,1,I)-P(J,4,I)+(8.D0*P(J,3,I)-3.5D0*P(J,2,I)+.8INSI-238
     1D0*P(J,5,I)-.1D0*P(J,6,I))/3.D0                                   INSI-239
      PD(J,3,I)=2.25D0*(P(J,4,I)-P(J,2,I))-.45D0*(P(J,5,I)-P(J,1,I))+.05INSI-240
     1D0*P(J,6,I)                                                       INSI-241
      FP=3.D0                                                           INSI-242
      DO 40 IS=4,IST                                                    INSI-243
      FP=FP+1.D0                                                        INSI-244
   40 PD(J,IS,I)=FP*(.75D0*(P(J,IS+1,I)-P(J,IS-1,I))-.15D0*(P(J,IS+2,I)-INSI-245
     1P(J,IS-2,I))+(P(J,IS+3,I)-P(J,IS-3,I))/60.D0)                     INSI-246
      FP=FP+1.D0                                                        INSI-247
      PD(J,ISM-2,I)=FP*(P(J,ISM-6,I)-8.D0*P(J,ISM-5,I)+30.D0*P(J,ISM-4,IINSI-248
     1)-80.D0*P(J,IST,I)+35.D0*P(J,ISM-2,I)+24.D0*P(J,ISM-1,I)-2.D0*P(J,INSI-249
     2ISM,I))/60.D0                                                     INSI-250
      FP=FP+1.D0                                                        INSI-251
      PD(J,ISM-1,I)=FP*(15.D0*P(J,ISM-5,I)-2.D0*P(J,ISM-6,I)-50.D0*P(J,IINSI-252
     1SM-4,I)+100.D0*P(J,IST,I)-150.D0*P(J,ISM-2,I)+77.D0*P(J,ISM-1,I)+1INSI-253
     20.D0*P(J,ISM,I))/60.D0                                            INSI-254
      FP=FP+1.D0                                                        INSI-255
      PD(J,ISM,I)=FP*(10.D0*P(J,ISM-6,I)-72.D0*P(J,ISM-5,I)+225.D0*P(J,IINSI-256
     1SM-4,I)-400.D0*P(J,IST,I)+450.D0*P(J,ISM-2,I)-360.D0*P(J,ISM-1,I)+INSI-257
     2147.D0*P(J,ISM,I))/60.D0                                          INSI-258
   41 CONTINUE                                                          INSI-259
      RETURN                                                            INSI-260
   42 IF (LO(22)) RETURN                                                INSI-261
      PAD1(1,KITER)=0.D0                                                INSI-262
      PAD1(2,KITER)=0.D0                                                INSI-263
      RETURN                                                            INSI-264
 1000 FORMAT (5X,I5,2D30.15,I10,10X,I5)                                 INSI-265
      END                                                               INSI-266
