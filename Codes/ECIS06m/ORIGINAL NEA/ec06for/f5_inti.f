C 08/03/07                                                      ECIS06  INTI-000
      SUBROUTINE INTI(FAM,X,FR,GR,WRE,PAD,IPE,ISM,KAB,W,ITERM,NC,V,NVI,MINTI-001
     1C,NAT,AT,LMD,AG,NCIN,NNI,ITERR,LO)                                INTI-002
C  E. C. I. S. METHOD: SCHROEDINGER EQUATION DRIVING ROUTINE.           INTI-003
C  INTI CALLS INSH  TO SOLVE THE SINGLE HOMOGENEOUS EQUATIONS,          INTI-004
C             INSI  TO SOLVE THE SINGLE INHOMOGENEOUS EQUATIONS.        INTI-005
C  PADE APPROXIMANTS OF TYPE I MAY BE USED TO ACCELERATE THE CONVERGENCEINTI-006
C INPUT:     FAM(*,I):MATCHING VALUES FOR I=1 TO 6, CONSTANTS OF EACH   INTI-007
C                     EQUATION FOR I=7 TO 12.                           INTI-008
C            ISM:     NUMBER OF RADIAL POINTS.                          INTI-009
C            KAB:     MAXIMUM NUMBER OF COUPLED CHANNELS.               INTI-010
C            ITERM:   MAXIMUM NUMBER OF ITERATIONS AND DIM. OF PAD.     INTI-011
C            NC:      NUMBER OF COUPLED EQUATIONS.                      INTI-012
C            V:       REAL, IMAGINARY POTENTIALS AND COUPLINGS.         INTI-013
C            NVI:     ADDRESSES OF COUPLINGS.                           INTI-014
C            MC:      NUCLEAR STATE NUMBER, ANGULAR MOMENTA....         INTI-015
C            NAT,AT:  TABLE OF COUPLING COEFFICIENTS.                   INTI-016
C            LMD:     FIRST DIMENSION OF TABLES NAT AND AT.             INTI-017
C            AG:      COULOMB INTEGRALS FOR COULOMB CORRECTIONS.        INTI-018
C            NCIN:    NUMBER OF SOLUTIONS.                              INTI-019
C            LO(I):   LOGICAL CONTROLS:                                 INTI-020
C               LO(22) =.TRUE. NO USE OF PADE APPROXIMANTS.             INTI-021
C               LO(23) =.TRUE. NO USE OF PADE AND SHIFT TO USUAL COUPLEDINTI-022
C                              EQUATIONS WHEN THERE IS NO CONVERGENCE.  INTI-023
C               LO(29) =.TRUE. NO DIAGONAL TERMS IN SECOND MEMBER.      INTI-024
C               LO(30) =.TRUE. PURE DWBA CALCULATION.                   INTI-025
C               LO(44) =.TRUE. COULOMB CORRECTIONS.                     INTI-026
C               LO(53) =.TRUE. OUTPUT OF THE NUMBER OF ITERATIONS.      INTI-027
C               LO(57) =.TRUE. OUTPUT PHASE-SHIFTS AT EACH ITERATION.   INTI-028
C               LO(104)=.TRUE. CONVERGENCE IS OBTAINED IN THE ITERATION.INTI-029
C               LO(106)=.TRUE. WHEN THE ITERATION IS NOT THE LAST ONE   INTI-030
C                              ALLOWED.                                 INTI-031
C               LO(107)=.TRUE. ALL THE COUPLINGS CALCULATED BEFOREHAND. INTI-032
C               LO(110)=.TRUE. DERIVATIVES ARE NEEDED.                  INTI-033
C               LO(121)=.TRUE. OPTICAL MODEL WITHOUT COUPLING.          INTI-034
C OUTPUT:    FAM(*,I):SCATTERING COEFFICIENTS FOR I>9.                  INTI-035
C            NNI:     NUMBER OF EQUATIONS WITH NEGLIGIBLE INHOMOGENEOUS INTI-036
C                     TERM.                                             INTI-037
C            ITERR:   LARGEST NUMBER OF ITERATIONS DONE FOR THIS J.     INTI-038
C WORKING AREAS:                                                        INTI-039
C            X:       USED IN INSH AND INSI.                            INTI-040
C            FR:      SOLUTIONS OF HOMOGENEOUS EQUATIONS.               INTI-041
C            GR:      SOLUTIONS OF COUPLED EQUATIONS.                   INTI-042
C            WRE:     REAL/IMAGINARY INHOMOGENEOUS TERM.                INTI-043
C            PAD:     PADE APPROXIMANTS.                                INTI-044
C            IPE(I):  FIRST NON NEGLIGIBLE POINT OF UNCOUPLED FUNCTIONS INTI-045
C                     FOR I=1 TO NC, OF COUPLED FUNCTIONS FOR I>NC.     INTI-046
C            W:       FREE PART OF THE STORAGE FOR COUPLINGS            INTI-047
C                                                                       INTI-048
C FOR THE COMMON  /CONVE/ SEE CALX.                                     INTI-049
C FOR THE COMMON  /POTE2/ SEE REDM.                                     INTI-050
C                                                                       INTI-051
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /CONVE/:                     INTI-052
C  BJM:       CONVERGENCE COEFFICIENT OF IMAGINARY POTENTIAL.           INTI-053
C   USED:     BJM.                                                      INTI-054
C                                                                       INTI-055
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /POTE2/:                     INTI-056
C  ITY(2):    STARTING ADDRESS OF IMAGINARY CENTRAL POTENTIAL.          INTI-057
C  ITY(9):    STARTING ADDRESS OF COULOMB CENTRAL POTENTIAL.            INTI-058
C   USED:     ITY(2),ITY(9).                                            INTI-059
C                                                                       INTI-060
C***********************************************************************INTI-061
      IMPLICIT REAL*8 (A-H,O-Z)                                         INTI-062
      LOGICAL LO(150)                                                   INTI-063
      DIMENSION FAM(KAB,*),X(*),FR(ISM,4,*),GR(2,ISM,*),WRE(2,*),PAD(2,IINTI-064
     1TERM,*),IPE(*),W(ISM,*),V(ISM,*),NVI(KAB,KAB,4),MC(KAB,6),NAT(2*LMINTI-065
     2D,*),AT(LMD,*),AG(KAB,KAB,4),Z(4),Y(2)                            INTI-066
      COMMON /CONVE/ H,BJM,EITER,ACONV,CONJ,HCONV                       INTI-067
      COMMON /INOUT/ MR,MW,MS                                           INTI-068
      COMMON /POTE2/ ITY(12),INVT,INTV,INSL,NPX                         INTI-069
C COMPUTE ALL COUPLING POTENTIALS AND SOLVE ALL HOMOGENEOUS EQUATIONS.  INTI-070
      DO 5 I=1,NC                                                       INTI-071
      I1=MC(I,4)                                                        INTI-072
      IF (I1.LT.0) GO TO 1                                              INTI-073
      CALL INSH(FR(1,1,I),I,I1,KT,FAM,X,KAB,ISM,V,NAT,AT,LMD,NVI,LO)    INTI-074
      IPE(I)=KT                                                         INTI-075
      GO TO 4                                                           INTI-076
    1 I1=-I1                                                            INTI-077
      IF (I.EQ.I1) GO TO 4                                              INTI-078
      DO 3 J=1,4                                                        INTI-079
      DO 2 IS=1,ISM                                                     INTI-080
    2 FR(IS,J,I)=FR(IS,J,I1)                                            INTI-081
    3 CONTINUE                                                          INTI-082
      IPE(I)=IPE(I1)                                                    INTI-083
    4 IF (LO(57)) WRITE (MW,1000) FAM(I,9),FAM(I,10),I,IPE(I)           INTI-084
    5 CONTINUE                                                          INTI-085
      IF ((.NOT.LO(107)).OR.LO(121)) GO TO 24                           INTI-086
      IVV=0                                                             INTI-087
C CALCULATION OF THE NON DIAGONAL COUPLING POTENTIALS.                  INTI-088
      DO 23 IC=1,NC                                                     INTI-089
      DO 19 IP=1,NC                                                     INTI-090
      K1=NVI(IP,IC,1)                                                   INTI-091
      K2=NVI(IP,IC,2)                                                   INTI-092
      K3=NVI(IP,IC,3)                                                   INTI-093
      IF (K2.EQ.K3.AND.IP.LT.IC) GO TO 17                               INTI-094
      DO 6 J=1,4                                                        INTI-095
    6 NVI(IP,IC,J)=0                                                    INTI-096
C NON DERIVATIVE COUPLING POTENTIALS FOR IK=1.                          INTI-097
C DERIVATIVE COUPLING POTENTIALS FOR IK=3.                              INTI-098
      DO 16 IK=1,3,2                                                    INTI-099
      IF (K1.GT.K2) GO TO 15                                            INTI-100
      IF (LO(29).AND.(IC.EQ.IP)) GO TO 15                               INTI-101
      IVV=IVV+1                                                         INTI-102
      NVI(IP,IC,IK)=IVV                                                 INTI-103
      NIMA=0                                                            INTI-104
      DO 7 IS=1,ISM                                                     INTI-105
    7 W(IS,IVV)=0.D0                                                    INTI-106
      DO 11 K=K1,K2                                                     INTI-107
      IF (AT(2,K).EQ.0.D0) GO TO 9                                      INTI-108
      KT=NAT(1,K)                                                       INTI-109
      IF (NAT(2,K).NE.0) NIMA=NIMA+1                                    INTI-110
      DO 8 IS=1,ISM                                                     INTI-111
    8 W(IS,IVV)=W(IS,IVV)+AT(2,K)*V(IS,KT)                              INTI-112
    9 IF (.NOT.LO(133)) GO TO 11                                        INTI-113
      IF (AT(3,K).EQ.0.D0) GO TO 11                                     INTI-114
      KT=NAT(1,K)+ITY(9)                                                INTI-115
      DO 10 IS=1,ISM                                                    INTI-116
   10 W(IS,IVV)=W(IS,IVV)+AT(3,K)*V(IS,KT)                              INTI-117
   11 CONTINUE                                                          INTI-118
      IF (NIMA.EQ.0) GO TO 15                                           INTI-119
      IVV=IVV+1                                                         INTI-120
      NVI(IP,IC,IK+1)=IVV                                               INTI-121
      DO 12 IS=1,ISM                                                    INTI-122
   12 W(IS,IVV)=0.D0                                                    INTI-123
      DO 14 K=K1,K2                                                     INTI-124
      IF (NAT(2,K).EQ.0) GO TO 14                                       INTI-125
      KT=NAT(2,K)+ITY(2)                                                INTI-126
      DO 13 IS=1,ISM                                                    INTI-127
   13 W(IS,IVV)=W(IS,IVV)+AT(2,K)*V(IS,KT)                              INTI-128
   14 CONTINUE                                                          INTI-129
C DERIVATIVE COUPLING POTENTIALS.                                       INTI-130
   15 K1=K2+1                                                           INTI-131
   16 K2=K3                                                             INTI-132
      GO TO 19                                                          INTI-133
C SYMMETRISATION OF THE TABLE.                                          INTI-134
   17 DO 18 K=1,4                                                       INTI-135
   18 NVI(IP,IC,K)=NVI(IC,IP,K)                                         INTI-136
   19 CONTINUE                                                          INTI-137
C CORRECTION TO AN INCREASE OF THE IMAGINARY POTENTIAL.                 INTI-138
      IF (BJM.EQ.0.D0) GO TO 23                                         INTI-139
      I1=MC(IC,4)+ITY(2)                                                INTI-140
      IF (NVI(IC,IC,2).EQ.0) GO TO 21                                   INTI-141
      K=NVI(IC,IC,2)                                                    INTI-142
      DO 20 IS=1,ISM                                                    INTI-143
   20 W(IS,K)=W(IS,K)-FAM(IC,7)*V(IS,I1)*BJM                            INTI-144
      GO TO 23                                                          INTI-145
   21 IVV=IVV+1                                                         INTI-146
      NVI(IC,IC,2)=IVV                                                  INTI-147
      DO 22 IS=1,ISM                                                    INTI-148
   22 W(IS,IVV)=-FAM(IC,7)*V(IS,I1)*BJM                                 INTI-149
   23 CONTINUE                                                          INTI-150
   24 NNI=0                                                             INTI-151
C LOOP ON SOLUTIONS.                                                    INTI-152
      DO 38 NCI=1,NCIN                                                  INTI-153
      NCI1=NCI+10                                                       INTI-154
      NCI2=NCI1+NCIN                                                    INTI-155
C SET THE ZERO'S ORDER SOLUTION.                                        INTI-156
      DO 25 IC=1,NC                                                     INTI-157
      IPE(IC+NC)=ISM+1                                                  INTI-158
      FAM(IC,NCI1)=0.D0                                                 INTI-159
   25 FAM(IC,NCI2)=0.D0                                                 INTI-160
      FAM(NCI,NCI1)=FAM(NCI,9)                                          INTI-161
      FAM(NCI,NCI2)=FAM(NCI,10)                                         INTI-162
      IF ((IPE(NCI).GT.ISM-3).OR.LO(121)) GO TO 38                      INTI-163
      DO 26 IS=1,ISM                                                    INTI-164
      GR(1,IS,NCI)=FR(IS,1,NCI)                                         INTI-165
   26 GR(2,IS,NCI)=FR(IS,2,NCI)                                         INTI-166
      IPE(NCI+NC)=IPE(NCI)                                              INTI-167
C DERIVE THE ZERO'S ORDER EQUATION.                                     INTI-168
      IF (LO(110)) CALL INSI(WRE,GR,GR(1,1,NC+1),FR,FAM,X,PAD,PAD,1,KAB,INTI-169
     1ISM,IPE,NCI,V,W,NAT,AT,LMD,NVI,NC,Y,Y,MC,LO,.TRUE.,Z)             INTI-170
      AMAX=0.D0                                                         INTI-171
C E. C. I. S. LOOP.                                                     INTI-172
      DO 36 KITER=1,ITERM                                               INTI-173
      ITERR=MAX0(ITERR,KITER)                                           INTI-174
      LO(104)=.TRUE.                                                    INTI-175
      LO(106)=((KITER.NE.ITERM).AND.(AMAX.LT.1.D10)).OR.LO(23)          INTI-176
      IF (AMAX.GT.1.D10) WRITE (MW,1001) AMAX,KITER                     INTI-177
      DO 35 IC=1,NC                                                     INTI-178
      I=MOD(IC+NCI-1,NC)+1                                              INTI-179
C COMPUTE SECOND MEMBERS AND SOLVE THE INHOMOGENEOUS EQUATIONS.         INTI-180
      DO 27 J=1,4                                                       INTI-181
   27 Z(J)=0.D0                                                         INTI-182
      IF (.NOT.LO(44)) GO TO 32                                         INTI-183
      DO 30 J=1,NC                                                      INTI-184
      IF ((LO(29).AND.(I.EQ.J)).OR.(IPE(J+NC).GE.ISM)) GO TO 30         INTI-185
C ORDER IN FG   F(EI)*F(EF),G(EI)*F(EF),F(EI)*G(EF),G(EI)*G(EF).        INTI-186
      IF (LO(22)) GO TO 28                                              INTI-187
      IJ=KITER                                                          INTI-188
      IF (1+MOD(J+NC-NCI-1,NC).GE.IC) IJ=IJ-1                           INTI-189
      IF (IJ.EQ.0) GO TO 28                                             INTI-190
      Y(1)=PAD(1,IJ,J)                                                  INTI-191
      Y(2)=PAD(2,IJ,J)                                                  INTI-192
      GO TO 29                                                          INTI-193
   28 Y(1)=FAM(J,NCI1)                                                  INTI-194
      Y(2)=FAM(J,NCI2)                                                  INTI-195
   29 Z(1)=Z(1)+Y(1)*AG(I,J,3)-Y(2)*AG(I,J,1)                           INTI-196
      Z(2)=Z(2)+Y(1)*AG(I,J,1)+Y(2)*AG(I,J,3)                           INTI-197
      Z(3)=Z(3)+Y(1)*AG(I,J,4)-Y(2)*AG(I,J,2)                           INTI-198
      Z(4)=Z(4)+Y(1)*AG(I,J,2)+Y(2)*AG(I,J,4)                           INTI-199
   30 CONTINUE                                                          INTI-200
      IF (I.EQ.NCI) GO TO 31                                            INTI-201
      Z(1)=Z(1)+AG(I,NCI,1)                                             INTI-202
      Z(3)=Z(3)+AG(I,NCI,2)                                             INTI-203
   31 Z(3)=Z(3)-Z(2)                                                    INTI-204
      Z(4)=Z(4)+Z(1)                                                    INTI-205
      Z(1)=Z(1)+FAM(I,9)*Z(3)-FAM(I,10)*Z(4)                            INTI-206
      Z(2)=Z(2)+Z(3)*FAM(I,10)+Z(4)*FAM(I,9)                            INTI-207
   32 IF (I.NE.NCI) GO TO 33                                            INTI-208
      Z(1)=Z(1)-FAM(I,9)                                                INTI-209
      Z(2)=Z(2)-FAM(I,10)                                               INTI-210
      Z(3)=Z(3)-1.D0                                                    INTI-211
   33 IF (LO(23).AND.(.NOT.LO(106))) LO(104)=.TRUE.                     INTI-212
      CALL INSI(WRE,GR,GR(1,1,NC+1),FR,FAM,X,PAD(1,1,I),PAD(1,1,NC+1),KIINTI-213
     1TER,KAB,ISM,IPE(NC+1),I,V,W,NAT,AT,LMD,NVI,NC,FAM(I,NCI1),FAM(I,NCINTI-214
     2I2),MC,LO,.FALSE.,Z)                                              INTI-215
      IF (IPE(I+NC).GE.ISM) GO TO 34                                    INTI-216
      AMAX=DMAX1(AMAX,DABS(FAM(I,NCI1))+DABS(FAM(I,NCI2)))              INTI-217
      IF (I.NE.NCI.AND.LO(30)) IPE(I+NC)=ISM+1                          INTI-218
      GO TO 35                                                          INTI-219
   34 FAM(I,NCI1)=-Z(1)                                                 INTI-220
      FAM(I,NCI2)=-Z(2)                                                 INTI-221
      IF (LO(57)) WRITE (MW,1002)                                       INTI-222
      IF (KITER.EQ.1.AND.(.NOT.LO(44))) NNI=NNI+1                       INTI-223
      IF (I.EQ.NCI) IPE(NCI+NC)=IPE(NCI)                                INTI-224
   35 CONTINUE                                                          INTI-225
      IF (LO(104).OR.(AMAX.GT.1.D10.AND.LO(23))) GO TO 37               INTI-226
   36 CONTINUE                                                          INTI-227
      KITER=MIN0(KITER,ITERM)                                           INTI-228
   37 IF (LO(53)) WRITE (MW,1003) KITER                                 INTI-229
   38 CONTINUE                                                          INTI-230
      NNI=NNI/NCIN                                                      INTI-231
      RETURN                                                            INTI-232
 1000 FORMAT (2D30.15,5X,2I5)                                           INTI-233
 1001 FORMAT (' MAXIMUM',D15.6,' OBTAINED IN PREVIOUS ITERATION. LAST ONINTI-234
     1E IS',I3)                                                         INTI-235
 1002 FORMAT (5X,'THE INHOMOGENEOUS TERM IS NEGLECTED.')                INTI-236
 1003 FORMAT (5X,I5,' ITERATIONS.')                                     INTI-237
      END                                                               INTI-238
