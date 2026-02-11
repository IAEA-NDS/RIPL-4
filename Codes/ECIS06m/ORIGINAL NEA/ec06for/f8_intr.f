C 27/06/06                                                      ECIS06  INTR-000
      SUBROUTINE INTR(FAM,X,FR,GR,WRE,PAD,IPE,ISM,KAB,W,ITERM,NC,V,VR,NVINTR-001
     1I,MC,NAT,AT,AG,NCIN,NNI,ITERR,LO)                                 INTR-002
C  E. C. I. S. METHOD: DIRAC EQUATION DRIVING ROUTINE.                  INTR-003
C             INRH  TO SOLVE THE SINGLE HOMOGENEOUS EQUATIONS.          INTR-004
C             INRI  TO SOLVE THE SINGLE INHOMOGENEOUS EQUATIONS.        INTR-005
C  PADE APPROXIMANTS OF TYPE I MAY BE USED TO ACCELERATE THE CONVERGENCEINTR-006
C INPUT:     FAM(*,I):MATCHING VALUES FOR I=1 TO 6, CONSTANTS OF EACH   INTR-007
C                     EQUATION FOR I=7 TO 10.                           INTR-008
C            ISM:     NUMBER OF RADIAL POINTS.                          INTR-009
C            KAB:     MAXIMUM NUMBER OF COUPLED CHANNELS.               INTR-010
C            ITERM:   MAXIMUM NUMBER OF ITERATIONS AND DIM. OF PAD.     INTR-011
C            NC:      NUMBER OF COUPLED EQUATIONS.                      INTR-012
C            V,VR:    REAL, IMAGINARY POTENTIALS AND COUPLINGS.         INTR-013
C            NVI:     ADDRESSES OF COUPLINGS.                           INTR-014
C            MC:      NUCLEAR STATE NUMBER, ANGULAR MOMENTA....         INTR-015
C            NAT,AT:  TABLE OF COUPLING COEFFICIENTS.                   INTR-016
C            AG:      COULOMB INTEGRALS FOR COULOMB CORRECTIONS.        INTR-017
C            NCIN:    NUMBER OF SOLUTIONS.                              INTR-018
C            LO(I):   LOGICAL CONTROLS:                                 INTR-019
C               LO(22) =.TRUE. NO USE OF PADE APPROXIMANTS.             INTR-020
C               LO(23) =.TRUE. NO USE OF PADE AND SHIFT TO USUAL COUPLEDINTR-021
C                              EQUATIONS WHEN THERE IS NO CONVERGENCE.  INTR-022
C               LO(30) =.TRUE. PURE DWBA CALCULATION.                   INTR-023
C               LO(44) =.TRUE. COULOMB CORRECTIONS.                     INTR-024
C               LO(53) =.TRUE. OUTPUT OF THE NUMBER OF ITERATIONS.      INTR-025
C               LO(57) =.TRUE. OUTPUT PHASE-SHIFTS AT EACH ITERATION.   INTR-026
C               LO(74) =.TRUE. OUTPUT OF TIME IN DIFFERENT STEPS.       INTR-027
C               LO(104)=.TRUE. CONVERGENCE IS OBTAINED IN THE ITERATION.INTR-028
C               LO(106)=.TRUE. WHEN THE ITERATION IS NOT THE LAST ONE   INTR-029
C                              ALLOWED.                                 INTR-030
C               LO(107)=.TRUE. ALL THE COUPLINGS CALCULATED BEFOREHAND. INTR-031
C               LO(121)=.TRUE. OPTICAL MODEL WITHOUT COUPLING.          INTR-032
C OUTPUT:    FAM(*,I):SCATTERING COEFFICIENTS FOR I>9.                  INTR-033
C            NNI:     NUMBER OF EQUATIONS WITH NEGLIGIBLE INHOMOGENEOUS INTR-034
C                     TERM.                                             INTR-035
C            ITERR:   LARGEST NUMBER OF ITERATIONS DONE FOR THIS J.     INTR-036
C WORKING AREAS:                                                        INTR-037
C            X:       USED IN INRH AND INRI.                            INTR-038
C            FR:      SOLUTIONS OF HOMOGENEOUS EQUATIONS.               INTR-039
C            GR:      SOLUTIONS OF COUPLED EQUATIONS.                   INTR-040
C            WRE:     REAL/IMAGINARY INHOMOGENEOUS TERM.                INTR-041
C            PAD:     PADE APPROXIMANTS.                                INTR-042
C            IPE(*,I):FIRST NON NEGLIGIBLE POINT OF UNCOUPLED FUNCTIONS INTR-043
C                     FOR I=1, OF COUPLED FUNCTIONS FOR I=2.            INTR-044
C            W:       FREE PART OF THE STORAGE FOR COUPLINGS            INTR-045
C***********************************************************************INTR-046
      IMPLICIT REAL*8 (A-H,O-Z)                                         INTR-047
      LOGICAL LO(150)                                                   INTR-048
      DIMENSION FAM(KAB,*),X(*),FR(ISM,8,*),GR(ISM,4,*),WRE(ISM,*),PAD(2INTR-049
     1,ITERM,*),IPE(NC,2),W(ISM,4,*),V(ISM,14,*),VR(ISM,4,*),NVI(KAB,KABINTR-050
     2,4),MC(KAB,6),NAT(6,*),AT(3,*),AG(KAB,KAB,4),Y(2),Z(4)            INTR-051
      COMMON /INOUT/ MR,MW,MS                                           INTR-052
      NNI=0                                                             INTR-053
C COMPUTE ALL COUPLING POTENTIALS AND SOLVE ALL HOMOGENEOUS EQUATIONS.  INTR-054
      DO 5 I=1,NC                                                       INTR-055
      I1=MC(I,4)                                                        INTR-056
      IF (I1.LT.0) GO TO 1                                              INTR-057
      CALL INRH(FR(1,1,I),I,K,FAM,X,KAB,ISM,DFLOAT(MC(I,6)+1),V(1,1,I1),INTR-058
     1LO)                                                               INTR-059
      IPE(I,1)=K                                                        INTR-060
      GO TO 4                                                           INTR-061
    1 I1=-I1                                                            INTR-062
      IF (I1.EQ.I) GO TO 4                                              INTR-063
      DO 3 J=1,8                                                        INTR-064
      DO 2 IS=1,ISM                                                     INTR-065
    2 FR(IS,J,I)=FR(IS,J,I1)                                            INTR-066
    3 CONTINUE                                                          INTR-067
      IPE(I,1)=IPE(I1,1)                                                INTR-068
    4 IF (LO(57)) WRITE (MW,1000) FAM(I,9),FAM(I,10),I,IPE(I,1)         INTR-069
    5 CONTINUE                                                          INTR-070
      IF (.NOT.LO(107).OR.LO(121)) GO TO 16                             INTR-071
      IW=0                                                              INTR-072
C CALCULATION OF THE NON DIAGONAL COUPLING POTENTIALS.                  INTR-073
      DO 15 IC=1,NC                                                     INTR-074
      DO 14 IP=1,IC                                                     INTR-075
      K1=NVI(IP,IC,1)                                                   INTR-076
      K2=NVI(IP,IC,2)                                                   INTR-077
      K3=NVI(IP,IC,3)                                                   INTR-078
      NVI(IP,IC,3)=0                                                    INTR-079
      NVI(IP,IC,4)=0                                                    INTR-080
      IF (K1.GT.K2) GO TO 9                                             INTR-081
      IW=IW+1                                                           INTR-082
      NVI(IP,IC,3)=IW                                                   INTR-083
      K=NAT(1,K1)                                                       INTR-084
      DO 6 IS=1,ISM                                                     INTR-085
      W(IS,1,IW)=AT(2,K1)*VR(IS,1,K)                                    INTR-086
      W(IS,2,IW)=AT(2,K1)*VR(IS,2,K)                                    INTR-087
      W(IS,3,IW)=AT(3,K1)*VR(IS,3,K)                                    INTR-088
    6 W(IS,4,IW)=AT(3,K1)*VR(IS,4,K)                                    INTR-089
      K4=K1+1                                                           INTR-090
      IF (K4.GT.K2) GO TO 9                                             INTR-091
      DO 8 K1=K4,K2                                                     INTR-092
      K=NAT(1,K1)                                                       INTR-093
      DO 7 IS=1,ISM                                                     INTR-094
      W(IS,1,IW)=W(IS,1,IW)+AT(2,K1)*VR(IS,1,K)                         INTR-095
      W(IS,2,IW)=W(IS,2,IW)+AT(2,K1)*VR(IS,2,K)                         INTR-096
      W(IS,3,IW)=W(IS,3,IW)+AT(3,K1)*VR(IS,3,K)                         INTR-097
    7 W(IS,4,IW)=W(IS,4,IW)+AT(3,K1)*VR(IS,4,K)                         INTR-098
    8 CONTINUE                                                          INTR-099
    9 IF (K2.GE.K3) GO TO 13                                            INTR-100
      IW=IW+1                                                           INTR-101
      NVI(IP,IC,4)=IW                                                   INTR-102
      K1=K2+1                                                           INTR-103
      K=NAT(1,K1)                                                       INTR-104
      DO 10 IS=1,ISM                                                    INTR-105
      W(IS,1,IW)=AT(2,K1)*VR(IS,1,K)+AT(3,K1)*VR(IS,3,K)                INTR-106
      W(IS,2,IW)=AT(2,K1)*VR(IS,2,K)+AT(3,K1)*VR(IS,4,K)                INTR-107
      W(IS,3,IW)=AT(2,K1)*VR(IS,1,K)-AT(3,K1)*VR(IS,3,K)                INTR-108
   10 W(IS,4,IW)=AT(2,K1)*VR(IS,2,K)-AT(3,K1)*VR(IS,4,K)                INTR-109
      K4=K1+1                                                           INTR-110
      IF (K4.GT.K3) GO TO 13                                            INTR-111
      DO 12 K1=K4,K3                                                    INTR-112
      K=NAT(1,K1)                                                       INTR-113
      DO 11 IS=1,ISM                                                    INTR-114
      W(IS,1,IW)=W(IS,1,IW)+AT(2,K1)*VR(IS,1,K)+AT(3,K1)*VR(IS,3,K)     INTR-115
      W(IS,2,IW)=W(IS,2,IW)+AT(2,K1)*VR(IS,2,K)+AT(3,K1)*VR(IS,4,K)     INTR-116
      W(IS,3,IW)=W(IS,3,IW)+AT(2,K1)*VR(IS,1,K)-AT(3,K1)*VR(IS,3,K)     INTR-117
   11 W(IS,4,IW)=W(IS,4,IW)+AT(2,K1)*VR(IS,2,K)-AT(3,K1)*VR(IS,4,K)     INTR-118
   12 CONTINUE                                                          INTR-119
   13 NVI(IC,IP,3)=NVI(IP,IC,3)                                         INTR-120
      NVI(IC,IP,4)=NVI(IP,IC,4)                                         INTR-121
   14 CONTINUE                                                          INTR-122
   15 CONTINUE                                                          INTR-123
   16 IF (LO(74)) CALL HORA                                             INTR-124
C LOOP ON SOLUTIONS.                                                    INTR-125
      DO 30 NCI=1,NCIN                                                  INTR-126
      NCI1=NCI+10                                                       INTR-127
      NCI2=NCI1+NCIN                                                    INTR-128
C SET THE ZERO'S ORDER SOLUTION.                                        INTR-129
      DO 17 IC=1,NC                                                     INTR-130
      IPE(IC,2)=ISM+1                                                   INTR-131
      FAM(IC,NCI1)=0.D0                                                 INTR-132
   17 FAM(IC,NCI2)=0.D0                                                 INTR-133
      FAM(NCI,NCI1)=FAM(NCI,9)                                          INTR-134
      FAM(NCI,NCI2)=FAM(NCI,10)                                         INTR-135
      IF ((IPE(NCI,1).GT.ISM-3).OR.LO(121)) GO TO 30                    INTR-136
      DO 18 IS=1,ISM                                                    INTR-137
      GR(IS,1,NCI)=FR(IS,1,NCI)                                         INTR-138
      GR(IS,2,NCI)=FR(IS,2,NCI)                                         INTR-139
      GR(IS,3,NCI)=FR(IS,5,NCI)                                         INTR-140
   18 GR(IS,4,NCI)=FR(IS,6,NCI)                                         INTR-141
      IPE(NCI,2)=IPE(NCI,1)                                             INTR-142
      AMAX=0.D0                                                         INTR-143
C E. C. I. S. LOOP.                                                     INTR-144
      DO 28 KITER=1,ITERM                                               INTR-145
      ITERR=MAX0(ITERR,KITER)                                           INTR-146
      LO(104)=.TRUE.                                                    INTR-147
      LO(106)=((KITER.NE.ITERM).AND.(AMAX.LT.1.D10)).OR.LO(23)          INTR-148
      IF (AMAX.GT.1.D10) WRITE (MW,1001) AMAX,KITER                     INTR-149
      DO 27 IC=1,NC                                                     INTR-150
      I=MOD(IC+NCI-1,NC)+1                                              INTR-151
      I1=MC(I,1)                                                        INTR-152
      DO 19 J=1,4                                                       INTR-153
   19 Z(J)=0.D0                                                         INTR-154
      IF (.NOT.LO(44)) GO TO 24                                         INTR-155
      DO 22 J=1,NC                                                      INTR-156
      IF (IPE(J,2).GE.ISM) GO TO 22                                     INTR-157
C ORDER IN FG   F(EI)*F(EF),G(EI)*F(EF),F(EI)*G(EF),G(EI)*G(EF).        INTR-158
      IF (LO(22)) GO TO 20                                              INTR-159
      IJ=KITER                                                          INTR-160
      IF (1+MOD(J+NC-NCI-1,NC).GE.IC) IJ=IJ-1                           INTR-161
      IF (IJ.EQ.0) GO TO 20                                             INTR-162
      Y(1)=PAD(1,IJ,J)                                                  INTR-163
      Y(2)=PAD(2,IJ,J)                                                  INTR-164
      GO TO 21                                                          INTR-165
   20 Y(1)=FAM(J,NCI1)                                                  INTR-166
      Y(2)=FAM(J,NCI2)                                                  INTR-167
   21 Z(1)=Z(1)+Y(1)*AG(I,J,3)-Y(2)*AG(I,J,1)                           INTR-168
      Z(2)=Z(2)+Y(1)*AG(I,J,1)+Y(2)*AG(I,J,3)                           INTR-169
      Z(3)=Z(3)+Y(1)*AG(I,J,4)-Y(2)*AG(I,J,2)                           INTR-170
      Z(4)=Z(4)+Y(1)*AG(I,J,2)+Y(2)*AG(I,J,4)                           INTR-171
   22 CONTINUE                                                          INTR-172
      IF (I.EQ.NCI) GO TO 23                                            INTR-173
      Z(1)=Z(1)+AG(I,NCI,1)                                             INTR-174
      Z(3)=Z(3)+AG(I,NCI,2)                                             INTR-175
   23 Z(3)=Z(3)-Z(2)                                                    INTR-176
      Z(4)=Z(4)+Z(1)                                                    INTR-177
      Z(1)=Z(1)+FAM(I,9)*Z(3)-FAM(I,10)*Z(4)                            INTR-178
      Z(2)=Z(2)+Z(3)*FAM(I,10)+Z(4)*FAM(I,9)                            INTR-179
   24 IF (I.NE.NCI) GO TO 25                                            INTR-180
      Z(1)=Z(1)-FAM(I,9)                                                INTR-181
      Z(2)=Z(2)-FAM(I,10)                                               INTR-182
      Z(3)=Z(3)-1.D0                                                    INTR-183
   25 IF (LO(23).AND.(.NOT.LO(106))) LO(104)=.TRUE.                     INTR-184
C COMPUTE SECOND MEMBERS AND SOLVE THE INHOMOGENEOUS EQUATIONS.         INTR-185
      CALL INRI(WRE,FR(1,1,I),GR,W,NVI,FAM(I,5),X,PAD(1,1,I),PAD(1,1,NC+INTR-186
     11),KITER,KAB,NC,ISM,IPE(1,2),I,NAT,AT,VR,FAM(I,NCI1),FAM(I,NCI2),DINTR-187
     2FLOAT(MC(I,6)+1),V(1,1,I1),LO,Z)                                  INTR-188
      IF (IPE(I,2).GE.ISM) GO TO 26                                     INTR-189
      IF (I.NE.NCI.AND.LO(30)) IPE(I,2)=ISM+1                           INTR-190
      AMAX=DMAX1(AMAX,DABS(FAM(I,NCI1))+DABS(FAM(I,NCI2)))              INTR-191
      GO TO 27                                                          INTR-192
   26 FAM(I,NCI1)=-Z(1)                                                 INTR-193
      FAM(I,NCI2)=-Z(2)                                                 INTR-194
      IF (LO(57)) WRITE (MW,1002)                                       INTR-195
      IF (KITER.EQ.1.AND.LO(44)) NNI=NNI+1                              INTR-196
      IF (I.EQ.NCI) IPE(NCI,2)=IPE(NCI,1)                               INTR-197
   27 CONTINUE                                                          INTR-198
      IF (LO(74)) CALL HORA                                             INTR-199
      IF (LO(104).OR.(AMAX.GT.1.D10.AND.LO(23))) GO TO 29               INTR-200
   28 CONTINUE                                                          INTR-201
      KITER=MIN0(KITER,ITERM)                                           INTR-202
   29 IF (LO(53)) WRITE (MW,1003) KITER                                 INTR-203
   30 CONTINUE                                                          INTR-204
      NNI=NNI/NCIN                                                      INTR-205
      RETURN                                                            INTR-206
 1000 FORMAT (2D30.15,5X,2I5)                                           INTR-207
 1001 FORMAT (' MAXIMUM',D15.6,' OBTAINED IN PREVIOUS ITERATION. LAST ONINTR-208
     1E IS',I3)                                                         INTR-209
 1002 FORMAT (5X,'THE INHOMOGENEOUS TERM IS NEGLECTED.')                INTR-210
 1003 FORMAT (5X,I5,' ITERATIONS.')                                     INTR-211
      END                                                               INTR-212
