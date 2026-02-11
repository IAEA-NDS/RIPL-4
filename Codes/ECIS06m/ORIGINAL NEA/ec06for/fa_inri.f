C 27/06/06                                                      ECIS06  INRI-000
      SUBROUTINE INRI(W,P,Q,WW,NVI,FAM,X,PAD1,PAD2,KITER,KAB,NC,ISM,IPI,INRI-001
     1J,NAT,AT,VR,FAR,FAI,CC,V,LO,Z)                                    INRI-002
C  E. C. I. S. METHOD: INTEGRATION OF A SINGLE INHOMOGENEOUS EQUATION   INRI-003
C  BY THE NUMEROV METHOD      - DIRAC EQUATION -                        INRI-004
C INPUT:     P:       UNCOUPLED SOLUTIONS.                              INRI-005
C            Q:       COUPLED SOLUTION.                                 INRI-006
C            WW:      COUPLING BETWEEN EQUATIONS.                       INRI-007
C            NVI:     TABLE OF ADDRESSES OF COUPLINGS.                  INRI-008
C            FAM:     WAVE NUMBER.                                      INRI-009
C            KITER:   CURRENT ITERATION NUMBER.                         INRI-010
C            KAB:     DIMENSION OF TABLE NVI.                           INRI-011
C            NC:      NUMBER OF COUPLED CHANNELS.                       INRI-012
C            ISM:     NUMBER OF RADIAL POINTS.                          INRI-013
C            IPI:     THE FUNCTION IS NEGLIGIBLE FOR R < IPI*H.         INRI-014
C            J:       CHANNEL NUMBER OF THE EQUATION.                   INRI-015
C            NAT,AT:  COEFFICIENTS AND ADDRESSES OF COUPLINGS.          INRI-016
C            VR:      COUPLING POTENTIALS.                              INRI-017
C            CC:      EIGENVALUE OF L.S + 1.                            INRI-018
C            V:       D(R) AND TENSOR POTENTIALS FOR H**4 CORRECTION.   INRI-019
C            Z:       COULOMB INTEGRAL FOR CORRECTIONS.                 INRI-020
C            LO(I):   LOGICAL CONTROLS:                                 INRI-021
C               LO(22) =.TRUE. NO USE OF PADE APPROXIMANTS.             INRI-022
C               LO(30) =.TRUE. PURE DWBA CALCULATION.                   INRI-023
C               LO(57) =.TRUE. OUTPUT PHASE-SHIFTS AT EACH ITERATION.   INRI-024
C               LO(104)=.TRUE. CONVERGENCE IS OBTAINED IN THE ITERATION.INRI-025
C               LO(105)=.TRUE. CONVERGENCE OBTAINED FOR THIS EQUATION.  INRI-026
C               LO(107)=.TRUE. ALL THE COUPLINGS CALCULATED BEFOREHAND. INRI-027
C OUTPUT:    FAR,FAI: REAL/IMAGINARY PART OF SCATTERING COEFFICIENTS.   INRI-028
C WORKING AREAS:                                                        INRI-029
C            W(ISM,4):SECOND MEMBERS.                                   INRI-030
C            X:       INTEGRAL OF REGULAR FUNCTION WITH SECOND MEMBER.  INRI-031
C            PAD1:    ITERATION RESULTS FOR THE CHANNEL J.              INRI-032
C            PAD2:    WORKING FIELD OF PADE,TWICE LONGER THAN PAD1.     INRI-033
C                                                                       INRI-034
C FOR THE COMMON  /CONVE/ SEE CALX.                                     INRI-035
C                                                                       INRI-036
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /CONVE/:                     INRI-037
C  H:         STEP SIZE FOR INTEGRATION.                                INRI-038
C  EITER:     CONVERGENCE CRITERION FOR S-MATRIX.                       INRI-039
C  ACONV:     CONVERGENCE CRITERION FOR POTENTIAL AND FUNCTION.         INRI-040
C  HCONV:     ACONV*H*H CONVERGENCE CRITERION FOR SECOND MEMBERS.       INRI-041
C   USED:     H,EITER,ACONV,HCONV.                                      INRI-042
C                                                                       INRI-043
C***********************************************************************INRI-044
      IMPLICIT REAL*8 (A-H,O-Z)                                         INRI-045
      LOGICAL LO(150)                                                   INRI-046
      DIMENSION W(ISM,4),P(ISM,8),Q(ISM,4,*),WW(ISM,4,*),NVI(KAB,KAB,4),INRI-047
     1FAM(*),X(2,*),PAD1(2,*),PAD2(*),IPI(*),NAT(6,*),AT(3,*),VR(ISM,4,*INRI-048
     2),FAR(*),FAI(*),V(ISM,14),Z(4),WX(4)                              INRI-049
      COMMON /CONVE/ H,BJM,EITER,ACONV,CONJ,HCONV                       INRI-050
      COMMON /INOUT/ MR,MW,MS                                           INRI-051
      ISM1=ISM+1                                                        INRI-052
C PUT ZERO IN THE INHOMOGENEOUS TERMS.                                  INRI-053
      DO 2 K=1,4                                                        INRI-054
      DO 1 IS=1,ISM                                                     INRI-055
    1 W(IS,K)=0.D0                                                      INRI-056
    2 CONTINUE                                                          INRI-057
      IDP=ISM1                                                          INRI-058
      DO 15 IC=1,NC                                                     INRI-059
      IF (IPI(IC).GT.ISM) GO TO 15                                      INRI-060
      K3=IPI(IC)                                                        INRI-061
      IDP=MIN0(IDP,K3)                                                  INRI-062
      IF (LO(107)) GO TO 11                                             INRI-063
C NO PREVIOUS CALCULATION OF COUPLING POTENTIALS.                       INRI-064
      L1=NVI(J,IC,1)                                                    INRI-065
      L2=NVI(J,IC,2)                                                    INRI-066
      L3=NVI(J,IC,3)                                                    INRI-067
      IF (L1.GT.L2) GO TO 6                                             INRI-068
C SCALAR AND VECTOR TERMS OF COUPLING POTENTIALS.                       INRI-069
      DO 5 IS=K3,ISM                                                    INRI-070
      DO 3 I=1,4                                                        INRI-071
    3 WX(I)=0.D0                                                        INRI-072
      DO 4 K1=L1,L2                                                     INRI-073
      K=NAT(1,K1)                                                       INRI-074
      WX(1)=WX(1)+AT(2,K1)*VR(IS,1,K)                                   INRI-075
      WX(2)=WX(2)+AT(2,K1)*VR(IS,2,K)                                   INRI-076
      WX(3)=WX(3)+AT(3,K1)*VR(IS,3,K)                                   INRI-077
    4 WX(4)=WX(4)+AT(3,K1)*VR(IS,4,K)                                   INRI-078
      W(IS,1)=W(IS,1)+WX(1)*Q(IS,1,IC)-WX(2)*Q(IS,2,IC)                 INRI-079
      W(IS,2)=W(IS,2)+WX(1)*Q(IS,2,IC)+WX(2)*Q(IS,1,IC)                 INRI-080
      W(IS,3)=W(IS,3)+WX(3)*Q(IS,3,IC)-WX(4)*Q(IS,4,IC)                 INRI-081
    5 W(IS,4)=W(IS,4)+WX(3)*Q(IS,4,IC)+WX(4)*Q(IS,3,IC)                 INRI-082
    6 IF (L2.GE.L3) GO TO 13                                            INRI-083
C TENSOR TERMS OF COUPLING POTENTIALS.                                  INRI-084
      L2=L2+1                                                           INRI-085
      DO 10 IS=K3,ISM                                                   INRI-086
      DO 7 I=1,4                                                        INRI-087
    7 WX(I)=0.D0                                                        INRI-088
      DO 8 K1=L2,L3                                                     INRI-089
      K=NAT(1,K1)                                                       INRI-090
      WX(1)=WX(1)+AT(2,K1)*VR(IS,1,K)                                   INRI-091
      WX(2)=WX(2)+AT(2,K1)*VR(IS,2,K)                                   INRI-092
      WX(3)=WX(3)+AT(3,K1)*VR(IS,3,K)                                   INRI-093
    8 WX(4)=WX(4)+AT(3,K1)*VR(IS,4,K)                                   INRI-094
      IF (J.GT.IC) GO TO 9                                              INRI-095
      WX(3)=-WX(3)                                                      INRI-096
      WX(4)=-WX(4)                                                      INRI-097
    9 W(IS,1)=W(IS,1)+(WX(1)+WX(3))*Q(IS,3,IC)-(WX(2)+WX(4))*Q(IS,4,IC) INRI-098
      W(IS,2)=W(IS,2)+(WX(1)+WX(3))*Q(IS,4,IC)+(WX(2)+WX(4))*Q(IS,3,IC) INRI-099
      W(IS,3)=W(IS,3)+(WX(1)-WX(3))*Q(IS,1,IC)-(WX(2)-WX(4))*Q(IS,2,IC) INRI-100
   10 W(IS,4)=W(IS,4)+(WX(1)-WX(3))*Q(IS,2,IC)+(WX(2)-WX(4))*Q(IS,1,IC) INRI-101
      GO TO 15                                                          INRI-102
C COUPLING POTENTIALS COMPUTED IN INTR.                                 INRI-103
   11 K1=NVI(IC,J,3)                                                    INRI-104
      IF (K1.EQ.0) GO TO 13                                             INRI-105
C SCALAR AND VECTOR TERMS OF COUPLING POTENTIALS.                       INRI-106
      DO 12 IS=K3,ISM                                                   INRI-107
      W(IS,1)=W(IS,1)+WW(IS,1,K1)*Q(IS,1,IC)-WW(IS,2,K1)*Q(IS,2,IC)     INRI-108
      W(IS,2)=W(IS,2)+WW(IS,1,K1)*Q(IS,2,IC)+WW(IS,2,K1)*Q(IS,1,IC)     INRI-109
      W(IS,3)=W(IS,3)+WW(IS,3,K1)*Q(IS,3,IC)-WW(IS,4,K1)*Q(IS,4,IC)     INRI-110
   12 W(IS,4)=W(IS,4)+WW(IS,3,K1)*Q(IS,4,IC)+WW(IS,4,K1)*Q(IS,3,IC)     INRI-111
   13 K1=NVI(IC,J,4)                                                    INRI-112
C TENSOR TERMS OF COUPLING POTENTIALS.                                  INRI-113
      IF (K1.EQ.0) GO TO 15                                             INRI-114
      L1=1                                                              INRI-115
      IF (J.LT.IC) L1=3                                                 INRI-116
      L2=4-L1                                                           INRI-117
      DO 14 IS=K3,ISM                                                   INRI-118
      W(IS,1)=W(IS,1)+WW(IS,L1,K1)*Q(IS,3,IC)-WW(IS,L1+1,K1)*Q(IS,4,IC) INRI-119
      W(IS,2)=W(IS,2)+WW(IS,L1,K1)*Q(IS,4,IC)+WW(IS,L1+1,K1)*Q(IS,3,IC) INRI-120
      W(IS,3)=W(IS,3)+WW(IS,L2,K1)*Q(IS,1,IC)-WW(IS,L2+1,K1)*Q(IS,2,IC) INRI-121
   14 W(IS,4)=W(IS,4)+WW(IS,L2,K1)*Q(IS,2,IC)+WW(IS,L2+1,K1)*Q(IS,1,IC) INRI-122
   15 CONTINUE                                                          INRI-123
      IPI(J)=IDP                                                        INRI-124
C SEARCH FOR THE FIRST NON NEGLIGIBLE VALUE.                            INRI-125
      IF (IDP.GT.ISM) GO TO 28                                          INRI-126
      DO 16 IS=IDP,ISM                                                  INRI-127
      IF ((DABS(W(IS,1))+DABS(W(IS,2))+DABS(W(IS,3))+DABS(W(IS,4))).GT.HINRI-128
     1CONV) GO TO 17                                                    INRI-129
   16 CONTINUE                                                          INRI-130
   17 IPI(J)=IS                                                         INRI-131
      IF (IPI(J).GE.ISM) GO TO 28                                       INRI-132
      KT=MAX0(IPI(J),1)                                                 INRI-133
C INTEGRAL OF THE REGULAR FUNCTION WITH THE SECOND MEMBER.              INRI-134
      X(1,KT)=P(KT,2)*W(KT,2)-P(KT,1)*W(KT,1)+P(KT,6)*W(KT,4)-P(KT,5)*W(INRI-135
     1KT,3)                                                             INRI-136
      X(2,KT)=-P(KT,1)*W(KT,2)-P(KT,2)*W(KT,1)-P(KT,5)*W(KT,4)-P(KT,6)*WINRI-137
     1(KT,3)                                                            INRI-138
      KKT=KT+1                                                          INRI-139
      IF (KKT.GT.ISM) GO TO 19                                          INRI-140
      DO 18 IS=KKT,ISM                                                  INRI-141
      X(1,IS)=X(1,IS-1)-P(IS,1)*W(IS,1)+P(IS,2)*W(IS,2)-P(IS,5)*W(IS,3)+INRI-142
     1P(IS,6)*W(IS,4)                                                   INRI-143
   18 X(2,IS)=X(2,IS-1)-P(IS,1)*W(IS,2)-P(IS,2)*W(IS,1)-P(IS,5)*W(IS,4)-INRI-144
     1P(IS,6)*W(IS,3)                                                   INRI-145
C SCATTERING COEFFICIENTS.                                              INRI-146
   19 HX=Z(3)*H                                                         INRI-147
      HY=Z(4)*H                                                         INRI-148
      A=X(1,ISM)/FAM(1)-Z(1)                                            INRI-149
      B=X(2,ISM)/FAM(1)-Z(2)                                            INRI-150
      IF (LO(30)) GO TO 27                                              INRI-151
      IST=ISM+KT                                                        INRI-152
      R=H*DFLOAT(ISM)                                                   INRI-153
      DO 20 IS=KT,ISM                                                   INRI-154
      JS=IST-IS                                                         INRI-155
      Q(JS,1,J)=(X(1,JS)*P(JS,3)-X(2,JS)*P(JS,4)-HX*P(JS,1)+HY*P(JS,2)-.INRI-156
     15D0*W(JS,3))/H+(V(JS,5)*W(JS,1)-V(JS,6)*W(JS,2)+W(JS,3)*(CC/R+V(JSINRI-157
     2,13))-W(JS,4)*V(JS,14))/12.D0                                     INRI-158
      Q(JS,2,J)=(X(1,JS)*P(JS,4)+X(2,JS)*P(JS,3)-HY*P(JS,1)-HX*P(JS,2)-.INRI-159
     15D0*W(JS,4))/H+(V(JS,5)*W(JS,2)+V(JS,6)*W(JS,1)+W(JS,4)*(CC/R+V(JSINRI-160
     2,13))+W(JS,3)*V(JS,14))/12.D0                                     INRI-161
      Q(JS,3,J)=(X(1,JS)*P(JS,7)-X(2,JS)*P(JS,8)-HX*P(JS,5)+HY*P(JS,6)+.INRI-162
     15D0*W(JS,1))/H+(V(JS,7)*W(JS,3)-V(JS,8)*W(JS,4)+W(JS,1)*(CC/R+V(JSINRI-163
     2,13))-W(JS,2)*V(JS,14))/12.D0                                     INRI-164
      Q(JS,4,J)=(X(1,JS)*P(JS,8)+X(2,JS)*P(JS,7)-HY*P(JS,5)-HX*P(JS,6)+.INRI-165
     15D0*W(JS,2))/H+(V(JS,7)*W(JS,4)+V(JS,8)*W(JS,3)+W(JS,2)*(CC/R+V(JSINRI-166
     2,13))+W(JS,1)*V(JS,14))/12.D0                                     INRI-167
      HX=HX+P(JS,3)*W(JS,1)-P(JS,4)*W(JS,2)+P(JS,7)*W(JS,3)-P(JS,8)*W(JSINRI-168
     1,4)                                                               INRI-169
      HY=HY+P(JS,3)*W(JS,2)+P(JS,4)*W(JS,1)+P(JS,7)*W(JS,4)+P(JS,8)*W(JSINRI-170
     1,3)                                                               INRI-171
   20 R=R-H                                                             INRI-172
C CORRECTIONS OF ORDER H**4.                                            INRI-173
      KT1=KT+3                                                          INRI-174
      ISM3=ISM-3                                                        INRI-175
      IF (KT1.GT.ISM3) GO TO 22                                         INRI-176
      HZ=H*720.D0                                                       INRI-177
      DO 21 IS=KT1,ISM3                                                 INRI-178
      Q(IS,1,J)=Q(IS,1,J)-(45.D0*(W(IS+1,3)-W(IS-1,3))-9.D0*(W(IS+2,3)-WINRI-179
     1(IS-2,3))+W(IS+3,3)-W(IS-3,3))/HZ                                 INRI-180
      Q(IS,2,J)=Q(IS,2,J)-(45.D0*(W(IS+1,4)-W(IS-1,4))-9.D0*(W(IS+2,4)-WINRI-181
     1(IS-2,4))+W(IS+3,4)-W(IS-3,4))/HZ                                 INRI-182
      Q(IS,3,J)=Q(IS,3,J)+(45.D0*(W(IS+1,1)-W(IS-1,1))-9.D0*(W(IS+2,1)-WINRI-183
     1(IS-2,1))+W(IS+3,1)-W(IS-3,1))/HZ                                 INRI-184
   21 Q(IS,4,J)=Q(IS,4,J)+(45.D0*(W(IS+1,2)-W(IS-1,2))-9.D0*(W(IS+2,2)-WINRI-185
     1(IS-2,2))+W(IS+3,2)-W(IS-3,2))/HZ                                 INRI-186
   22 IF (KT.EQ.1) GO TO 24                                             INRI-187
      K=KT-1                                                            INRI-188
      DO 23 IS=1,K                                                      INRI-189
      Q(IS,1,J)=-(HX*P(IS,1)-HY*P(IS,2))/H                              INRI-190
      Q(IS,2,J)=-(HY*P(IS,1)+HX*P(IS,2))/H                              INRI-191
      Q(IS,3,J)=-(HX*P(IS,5)-HY*P(IS,6))/H                              INRI-192
   23 Q(IS,4,J)=-(HY*P(IS,5)+HX*P(IS,6))/H                              INRI-193
   24 DO 25 KT=1,ISM                                                    INRI-194
      IF ((DABS(Q(KT,1,J))+DABS(Q(KT,2,J))+DABS(Q(KT,3,J))+DABS(Q(KT,4,JINRI-195
     1))).GT.ACONV) GO TO 26                                            INRI-196
   25 CONTINUE                                                          INRI-197
   26 IPI(J)=KT                                                         INRI-198
      IF (LO(22)) GO TO 27                                              INRI-199
      PAD1(1,KITER)=A                                                   INRI-200
      PAD1(2,KITER)=B                                                   INRI-201
C  TEST OF CONVERGENCE.                                                 INRI-202
   27 LO(105)=(DABS(B-FAI(1)).LE.EITER.AND.(DABS(A-FAR(1)).LE.EITER))   INRI-203
      IF ((.NOT.LO(22)).AND.LO(104).AND.KITER.GT.3.AND.(.NOT.LO(105))) CINRI-204
     1ALL PADE(PAD1,PAD2,KITER,A,B,EITER,1.D0,0.D0,LO)                  INRI-205
      LO(104)=LO(104).AND.LO(105)                                       INRI-206
      FAR(1)=A                                                          INRI-207
      FAI(1)=B                                                          INRI-208
      IF (LO(57)) WRITE (MW,1000) J,FAR(1),FAI(1),KITER,KT              INRI-209
      RETURN                                                            INRI-210
   28 IF (LO(22)) RETURN                                                INRI-211
      PAD1(1,KITER)=0.D0                                                INRI-212
      PAD1(2,KITER)=0.D0                                                INRI-213
      RETURN                                                            INRI-214
 1000 FORMAT (5X,I5,2D30.15,I10,10X,I5)                                 INRI-215
      END                                                               INRI-216
