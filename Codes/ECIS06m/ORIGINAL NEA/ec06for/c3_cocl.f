C 31/08/06                                                      ECIS06  COCL-000
      SUBROUTINE COCL(ETA,R,L,G,GD,F,FD,SIGMA)                          COCL-001
C COMPUTATION OF LOGARITHMIC DERIVATIVE AT THE MATCHING RADIUS OF       COCL-002
C INCREASING AND DECREASING SOLUTIONS WITH WRONSKIAN EQUAL TO UNITY     COCL-003
C FOR CLOSED CHANNELS..                                                 COCL-004
C INPUT:     ETA:     COULOMB PARAMETER.                                COCL-005
C            R:       |K|*R VALUE.                                      COCL-006
C            L:       MAXIMUM L VALUE.                                  COCL-007
C OUTPUT:    F:       DECREASING SOLUTION AT (ETA,RHO) FOR I = 1 TO L+1.COCL-008
C            FD:      DERIVATIVE OF F(I) FOR I = 1 TO L+1.              COCL-009
C            G:       INCREASING SOLUTION FOR I = 1 TO L+1.             COCL-010
C            GD:      DERIVATIVE OF F(I) FOR I = 1 TO L+1.              COCL-011
C            SIGMA:   RETURNS 0 AT THE SAME NUMBER.                     COCL-012
C                                                                       COCL-013
C FOR NEGATIVE VALUES OF ETA, THE LOGARITHMIC DERIVATIVE OF THE         COCL-014
C DECREASING SOLUTION IS COMPUTED, THE LOGARITHMIC DERIVATIVE OF THE    COCL-015
C INCREASING FUNCTION IS ASSUMED OPPOSITE. IF A NEGATIVE ETA DIFFERS    COCL-016
C BY LESS THAN 1/4 PERCENT OF AN INTEGER THE NEXT INTEGER VALUE IS      COCL-017
C USED (LAGUERRE POLYNOMIAL); IN OTHER CASES, BACKWARDS INTEGRATION     COCL-018
C STARTING BEYOND THE OSCILLATIONS ESTIMATED UP TO 2*|ETA| IS USED.     COCL-019
C                                                                       COCL-020
C THE FUNCTIONS ARE RENORMALISED TO F=G AND SUCH THAT  F*GD-G*FD=1      COCL-021
C***********************************************************************COCL-022
      IMPLICIT REAL*8 (A-H,O-Z)                                         COCL-023
      DIMENSION G(*),GD(*),F(*),FD(*),SIGMA(*),S(7)                     COCL-024
      DATA GAMA /0.577215664901533D0/                                   COCL-025
C FOR L=0.                                                              COCL-026
      IF (ETA.LE.-1.D-6) GO TO 17                                       COCL-027
      IF (ETA.GT.1.D-6) GO TO 1                                         COCL-028
C NO COULOMB POTENTIAL.                                                 COCL-029
      FP=-1.D0                                                          COCL-030
      GP=1.D0                                                           COCL-031
      GO TO 14                                                          COCL-032
C POSITIVE VALUES OF ETA.                                               COCL-033
    1 IF ((ETA+1.D0)*R.GT.8.D0) GO TO 9                                 COCL-034
C SERIES EXPANSION.                                                     COCL-035
      R2=-R*R                                                           COCL-036
      ET2=ETA+ETA                                                       COCL-037
      U0=0.D0                                                           COCL-038
      U1=R                                                              COCL-039
      V0=1.D0                                                           COCL-040
      V1=0.D0                                                           COCL-041
      U=U0+U1                                                           COCL-042
      V=V0+V1                                                           COCL-043
      UP=1.D0                                                           COCL-044
      VP=0.D0                                                           COCL-045
      DO 3 N=2,10000                                                    COCL-046
      XN=DFLOAT(N)                                                      COCL-047
      XN1=XN*(XN-1.D0)                                                  COCL-048
      U2=(ET2*R*U1-R2*U0)/XN1                                           COCL-049
      U=U+U2                                                            COCL-050
      V2=(ET2*R*V1-R2*V0-ET2*(2.D0*XN-1.D0)*U2)/XN1                     COCL-051
      V=V+V2                                                            COCL-052
      UP=UP+XN*U2/R                                                     COCL-053
      VP=VP+XN*V2/R                                                     COCL-054
      IF (DABS(U2).GT.1.D-16*DABS(U)) GO TO 2                           COCL-055
      IF (DABS(V2).LE.1.D-16*DABS(V)) GO TO 4                           COCL-056
    2 U0=U1                                                             COCL-057
      U1=U2                                                             COCL-058
      V0=V1                                                             COCL-059
    3 V1=V2                                                             COCL-060
    4 XX=DABS(ETA)                                                      COCL-061
      IF (XX.LT.1.D-8) GO TO 7                                          COCL-062
      Y=XX                                                              COCL-063
      K=0                                                               COCL-064
      IF (XX.LE.7.5D0) K=IDINT(8.5D0-XX)                                COCL-065
      X=1.D0+XX+DFLOAT(K)                                               COCL-066
      UU=1.D0/X**2                                                      COCL-067
      PSR=DLOG(X)-.5D0/X-UU/12.D0+UU**2/120.D0-UU**3/252.D0+UU**4/240.D0COCL-068
     1-UU**5/132.D0+UU**6*691.D0/32760.D0                               COCL-069
      IF (K.EQ.0) GO TO 6                                               COCL-070
      DO 5 I=1,K                                                        COCL-071
    5 PSR=PSR-1.D0/(X-DFLOAT(I))                                        COCL-072
    6 PSR=PSR-.5D0/Y+2.D0*GAMA-1.D0                                     COCL-073
      GO TO 8                                                           COCL-074
    7 PSR=GAMA-1.D0                                                     COCL-075
    8 CE=ET2*(PSR+DLOG(2.D0*R))                                         COCL-076
      FP=(VP+UP*CE+ET2*U/R)/(V+U*CE)                                    COCL-077
      GO TO 14                                                          COCL-078
    9 IF (R.LT.10.D0*(ETA+1.D0)) GO TO 11                               COCL-079
C ASYMPTOTIC EXPANSION.                                                 COCL-080
      C=1.D0/R                                                          COCL-081
      A=1.D0                                                            COCL-082
      B=A                                                               COCL-083
      D=0.D0                                                            COCL-084
      DO 10 M=1,26                                                      COCL-085
      AM=DFLOAT(M)                                                      COCL-086
      A=-A*0.5D0*(ETA+AM-1.D0)*(ETA+AM)*C/AM                            COCL-087
      B=B+A                                                             COCL-088
   10 D=D-A*AM*C                                                        COCL-089
      FP=D/B-1.D0-ETA/R                                                 COCL-090
      GO TO 14                                                          COCL-091
C LONG RANGE INTEGRATION.                                               COCL-092
   11 H=DMIN1(.001953125D0,.25D0*R)                                     COCL-093
      N=1+IDINT(10.D0/H)                                                COCL-094
      S(6)=DEXP(-H)                                                     COCL-095
      S(7)=1.D0                                                         COCL-096
      V2=H**2*(1.D0+2.D0*ETA/(R+H*DFLOAT(N)))/12.D0                     COCL-097
      V3=H**2*(1.D0+2.D0*ETA/(R+H*DFLOAT(N-1)))/12.D0                   COCL-098
      DO 13 I=N-2,-3,-1                                                 COCL-099
      DO 12 J=1,6                                                       COCL-100
   12 S(J)=S(J+1)/S(7)                                                  COCL-101
      V1=V2                                                             COCL-102
      V2=V3                                                             COCL-103
      V3=H**2*(1.D0+2.D0*ETA/(R+H*DFLOAT(I)))/12.D0                     COCL-104
   13 S(7)=(S(6)*(2.D0+10.D0*V2)-S(5)*(1.D0-V1))/(1.D0-V3)              COCL-105
      FP=((S(1)-S(7))/60.D0+.15D0*(S(6)-S(2))+.75D0*(S(3)-S(5)))/H/S(4) COCL-106
   14 XX=1.D0                                                           COCL-107
      MP=L+25+IDINT(5.D0*DABS(ETA))                                     COCL-108
      DO 15 M=MP,1,-1                                                   COCL-109
      XM=DFLOAT(M)                                                      COCL-110
      A=ETA/XM                                                          COCL-111
      B=A+XM/R                                                          COCL-112
      XX=-(A*A-1.D0)/(B+XX)+B                                           COCL-113
      IF (M.LE.L) GD(M)=XX                                              COCL-114
   15 CONTINUE                                                          COCL-115
      DO 16 M=1,L+1                                                     COCL-116
      X=DSQRT(DABS(GD(M)-FP))                                           COCL-117
      XM=DFLOAT(M)                                                      COCL-118
      F(M)=1.D0/X                                                       COCL-119
      FD(M)=FP/X                                                        COCL-120
      G(M)=1.D0/X                                                       COCL-121
      GD(M)=GD(M)/X                                                     COCL-122
      A=ETA/XM                                                          COCL-123
      B=A+XM/R                                                          COCL-124
   16 FP=(A*A-1.D0)/(B-FD(M)/F(M))-B                                    COCL-125
      GO TO 25                                                          COCL-126
C NEGATIVE VALUES OF ETA.                                               COCL-127
   17 LP=IDINT(-ETA+.5D0)                                               COCL-128
      IF (DABS((DFLOAT(LP)+ETA)/ETA).LT..0025D0) GO TO 22               COCL-129
      LL=MIN0(L,IDINT(-ETA))                                            COCL-130
      RP=DMAX1(-2.D0*ETA-R,0.D0)                                        COCL-131
      DO 20 M=0,LL                                                      COCL-132
      AL=M*(M+1)                                                        COCL-133
      H=DMIN1(.001953125D0,R/DFLOAT(4+M))                               COCL-134
      N=IDINT(1.D0+(40.D0+RP)/H)                                        COCL-135
      RR=R+DFLOAT(N)*H                                                  COCL-136
      V2=H**2*(1.D0+AL/RR**2+2.D0*ETA/RR)/12.D0                         COCL-137
      RR=R+DFLOAT(N-1)*H                                                COCL-138
      V3=H**2*(1.D0+AL/RR**2+2.D0*ETA/RR)/12.D0                         COCL-139
      S(6)=DEXP(-H)                                                     COCL-140
      S(7)=1.D0                                                         COCL-141
      DO 19 I=N-2,-3,-1                                                 COCL-142
      DO 18 J=1,6                                                       COCL-143
   18 S(J)=S(J+1)/S(7)                                                  COCL-144
      V1=V2                                                             COCL-145
      V2=V3                                                             COCL-146
      RR=R+DFLOAT(I)*H                                                  COCL-147
      V3=H**2*(1.D0+AL/RR**2+2.D0*ETA/RR)/12.D0                         COCL-148
   19 S(7)=(S(6)*(2.D0+10.D0*V2)-S(5)*(1.D0-V1))/(1.D0-V3)              COCL-149
      FP=((S(1)-S(7))/60.D0+.15D0*(S(6)-S(2))+.75D0*(S(3)-S(5)))/H      COCL-150
      FF=S(4)                                                           COCL-151
      GP=1.D0                                                           COCL-152
      IF (FP.NE.0.D0) GP=-FP                                            COCL-153
      GG=1.D0                                                           COCL-154
      IF (FF.NE.0.D0) GG=FF                                             COCL-155
      X=DSQRT(DABS(GP*FF-GG*FP))                                        COCL-156
      F(M+1)=FF/X                                                       COCL-157
      FD(M+1)=FP/X                                                      COCL-158
      G(M+1)=GG/X                                                       COCL-159
      GD(M+1)=GP/X                                                      COCL-160
   20 CONTINUE                                                          COCL-161
      IF (LL.EQ.L) GO TO 25                                             COCL-162
      FP=FD(LL+1)/F(LL+1)                                               COCL-163
      DO 21 M=LL+1,L                                                    COCL-164
      XM=DFLOAT(M)                                                      COCL-165
      A=ETA/XM                                                          COCL-166
      B=A+XM/R                                                          COCL-167
      FP=(A*A-1.D0)/(B-FD(M)/F(M))-B                                    COCL-168
      X=DSQRT(DABS(2.D0*FP))                                            COCL-169
      F(M+1)=1.D0/X                                                     COCL-170
      FD(M+1)=FP/X                                                      COCL-171
      G(M+1)=1.D0/X                                                     COCL-172
   21 GD(M+1)=-FP/X                                                     COCL-173
      GO TO 25                                                          COCL-174
   22 FF=1.D0                                                           COCL-175
      FP=1.D0                                                           COCL-176
      MP=L+50-IDINT(5.D0*ETA)                                           COCL-177
      DO 24 M=MP,1,-1                                                   COCL-178
      FFM=DFLOAT(M)                                                     COCL-179
      A=ETA/FFM                                                         COCL-180
      B=A+FFM/R                                                         COCL-181
      Z=FP+B*FF                                                         COCL-182
      FP=B*Z-(A*A-1.D0)*FF                                              COCL-183
      FF=Z                                                              COCL-184
      IF (DABS(FF).LT.1.D0) GO TO 23                                    COCL-185
      FP=FP/FF                                                          COCL-186
      FF=1.D0                                                           COCL-187
   23 IF (M.GT.L+1) GO TO 24                                            COCL-188
      GP=1.D0                                                           COCL-189
      IF (FP.NE.0.D0) GP=-FP                                            COCL-190
      GG=1.D0                                                           COCL-191
      IF (FF.NE.0.D0) GG=FF                                             COCL-192
      X=DSQRT(DABS(GP*FF-GG*FP))                                        COCL-193
      F(M)=FF/X                                                         COCL-194
      FD(M)=FP/X                                                        COCL-195
      G(M)=GG/X                                                         COCL-196
      GD(M)=GP/X                                                        COCL-197
   24 CONTINUE                                                          COCL-198
   25 DO 26 M=0,L                                                       COCL-199
   26 SIGMA(M+1)=0.D0                                                   COCL-200
      RETURN                                                            COCL-201
      END                                                               COCL-202
