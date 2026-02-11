C 02/02/06                                                      ECIS06  CORA-000
      SUBROUTINE CORA(LQ,L1,L2,EI,FI,VI,WI,B,C,LT)                      CORA-001
C COEFFICIENTS FOR COULOMB INTEGRALS                                    CORA-002
C M(L1,L2,LQ,R)=SUM FROM R TO INFINITY OF H(L1)*K(L2)/R**(LQ+1)         CORA-003
C EXPRESSED AS SUM ON I OF B(I)*M(L+I-2,L+I-2,1) + C1(R)*H(L)*K(L)      CORA-004
C + C2(R)*H(L)*K(L+1) + C3(R)*H(L+1)*K(L) + C4(R)*H(L+1)*K(L+1)         CORA-005
C WHERE H AND K ARE REGULAR OR IRREGULAR COULOMB FUNCTIONS              CORA-006
C L=INTEGER PART OF (L1+L2-LQ+3)/2 AND I RANGES FROM 1 TO 4.            CORA-007
C INPUT:     LQ:      LQ IN THE EXPRESSIONS ABOVE.                      CORA-008
C            L1,L2:   L1, L2 IN THE EXPRESSIONS ABOVE.                  CORA-009
C            EI,EF:   COULOMB PARAMETERS.                               CORA-010
C            VI,WI:   PRODUCT OF R WITH THE WAVE NUMBERS.               CORA-011
C            LT:      .TRUE. FOR ONLY INTEGRALS FROM 0 TO INFINITY.     CORA-012
C OUTPUT:    B(4):    COEFFICIENTS B IN THE EXPRESSIONS ABOVE.          CORA-013
C            C(4):    COEFFICIENTS C.                                   CORA-014
C                                                                       CORA-015
C THESE RESULTS ARE OBTAINED BY EXPRESSING                              CORA-016
C H(L1)*K(L2)/R**(LQ+1)- SUM ON I OF B(I)*H(L+I-2)*K(L+I-2)/R**2        CORA-017
C AS P1*H(L)*K(L) + P2*H(L)*K(L+1) + P3*H(L+1)*K(L) + P4*H(L+1)*Q(L+1)  CORA-018
C WHERE P1,P2,P3,P4 ARE POLYNOMIAL IN 1/R AND IDENTIFICATION TERM BY    CORA-019
C TERM STARTING FROM THE HIGHEST DEGREE IN 1/R WITH THE DERIVATIVE OF   CORA-020
C  Q1*H(L)*K(L) + Q2*H(L)*K(L+1) + Q3*H(L+1)*K(L) + Q4*H(L+1)*Q(L+1)    CORA-021
C THIS AS BEEN DONE USING AMP ( J.-M. DROUFFE, AMP LANGUAGE REFERENCE   CORA-022
C MANUAL - VERSION 6 - NOTE CEA-N-2297 1982).                           CORA-023
C                                                                       CORA-024
C THE NON RELATIVISTIC RESULTS FOR LQ=6 L2-L1=0,2,4 AND ALL THE         CORA-025
C COEFFICIENTS UP TO LQ=5 HAVE BEEN OBTAINED. THIS ALLOW TRANSFER OF    CORA-026
C ANGULAR MOMENTUM UP TO 5 FOR CENTRAL TERM AND 4 FOR SPIN-ORBIT.       CORA-027
C                                                                       CORA-028
C***********************************************************************CORA-029
      IMPLICIT REAL*8 (A-H,O-Z)                                         CORA-030
      DIMENSION B(4),C(4)                                               CORA-031
      LOGICAL LT                                                        CORA-032
      COMMON /INOUT/ MR,MW,MS                                           CORA-033
      DO 1 I=1,4                                                        CORA-034
      B(I)=0.D0                                                         CORA-035
    1 C(I)=0.D0                                                         CORA-036
      IF (L1.LE.L2) GO TO 2                                             CORA-037
      LL=L1-L2+1                                                        CORA-038
      E=FI                                                              CORA-039
      F=EI                                                              CORA-040
      V=WI                                                              CORA-041
      W=VI                                                              CORA-042
      GO TO 3                                                           CORA-043
    2 LL=L2-L1+1                                                        CORA-044
      E=EI                                                              CORA-045
      F=FI                                                              CORA-046
      V=VI                                                              CORA-047
      W=WI                                                              CORA-048
    3 IF ((LQ.GT.6).OR.(LL*(LQ-LL+2).LE.0)) GO TO 71                    CORA-049
      A=DFLOAT((L1+L2-LQ+5)/2)                                          CORA-050
      S=V/W                                                             CORA-051
      X=E*S-F                                                           CORA-052
      IF (DABS(X).LT.1.D-14) X=0.D0                                     CORA-053
      U=S**2-1.D0                                                       CORA-054
      IF (DABS(U).LT.1.D-14) U=0.D0                                     CORA-055
      IF ((X.EQ.0.D0).AND.(E.EQ.0.D0)) GO TO 49                         CORA-056
      P=S**2*(A**2+E**2)                                                CORA-057
      Q=A**2+F**2                                                       CORA-058
      Y=E*S+F                                                           CORA-059
      Z=E**2-F**2                                                       CORA-060
      IF (DABS(Z).LT.1.D-14) Z=0.D0                                     CORA-061
      IF (X.EQ.0.D0.AND.Z.EQ.0.D0) Z=1.D0                               CORA-062
      DEN=2.D0*X**2*(P-Q)+8.D0*Z*S**2                                   CORA-063
      Z=U                                                               CORA-064
      IF (X.EQ.0.D0.AND.Z.EQ.0.D0) Z=-1.D0/(E*F)                        CORA-065
      GO TO ( 4 , 9 , 13 , 20 , 28 , 41 ) , LQ                          CORA-066
    4 GO TO ( 5 , 6 ) , LL                                              CORA-067
C MULTIPOLE LQ=1 L1=L2.                                                 CORA-068
    5 B(1)=1.D0                                                         CORA-069
      RETURN                                                            CORA-070
    6 IF ((E.NE.F).OR.(V.NE.W)) GO TO 7                                 CORA-071
C MULTIPOLE LQ=1 L1=L2-1 WITH E=F AND V=W.                              CORA-072
      B1=DSQRT((A-1.D0)**2+E**2)                                        CORA-073
      B(2)=-(A-.5D0)*Q/(E*B1)                                           CORA-074
      B(3)=(A+.5D0)*Q/(E*B1)                                            CORA-075
      IF (LT) RETURN                                                    CORA-076
      C(1)=(Q+(2.D0*A-1.D0)*E/V)/(2.D0*E*V*B1)                          CORA-077
      C(3)=-DSQRT(Q)/(2.D0*V*A*B1)                                      CORA-078
      C(2)=C(3)*(2.D0*A-1.D0)                                           CORA-079
      C(4)=Q/(2.D0*E*V*B1)                                              CORA-080
      GO TO 70                                                          CORA-081
C MULTIPOLE LQ=1 L1=L2-1.                                               CORA-082
    7 A1=-2.D0*F*(2.D0*(A-1.D0)*Q-A**2*U)                               CORA-083
      A2=F*(4.D0*Q*(A-1.D0)*((6.D0*A+1.D0)*Q-4.D0*A**2+1.D0)+U*((((16.D0CORA-084
     1*A-24.D0)*A+4.D0)*A**2-2.D0)*Q+A**2*(2.D0*A-1.D0-U*(A-1.D0)**2)*(4CORA-085
     2.D0*A+2.D0)))                                                     CORA-086
      A3=-F*(A-1.D0)*((24.D0*A-4.D0)*Q+16.D0*A**2-4.D0+U*(((8.D0*A+8.D0)CORA-087
     1*A+2.D0)*A-2.D0))                                                 CORA-088
      A4=4.D0*(A-1.D0)*F                                                CORA-089
      IF (X.EQ.0.D0) GO TO 8                                            CORA-090
      A1=A1-X*(Q*(A-1.D0)*(2.D0*A+6.D0)+A**2*(4.D0+A*U)+X*((2.D0*A-2.D0)CORA-091
     1*F+A*X))                                                          CORA-092
      A2=A2+X*(Q*((A-1.D0)*(A+5.D0)*(12.D0*A+2.D0)*Q-(((40.D0*A-16.D0)*ACORA-093
     1-34.D0)*A-8.D0)*A-6.D0+(((((8.D0*A+12.D0)*A-38.D0)*A+8.D0)*A+5.D0)CORA-094
     2*A-4.D0)*U+(((16.D0*A+40.D0)*A-44.D0)*A-12.D0)*F*X+(((4.D0*A+20.D0CORA-095
     3)*A-13.D0)*A-8.D0)*X**2)+(2.D0*A+1.D0)*(A**2*(4.D0-8.D0*A+((6.D0*ACORA-096
     4-7.D0)*A+4.D0)*U+A*(A-1.D0)**2*U**2)+X*(((4.D0*A+6.D0)*A-2.D0+(A-1CORA-097
     5.D0)*((4.D0*A-4.D0)*A+2.D0)*U)*F+X*(A*(((2.D0*A-2.D0)*A+1.D0)*U-(4CORA-098
     6.D0*A-6.D0)*A+1.D0)+(4.D0*A-2.D0)*F*X+A*X**2))))                  CORA-099
      A3=A3+X*(A-1.D0)*(((8.D0*A-8.D0)*A+6.D0)*A+2.D0-(2.D0*A+6.D0)*(6.DCORA-100
     10*A-1.D0)*Q-((((4.D0*A+12.D0)*A+5.D0)*A-3.D0)*A-2.D0)*U-2.D0*X*(((CORA-101
     24.D0*A+8.D0)*A-3.D0)*F-X))                                        CORA-102
      A4=A4+X*2.D0*(A**2-1.D0)                                          CORA-103
    8 BD=.5D0*DEN*DSQRT((A-1.D0)**2+E**2)                               CORA-104
      IF (.NOT.LT) C(1)=.5D0-A                                          CORA-105
      GO TO 48                                                          CORA-106
    9 GO TO ( 10 , 11 , 12 ) , LL                                       CORA-107
C MULTIPOLE LQ=2 L1=L2.                                                 CORA-108
   10 A1=2.D0*Z*A**2*Y                                                  CORA-109
      A2=-Z*Y*(2.D0*A+1.D0)*(((A-1.D0)*A-.5D0)*Y**2+2.D0*A**2*(A-1.D0)**CORA-110
     12*(S**2+1.D0))                                                    CORA-111
      A3=2.D0*Z*(2.D0*A+1.D0)*(A-1.D0)**2*Y                             CORA-112
      A4=0.D0                                                           CORA-113
      BD=DEN*(A-1.D0)/DSQRT(S)                                          CORA-114
      IF (.NOT.LT) C(1)=-1.D0                                           CORA-115
      IF (X.EQ.0) GO TO 48                                              CORA-116
      A1=A1+X*(4.D0*F*E*S-2.D0*A*(A+Q+A*S**2+P))                        CORA-117
      A2=A2+X*((A+.5D0)*((((4.D0*A-6.D0)*A+2.D0)*A-1.D0)*(Q+P*S**2)+4.D0CORA-118
     1*(2.D0*A-3.D0)*F*E*S*(Q+P+S**2+1.D0)+(((4.D0*A+2.D0)*A-14.D0)*A+1.CORA-119
     2D0)*(Q*S**2+P)+A**2*((S**4+1.D0)*((2.D0*A-6.D0)*A+5.D0)-S**2*(4.D0CORA-120
     3*(A+1.D0)*A-14.D0)))+P*Q*((8.D0*A-4.D0)*A+8.D0))                  CORA-121
      A3=A3-2.D0*X*(A-1.D0)*((2.D0*A-3.D0)*(F**2+E**2*S**2)+4.D0*(2.D0*ACORA-122
     1+1.D0)*F*E*S+(S**2+1.D0)*(A+1.D0)*((6.D0*A+1.D0)*A-3.D0))         CORA-123
      A4=4.D0*(A-1.D0)*X                                                CORA-124
      GO TO 48                                                          CORA-125
C MULTIPOLE LQ=2 L1=L2-1.                                               CORA-126
   11 B1=X*(A-1.D0)-E*S                                                 CORA-127
      B3=-(A-1.D0)**2                                                   CORA-128
      B2=B3/(A-.5D0)                                                    CORA-129
      BD=-DEN*A*(A-1.D0)*(2.D0*A-3.D0)*DSQRT(S*((A-1.D0)**2+E**2))/(2.D0CORA-130
     1*A-1.D0)                                                          CORA-131
      IF (LT) GO TO 15                                                  CORA-132
      C(1)=B1+2.D0*A*(A-1.D0)**2/W                                      CORA-133
      C(2)=B3                                                           CORA-134
      C(3)=B2                                                           CORA-135
      GO TO 15                                                          CORA-136
C MULTIPOLE LQ=2 L1=L2-2.                                               CORA-137
   12 B1=Q                                                              CORA-138
      B3=-3.D0*E*S-X*(A-1.D0)                                           CORA-139
      B4=(A-1.D0)/(A-.5D0)                                              CORA-140
      BD=-3.D0*DEN*A*DSQRT(S*Q*((A-1.D0)**2+E**2))/(2.D0*A-1.D0)        CORA-141
      IF (LT) GO TO 17                                                  CORA-142
      C(1)=B1                                                           CORA-143
      C(2)=B3-2.D0*A*(A-1.D0)/W                                         CORA-144
      C(4)=B4                                                           CORA-145
      GO TO 17                                                          CORA-146
   13 GO TO ( 14 , 16 , 18 , 19 ) , LL                                  CORA-147
C MULTIPOLE LQ=3 L1=L2.                                                 CORA-148
   14 B1=Y/(A-1.D0)                                                     CORA-149
      B2=1.D0                                                           CORA-150
      B3=1.D0                                                           CORA-151
      BD=DEN*A*(2.D0*A-3.D0)                                            CORA-152
      IF (LT) GO TO 15                                                  CORA-153
      C(1)=B1-2.D0*A/W                                                  CORA-154
      C(2)=B3                                                           CORA-155
      C(3)=B2                                                           CORA-156
   15 A1=-2.D0*A**2*Y*B1*Z                                              CORA-157
      A2=Y*((2.D0*A+1.D0)*((2.D0*(A**2-A)**2*(1.D0+S**2)+Y**2*(A**2-A-.5CORA-158
     1D0))*B1-Y*(P*B2+Q*B3)))*Z                                         CORA-159
      A3=-Y*((2.D0*A-1.D0)*(-Y*(B2+B3))+(4.D0*A+2.D0)*(A-1.D0)**2*B1)*Z CORA-160
      A4=0.D0                                                           CORA-161
      IF (X.EQ.0.D0) GO TO 48                                           CORA-162
      A1=A1+X*2.D0*P*(A*X-2.D0*F)*B2                                    CORA-163
      A2=A2+X*P*(2.D0*Q*((2.D0*A+3.D0)*(A+1.D0)*F+((2.D0*A+3.D0)*A-4.D0)CORA-164
     1*E*S)-2.D0*P*(2.D0*A-3.D0)*(2.D0*A+1.D0)*F+(2.D0*A+1.D0)*((((2.D0*CORA-165
     2A-1.D0)*A-8.D0)*A+5.D0)*F-S*((A+1.D0)*((6.D0*A-9.D0)*A+1.D0)*E-S*(CORA-166
     3(6.D0*A-7.D0)*(A**2-1.D0)*F-S*((((2.D0*A-5.D0)*A+2.D0)*A-1.D0)*E))CORA-167
     4)))*B2                                                            CORA-168
      A3=A3+X*(2.D0*P*((2.D0*A-3.D0)*(A-1.D0)*E*S+F*((2.D0*A-3.D0)*A-4.DCORA-169
     10))-2.D0*Q*E*S*(2.D0*A+3.D0)*(2.D0*A-1.D0)-(2.D0*A-1.D0)*((((2.D0*CORA-170
     2A+5.D0)*A+2.D0)*A+1.D0)*F-S*((6.D0*A+7.D0)*(A**2-1.D0)*E-S*((A-1.DCORA-171
     30)*((6.D0*A+9.D0)*A+1.D0)*F-S*((((2.D0*A+1.D0)*A-8.D0)*A-5.D0)*E))CORA-172
     4)))*B2                                                            CORA-173
      IF (LQ.EQ.3) GO TO 47                                             CORA-174
      A4=A4+X*2.D0*(2.D0*E*S-A*X)*B2                                    CORA-175
      GO TO 46                                                          CORA-176
C MULTIPOLE LQ=3 L1=L2-1.                                               CORA-177
   16 B1=Q/(A-1.D0)                                                     CORA-178
      B3=-X                                                             CORA-179
      B4=1.D0/(A+1.D0)                                                  CORA-180
      BD=-3.D0*DEN*A*DSQRT(Q)                                           CORA-181
      IF (LT) GO TO 17                                                  CORA-182
      C(1)=B1                                                           CORA-183
      C(2)=B3-2.D0*A/W                                                  CORA-184
      C(4)=B4                                                           CORA-185
   17 A1=-2.D0*A**2*Y*B1*Z                                              CORA-186
      A2=Y*((2.D0*A+1.D0)*((2.D0*(A**2-A)**2*(1.D0+S**2)+Y**2*(A**2-A-.5CORA-187
     1D0))*B1-Y*Q*B3)+(4.D0*A-2.D0)*(A+1.D0)**2*P*Q*B4)*Z               CORA-188
      A3=-Y*((2.D0*A-1.D0)*((2.D0*(A**2+A)**2*(1.D0+S**2)+Y**2*(A**2+A-.CORA-189
     15D0))*B4-Y*B3)+(4.D0*A+2.D0)*(A-1.D0)**2*B1)*Z                    CORA-190
      A4=2.D0*A**2*Y*B4*Z                                               CORA-191
      IF (X.EQ.0.D0) GO TO 48                                           CORA-192
      A1=A1-X*4.D0*P*Q*(A+1.D0)*B4                                      CORA-193
      A2=A2+X*2.D0*P*Q*(A+1.D0)*((2.D0*A+3.D0)*(Q+P)+(2.D0*A-1.D0)*(4.D0CORA-194
     1*E*F*S+(1.D0+S**2)*((2.D0*A-4.D0)*A-3.D0)))*B4                    CORA-195
      A3=A3+X*(-4.D0*P*Q*((2.D0*A+1.D0)*A+2.D0)-(2.D0*A-1.D0)*(2.D0*S*E*CORA-196
     1F*(2.D0*A+3.D0)*(Q+P+S**2+1.D0)-A**2*(((A+3.D0)*A+2.5D0)*(1.D0+S**CORA-197
     24)-S**2*((2.D0*A-2.D0)*A-7.D0))+A*(A-1.D0)*(2.D0*A+3.D0)*(1.D0+S**CORA-198
     32)*(P+Q)-(S**2-1.D0)*((2.D0*A+4.D0)*A+.5D0)*(Q-P)))*B4            CORA-199
      A4=A4+X*2.D0*((A-1.D0)*(P+Q)+Y**2)*B4                             CORA-200
      GO TO 46                                                          CORA-201
C MULTIPOLE LQ=3 L1=L2-2.                                               CORA-202
   18 B1=-Q*(2.D0*A-1.D0)*((2.D0*A-3.D0)*(A+1.D0)*X+(A-3.D0)*Y)/(A-1.D0)CORA-203
      B2=6.D0*Q*(A-1.D0)                                                CORA-204
      B3=((A-1.D0)*(6.D0*A-3.D0)*Q+(2.D0*A-1.D0)*((A+1.D0)*X+2.D0*Y)*(2.CORA-205
     1D0*A-3.D0)*X-(A-1.D0)*(6.D0*A-9.D0)*P)                            CORA-206
      B4=-((2.D0*A-1.D0)*X+Y)*(2.D0*A-3.D0)                             CORA-207
      BD=12.D0*A**2*(2.D0*A-3.D0)*DEN*S*DSQRT(Q*((A-1.D0)**2+E**2))     CORA-208
      IF (LT) GO TO 45                                                  CORA-209
      C(1)=B1-6.D0*Q*A*(A-1.D0)*(2.D0*A-1.D0)/W                         CORA-210
      C(2)=B3+A*(4.D0*A-6.D0)*(2.D0*A-1.D0)*((A+1.D0)*X+2.D0*Y+3.D0*A*(ACORA-211
     1-1.D0)/W)/W                                                       CORA-212
      C(3)=B2                                                           CORA-213
      C(4)=B4-6.D0*A*(A-1.D0)*(2.D0*A-3.D0)/W                           CORA-214
      GO TO 45                                                          CORA-215
C MULTIPOLE LQ=3 L1=L2-3.                                               CORA-216
   19 B1=-Q*(3.D0*E*S*Y+A*(3.D0*X*Y+2.D0*(A+1.D0)*(X**2+2.D0*A-1.D0-2.D0CORA-217
     1*(A-1.D0)*S**2)))                                                 CORA-218
      B2=3.D0*Q*(A*X+Y)                                                 CORA-219
      B4=-(3.D0*F*Y+A*(3.D0*X*Y+2.D0*(A-1.D0)*(X**2+2.D0*A+2.D0-(2.D0*A+CORA-220
     11.D0)*S**2)))                                                     CORA-221
      B3=-X*B4*(A+1.D0)+3.D0*Q*(Y+A*X)+A*(3.D0*A+1.5D0)*(1.D0-S**2)*((3.CORA-222
     1D0*A-2.D0)*Y+((A-2.D0)*A+2.D0)*X)                                 CORA-223
      BD=30.D0*A**2*DEN*S*DSQRT(Q*((A-1.D0)**2+E**2)*((A+1.D0)**2+F**2))CORA-224
      IF (LT) GO TO 45                                                  CORA-225
      C(1)=B1-3.D0*Q*(2.D0*A-1.D0)*((6.D0*A+4.D0)*F+(A+1.D0)*(A+4.D0)*X+CORA-226
     14.D0*A*(A**2-1.D0)/W)/W                                           CORA-227
      C(2)=B3+((4.D0*A**2*(A**2-1.D0)*(2.D0*A-1.D0-S**2*(2.D0*A+1.D0))+(CORA-228
     118.D0*A**2-3.D0)*Y**2+A*(12.D0*A**2+3.D0)*X*Y+(A**2-1.D0)*(4.D0*A*CORA-229
     2*2-3.D0)*X**2)+3.D0*A*(4.D0*A**2-1.D0)*(5.D0*A*Y+(A**2+4.D0)*X+4.DCORA-230
     30*A*(A**2-1.D0)/W)/W)/W                                           CORA-231
      C(3)=B2+12.D0*Q*(A**2-1.D0)/W                                     CORA-232
      C(4)=B4-3.D0*(2.D0*A+1.D0)*((6.D0*A-4.D0)*F+A*(A+1.D0)*X+4.D0*A*(ACORA-233
     1**2-1.D0)/W)/W                                                    CORA-234
      GO TO 45                                                          CORA-235
   20 GO TO ( 21 , 22 , 23 , 24 , 26 ) , LL                             CORA-236
C MULTIPOLE LQ=4 L1=L2.                                                 CORA-237
   21 B1=-P*Q*(2.D0*A+3.D0)/(A-1.D0)                                    CORA-238
      B2=-Q*(3.D0*E*S+A*X)                                              CORA-239
      B3=-P*(3.D0*F-A*X)                                                CORA-240
      B4=-((P+Q)*(A+1.5D0)-1.5D0*Y**2)/(A+1.D0)                         CORA-241
      BD=3.D0*DEN*A**2*(2.D0*A+3.D0)*(A+2.D0)*DSQRT(S*P*Q)              CORA-242
      IF (LT) GO TO 45                                                  CORA-243
      C(1)=B1                                                           CORA-244
      C(2)=(B3+P*A*(2.D0*A+3.D0)/W)                                     CORA-245
      C(3)=(B2+Q*A*(2.D0*A+3.D0)/W)                                     CORA-246
      C(4)=(B4-3.D0*A*(Y-A*(2.D0*A+3.D0)/W)/W)                          CORA-247
      GO TO 45                                                          CORA-248
C MULTIPOLE LQ=4 L1=L2-1.                                               CORA-249
   22 B1=Q*(3.D0*E*S-A*X)*(4.D0*A+6.D0)/(A-1.D0)                        CORA-250
      B2=18.D0*Q                                                        CORA-251
      B3=((2.D0*A+3.D0)*(3.D0*Q+(2.D0*A-3.D0)*X**2)-(6.D0*A-9.D0)*P)    CORA-252
      B4=(3.D0*F-A*X)*(4.D0*A-6.D0)/(A+1.D0)                            CORA-253
      BD=12.D0*DEN*A**2*DSQRT(S*Q)*(4.D0*A**2-9.D0)                     CORA-254
      IF (LT) GO TO 45                                                  CORA-255
      C(1)=B1-6.D0*Q*A*(2.D0*A+3.D0)/W                                  CORA-256
      C(2)=B3+2.D0*A*(4.D0*A**2-9.D0)*(X+3.D0*A/W)/W                    CORA-257
      C(3)=B2                                                           CORA-258
      C(4)=B4-6.D0*A*(2.D0*A-3.D0)/W                                    CORA-259
      GO TO 45                                                          CORA-260
C MULTIPOLE LQ=4 L1=L2-2.                                               CORA-261
   23 B1=((2.D0*A-E**2)*S**2+(A+2.D0)*(2.D0*E*F*S-A*(2.D0*A+2.D0+F**2)/(CORA-262
     1A-1.D0))/(A+1.D0))                                                CORA-263
      B2=6.D0*Q*((A+1.D0)*X-3.D0*F)                                     CORA-264
      B3=(-X*(4.D0*A+6.D0)*B1*(A**2-1.D0)+3.D0*(P*((2.D0*A+5.D0)*X-3.D0*CORA-265
     1Y)+(1.D0-S**2)*A*(A+2.D0)*(2.D0*A+3.D0)*(A*X+E*S)))               CORA-266
      B4=((4.D0*A+6.D0)*(B1*(A-1.D0)+3.D0*A*(S**2-1.D0))+3.D0*Y*(3.D0*Y-CORA-267
     1(2.D0*A+5.D0)*X)/(A+1.D0))                                        CORA-268
      B1=Q*(A+1.D0)*(4.D0*A+6.D0)*B1                                    CORA-269
      BD=60.D0*DEN*A**2*(A+2.D0)*(2.D0*A+3.D0)*DSQRT(S*Q*((A+1.D0)**2+F*CORA-270
     1*2))                                                              CORA-271
      IF (LT) GO TO 45                                                  CORA-272
      C(1)=B1-6.D0*Q*(A+2.D0)*(2.D0*A+3.D0)*(E*S+A*(X+4.D0*(A+1.D0)/W))/CORA-273
     1W                                                                 CORA-274
      C(2)=B3+2.D0*(2.D0*A+3.D0)*(6.D0*E*F*S+A*(F*((2.D0*A+13.D0)*X+6.D0CORA-275
     1*F)+(A+1.D0)*(2.D0*A*(2.D0*A+4.D0-(2.D0*A+1.D0)*S**2)+(2.D0*A+1.D0CORA-276
     2)*X**2))+3.D0*A*(A+2.D0)*(2.D0*A+1.D0)*((A+1.D0)*X+5.D0*F+A*(4.D0*CORA-277
     2A+4.D0)/W)/W)/W                                                   CORA-278
      C(3)=B2+12.D0*Q*(A+1.D0)*(2.D0*A+3.D0)/W                          CORA-279
      C(4)=B4-(12.D0*A+6.D0)*((A+6.D0)*F+A*(A+1.D0)*(X+(4.D0*A+6.D0)/W))CORA-280
     1/W                                                                CORA-281
      GO TO 45                                                          CORA-282
C MULTIPOLE LQ=4 L1=L2-3.                                               CORA-283
   24 B1=F*Q*(12.D0*(Q-(2.D0*A-1.D0)*(2.D0*A-5.D0))/((2.D0*A-3.D0)*(A-1.CORA-284
     1D0))+(26.D0*A-10.D0)*U)                                           CORA-285
      B2=Q*(36.D0*(Q+4.D0*A**2-5.D0)/(2.D0*A-3.D0)-30.D0*(A**2-1.D0)*U)/CORA-286
     1(2.D0*A+3.D0)                                                     CORA-287
      B3=B2+U*(A**2*(6.D0*A+3.D0)*((4.D0*A+10.D0)+5.D0*(A**2-1.D0)*U)-18CORA-288
     1.D0*Q*A*(3.D0*A+4.D0))/(2.D0*A+3.D0)                              CORA-289
      B4=F*(12.D0*(Q-(2.D0*A+5.D0)*(2.D0*A+1.D0))+U*(4.D0*A+2.D0)*((13.DCORA-290
     10*A+22.D0)*A-15.D0))/((A+1.D0)*(2.D0*A+3.D0))                     CORA-291
      IF (X.EQ.0.D0) GO TO 25                                           CORA-292
      B1=B1-2.D0*Q*X*((6.D0*Q*(A-2.D0)+(((28.D0*A-30.D0)*A-46.D0)*A+30.DCORA-293
     10))/((2.D0*A-3.D0)*(A-1.D0))-(A+1.D0)*((9.D0*A-5.D0)*U-2.D0*X**2))CORA-294
      B2=B2-12.D0*Q*X*(F-X*(A+1.D0))/(2.D0*A+3.D0)                      CORA-295
      B3=B3+X*(24.D0*Q*(A*X-F)+(2.D0*A+1.D0)*(12.D0*(2.D0*A+5.D0)*F-((44CORA-296
     1.D0*A+80.D0)*A-30.D0)*F*U+X*A*(28.D0*A+52.D0-U*(24.D0*A+2.D0)*(A+1CORA-297
     2.D0))+12.D0*F*X**2+4.D0*A*X**3*(A+1.D0)))/(2.D0*A+3.D0)           CORA-298
      B4=B4-X*(12.D0*(A-1.D0)*Q+(4.D0*A+2.D0)*(6.D0*F*X+A*(14.D0*A+26.D0CORA-299
     1-(A+1.D0)*((9.D0*A+1.D0)*U-2.D0*X**2))))/((2.D0*A+3.D0)*(A+1.D0)) CORA-300
   25 BD=360.D0*DEN*A**2*S*DSQRT(S*Q*((A+1.D0)**2+F**2)*((A-1.D0)**2+E**CORA-301
     12))                                                               CORA-302
      IF (LT) GO TO 45                                                  CORA-303
      C(1)=B1+Q*(12.D0*(A*(2.D0*A-1.D0)*(2.D0*A-5.D0)-Q*(9.D0*A-12.D0))/CORA-304
     1(2.D0*A-3.D0)+30.D0*A*U*(A**2-1.D0)-24.D0*(2.D0*A-1.D0)*(F*(7.D0*ACORA-305
     2+5.D0)+5.D0*A*(A**2-1.D0)/W)/W)/W                                 CORA-306
      C(2)=B3+(F*(48.D0*Q-8.D0*A*U*(2.D0*A+1.D0)*(7.D0*A-5.D0))+6.D0*(10CORA-307
     1.D0*Q*(11.D0*A**2-2.D0)-A**2*(2.D0*A+1.D0)*(30.D0*(2.D0*A-1.D0)+5.CORA-308
     2D0*U*(A**2-1.D0)-(8.D0*A-4.D0)*(12.D0*F+5.D0*(A**2-1.D0)/W)/W))/W)CORA-309
     3/W                                                                CORA-310
      C(3)=B2+24.D0*Q*(2.D0*F+5.D0*(A**2-1.D0)/W)/W                     CORA-311
      C(4)=B4+((12.D0*A*(2.D0*A+1.D0)*(2.D0*A+5.D0)-36.D0*Q*(3.D0*A+4.D0CORA-312
     1)+30.D0*A*U*(2.D0*A+1.D0)*(A**2-1.D0))/(2.D0*A+3.D0)-24.D0*(2.D0*ACORA-313
     2+1.D0)*(F*(7.D0*A-5.D0)+5.D0*A*(A**2-1.D0)/W)/W)/W                CORA-314
      IF (X.EQ.0.D0) GO TO 45                                           CORA-315
      C(1)=C(1)-12.D0*X*Q*((A+2.D0)*(3.D0*F+(A+1.D0)*X)+(4.D0*A-2.D0)*(ACORA-316
     1+1.D0)*(A+5.D0)/W)/W                                              CORA-317
      C(2)=C(2)+X*((12.D0*Q*(5.D0*A+6.D0)+(2.D0*A+1.D0)*(12.D0*F*X*(A+2.CORA-318
     1D0)-2.D0*A*(4.D0*A+10.D0+(A+1.D0)*(U*(9.D0*A-5.D0)-2.D0*X**2))))+1CORA-319
     22.D0*(2.D0*A+1.D0)*(F*((7.D0*A+24.D0)*A-10.D0)+A*(A+1.D0)*((A+2.D0CORA-320
     3)*X+(4.D0*A-2.D0)*(A+5.D0)/W))/W)/W                               CORA-321
      C(3)=C(3)+24.D0*X*Q*(A+1.D0)/W                                    CORA-322
      C(4)=C(4)-12.D0*X*(2.D0*A+1.D0)*((3.D0*F*(A+2.D0)+A*X*(A+1.D0))/(2CORA-323
     1.D0*A+3.D0)+2.D0*A*(A+1.D0)/W)/W                                  CORA-324
      GO TO 45                                                          CORA-325
C MULTIPOLE LQ=4 L1=L2-4.                                               CORA-326
   26 B1=6.D0*Q*(Q*((42.D0*A+3.D0)*U-36.D0)+A*(24.D0*A-12.D0-U*((34.D0*ACORA-327
     1+3.D0)*A+8.D0)-(8.D0*A+12.D0)*(A**2-1.D0)*U**2))                  CORA-328
      B2=18.D0*Q*(12.D0-U*(9.D0*A+1.D0))*F                              CORA-329
      B3=9.D0*F*(Q*(24.D0-U*(38.D0*A+2.D0))+U*A*(2.D0*A+1.D0)*(A+1.D0)*(CORA-330
     14.D0+(17.D0*A-12.D0)*U))                                          CORA-331
      B4=6.D0*(Q*((42.D0*A+3.D0)*U-36.D0)+A*(2.D0*A+1.D0)*(12.D0-(17.D0*CORA-332
     1A-2.D0)*U-U**2*(4.D0*A+6.D0)*(A-1.D0)))                           CORA-333
      IF (X.EQ.0.D0) GO TO 27                                           CORA-334
      B1=B1-Q*X*(X*(24.D0*Q+(56.D0*A+188.D0)*A+108.D0)+108.D0*(2.D0*A+3.CORA-335
     1D0)*F-2.D0*(2.D0*A+3.D0)*(58.D0*A+7.D0)*F*U-2.D0*(A+1.D0)*(2.D0*A+CORA-336
     23.D0)*(15.D0*A+4.D0)*U*X+12.D0*(2.D0*A+3.D0)*F*X**2+4.D0*(A+1.D0)*CORA-337
     3(2.D0*A+3.D0)*X**3)                                               CORA-338
      B2=B2+Q*X*(24.D0*F*X+(A+1.D0)*(108.D0-(66.D0*A+24.D0)*U+12.D0*X**2CORA-339
     1))                                                                CORA-340
      B3=B3+X*(24.D0*F*Q*X+(A+1.D0)*(Q*(324.D0+36.D0*X**2-(498.D0*A+42.DCORA-341
     10)*U)+(2.D0*A+1.D0)*(F*X*(108.D0-(152.D0*A+24.D0)*U+12.D0*X**2)-A*CORA-342
     2(72.D0-(168.D0*A-12.D0)*U-((57.D0*A+24.D0)*A-36.D0)*U**2-X**2*(28.CORA-343
     3D0-(36.D0*A+20.D0)*U+4.D0*X**2)))))                               CORA-344
      B4=B4+X*(X*(A*(2.D0*A+1.D0)*(30.D0*A+20.D0)*U-24.D0*Q)-(8.D0*A+4.DCORA-345
     10)*(F*(27.D0-(29.D0*A+6.D0)*U)+X*(7.D0*A+3.D0*F*X+A*X**2)))       CORA-346
   27 BD=2520.D0*DEN*A**2*S*DSQRT(S*Q*((A+1.D0)**2+F**2)*((A-1.D0)**2+E*CORA-347
     1*2)*((A+2.D0)**2+F**2))                                           CORA-348
      IF (LT) GO TO 45                                                  CORA-349
      C(1)=B1+18.D0*Q*(F*((A-1.D0)*(8.D0*A-4.D0)-20.D0*Q+(2.D0*A+3.D0)*(CORA-350
     1(17.D0*A-3.D0)*A-4.D0)*U)-4.D0*(15.D0*Q*((4.D0*A+4.D0)*A-1.D0)-A*(CORA-351
     22.D0*A-1.D0)*((34.D0*A+55.D0)*A+16.D0)-(2.D0*A+3.D0)*(2.D0*A*U*(A+CORA-352
     31.D0)*(A-1.D0)*(A+2.D0)-5.D0*(2.D0*A-1.D0)*(F*((5.D0*A+5.D0)*A+4.DCORA-353
     40)+2.D0*A*(A**2-1.D0)*(A+2.D0)/W)/W))/W)/W                        CORA-354
      C(2)=B3+(360.D0*Q**2-72.D0*Q*(4.D0*A**2+1.D0)-24.D0*(A+1.D0)*U*(Q*CORA-355
     1((66.D0*A-12.D0)*A-9.D0)-A**2*(2.D0*A+1.D0)*(34.D0*A-19.D0+U*(A-1.CORA-356
     2D0)*(2.D0*A+3.D0)))+18.D0*(20.D0*Q*F*((20.D0*A+14.D0)*A-3.D0)-(A+1CORA-357
     3.D0)*(A*F*((4.D0*A**2-1.D0)*140.D0+U*(10.D0*A+15.D0)*(2.D0*A+1.D0)CORA-358
     4*(5.D0*A-4.D0))-4.D0*(10.D0*Q*(((30.D0*A+19.D0)*A-1.D0)*A-6.D0)-A*CORA-359
     5*2*(2.D0*A+1.D0)*((2.D0*A-1.D0)*(84.D0*A+56.D0)+(2.D0*A+3.D0)*((2.CORA-360
     6D0*A-2.D0)*(A+2.D0)*U-(10.D0*A-5.D0)*(7.D0*F+2.D0*(A-1.D0)*(A+2.D0CORA-361
     7)/W)/W)))/W))/W)/W                                                CORA-362
      C(3)=B2+72.D0*Q*(5.D0*Q-4.D0*A**2-1.D0-(2.D0*A+3.D0)*((A**2-1.D0)*CORA-363
     1U-5.D0*(F*(3.D0*A-1.D0)+2.D0*(A**2-1.D0)*(A+2.D0)/W)/W))/W        CORA-364
      C(4)=B4-18.D0*(20.D0*Q*F-(A+1.D0)*((2.D0*A+1.D0)*F*(4.D0+U*(17.D0*CORA-365
     1A-12.D0))-4.D0*(15.D0*Q*(4.D0*A-1.D0)-(2.D0*A+1.D0)*(A*(34.D0*A-19CORA-366
     2.D0)+(2.D0*A+3.D0)*(A*(A-1.D0)*U-5.D0*(F*(5.D0*A-4.D0)+2.D0*A*(A-1CORA-367
     3.D0)*(A+2.D0)/W)/W)))/W))/W                                       CORA-368
      IF (X.EQ.0.D0) GO TO 45                                           CORA-369
      C(1)=C(1)-6.D0*Q*X*((A+2.D0)*(60.D0*Q-(44.D0*A-10.D0)*A-6.D0-(2.D0CORA-370
     1*A+3.D0)*((A+1.D0)*((11.D0*A-6.D0)*U-2.D0*X**2)-12.D0*F*X))+(8.D0*CORA-371
     2A+12.D0)*(((20.D0*A+90.D0)*A-5.D0)*F+(A+1.D0)*(A+2.D0)*((2.D0*A+5.CORA-372
     3D0)*X+(10.D0*A-5.D0)*(A+6.D0)/W))/W)/W                            CORA-373
      C(2)=C(2)+X*((2.D0*A+3.D0)*(240.D0*F*Q+(A+1.D0)*(144.D0*Q*X-(2.D0*CORA-374
     1A+1.D0)*((24.D0-24.D0*X**2+(182.D0*A-72.D0)*U)*F+A*X*(44.D0+(30.D0CORA-375
     2*A+8.D0)*U-4.D0*X**2))))+(A+1.D0)*(Q*((2400.D0*A+8400.D0)*A-360.D0CORA-376
     3)-(24.D0*A+12.D0)*(A*((102.D0*A+357.D0)*A-114.D0)-(2.D0*A+3.D0)*((CORA-377
     4A+2.D0)*(10.D0*F*X+A*(X**2-(5.5D0*A-3.D0)*U))+(F*((60.D0*A+270.D0)CORA-378
     5*A-120.D0)+A*(A+2.D0)*(X*(4.D0*A+10.D0)+(20.D0*A-10.D0)*(A+6.D0)/WCORA-379
     6))/W)))/W)/W                                                      CORA-380
      C(3)=C(3)+24.D0*X*Q*(2.D0*A+3.D0)*(5.D0*F+(A+1.D0)*(X+5.D0*(A+2.D0CORA-381
     1)/W))/W                                                           CORA-382
      C(4)=C(4)-6.D0*X*(A+1.D0)*(60.D0*Q-(2.D0*A+1.D0)*(22.D0*A+(11.D0*ACORA-383
     1+4.D0)*A*U-12.D0*F*X-2.D0*A*X**2-(8.D0*A+12.D0)*(10.D0*F+A*(X+5.D0CORA-384
     2*(A+2.D0)/W))/W))/W                                               CORA-385
      GO TO 45                                                          CORA-386
   28 GO TO ( 29 , 31 , 33 , 35 , 37 , 39 ) , LL                        CORA-387
C MULTIPOLE LQ=5 L1=L2.                                                 CORA-388
   29 BD=12.D0*DEN*A**2*(A+2.D0)*(2.D0*A+5.D0)*(4.D0*A**2-9.D0)*S*DSQRT(CORA-389
     1P*Q)                                                              CORA-390
      B1=28.D0*P*Q*F*(2.D0*A+3.D0)/(A-1.D0)                             CORA-391
      B2=Q*(84.D0*Q-A*(48.D0*A-72.D0-U*(A+2.D0)*(6.D0*A+27.D0)))        CORA-392
      B3=P*(84.D0*Q-A*(6.D0*A-9.D0)*(8.D0+(A+2.D0)*U))                  CORA-393
      B4=(4.D0*A-6.D0)*(14.D0*Q+A*(24.D0+(7.D0*A+12.D0)*U))*F/(A+1.D0)  CORA-394
      IF (X.EQ.0.D0) GO TO 30                                           CORA-395
      B1=B1+(28.D0*A+42.D0)*X*P*Q/(A-1.D0)                              CORA-396
      B2=B2+Q*X*((28.D0*A+126.D0)*F+((4.D0*A+22.D0)*A+42.D0)*X)         CORA-397
      B3=B3-P*X*(2.D0*A-3.D0)*(14.D0*F-2.D0*A*X)                        CORA-398
      B4=B4+X*(4.D0*A-6.D0)*(21.D0*Q-A*(14.D0*A-12.D0+(A**2-3.D0)*U)-X*(CORA-399
     1(2.D0*A-7.D0)*F+A*X))/(A+1.D0)                                    CORA-400
   30 IF (LT) GO TO 45                                                  CORA-401
      C(1)=B1-P*Q*(A+2.D0)*(24.D0*A+36.D0)/W                            CORA-402
      C(2)=B3-P*(4.D0*A**2-9.D0)*(8.D0*F-A*(6.D0*A+12.D0)/W)/W          CORA-403
      C(3)=B2-Q*(4.D0*A**2-9.D0)*(8.D0*F-A*(6.D0*A+12.D0)/W)/W          CORA-404
      C(4)=B4-(12.D0*A-18.D0)*((2.D0*A-4.D0)*Q+A*(A*(8.D0+(A+2.D0)*U)+(8CORA-405
     1.D0*A+12.D0)*(F-(A+2.D0)*A/W)/W))/W                               CORA-406
      IF (X.EQ.0.D0) GO TO 45                                           CORA-407
      C(3)=C(3)-X*Q*(8.D0*A**2-18.D0)*(A+4.D0)/W                        CORA-408
      C(2)=C(2)+X*P*A*(8.D0*A**2-18.D0)/W                               CORA-409
      C(4)=C(4)-X*(12.D0*A-18.D0)*((2.D0*A-4.D0)*F+A*(X+(4.D0*A+6.D0)/W)CORA-410
     1)/W                                                               CORA-411
      GO TO 45                                                          CORA-412
C MULTIPOLE LQ=5 L1=L2-1.                                               CORA-413
   31 BD=-60.D0*A**2*(2.D0*A+5.D0)*(2.D0*A+3.D0)*(A+3.D0)*(A+2.D0)*DEN*SCORA-414
     1*DSQRT(P*Q*((A+1.D0)**2+F**2))                                    CORA-415
      B1=P*Q*(4.D0*A+6.D0)*((28.D0*Q-(12.D0*A-56.D0)*A+40.D0)/(A-1.D0)-(CORA-416
     12.D0*A+5.D0)*(A+1.D0)*U)                                          CORA-417
      B2=Q*F*(168.D0*Q+480.D0*A+240.D0+((48.D0*A+120.D0)*A+30.D0)*U)    CORA-418
      B3=P*F*(168.D0*Q+480.D0*A+240.D0-(12.D0*A+6.D0)*((2.D0*A+2.D0)*A-5CORA-419
     1.D0)*U)                                                           CORA-420
      B4=((2.D0*A-3.D0)*(56.D0*Q-(24.D0*A-208.D0)*A+80.D0-(((4.D0*A-18.DCORA-421
     10)*A-52.D0)*A-10.D0)*U)*Q+A**2*(2.D0*A+1.D0)*(480.D0-((24.D0*A-56.CORA-422
     2D0)*A-180.D0)*U-(2.D0*A+5.D0)*(2.D0*A+3.D0)*(A+1.D0)*U**2))/(A+1.DCORA-423
     30)                                                                CORA-424
      IF (X.EQ.0.D0) GO TO 32                                           CORA-425
      B1=B1-X*P*Q*(4.D0*A+6.D0)*(7.D0*F-X*(A+1.D0))                     CORA-426
      B2=B2+X*Q*((14.D0*A+210.D0)*Q-(56.D0*A-380.D0)*A+240.D0-(6.D0*A+6.CORA-427
     1D0)*(A**2-5.D0)*U-(12.D0*A-36.D0)*F*X-(6.D0*A+6.D0)*X**2)         CORA-428
      B3=B3-X*P*((98.D0*A-42.D0)*Q-(2.D0*A+1.D0)*(A*(24.D0*A-116.D0+(A+1CORA-429
     1.D0)*((7.D0*A+16.D0)*U-2.D0*X**2))+(14.D0*A-6.D0)*F*X))           CORA-430
      B4=B4-X*((2.D0*A-3.D0)*(((2.D0*A+1.D0)*(24.D0*A-80)+(14.D0*A-70.D0CORA-431
     1)*Q)*F-((2.D0*A-22.D0)*A+12.D0)*Q*X)+(4.D0*A+2.D0)*((((11.D0*A+17.CORA-432
     2D0)*A+1.D0)*A+15.D0)*F*U-X*((A-1.D0)*(2.D0*A-3.D0)*F*X+A*((14.D0*ACORA-433
     3-18.D0)*A+58.D0+(A+1.D0)*((A+2.D0)*(A-4.D0)*U+X**2)))))/(A+1.D0)  CORA-434
   32 IF (LT) GO TO 45                                                  CORA-435
      C(1)=B1-12.D0*P*Q*(2.D0*A+3.D0)*(A+2.D0)*(2.D0*F-(2.D0*A+5.D0)*(A+CORA-436
     11.D0)/W)/W                                                        CORA-437
      C(2)=B3-P*(4.D0*A+6.D0)*((16.D0*A-24.D0)*Q+(2.D0*A+1.D0)*(40.D0*A-CORA-438
     1A*(2.D0*A+5.D0)*(A+1.D0)*(U-(6.D0*A+12.D0)/W**2)+(30.D0*A+60.D0)*FCORA-439
     2/W))/W                                                            CORA-440
      C(3)=B2-Q*(4.D0*A+6.D0)*((16.D0*A-24.D0)*Q+(2.D0*A+1.D0)*(40.D0*A+CORA-441
     1(30.D0*A+60.D0)*F/W)+A*(A+1.D0)*(6.D0*A+15.D0)*(U-(10.D0*A+20.D0)/CORA-442
     2W**2))/W                                                          CORA-443
      C(4)=B4-6.D0*(F*(((8.D0*A-28.D0)*A+24.D0)*Q+A*(4.D0*A+2.D0)*(((2.DCORA-444
     10*A+2.D0)*A-5.D0)*U-40.D0))-(2.D0*A+3.D0)*((((4.D0*A+6.D0)*A+62.D0CORA-445
     2)*A+20.D0)*Q+A*(2.D0*A+1.D0)*((2.D0*A+5.D0)*(A+1.D0)*U*A-40.D0*A-(CORA-446
     310.D0*A+20.D0)*(4.D0*F+A*(A+1.D0)*(2.D0*A+5.D0)/W)/W))/W)/W       CORA-447
      IF (X.EQ.0.D0) GO TO 45                                           CORA-448
      C(1)=C(1)+X*P*Q*(12.D0*A+18.D0)*(A+1.D0)*(A+2.D0)/W               CORA-449
      C(2)=C(2)+X*P*(4.D0*A+6.D0)*(((8.D0*A-8.D0)*A-6.D0)*F-A*(A+1.D0)*(CORA-450
     12.D0*A+1.D0)*(X+(3.D0*A+6.D0)/W))/W                               CORA-451
      C(3)=C(3)-X*Q*(24.D0*A+36.D0)*((2.D0*A-3.D0)*F+(A+1.D0)*(X+(5.D0*ACORA-452
     1+10.D0)/W))/W                                                     CORA-453
      C(4)=C(4)+6.D0*X*((A-2.D0)*(A-3.D0)*(2.D0*A-3.D0)*Q+(2.D0*A+1.D0)*CORA-454
     1(((8.D0*A-12.D0)*A+40.D0)*A+((2.D0*A-2.D0)*A+6.D0)*F*X+A*(A+1.D0)*CORA-455
     2((A**2-5.D0)*U+X**2)+(4.D0*A+6.D0)*(((2.D0*A+2.D0)*A+10.D0)*F+A*(ACORA-456
     3+1.D0)*(X+(5.D0*A+10.D0)/W))/W))/W                                CORA-457
      GO TO 45                                                          CORA-458
C MULTIPOLE LQ=5 L1=L2-2.                                               CORA-459
   33 DN=720.D0*A**2*(4.D0*A**2-9.D0)*(2.D0*A+5.D0)*(A+2.D0)            CORA-460
      BD=DEN*DN*S*DSQRT(Q*((A+1.D0)**2+F**2))/2.D0                      CORA-461
      B1=2.D0*(2.D0*A+3.D0)*(24.D0*((2.D0*A+5.D0)*(A+2.D0)+7.D0*Q)/(A-1.CORA-462
     1D0)-(2.D0*A-3.D0)*(11.D0*A-10.D0)*U)*F                            CORA-463
      B2=-6.D0*(24.D0*(2.D0*A**2.D0-15.D0*A-10.D0)-168.D0*Q+5.D0*(A+1.D0CORA-464
     1)*(2.D0*A**2+A-6.D0)*U)                                           CORA-465
      B3=((1008.D0*Q-144.D0*(2.D0*A**2-15.D0*A-10.D0))*Q-(2.D0*A-3.D0)*UCORA-466
     1*(6.D0*(17.D0*A**2-9.D0*A+10.D0)*Q+A**2*(2.D0*A+1.D0)*(24.D0*(A+10CORA-467
     2.D0)-15.D0*(A+1.D0)*(A+2.D0)*U)))                                 CORA-468
      B4=2.D0*(2.D0*A-3.D0)*(168.D0*Q+(2.D0*A+1.D0)*(24.D0*(A+10.D0)-(11CORA-469
     1.D0*A**2-A-30.D0)*U))*F/(A+1.D0)                                  CORA-470
      IF (X.EQ.0.D0) GO TO 34                                           CORA-471
      B1=B1+X*2.D0*(2.D0*A+3.D0)*((20.D0*A**3-198.D0*A**2+316.D0*A+240.DCORA-472
     10-42.D0*(2.D0*A-5.D0)*Q)/(A-1.D0)+(2.D0*A-3.D0)*(12.D0*F*X+(A+1.D0CORA-473
     2)*((9.D0*A+10.D0)*U-2.D0*X**2)))                                  CORA-474
      B2=B2-X*12.D0*(7.D0*(2.D0*A-3.D0)*F-(A+1.D0)*(2.D0*A-3.D0)*X)     CORA-475
      B3=B3-X*(2.D0*A-3.D0)*((420.D0*F-36.D0*(5.D0*A-2.D0)*X)*Q+(2.D0*A+CORA-476
     11.D0)*(48.D0*(A+10.D0)*F-(40.D0*A**2-38.D0*A-60.D0)*F*U+4.D0*(5.D0CORA-477
     2*A-58.D0)*A*X+12.D0*(2.D0*A-1.D0)*F*X**2+A*(A+1.D0)*X*(8.D0*(3.D0*CORA-478
     3A+4.D0)*U-4.D0*X**2)))                                            CORA-479
      B4=B4-X*2.D0*(2.D0*A-3.D0)*(42.D0*(2.D0*A-1.D0)*Q-(2.D0*A+1.D0)*(6CORA-480
     1.D0*(2.D0*A-1.D0)*F*X+A*((10.D0*A-116.D0)+(A+1.D0)*((9.D0*A+16.D0)CORA-481
     2*U-2.D0*X**2))))/(A+1.D0)                                         CORA-482
   34 IF (LT) GO TO 45                                                  CORA-483
      B1=B1*Q                                                           CORA-484
      B2=B2*Q                                                           CORA-485
      C(1)=B1-6.D0*Q*(A+2.D0)*(2.D0*A+3.D0)*(8.D0*A*(2.D0*A+5.D0)+24.D0*CORA-486
     1Q-5.D0*A*(A+1.D0)*(2.D0*A-3.D0)*U+4.D0*(2.D0*A-3.D0)*(2.D0*A+5.D0)CORA-487
     2*(F+5.D0*A*(A+1.D0)/W)/W)/W                                       CORA-488
      C(2)=B3-(4.D0*A**2-9.D0)*(8.D0*(12.D0*Q+A*(A+10.D0)*(2.D0*A+1.D0)*CORA-489
     1U)*F-(6.D0*A+12.D0)*((A+1.D0)*(20.D0*Q-A**2*(10.D0*A+5.D0)*U)+(8.DCORA-490
     20*A+4.D0)*A*(2.D0*A+5.D0)*(6.D0*F+5.D0*A*(A+1.D0)/W)/W)/W)/W      CORA-491
      C(3)=B2-Q*24.D0*(4.D0*A**2-9.D0)*(4.D0*F-5.D0*(A+1.D0)*(A+2.D0)/W)CORA-492
     1/W                                                                CORA-493
      C(4)=B4-6.D0*(2.D0*A-3.D0)*(24.D0*(A-2.D0)*Q+A*(2.D0*A+1.D0)*(8.D0CORA-494
     1*(A+10.D0)-5.D0*(A+1.D0)*(A+2.D0)*U)+4.D0*(2.D0*A+3.D0)*(2.D0*A+1.CORA-495
     2D0)*((A+10.D0)*F+5.D0*A*(A+1.D0)*(A+2.D0)/W)/W)/W                 CORA-496
      IF (X.EQ.0.D0) GO TO 45                                           CORA-497
      C(1)=C(1)+X*Q*(12.D0*A+24.D0)*(4.D0*A**2-9.D0)*(3.D0*F-(A+1.D0)*(XCORA-498
     1+(4.D0*A+10.D0)/W))/W                                             CORA-499
      C(2)=C(2)+X*(8.D0*A**2-18.D0)*(36.D0*(A-1.D0)*Q-(2.D0*A+1.D0)*(6.DCORA-500
     10*(A-2.D0)*F*X-A*(8.D0*(A+10.D0)-(A+1.D0)*((9.D0*A+10.D0)*U-2.D0*XCORA-501
     2**2))-6.D0*(A+2.D0)*((A+10.D0)*F+A*(A+1.D0)*(X+(4.D0*A+10.D0)/W))/CORA-502
     3W))/W                                                             CORA-503
      C(3)=C(3)+X*Q*24.D0*(4.D0*A**2-9.D0)*(A+1.D0)/W                   CORA-504
      C(4)=C(4)+X*12.D0*(2.D0*A+1.D0)*(2.D0*A-3.D0)*(3.D0*(A-2.D0)*F-A*(CORA-505
     1A+1.D0)*(X+(4.D0*A+6.D0)/W))/W                                    CORA-506
      GO TO 45                                                          CORA-507
C MULTIPOLE LQ=5 L1=L2-3.                                               CORA-508
   35 BD=-2520.D0*DEN*A**2*(A+2.D0)*(A+3.D0)*(2.D0*A+3.D0)*(2.D0*A+5.D0)CORA-509
     1*S*DSQRT(Q*((A+1.D0)**2+F**2)*((A+2.D0)**2+F**2))                 CORA-510
      B1=Q*(4.D0*A+6.D0)*(((24.D0*((A+36.D0)*A+40.D0)+168.D0*Q)*Q+48.D0*CORA-511
     1A*(A+2.D0)*(A+3.D0)*(2.D0*A+5.D0))/(A-1.D0)-(A*(((80.D0*A+620.D0)*CORA-512
     2A+966.D0)*A+540.D0)+((28.D0*A-26.D0)*A+60.D0)*Q-6.D0*A*(A+1.D0)*(ACORA-513
     3+2.D0)*(2.D0*A+3.D0)*(2.D0*A+5.D0)*U)*U)                          CORA-514
      B2=6.D0*F*Q*(48.D0*((2.D0*A+21.D0)*A+20.D0)+168.D0*Q-(((16.D0*A+36CORA-515
     1.D0)*A-22.D0)*A-60.D0)*U)                                         CORA-516
      B3=6.D0*F*((48.D0*((2.D0*A+21.D0)*A+20.D0)+168.D0*Q-(10.D0*A-12.D0CORA-517
     1)*(4.D0*A**2+5.D0)*U)*Q+A*(A+1.D0)*(2.D0*A+1.D0)*(((16.D0*A+168.D0CORA-518
     2)*A+720.D0)-(2.D0*A+3.D0)*(A+2.D0)*(8.D0*A+45.D0)*U)*U)           CORA-519
      B4=((2.D0*A-3.D0)*(48.D0*((A+48.D0)*A+40.D0)+336.D0*Q)*Q/(A+1.D0)-CORA-520
     1(4.D0*A**2-9.D0)*(28.D0*A-40.D0)*Q*U+A*(2.D0*A+1.D0)*(96.D0*((2.D0CORA-521
     2*A+21.D0)*A+90.D0)-(4.D0*A+6.D0)*(((40.D0*A+282.D0)*A+180.D0)-(A+2CORA-522
     3.D0)*(2.D0*A+3.D0)*(6.D0*A+15.D0)*U)*U))                          CORA-523
      IF (X.EQ.0.D0) GO TO 36                                           CORA-524
      B1=B1+X*Q*(4.D0*A+6.D0)*((((28.D0*A+402.D0)*A+704.D0)*A+360.D0)*X-CORA-525
     1((88.D0*A+700.D0)*A+600.D0)*F-126.D0*F*Q+(54.D0*A+24.D0)*Q*X+(2.D0CORA-526
     2*A+3.D0)*(((5.D0*A-53.D0)*A-70.D0)*F*U-(8.D0*A+2.D0)*F*X**2-(A+1.DCORA-527
     30)*(A+2.D0)*((15.D0*A+25.D0)*U-2.D0*X**2)*X))                     CORA-528
      B2=B2-X*6.D0*Q*((((8.D0*A+372.D0)*A+696.D0)*A+360.D0)+(70.D0*A+42.CORA-529
     1D0)*Q-(2.D0*A+3.D0)*((10.D0*A+6.D0)*F*X+(A+1.D0)*(A+2.D0)*((11.D0*CORA-530
     2A+25.D0)*U-2.D0*X**2)))                                           CORA-531
      B3=B3+X*((48.D0*((13.D0*A-7.D0)*A-3.D0)*F*X-(((144.D0*A+6696.D0)*ACORA-532
     1+1104.D0)*A-3600.D0)-84.D0*(13.D0*A-9.D0)*Q)*Q+(A+1.D0)*(((((352.DCORA-533
     20*A+224.D0)*A+1518.D0)*A+1260.D0)*U-12.D0*((20.D0*A-14.D0)*A-3.D0)CORA-534
     3*X**2)*Q+(2.D0*A+1.D0)*((2.D0*A-3.D0)*((88.D0*A+720.D0)+(16.D0*A+2CORA-535
     44.D0)*X**2)*F*X+A*((((280.D0*A+2292.D0)*A+2772.D0)*A+1080.D0)*U-(1CORA-536
     592.D0*A+2016.D0)*A-8640.D0-28.D0*((2.D0*A+31.D0)*A+24.D0)*X**2)-(2CORA-537
     6.D0*A+3.D0)*(((16.D0*A-230.D0)*A-300.D0)*F*U*X+A*(A+2.D0)*(((57.D0CORA-538
     7*A+171.D0)*A+90.D0)*U**2-(36.D0*A+62.D0)*U*X**2+4.D0*X**4)))))    CORA-539
      B4=B4-X*((4.D0*A-6.D0)*(126.D0*A+42.D0)*F*Q/(A+1.D0)-36.D0*(2.D0*ACORA-540
     1-3.D0)*(3.D0*A+1.D0)*Q*X-(4.D0*A+2.D0)*(A*((28.D0*A+434.D0)*A+336.CORA-541
     2D0)*X-(2.D0*A-3.D0)*(44.D0*A+360.D0)*F+(2.D0*A+3.D0)*(((5.D0*A-79.CORA-542
     3D0)*A-150.D0)*F*U-(8.D0*A-12.D0)*F*X**2-A*(A+2.D0)*((15.D0*A+31.D0CORA-543
     4)*U-2.D0*X**2)*X)))                                               CORA-544
   36 IF (LT) GO TO 45                                                  CORA-545
      C(1)=B1-12.D0*Q*(2.D0*A+3.D0)*(A+2.D0)*((12.D0*Q+(2.D0*A+3.D0)*((8CORA-546
     1.D0*A+41.D0)*A+15.D0)*U)*F-(4.D0*A+10.D0)*((A+3.D0)*4.D0*F+((33.D0CORA-547
     2+13.D0*A)*Q-A*(A+3.D0)*(4.D0*A-6.D0)-(2.D0*A+3.D0)*(3.D0*A*(A+1.D0CORA-548
     3)*(A+2.D0)*U-10.D0*(A+3.D0)*((4.D0*A+1.D0)*F+3.D0*A*(A+1.D0)*(A+2.CORA-549
     4D0)/W)/W))/W))/W                                                  CORA-550
      C(3)=B2-12.D0*Q*(2.D0*A+3.D0)*((((32.D0*A+288.D0)*A+456.D0)*A+240.CORA-551
     1D0)+8.D0*(2.D0*A-3.D0)*Q-(A+2.D0)*(2.D0*A+3.D0)*(3.D0*(A+1.D0)*(2.CORA-552
     2D0*A+5.D0)*U-10.D0*((2.D0*A+11.D0)*F+3.D0*(2.D0*A+5.D0)*(A+1.D0)*(CORA-553
     3A+2.D0)/W)/W))/W                                                  CORA-554
      C(2)=B3-(4.D0*A+6.D0)*((((192.D0*A+1728.D0)*A+2736.D0)*A+1440.D0+4CORA-555
     18.D0*(2.D0*A-3.D0)*Q)*Q-(A+1.D0)*((((176.D0*A+964.D0)*A+2166.D0)*ACORA-556
     2+540.D0)*Q-A**2*(2.D0*A+1.D0)*(((16.D0*A+168.D0)*A+720.D0)+6.D0*(ACORA-557
     3+2.D0)*(2.D0*A+3.D0)*(2.D0*A+5.D0)*U))*U-6.D0*(A+2.D0)*((2.D0*A+3.CORA-558
     4D0)*((A+1.D0)*(2.D0*A+1.D0)*(20.D0*A+75.D0)*A*U-(20.D0*A+110.D0)*QCORA-559
     5)*F+(A+1.D0)*(2.D0*A+5.D0)*(6.D0*(28.D0*A**2*(2.D0*A+1.D0)*(A+3.D0CORA-560
     6)-((52.D0*A+161.D0)*A+30.D0)*Q+A**2*(A+2.D0)*(2.D0*A+1.D0)*(2.D0*ACORA-561
     7+3.D0)*U)-10.D0*A*(A+3.D0)*(2.D0*A+1.D0)*(2.D0*A+3.D0)*(14.D0*F+6.CORA-562
     8D0*A*(A+2.D0)/W)/W)/W)/W)/W                                       CORA-563
      C(4)=B4+12.D0*(((2.D0*A+1.D0)*(A+1.D0)*(((16.D0*A+168.D0)*A+720.D0CORA-564
     1)-(A+2.D0)*(2.D0*A+3.D0)*(8.D0*A+45.D0)*U)-12.D0*(2.D0*A-3.D0)*(A-CORA-565
     22.D0)*Q)*F+(2.D0*A+3.D0)*(A+1.D0)*(((52.D0*A+266.D0)*A+660.D0)*Q-ACORA-566
     3*(2.D0*A+1.D0)*((8.D0*A+84.D0)*A+360.D0)-(A+2.D0)*(2.D0*A+1.D0)*(2CORA-567
     4.D0*A+3.D0)*(3.D0*A*(2.D0*A+5.D0)*U-10.D0*((8.D0*A+30.D0)*F+3.D0*ACORA-568
     5*(2.D0*A+5.D0)*(A+2.D0)/W)/W))/W)/W                               CORA-569
      IF (X.EQ.0.D0) GO TO 45                                           CORA-570
      C(1)=C(1)+6.D0*Q*X*(2.D0*A+3.D0)*(A+2.D0)*(((40.D0*A+300.D0)*A+476CORA-571
     1.D0)*A+240.D0+(18.D0*A-6.D0)*Q-(2.D0*A+3.D0)*((2.D0*A-10.D0)*F*X+(CORA-572
     2A+1.D0)*(A+2.D0)*((11.D0*A+15.D0)*U-2.D0*X**2)-(8.D0*A+20.D0)*((3.CORA-573
     3D0*A+13.D0)*F+(A+1.D0)*(A+2.D0)*(X+(5.D0*A+15.D0)/W))/W))/W       CORA-574
      C(2)=C(2)+X*(4.D0*A+6.D0)*(12.D0*(2.D0*A-3.D0)*(4.D0*A-1.D0)*F*Q-(CORA-575
     1A+1.D0)*(((60.D0*A-42.D0)*A+180.D0)*Q*X+(2.D0*A+1.D0)*(((32.D0*A+3CORA-576
     236.D0)*A+1440.D0)*F-(2.D0*A+3.D0)*((28.D0*A+202.D0)*A+180.D0)*F*U+CORA-577
     3A*(((40.D0*A+308.D0)*A+120.D0)-(A+2.D0)*(2.D0*A+3.D0)*((15.D0*A+25CORA-578
     4.D0)*U-2.D0*X**2))*X-(2.D0*A+3.D0)*(2.D0*A-24.D0)*F*X**2)+(3.D0*A+CORA-579
     56.D0)*(((92.D0*A+406.D0)*A+780.D0)*Q-(2.D0*A+1.D0)*(A*((8.D0*A+84.CORA-580
     6D0)*A+360.D0)+(2.D0*A+3.D0)*(A*(A+2.D0)*((11.D0*A+15.D0)*U-2.D0*X*CORA-581
     7*2)-(6.D0*A+40.D0)*F*X-(8.D0*A+20.D0)*((30.D0+8.D0*A)*F+A*(A+2.D0)CORA-582
     8*(X+(5.D0*A+15.D0)/W))/W)))/W))/W                                 CORA-583
      C(3)=C(3)+24.D0*X*Q*(2.D0*A+3.D0)**2*((2.D0*A-3.D0)*F-(A+1.D0)*(A+CORA-584
     12.D0)*(X+5.D0*(A+2.D0)/W))/W                                      CORA-585
      C(4)=C(4)+X*6.D0*(A+1.D0)*(18.D0*(A-2.D0)*(2.D0*A-3.D0)*Q+(2.D0*A+CORA-586
     11.D0)*(A*((40.D0*A+308.D0)*A+120.D0)-(2.D0*A+3.D0)*((2.D0*A-24.D0)CORA-587
     2*F*X+A*(A+2.D0)*((11.D0*A+25.D0)*U-2.D0*X**2)-4.D0*(2.D0*A+3.D0)*(CORA-588
     3(20.D0+3.D0*A)*F+A*(A+2.D0)*(X+5.D0*(A+2.D0)/W))/W)))/W           CORA-589
      GO TO 45                                                          CORA-590
C MULTIPOLE LQ=5 L1=L2-4.                                               CORA-591
   37 BD=10080.D0*DEN*A**2*(2.D0*A+5.D0)*(A+2.D0)*S**2*DSQRT(Q*((A+1.D0)CORA-592
     1**2+F**2)*((A+2.D0)**2+F**2)*((A-1.D0)**2+E**2))                  CORA-593
      B1=Q*(((336.D0*Q+48.D0*((8.D0*A+22.D0)*A+47.D0))*Q-96.D0*(2.D0*A-1CORA-594
     1.D0)*(A+2.D0)*(2.D0*A+5.D0)*(3.D0*A-7.D0))/((A-1.D0)*(2.D0*A-3.D0)CORA-595
     2)+(4.D0*(((162.D0*A+683.D0)*A+538.D0)*A-280.D0)-4.D0*(7.D0*A-38.D0CORA-596
     3)*Q-(2.D0*A+3.D0)*(((163.D0*A+497.D0)*A-16.D0)*A-140.D0)*U)*U)*F  CORA-597
      B2=3.D0*Q*(((96.D0*((((12.D0*A+56.D0)*A+29.D0)*A-70.D0)*A-70.D0)+(CORA-598
     148.D0*((4.D0*A+28.D0)*A+47.D0)+336.D0*Q)*Q)/(2.D0*A-3.D0)-4.D0*(((CORA-599
     2((66.D0*A+307.D0)*A+172.D0)*A-280.D0)*A-280.D0)+(A-2.D0)*(4.D0*A-1CORA-600
     39.D0)*Q)*U)/(2.D0*A+3.D0)+35.D0*(A**2-1.D0)*(A+2.D0)**2*U**2)     CORA-601
      B3=(((1008.D0*Q+144.D0*((4.D0*A+28.D0)*A+47.D0))*Q+288.D0*((((12.DCORA-602
     10*A+56.D0)*A+29.D0)*A-70.D0)*A-70.D0))*Q/(2.D0*A-3.D0)-12.D0*(((((CORA-603
     2178.D0*A+967.D0)*A+1410.D0)*A+428.D0)*A-280.D0)+(A-2.D0)*(10.D0*A-CORA-604
     319.D0)*Q)*Q*U+(3.D0*A+3.D0)*(((((458.D0*A+2185.D0)*A+2899.D0)*A-35CORA-605
     42.D0)*A-420.D0)*Q*U+(4.D0*A+2.D0)*A**2*(8.D0*((6.D0*A+37.D0)*A+70.CORA-606
     5D0)-(((62.D0*A+349.D0)*A+524.D0)*A-140.D0)*U))*U)/(2.D0*A+3.D0)-52CORA-607
     6.5D0*A**2*(A**2-1.D0)*(A+2.D0)**2*(2.D0*A+1.D0)*U**3              CORA-608
      B4=(((336.D0*Q+48.D0*((8.D0*A+34.D0)*A+47.D0)-(4.D0*A-8.D0)*((14.DCORA-609
     10*A-41.D0)*A+57.D0)*U)*Q/(A+1.D0)+(2.D0*A+1.D0)*(4.D0*(((162.D0*A+CORA-610
     2863.D0)*A+1070.D0)*A-840.D0)*U-96.D0*((6.D0*A+37.D0)*A+70.D0)))/(2CORA-611
     3.D0*A+3.D0)-(2.D0*A+1.D0)*(((163.D0*A+631.D0)*A+256.D0)*A-420.D0)*CORA-612
     4U**2)*F                                                           CORA-613
      IF (X.EQ.0.D0) GO TO 38                                           CORA-614
      B1=B1+X*Q*(8.D0*((((((180.D0*A+604.D0)*A-425.D0)*A-1813.D0)*A-394.CORA-615
     1D0)*A+840.D0)+(((34.D0*A+147.D0)*A-193.D0)*A-450.D0)*Q+42.D0*(A-2.CORA-616
     2D0)*Q**2)/((1.D0-A)*(2.D0*A-3.D0))+96.D0*F*Q*X+16.D0*((((26.D0*A+1CORA-617
     326.D0)*A+133.D0)*A-13.D0)*A-70.D0)*U+4.D0*((41.D0*A+115.D0)*A+170.CORA-618
     4D0)*Q*U+4.D0*((2.D0*A-1.D0)*A-54.D0)*F*X-8.D0*(5.D0*A+2.D0)*Q*X**2CORA-619
     5-(((152.D0*A+892.D0)*A+1380.D0)*A+664.D0)*X**2+(2.D0*A+3.D0)*(4.D0CORA-620
     6*((20.D0*A+97.D0)*A+62.D0)*F*U*X+4.D0*(A-2.D0)*F*X**3+(A+1.D0)*(A+CORA-621
     72.D0)*(((44.D0*A+36.D0)*U-4.D0*X**2)*X**2-((87.D0*A+67.D0)*A-70.D0CORA-622
     8)*U**2)))                                                         CORA-623
      B2=B2-12.D0*X*Q*((28.D0*F*Q+2.D0*((18.D0*A+83.D0)*A+56.D0)*F-(((34CORA-624
     1.D0*A+217.D0)*A+341.D0)*A+166.D0)*X-4.D0*(3.D0*A+1.D0)*Q*X)/(2.D0*CORA-625
     2A+3.D0)+(((9.D0*A+58.D0)*A+44.D0)*U+2.D0*A*X**2)*F+(A+1.D0)*(A+2.DCORA-626
     30)*(9.D0*(A+1.D0)*U-X**2)*X)                                      CORA-627
      B3=B3+X*(((96.D0*(5.D0*A-3.D0)*Q+4.D0*(((170.D0*A+1101.D0)*A+901.DCORA-628
     10)*A+162.D0)-24.D0*((10.D0*A-5.D0)*A-2.D0)*F*X)*Q*X-4.D0*((((22.D0CORA-629
     2*A+543.D0)*A+629.D0)*A+510.D0)*U+168.D0*Q)*F*Q+(A+1.D0)*(4.D0*(((2CORA-630
     36.D0*A-5.D0)*A+18.D0)*X**2-(((166.D0*A+627.D0)*A+1415.D0)*A+558.D0CORA-631
     4)*U)*Q*X+(2.D0*A+1.D0)*((96.D0*((6.D0*A+37.D0)*A+70.D0)-4.D0*((2.DCORA-632
     50*A-105.D0)*A-498.D0)*X**2-8.D0*(((138.D0*A+751.D0)*A+1012.D0)*A-4CORA-633
     620.D0)*U)*F+A*(4.D0*((38.D0*A+197.D0)*A+154.D0)*X**2-(((572.D0*A+2CORA-634
     7362.D0)*A+960.D0)*A-40.D0)*U+8.D0*((90.D0*A+491.D0)*A+674.D0))*X))CORA-635
     8)/(2.D0*A+3.D0)-24.D0*(17.D0*A+50.D0)*F*Q+(A+1.D0)*(2.D0*A+1.D0)*(CORA-636
     9(((373.D0*A+1411.D0)*A+544.D0)*A-420.D0)*F*U**2-4.D0*((23.D0*A+125CORA-637
     A.D0)*A+54.D0)*F*U*X**2-4.D0*(A-6.D0)*F*X**4+A*(A+2.D0)*(((141.D0*ACORA-638
     B+199.D0)*A-22.D0)*U**2-(50.D0*A+48.D0)*U*X**2+4.D0*X**4)*X))      CORA-639
      B4=B4+X*(((48.D0*((4.D0*A-4.D0)*A-1.D0)*F*X-8.D0*(((34.D0*A+225.D0CORA-640
     1)*A-61.D0)*A-168.D0)-336.D0*(A-1.D0)*Q)*Q/(A+1.D0)+4.D0*((((82.D0*CORA-641
     2A+315.D0)*A+902.D0)*A+396.D0)*U-A*(20.D0*A-26.D0)*X**2)*Q+4.D0*(2.CORA-642
     3D0*A+1.D0)*(((2.D0*A-105.D0)*A-498.D0)*F*X+A*((((104.D0*A+430.D0)*CORA-643
     4A+171.D0)*A-10.D0)*U-2.D0*((90.D0*A+491.D0)*A+674.D0)-((38.D0*A+19CORA-644
     57.D0)*A+154.D0)*X**2)))/(2.D0*A+3.D0)+(2.D0*A+1.D0)*(4.D0*((20.D0*CORA-645
     6A+107.D0)*A+54.D0)*F*U*X+4.D0*(A-6.D0)*F*X**3+A*(A+2.D0)*(4.D0*(11CORA-646
     7.D0*A+12.D0)*U*X**2-((87.D0*A+145.D0)*A-22.D0)*U**2-4.D0*X**4)))  CORA-647
   38 IF (LT) GO TO 45                                                  CORA-648
      C(1)=B1-3.D0*Q*(A+2.D0)*(((16.D0*((56.D0*A+62.D0)*A-177.D0)+48.D0*CORA-649
     1Q)*Q-32.D0*A*(2.D0*A-1.D0)*(2.D0*A+5.D0)*(3.D0*A-7.D0))/(2.D0*A-3.CORA-650
     2D0)+(4.D0*A*(((62.D0*A+159.D0)*A+29.D0)*A+70.D0)-4.D0*((97.D0*A+22CORA-651
     39.D0)*A-6.D0)*Q+35.D0*A*(A**2-1.D0)*(2.D0*A+3.D0)*(A+2.D0)*U)*U-4.CORA-652
     4D0*(2.D0*A+5.D0)*((4.D0*(2.D0*A-1.D0)*(3.D0*A-7.D0)-92.D0*Q+(2.D0*CORA-653
     5A+3.D0)*((67.D0*A-17.D0)*A-14.D0)*U)*F-5.D0*(4.D0*((59.D0*A+59.D0)CORA-654
     6*A-18.D0)*Q-4.D0*A*(2.D0*A-1.D0)*(3.D0*A+1.D0)*(11.D0*A+14.D0)-(2*CORA-655
     7A+3.D0)*(7.D0*A*(A**2-1.D0)*(A+2.D0)*U-6.D0*(2.D0*A-1.D0)*(((17.D0CORA-656
     8*A+17.D0)*A+14.D0)*F+7.D0*A*(A**2-1.D0)*(A+2.D0)/W)/W))/W)/W)/W   CORA-657
      C(2)=B3+(((96.D0*((12.D0*A+56.D0)*A+59.D0-Q)-8.D0*(((188.D0*A+864.CORA-658
     1D0)*A+889.D0)*A-18.D0)*U)*Q+A*(A+1.D0)*(2.D0*A+1.D0)*(16.D0*((6.D0CORA-659
     2*A+37.D0)*A+70.D0)*U+4.D0*(2.D0*A+3.D0)*((67.D0*A+136.D0)*A-140.D0CORA-660
     3)*U**2))*F+3.D0*(A+2.D0)*(40.D0*((19.D0*A+46.D0)*Q-((12.D0*A+32.D0CORA-661
     4)*A+7.D0)*A-14.D0)*Q-(A+1.D0)*(20.D0*(((132.D0*A+292.D0)*A-97.D0)*CORA-662
     5A-42.D0)*Q-A**2*(2.D0*A+1.D0)*(20.D0*((66.D0*A+127.D0)*A-98.D0)+35CORA-663
     6.D0*(A-1.D0)*(2.D0*A+3.D0)*(A+2.D0)*U))*U+4.D0*(2.D0*A+5.D0)*(36.DCORA-664
     70*((54.D0*A+37.D0)*A-10.D0)*Q*F-(A+1.D0)*(A*(2.D0*A+1.D0)*(672.D0*CORA-665
     8(2.D0*A-1.D0)+6.D0*(2.D0*A+3.D0)*(17.D0*A-14.D0)*U)*F-5.D0*(14.D0*CORA-666
     9(((86.D0*A+55.D0)*A-3.D0)*A-18.D0)*Q-A**2*(2.D0*A+1.D0)*(112.D0*(2CORA-667
     A.D0*A-1.D0)*(3.D0*A+2.D0)+(2.D0*A+3.D0)*(7.D0*(A-1.D0)*(A+2.D0)*U-CORA-668
     B6.D0*(2.D0*A-1.D0)*(24.D0*F+7.D0*(A-1.D0)*(A+2.D0)/W)/W)))/W))/W)/CORA-669
     CW)/W                                                              CORA-670
      C(3)=B2+12.D0*Q*((8.D0*((12.D0*A+56.D0)*A+59.D0)-8.D0*Q-(2.D0*A+3.CORA-671
     1D0)*((32.D0*A+98.D0)*A-4.D0)*U)*F+5.D0*(A+2.D0)*((38.D0*A+92.D0)*QCORA-672
     2-(((24.D0*A+64.D0)*A+14.D0)*A+28.D0)-(2.D0*A+3.D0)*(7.D0*(A**2-1.DCORA-673
     30)*(A+2.D0)*U-6.D0*(2.D0*A+5.D0)*((10.D0*A-4.D0)*F+7.D0*(A**2-1.D0CORA-674
     4)*(A+2.D0)/W)/W))/W)/W                                            CORA-675
      C(4)=B4-3.D0*(((16.D0*(((56.D0*A+330.D0)*A+619.D0)*A+354.D0+(3.D0*CORA-676
     1A-6.D0)*Q)*Q-(A+1.D0)*(4.D0*(((194.D0*A+935.D0)*A+1292.D0)*A-36.D0CORA-677
     2)*Q*U+A*(2.D0*A+1.D0)*(32.D0*((6.D0*A+37.D0)*A+70.D0)-4.D0*(((62.DCORA-678
     30*A+349.D0)*A+524.D0)*A-140.D0)*U-35.D0*(A-1.D0)*(2.D0*A+3.D0)*(A+CORA-679
     42.D0)**2*U**2))))/(2.D0*A+3.D0)+4.D0*(4.D0*((46.D0*A+213.D0)*A+230CORA-680
     5.D0)*Q*F-(A+1.D0)*((2.D0*A+1.D0)*((24.D0*A+148.D0)*A+280.D0+(2.D0*CORA-681
     6A+3.D0)*((67.D0*A+136.D0)*A-140.D0)*U)*F-5.D0*(A+2.D0)*(4.D0*((118CORA-682
     7.D0*A+257.D0)*A-90.D0)*Q-4.D0*A*(2.D0*A+1.D0)*((66.D0*A+127.D0)*A-CORA-683
     898.D0)-(2.D0*A+1.D0)*(2.D0*A+3.D0)*(7.D0*A*(A-1.D0)*(A+2.D0)*U-6.DCORA-684
     90*(2.D0*A+5.D0)*((17.D0*A-14.D0)*F+7.D0*A*(A+2.D0)*(A-1.D0)/W)/W))CORA-685
     A/W))/W)/W                                                         CORA-686
      IF (X.EQ.0.D0) GO TO 45                                           CORA-687
      C(1)=C(1)-Q*X*(3.D0*A+6.D0)*(8.D0*((38.D0*A+167.D0)*A+171.D0-3.D0*CORA-688
     1Q)*F+4.D0*(((26.D0*A+145.D0)*A+231.D0)*A+106.D0)*X+(40.D0*A+64.D0)CORA-689
     2*Q*X+(2.D0*A+3.D0)*((8.D0*A+32.D0)*F*X**2-4.D0*((35.D0*A+101.D0)*ACORA-690
     3+6.D0)*F*U+4.D0*(A+1.D0)*(A+2.D0)*(X**3-(9.D0*A+4.D0)*U*X))-4.D0*(CORA-691
     42.D0*A+5.D0)*((((44.D0*A+94.D0)*A-22.D0)*A+28.D0)-4.D0*(19.D0*A+44CORA-692
     5.D0)*Q-(2.D0*A+3.D0)*(4.D0*(3.D0*A+8.D0)*F*X+(A+1.D0)*(A+2.D0)*(2.CORA-693
     6D0*X**2-(13.D0*A-7.D0)*U)+10.D0*(((11.D0*A+59.D0)*A-6.D0)*F+(A+1.DCORA-694
     70)*(A+2.D0)*((A+3.D0)*X+(3.D0*A+21.D0)*(2.D0*A-1.D0)/W))/W))/W)/W CORA-695
      C(2)=C(2)+X*(8.D0*(((196.D0*A+1104.D0)*A+1901.D0)*A+1026.D0+(15.D0CORA-696
     1*A-18.D0)*Q+12.D0*(5.D0*A+4.D0)*F*X)*Q+(A+1.D0)*(32.D0*((4.D0*A+11CORA-697
     2.D0)*A+18.D0)*Q*X**2-8.D0*(((212.D0*A+966.D0)*A+1267.D0)*A+54.D0)*CORA-698
     3Q*U+(2.D0*A+1.D0)*(8.D0*((38.D0*A+213.D0)*A+318.D0)*F*X+4.D0*A*(((CORA-699
     4(106.D0*A+599.D0)*A+912.D0)*A-140.D0)*U+((26.D0*A+107.D0)*A+46.D0)CORA-700
     5*X**2-8.D0*((6.D0*A+37.D0)*A+70.D0))+(2.D0*A+3.D0)*(8.D0*(A+6.D0)*CORA-701
     6F*X**3-4.D0*((47.D0*A+170.D0)*A+48.D0)*F*U*X+A*(A+2.D0)*(((87.D0*ACORA-702
     7+67.D0)*A-70.D0)*U**2-4.D0*(11.D0*A+9.D0)*U*X**2+4.D0*X**4))))+3.DCORA-703
     80*(A+2.D0)*(32.D0*((24.D0*A+107.D0)*A+110.D0)*F*Q+(A+1.D0)*(96.D0*CORA-704
     9((4.D0*A+17.D0)*A+20.D0)*Q*X-(8.D0*A+4.D0)*(2.D0*((6.D0*A+37.D0)*ACORA-705
     A+70.D0)*F+A*((22.D0*A+125.D0)*A+194.D0)*X+(2.D0*A+3.D0)*(((61.D0*ACORA-706
     B+143.D0)*A-70.D0)*F*U+A*(A+2.D0)*((9.D0*A+4.D0)*U-X**2)*X-(6.D0*A+CORA-707
     C20.D0)*F*X**2))+(8.D0*A+20.D0)*(((532.D0*A+2246.D0)*A-180.D0)*Q-(2CORA-708
     D.D0*A+1.D0)*(4.D0*A*((66.D0*A+295.D0)*A-98.D0)-(2.D0*A+3.D0)*((22.CORA-709
     ED0*A+60.D0)*F*X+A*(A+2.D0)*(2.D0*X**2-(13.D0*A-7.D0)*U)+10.D0*(((1CORA-710
     F7.D0*A+93.D0)*A-42.D0)*F+A*(A+2.D0)*((A+3.D0)*X+3.D0*(A+7.D0)*(2.DCORA-711
     G0*A-1.D0)/W))/W)))/W))/W)/W                                       CORA-712
      C(3)=C(3)+12.D0*Q*X*(2.D0*(((28.D0*A+148.D0)*A+221.D0)*A+106.D0)+(CORA-713
     16.D0*A-4.D0)*Q+(2.D0*A+3.D0)*((2.D0*A+12.D0)*F*X+(A+2.D0)*((A+1.D0CORA-714
     2)*(2.D0*X**2-(13.D0*A+8.D0)*U)+5.D0*((10.D0*A+28.D0)*F+2.D0*(A+1.DCORA-715
     30)*(A+2.D0)*(X+(6.D0*A+15.D0)/W))/W)))/W                          CORA-716
      C(4)=C(4)+3.D0*X*((24.D0*(2.D0*A-1.D0)*(A-2.D0)*F*Q-(4.D0*A+4.D0)*CORA-717
     1(((20.D0*A+46.D0)*A+108.D0)*Q*X+(2.D0*A+1.D0)*(((76.D0*A+426.D0)*ACORA-718
     2+636.D0)*F+A*((26.D0*A+107.D0)*A+46.D0)*X)))/(2.D0*A+3.D0)+(4.D0*ACORA-719
     3+4.D0)*((2.D0*A+1.D0)*((35.D0*A+130.D0)*A+48.D0)*F*U+(2.D0*A+1.D0)CORA-720
     4*(A*(A+2.D0)*(9.D0*(A+1.D0)*U-X**2)*X-(2.D0*A+12.D0)*F*X**2)-(4.D0CORA-721
     5*((38.D0*A+169.D0)*A+210.D0)*Q-(2.D0*A+1.D0)*(A*((44.D0*A+250.D0)*CORA-722
     6A+388.D0)-(2.D0*A+3.D0)*((12.D0*A+40.D0)*F*X-(A+2.D0)*(A*((13.D0*ACORA-723
     7+8.D0)*U-2.D0*X**2)-10.D0*((11.D0*A+30.D0)*F+A*(A+2.D0)*(X+(6.D0*ACORA-724
     8+15.D0)/W))/W))))/W))/W                                           CORA-725
      GO TO 45                                                          CORA-726
C MULTIPOLE LQ=5 L1=L2-5.                                               CORA-727
   39 BD=181440.D0*DEN*S**2*DSQRT(Q*((A+1.D0)**2+F**2)*((A+2.D0)**2+F**2CORA-728
     1)*((A-1.D0)**2+E**2)*((A+3.D0)**2+F**2))*A**2                     CORA-729
      B1=Q*(2304.D0*A*(2.D0*A-1.D0)-6912.D0*Q+864.D0*(9.D0*A+1.D0)*Q*U-2CORA-730
     188.D0*A*((22.D0*A+3.D0)*A+5.D0)*U+8.D0*A*(((530.D0*A+611.D0)*A+12.CORA-731
     2D0)*A+72.D0)*U**2-8.D0*((602.D0*A+683.D0)*A-60.D0)*Q*U**2+96.D0*A*CORA-732
     3(A**2-1.D0)*(2.D0*A+5.D0)*(2.D0*A+3.D0)*U**3)                     CORA-733
      B2=24.D0*Q*F*(288.D0-(204.D0*A+36.D0)*U+((89.D0*A+141.D0)*A-20.D0)CORA-734
     1*U**2)                                                            CORA-735
      B3=12.D0*F*(Q*(576.D0-24.D0*(37.D0*A+3.D0)*U+((718.D0*A+582.D0)*A-CORA-736
     140.D0)*U**2)+A*(A+1.D0)*(2.D0*A+1.D0)*U*(96.D0-(148.D0*A+12.D0)*U-CORA-737
     2((137.D0*A+148.D0)*A-180.D0)*U**2))                               CORA-738
      B4=Q*(864.D0*(9.D0*A+1.D0)*U-6912.D0-8.D0*(602.D0*A**2+613.D0*A-60CORA-739
     1.D0)*U**2)+A*(2.D0*A+1.D0)*(2304.D0-288.D0*(11.D0*A-1.D0)*U+8.D0*(CORA-740
     2(265.D0*A+174.D0)*A-54.D0)*U**2+48.D0*(2.D0*A+3.D0)*(2.D0*A+5.D0)*CORA-741
     3(A-1.D0)*U**3)                                                    CORA-742
      IF (X.EQ.0.D0) GO TO 40                                           CORA-743
      B1=B1+X*Q*((1176.D0*F*U+((2040.D0*A+2616.D0)*U-1200.D0-48.D0*X**2)CORA-744
     1*X)*Q+(((2944.D0*A+4432.D0)*A+912.D0)*X+8.D0*((794.D0*A+1535.D0)*ACORA-745
     2+222.D0)*F)*U-(16.D0*((94.D0*A+385.D0)*A+216.D0)+((368.D0*A+1016.DCORA-746
     30)*A+600.D0)*X**2)*X+(2*A+3.D0)*((((608.D0*A+632.D0)*U-600.D0-24.DCORA-747
     40*X**2)*X**2-3456.D0-((1950.D0*A+2614.D0)*A-140.D0)*U**2)*F+(A+1.DCORA-748
     50)*((120.D0*A+152.D0)*U*X**3-((414.D0*A+798.D0)*A+20.D0)*U**2*X-8.CORA-749
     6D0*X**5)))                                                        CORA-750
      B2=B2+X*Q*(1200.D0*F*X-1176.D0*Q*U-24.D0*((33.D0*A+120.D0)*A+38.D0CORA-751
     1)*U-24.D0*(53.D0*A+60.D0)*F*U*X+48.D0*F*X**3+(A+1.D0)*(3456.D0+((7CORA-752
     238.D0*A+1722.D0)*A+60.D0)*U**2+24.D0*(25.D0-(13.D0*A+19.D0)*U+X**2CORA-753
     3)*X**2))                                                          CORA-754
      B3=B3+X*((1200.D0*F*X-24.D0*((597.D0*A+720.D0)*A+74.D0)*U-1176.D0*CORA-755
     1Q*U-(2808.D0*A+2616.D0)*F*U*X+48.D0*F*X**3)*Q+(A+1.D0)*((10368.D0-CORA-756
     224.D0*(113.D0*A+79.D0)*U*X**2+(72.D0*X**2+1800.D0)*X**2+((11782.D0CORA-757
     3*A+10730.D0)*A-420.D0)*U**2)*Q+(2*A+1.D0)*((((3090.D0*A+3898.D0)*ACORA-758
     4+60.D0)*U**2+3456.D0+600.D0*X**2-304.D0*(14.D0*A+3.D0)*U-(680.D0*ACORA-759
     5+456.D0)*U*X**2+24.D0*X**4)*F*X+A*(((8.D0*X**2+184.D0)*X**2+752.D0CORA-760
     6)*X**2-2304.D0-((3.D0*(((187.D0*A+479.D0)*A-6.D0)*A-240.D0)*U+4.D0CORA-761
     7*((1175.D0*A+624.D0)*A-108.D0))*U-96.D0*(55.D0*A-3.D0))*U+(((570.DCORA-762
     80*A+1278.D0)*A+356.D0)*U-60.D0*A-440.D0-44.D0*(3.D0*A+4.D0)*X**2)*CORA-763
     9U*X**2))))                                                        CORA-764
      B4=B4+X*((1176.D0*F*U-1200.D0*X+120.D0*(17.D0*A+12.D0)*U*X-48.D0*XCORA-765
     1**3)*Q+(2.D0*A+1.D0)*((8.D0*(397.D0*A+114.D0)*U-3456.D0-10.D0*((19CORA-766
     25.D0*A+289.D0)*A+6.D0)*U**2-600.D0*X**2+8.D0*(76.D0*A+57.D0)*U*X**CORA-767
     32-24.D0*X**4)*F+A*X*(440.D0*U-752.D0-((414.D0*A+1050.D0)*A+356.D0)CORA-768
     4*U**2-184.D0*X**2+8.D0*(15.D0*A+22.D0)*U*X**2-8.D0*X**4)))        CORA-769
   40 IF (LT) GO TO 45                                                  CORA-770
      C(1)=B1+24.D0*Q*((60.D0*((9.D0*A+5.D0)*U-8.D0)*Q+96.D0*(A-1.D0)*(2CORA-771
     1.D0*A-1.D0)-4.D0*(((74.D0*A+15.D0)*A+28.D0)*A+3.D0)*U-(2.D0*A+3.D0CORA-772
     2)*(A+2.D0)*((137.D0*A-37.D0)*A-30.D0)*U**2)*F+4.D0*(30.D0*((8.D0*ACORA-773
     3-6.D0)*A+5.D0-7.D0*Q)*Q+(A+2.D0)*(5.D0*((128.D0*A+80.D0)*A-33.D0)*CORA-774
     4Q-A*(((676.D0*A+544.D0)*A-201.D0)*A-144.D0))*U-(6.D0*A-6.D0)*((2.DCORA-775
     50*A-1.D0)*(4.D0*A-6.D0)*A+(A+1.D0)*(2.D0*A+3.D0)*(2.D0*A+5.D0)*A*(CORA-776
     6A+2.D0)*U**2)+5.D0*((12.D0*(2.D0*A-1.D0)*(((50.D0*A+153.D0)*A+136.CORA-777
     7D0)*A+39.D0)-84.D0*((10.D0*A+19.D0)*A+3.D0)*Q)*F+(A+2.D0)*((2.D0*ACORA-778
     8+3.D0)*(2.D0*A+5.D0)*((37.D0*A-17.D0)*A-6.D0)*U*F+6.D0*(A*(2.D0*A-CORA-779
     91.D0)*(((236.D0*A+648.D0)*A+625.D0)*A+66.D0)-3.D0*(((140.D0*A+294.CORA-780
     AD0)*A+196.D0)*A-105.D0)*Q+(2.D0*A+3.D0)*(2.D0*A+5.D0)*(A*(A**2-1.DCORA-781
     B0)*(2.D0*A+6.D0)*U-(14.D0*A-7.D0)*(((7.D0*A+5.D0)*A+6.D0)*F+A*(A**CORA-782
     C2-1.D0)*(2.D0*A+6.D0)/W)/W))/W))/W)/W)/W                          CORA-783
      C(2)=B3+8.D0*((1440.D0*Q-288.D0*(4.D0*A**2+1.D0)+12.D0*(((168.D0*ACORA-784
     1+132.D0)*A+62.D0)*A+3.D0)*U-30.D0*(68.D0*A+30.D0)*Q*U)*Q+(A+1.D0)*CORA-785
     2(8.D0*A+12.D0)*(((338.D0*A-83.D0)*A-45.D0)*Q*U-(2.D0*A+1.D0)*A**2*CORA-786
     3(12.D0+(169.D0*A-99.D0)*U+(A-1.D0)*(6.D0*A+15.D0)*U**2))*U+3.D0*((CORA-787
     4(840.D0*Q-120.D0*(4.D0*A**2+5.D0))*Q-(2.D0*A+3.D0)*(40.D0*((65.D0*CORA-788
     5A+30.D0)*A-11.D0)*Q-5.D0*A*(2.D0*A+1.D0)*(A+1.D0)*(12.D0*(25.D0*A-CORA-789
     613.D0)+(A+2.D0)*(37.D0*A-30.D0)*U))*U)*F+4.D0*((630.D0*((10.D0*A+1CORA-790
     75.D0)*A+2.D0)*Q-90.D0*((((200.D0*A+420.D0)*A+162.D0)*A-63.D0)*A-26CORA-791
     8.D0))*Q+(2.D0*A+3.D0)*((A+1.D0)*(1512.D0*(4.D0*A**2-1.D0)*A**2-(A+CORA-792
     92.D0)*(30.D0*((44.D0*A-18.D0)*A-5.D0)*Q-6.D0*A**2*(2.D0*A+1.D0)*(1CORA-793
     A18.D0*A-83.D0+(A-1.D0)*(2.D0*A+5.D0)*U))*U)+5.D0*(42.D0*(((70.D0*ACORA-794
     B+129.D0)*A+53.D0)*A-30.D0)*Q*F-(A+1.D0)*(A*(2.D0*A+1.D0)*(252.D0*(CORA-795
     C2.D0*A-1.D0)*(4.D0*A+5.D0)+7.D0*(A+2.D0)*(2.D0*A+5.D0)*(7.D0*A-6.DCORA-796
     D0)*U)*F-6.D0*(A+2.D0)*(14.D0*(((56.D0*A+36.D0)*A+4.D0)*A-15.D0)*Q-CORA-797
     EA**2*(2.D0*A+1.D0)*((4.D0*A+3.D0)*54.D0*(2.D0*A-1.D0)+(2.D0*A+5.D0CORA-798
     F)*((A-1.D0)*(2.D0*A+6.D0)*U-7.D0*(2.D0*A-1.D0)*(9.D0*F+(A-1.D0)*(2CORA-799
     G.D0*A+6.D0)/W)/W)))/W))/W))/W)/W)/W                               CORA-800
      C(3)=B2+96.D0*Q*((120.D0-25.D0*(4.D0*A+3.D0)*U)*Q-24.D0*(4.D0*A**2CORA-801
     1+1.D0)+((((88.D0*A+72.D0)*A+12.D0)*A+3.D0)+3.D0*(A**2-1.D0)*(2.D0*CORA-802
     2A+3.D0)*(2.D0*A+5.D0)*U)*U+5.D0*((42.D0*Q-(24.D0*A**2+30.D0))*F-(ACORA-803
     3+2.D0)*((2.D0*A+3.D0)*(25.D0*A-11.D0)*U*F-6.D0*(21.D0*(4.D0*A+1.D0CORA-804
     4)*Q-3.D0*(((32.D0*A+20.D0)*A-4.D0)*A-13.D0)-(2.D0*A+3.D0)*(2.D0*A+CORA-805
     55.D0)*((A**2-1.D0)*U-7.D0*((5.D0*A-3.D0)*F+(A**2-1.D0)*(2.D0*A+6.DCORA-806
     60)/W)/W))/W))/W)/W                                                CORA-807
      C(4)=B4+24.D0*((60.D0*((9.D0*A+5.D0)*U-8.D0)*Q+(A+1.D0)*(2.D0*A+1.CORA-808
     1D0)*(96.D0-4.D0*(37.D0*A+3.D0)*U-((137.D0*A+148.D0)*A-180.D0)*U**2CORA-809
     2))*F+4.D0*((30.D0*((8.D0*A+6.D0)*A+5.D0)-210.D0*Q)*Q+(2.D0*A+3.D0)CORA-810
     3*((A+1.D0)*(10.D0*(32.D0*A-11.D0)*Q*U-(2.D0*A+1.D0)*(A*(12.D0+(169CORA-811
     4.D0*A-99.D0)*U)+3.D0*(A-1.D0)*(2.D0*A+5.D0)*A*U**2))-5.D0*(84.D0*(CORA-812
     55.D0*A+1.D0)*Q*F-(A+1.D0)*((2.D0*A+1.D0)*(12.D0*(25.D0*A-13.D0)+(ACORA-813
     6+2.D0)*(37.D0*A-30.D0)*U)*F-6.D0*(A+2.D0)*(105.D0*(2.D0*A-1.D0)*Q-CORA-814
     7(2.D0*A+1.D0)*((118.D0*A-83.D0)*A+(2.D0*A+5.D0)*(A*(A-1.D0)*U-7.D0CORA-815
     8*((7.D0*A-6.D0)*F+A*(A-1.D0)*(2.D0*A+6.D0)/W)/W)))/W))/W))/W)/W   CORA-816
      IF (X.EQ.0.D0) GO TO 45                                           CORA-817
      C(1)=C(1)-Q*X*(720.D0*(((6.D0*A+23.D0)*A+19.D0)+Q)*F*X+(A+2.D0)*(7CORA-818
     120.D0*Q*X**2-384.D0*((22.D0*A-5.D0)*A+3.D0)+11520.D0*Q+24.D0*(((43CORA-819
     20.D0*A+255.D0)*A+29.D0)*A+6.D0)*U-72.D0*(173.D0*A+107.D0)*Q*U+120.CORA-820
     3D0*((2.D0*A+17.D0)*A+9.D0)*X**2+(2.D0*A+3.D0)*(144.D0*F*X**3-24.D0CORA-821
     4*(95.D0*A+71.D0)*F*U*X+(A+1.D0)*(6.D0*((123.D0*A+137.D0)*A-120.D0)CORA-822
     5*U**2-24.D0*(13.D0*A+14.D0)*U*X**2+24.D0*X**4)))-24.D0*((2.D0*A+5.CORA-823
     6D0)*((20.D0*((10.D0*A-7.D0)*A+9.D0)-420.D0*Q)*F+(A+2.D0)*((10.D0*(CORA-824
     7(14.D0*A-1.D0)*A+3.D0)-180.D0*Q)*X+(2.D0*A+3.D0)*((199.D0*A-31.D0)CORA-825
     8*F*U-20.D0*F*X**2+(A+1.D0)*((21.D0*A+13.D0)*U-2.D0*X**2)*X)))+5.D0CORA-826
     9*(A+2.D0)*(4*((((220.D0*A+1340.D0)*A+973.D0)*A+123.D0)*A-234.D0)-8CORA-827
     A4.D0*((10.D0*A+55.D0)*A+21.D0)*Q+(2.D0*A+3.D0)*(2.D0*A+5.D0)*((A+1CORA-828
     B.D0)*(A+3.D0)*((15.D0*A-8.D0)*U-2.D0*X**2)-30.D0*(A+3.D0)*F*X-6.D0CORA-829
     C*(21.D0*((2.D0*A+13.D0)*A-3.D0)*F+(A+1.D0)*(A+3.D0)*((2.D0*A+7.D0)CORA-830
     D*X+(14.D0*A-7.D0)*(A+8.D0)/W))/W))/W)/W)/W                        CORA-831
      C(2)=C(2)+X*(720.D0*((12.D0*A+32.D0)*A+19.D0+Q)*Q*X+(2.D0*A+3.D0)*CORA-832
     1((7680.D0+480.D0*X**2-48.D0*(222.D0*A+107.D0)*U)*F*Q+(A+1.D0)*((28CORA-833
     28.D0*X**2-48.D0*(127.D0*A+71.D0)*U)*Q*X+(2.D0*A+1.D0)*(((48.D0*X**CORA-834
     32+720.D0)*X**2-768.D0+(32.D0*(62.D0*A+3.D0)+16.D0*((168.D0*A+139.DCORA-835
     40)*A-90.D0)*U-(920.D0*A+672.D0)*X**2)*U)*F+A*X*((8.D0*X**2+40.D0-(CORA-836
     5120.D0*A+152.D0)*U)*X**2+((414.D0*A+798.D0)*A+20.D0)*U**2-1408.D0+CORA-837
     6(2280.D0*A+736.D0)*U))))+6.D0*((A+2.D0)*((4200.D0*Q-200.D0*((20.D0CORA-838
     7*A+8.D0)*A+9.D0))*Q+(2.D0*A+3.D0)*(1200.D0*F*Q*X+(A+1.D0)*((400.D0CORA-839
     8*X**2-(4960.D0*A-620.D0)*U)*Q+(2.D0*A+1.D0)*((40.D0*X**2-200.D0-(5CORA-840
     948.D0*A+260.D0)*U)*F*X+A*((((123.D0*A+137.D0)*A-120.D0)*U+(2452.D0CORA-841
     A*A-792.D0))*U+96.D0+4.D0*X**2*(X**2-35.D0-(13.D0*A+14.D0)*U))))))+CORA-842
     B(8.D0*A+12.D0)*(420.D0*((10.D0*A+45.D0)*A+14.D0)*F*Q-(A+1.D0)*(240CORA-843
     C.D0*(2.D0*A+1.D0)*((10.D0*A+45.D0)*A-13.D0)*F-(A+2.D0)*(450.D0*(2.CORA-844
     DD0*A+5.D0)*Q*X+(2.D0*A+1.D0)*(2.D0*A+5.D0)*(30.D0*F*X**2-(274.D0*ACORA-845
     E-120.D0)*F*U+A*X*(2.D0*X**2-220.D0-(21.D0*A+13.D0)*U))+5.D0*(42.D0CORA-846
     F*((42.D0*A+231.D0)*A-45.D0)*Q-(2.D0*A+1.D0)*(8.D0*A*((118.D0*A+649CORA-847
     G.D0)*A-249.D0)+(A+3.D0)*(2.D0*A+5.D0)*(A*(15.D0*A-8.D0)*U-42.D0*F*CORA-848
     HX-2.D0*A*X**2))+(12.D0*A+6.D0)*(2.D0*A+5.D0)*(28.D0*((2.D0*A+13.D0CORA-849
     I)*A-6.D0)*F+A*(A+3.D0)*((2.D0*A+7.D0)*X+(14.D0*A-7.D0)*(A+8.D0)/W)CORA-850
     J)/W)/W)))/W)/W)/W                                                 CORA-851
      C(3)=C(3)+24.D0*X*Q*((30.D0*Q+10.D0*((4.D0*A+16.D0)*A+9.D0))*X+(2.CORA-852
     1D0*A+3.D0)*((160.D0-(124.D0*A+114.D0)*U+10.D0*X**2)*F+(A+1.D0)*(2.CORA-853
     2D0*X**2-(21.D0*A+28.D0)*U)*X)+(5.D0*A+10.D0)*(126.D0*Q-2.D0*((52.DCORA-854
     30*A+4.D0)*A+15.D0)+(2.D0*A+3.D0)*(18.D0*F*X+(A+1.D0)*(2.D0*X**2-(1CORA-855
     45.D0*A+13.D0)*U))+(12.D0*A+18.D0)*(2.D0*A+5.D0)*(14.D0*F+(A+1.D0)*CORA-856
     5(X+(7.D0*A+21.D0)/W))/W)/W)/W                                     CORA-857
      C(4)=C(4)-X*(720.D0*F*Q*X+(A+1.D0)*((11520.D0+720.D0*X**2-72.D0*(1CORA-858
     173.D0*A+114.D0)*U)*Q-(2.D0*A+1.D0)*((24.D0*(95.D0*A+84.D0)*U-2160.CORA-859
     2D0-144.D0*X**2)*F*X+A*(4224.D0-24.D0*(215.D0*A+92.D0)*U-6.D0*((123CORA-860
     3.D0*A+287.D0)*A+10.D0)*U**2-120.D0*X**2+24.D0*(13.D0*A+19.D0)*U*X*CORA-861
     4*2-24.D0*X**4)))+24.D0*(2.D0*A+3.D0)*(420.D0*F*Q+(A+1.D0)*(180.D0*CORA-862
     5Q*X-7.D0*A*(2.D0*A+1.D0)*(3.D0*A+4.D0)*U*X-(2.D0*A+1.D0)*((100.D0+CORA-863
     6(199.D0*A+130.D0)*U-20.D0*X**2)*F+A*X*(70.D0-2.D0*X**2))+5.D0*(A+2CORA-864
     7.D0)*(420.D0*Q+(2.D0*A+1.D0)*(30.D0*F*X+A*(2.D0*X**2-220.D0-(15.D0CORA-865
     8*A+13.D0)*U)+6.D0*(2.D0*A+5.D0)*(21.D0*F+A*(X+(7.D0*A+21.D0)/W))/WCORA-866
     9))/W))/W)/W                                                       CORA-867
      GO TO 45                                                          CORA-868
   41 IF ((X.NE.0.D0).OR.(MOD(LL,2).NE.1).OR.(LL.EQ.7)) GO TO 72        CORA-869
      GO TO ( 42 , 71 , 43 , 71 , 44 , 71 , 71 ) , LL                   CORA-870
C MULTIPOLE LQ=6 L1=L2.                                                 CORA-871
   42 BD=120.D0*DEN*A**2*(2.D0*A+7.D0)*(2.D0*A+5.D0)*(A+2.D0)*(A+3.D0)*(CORA-872
     1A+4.D0)*DSQRT(S*Q*P*((A+1.D0)**2+F**2)*((A+1.D0)**2+E**2))*S**2   CORA-873
      B1=P*Q*((8.D0*(((154.D0*A-205.D0)*A-376.D0)*A-140.D0+((47.D0*A-282CORA-874
     1.D0)*A-143.D0-63.D0*Q)*Q)+4.D0*((((47.D0*A+26.D0)*A-553.D0)*A-752.CORA-875
     2D0)*A-280.D0-((79.D0*A+282.D0)*A+143.D0)*Q)*U)/(A-1.D0)+2.D0*(2.D0CORA-876
     3*A+7.D0)*(2.D0*A+5.D0)*(A+1.D0)**2*U**2)                          CORA-877
      B2=(240.D0*(2.D0*A+1.D0)**2*(A-14.D0)+24.D0*((20.D0*A-336.D0)*A-14CORA-878
     13.D0)*Q-1512.D0*Q**2-6.D0*((((36.D0*A+836.D0)*A+2797.D0)*A+2374.D0CORA-879
     2)*A+560.D0+(((18.D0*A+275.D0)*A+726.D0)*A+286.D0)*Q+(A+1.D0)*(((26CORA-880
     3.D0*A+146.D0)*A+210.D0)*A+35.D0)*U)*U)/(2.D0*A+3.D0)              CORA-881
      B3=P*F*(B2+A*(2.D0*A+4.D0)*U*(54.D0*Q+6.D0*(38.D0*A+29.D0)+3.D0*(ACORA-882
     1+1.D0)*(9.D0*A+29.D0)*U))                                         CORA-883
      B2=Q*F*B2                                                         CORA-884
      B4=(((((((32.D0*A-216.D0)*A-1692.D0)*A-2272.D0)*A+454.D0)*A+1402.DCORA-885
     10)*A+210.D0)*Q*U**2+(2.D0*A-3.D0)*((16.D0*A+8.D0)*((107.D0*A-444.DCORA-886
     20)*A-140.D0)+8.D0*((((47.D0*A-176.D0)*A-924.D0)*A-724.D0)*A-140.D0CORA-887
     3)*U+(8.D0*((47.D0*A-390.D0)*A-143.D0)-4.D0*((142.D0*A+390.D0)*A+14CORA-888
     43.D0)*U-504.D0*Q)*Q)*Q+A**2*(2.D0*A+1.D0)*(480.D0*(2.D0*A+1.D0)*(ACORA-889
     5-14.D0)+(4.D0*(((214.D0*A-849.D0)*A-3808.D0)*A-2100.D0)+((((204.D0CORA-890
     6*A+752.D0)*A-232.D0)*A-2550.D0)*A-1680.D0)*U)*U))/((A+1.D0)*(2.D0*CORA-891
     7A+3.D0))+A**2*(2.D0*A+5.D0)*(2.D0*A+7.D0)*(2.D0*A+1.D0)*(A+1.D0)*UCORA-892
     8**3                                                               CORA-893
      IF (LT) GO TO 45                                                  CORA-894
      C(1)=B1+12.D0*P*Q*(A+2.D0)*((58.D0+76.D0*A+18.D0*Q+(A+1.D0)*(9.D0*CORA-895
     1A+29.D0)*U)*F-(2.D0*A+5.D0)*((28.D0*A+32.D0)*A+14.D0+(4.D0*A-6.D0)CORA-896
     2*Q+(A+1.D0)*((A+1.D0)*(2.D0*A+7.D0)*U+10.D0*(A+3.D0)*(2.D0*F-(A+1.CORA-897
     3D0)*(2.D0*A+7.D0)/W)/W))/W)/W                                     CORA-898
      C(2)=B3-2.D0*P*(2.D0*A*((((48.D0*A+16.D0)*A-570.D0)*A-703.D0)*A-21CORA-899
     10.D0)*U+(2.D0*A-3.D0)*((((4.D0*A-34.D0)*A-172.D0)*A-116.D0)*U+8.D0CORA-900
     2*((5.D0*A-48.D0)*A-29.D0)-72.D0*Q)*Q-A*(2.D0*A+1.D0)*(40.D0*(2.D0*CORA-901
     3A+1.D0)*(14.D0-A)-(2.D0*A+7.D0)*(2.D0*A+5.D0)*(A+1.D0)**2*U**2)-6.CORA-902
     4D0*(A+2.D0)*(5.D0*((14.D0-A)*(2.D0*A+1.D0)**2-(A-2.D0)*(2.D0*A-3.DCORA-903
     50)*Q-(A+1.D0)*(2.D0*A+1.D0)*((A+1.D0)*A-7.D0)*U)*F+(2.D0*A+5.D0)*(CORA-904
     6(((4.D0*A+6.D0)*A+92.D0)*A+60.D0)*Q+(2.D0*A+1.D0)*((A-14.D0)*A*(4.CORA-905
     7D0*A+2.D0)+(A+1.D0)*(A*(A+1.D0)*(2.D0*A+7.D0)*U-10.D0*(A+3.D0)*(7.CORA-906
     8D0*F+(2.D0*A+7.D0)*A*(A+1.D0)/W)/W)))/W)/W)/W                     CORA-907
      C(3)=B2+2.D0*Q*(A*((((16.D0*A+600.D0)*A+2332.D0)*A+2298.D0)*A+700.CORA-908
     1D0)*U+(2.D0*A-3.D0)*((((4.D0*A+70.D0)*A+212.D0)*A+116.D0)*U-8.D0*(CORA-909
     2(5.D0*A-48.D0)*A-29.D0)+72.D0*Q)*Q+40.D0*(2.D0*A+1.D0)**2*A*(14.D0CORA-910
     3-A)+3.D0*A*(2.D0*A+5.D0)*(2.D0*A+7.D0)*(A+1.D0)**2*U**2+6.D0*(A+2.CORA-911
     4D0)*(5.D0*((14.D0-A)*(2.D0*A+1.D0)**2-(A-2.D0)*(2.D0*A-3.D0)*Q+(A+CORA-912
     51.D0)*((8.D0*A+28.D0)*A+7.D0)*U)*F+(2.D0*A+5.D0)*((((4.D0*A+6.D0)*CORA-913
     6A+92.D0)*A+60.D0)*Q+2.D0*(A-14.D0)*A*(2.D0*A+1.D0)**2-5.D0*A*(A+1.CORA-914
     7D0)**2*(2.D0*A+7.D0)*U-10.D0*(A+3.D0)*(A+1.D0)*(2.D0*A+1.D0)*(7.D0CORA-915
     8*F+(2.D0*A+7.D0)*A*(A+1.D0)/W)/W)/W)/W)/W                         CORA-916
      C(4)=B4+6.D0*(2.D0*F*((2.D0*A-3.D0)*(A-2.D0)*(((18.D0*A+58.D0)*A+2CORA-917
     19.D0)*U+58.D0*(2.D0*A+1.D0)+18.D0*Q)*Q+A*(2.D0*A+1.D0)*(40.D0*(2.DCORA-918
     20*A+1.D0)*(A-14.D0)+(((58.D0*A-123.D0)*A-906.D0)*A-560.D0)*U+(A+1.CORA-919
     3D0)*(((9.D0*A+43.D0)*A+34.D0)*A-35.D0)*U**2))/(2.D0*A+3.D0)-((((((CORA-920
     416.D0*A**2+60.D0)*A+820.D0)*A+794.D0)*A+140.D0)*U+4.D0*(2.D0*A-3.DCORA-921
     50)*(2.D0*A-5.D0)*(A-2.D0)*Q)*Q-(2.D0*A+1.D0)*(A**2*(40.D0*(14.D0-ACORA-922
     6)*(2.D0*A+1.D0)-(A+1.D0)**2*(2.D0*A+5.D0)*(2.D0*A+7.D0)*U**2)-4.D0CORA-923
     7*(((24.D0*A-44.D0)*A+257.D0)*A+70.D0)*Q-(((48.D0*A-8.D0)*A-566.D0)CORA-924
     8*A-420.D0)*A**2*U+10.D0*(A+2.D0)*(4.D0*((14.D0-A)*(2.D0*A+1.D0)*A-CORA-925
     9((2.D0*A+2.D0)*A+15.D0)*Q-A*(A+1.D0)*((A+1.D0)*A-7.D0)*U)*F+(2.D0*CORA-926
     AA+5.D0)*(2.D0*(A-14.D0)*(2.D0*A+1.D0)*A**2+2.D0*(((2.D0*A+9.D0)*A+CORA-927
     B43.D0)*A+21.D0)*Q+A*(A+1.D0)*(A*(A+1.D0)*(2.D0*A+7.D0)*U+(4.D0*A+2CORA-928
     C.D0)*(A+3.D0)*((A+7.D0)*2.D0*F+(2.D0*A+7.D0)*A*(A+1.D0)/W)/W))/W)/CORA-929
     DW))/W)/W                                                          CORA-930
      GO TO 45                                                          CORA-931
C MULTIPOLE LQ=6 L1=L2-2.                                               CORA-932
   43 BD=2520.D0*DEN*A**2*(2.D0*A+5.D0)*(2.D0*A+7.D0)*(A+2.D0)*(A+3.D0)*CORA-933
     1(A+4.D0)*S*DSQRT(S*Q*P*((A+1.D0)**2+F**2)*((A+2.D0)**2+F**2))     CORA-934
      B1=2*P*Q*(120.D0*((((2.D0*A+57.D0)*A-117.D0)*A-306.D0)*A-140.D0+((CORA-935
     111.D0*A-129.D0)*A-113.D0-21.D0*Q)*Q)/(A-1.D0)-((((74.D0*A-865.D0)*CORA-936
     2A-5208.D0)*A-7344.D0)*A-3360.D0)*U+((266.D0*A+863.D0)*A+312.D0)*Q*CORA-937
     3U-3.D0*(A+1.D0)*(A+2.D0)*(2.D0*A+3.D0)*(2.D0*A+5.D0)*(2.D0*A+7.D0)CORA-938
     4*U**2)                                                            CORA-939
      B2=6.D0*Q*F*((120.D0*((2.D0*A-147.D0)*A-113.D0-21.D0*Q)*Q+60.D0*((CORA-940
     1((4.D0*A-56.D0)*A-261.D0)*A-259.D0)*A-56.D0)*U+(((86.D0*A-573.D0)*CORA-941
     2A-2171.D0)*A-312.D0)*Q*U-16800.D0*(A+1.D0)*(2.D0*A+1.D0))/(2.D0*A+CORA-942
     33.D0)+(A+2.D0)*(((43.D0*A+193.D0)*A+210.D0)*A+105.D0)*U**2)       CORA-943
      B3=3*P*((240.D0*((2.D0*A-147.D0)*A-113.D0-21.D0*Q)*Q+40.D0*(A+1.D0CORA-944
     1)*(2.D0*A+1.D0)*(((2.D0*A+63.D0)*A+69.D0)*A-168.D0)*U+(((892.D0*A+CORA-945
     21374.D0)*A-2182.D0)*A-624.D0)*Q*U-33600.D0*(A+1.D0)*(2.D0*A+1.D0))CORA-946
     3/(2.D0*A+3.D0)-(A+1.D0)*(A+2.D0)*(2.D0*A+1.D0)*((19.D0*A-6.D0)*A-2CORA-947
     410.D0)*U**2)*F                                                    CORA-948
      B4=(((2.D0*A-3.D0)*(240.D0*((((2.D0*A+69.D0)*A-369.D0)*A-534.D0)*ACORA-949
     1-140.D0+((11.D0*A-165.D0)*A-113.D0-21.D0*Q)*Q)*Q+((((532.D0*A-1326CORA-950
     2.D0)*A-5422.D0)*A-624.D0)*Q-((((148.D0*A-4110.D0)*A+4130.D0)*A+451CORA-951
     320.D0)*A+40512.D0)*A-6720.D0)*Q*U)/(A+1.D0)+240.D0*A**2*(2.D0*A+1.CORA-952
     4D0)*((A+2.D0)*((2.D0*A+59.D0)*A-189.D0)*U-840.D0))/(2.D0*A+3.D0)-(CORA-953
     52.D0*(((((24.D0*A-26.D0)*A+317.D0)*A+2298.D0)*A+2001.D0)*A+630.D0)CORA-954
     6*Q+A**2*(2.D0*A+1.D0)*(((74.D0*A-1068.D0)*A-4812.D0)*A-3360.D0+3*(CORA-955
     7A+2.D0)*(2.D0*A+3.D0)*(2.D0*A+5.D0)*(2.D0*A+7.D0)*U))*U**2)       CORA-956
      IF (LT) GO TO 45                                                  CORA-957
      C(1)=B1+6.D0*P*Q*(A+2.D0)*((360.D0*Q+40.D0*((2.D0*A+53.D0)*A+57.D0CORA-958
     1)-(2.D0*A+3.D0)*((19.D0*A+37.D0)*A-72.D0)*U)*F-(4.D0*A+10.D0)*(20.CORA-959
     2D0*(2.D0*A-3.D0)*Q+10.D0*(((2.D0*A+45.D0)*A+75.D0)*A+42.D0)-(A+1.DCORA-960
     30)*(A+2.D0)*(6.D0*A+9.D0)*(2.D0*A+7.D0)*U+(A+3.D0)*(20.D0*A+30.D0)CORA-961
     4*((A+16.D0)*F+(A+2.D0)*(6.D0*A+21.D0)*(A+1.D0)/W)/W)/W)/W         CORA-962
      C(2)=B3+2.D0*P*((2.D0*A-3.D0)*(240.D0*(3.D0*Q-(A-21.D0)*A+19.D0)-(CORA-963
     1A+1.D0)*((116.D0*A+166.D0)*A-432.D0)*U)*Q+A*(A+1.D0)*(2.D0*A+1.D0)CORA-964
     2*(16800.D0-(20.D0*(((2.D0*A+63.D0)*A+174.D0)*A+42.D0)-(A+2.D0)*(6.CORA-965
     3D0*A+9.D0)*(2.D0*A+5.D0)*(2.D0*A+7.D0)*U)*U)-(3.D0*A+6.D0)*((20.D0CORA-966
     4*(A-2.D0)*(2.D0*A-3.D0)*Q-(A+1.D0)*(2.D0*A+1.D0)*(840.D0-(2.D0*A+3CORA-967
     5.D0)*((A+30.D0)*A+84.D0)*U))*5.D0*F-(4.D0*A+10.D0)*(A+1.D0)*(30.D0CORA-968
     6*((2.D0*A+7.D0)*A+48.D0)*Q-(2.D0*A+1.D0)*(A*(840.D0+(3.D0*A+6.D0)*CORA-969
     7(2.D0*A+3.D0)*(2.D0*A+7.D0)*U)-(A+3.D0)*(20.D0*A+30.D0)*((7.D0*A+4CORA-970
     82.D0)*F+A*(A+2.D0)*(6.D0*A+21.D0)/W)/W))/W)/W)/W                  CORA-971
      C(3)=B2+6.D0*Q*((2.D0*A-3.D0)*(80.D0*(3.D0*Q-(A-21.D0)*A+19.D0)-((CORA-972
     1(12.D0*A-26.D0)*A-222.D0)*A-144.D0)*U)*Q-40.D0*A*((((4.D0*A-6.D0)*CORA-973
     2A-90.D0)*A-138.D0)*A-63.D0)*U+A*(A+1.D0)*(5600.D0*(2.D0*A+1.D0)-3.CORA-974
     3D0*(A+2.D0)*(2.D0*A+3.D0)*(2.D0*A+5.D0)*(2.D0*A+7.D0)*U**2)+10.D0*CORA-975
     4(A+2.D0)*((420.D0*(A+1.D0)*(2.D0*A+1.D0)-10.D0*(A-2.D0)*(2.D0*A-3.CORA-976
     5D0)*Q-(2.D0*A+3.D0)*(((10.D0*A+34.D0)*A+21.D0)*A+42.D0)*U)*F+(6.D0CORA-977
     6*A+15.D0)*((A+1.D0)*(((4.D0*A+14.D0)*A+96.D0)*Q+A*((A+2.D0)*(2.D0*CORA-978
     7A+3.D0)*(2.D0*A+7.D0)*U-56.D0*(2.D0*A+1.D0)))-14.D0*(A+3.D0)*(2.D0CORA-979
     7*A+3.D0)*((3.D0*A-2.D0)*F+A*(A+2.D0)*(2.D0*A+7.D0)*(A+1.D0)/W)/W)/CORA-980
     8W)/W)/W                                                           CORA-981
      C(4)=B4+6.D0*(((A-2.D0)*(2.D0*A-3.D0)*(40.D0*(9.D0*Q+(2.D0*A+73.D0CORA-982
     1)*A+57.D0)-(((38.D0*A-49.D0)*A-433.D0)*A-216.D0)*U)*Q+A*(A+1.D0)*(CORA-983
     22.D0*A+1.D0)*(40.D0*(((2.D0*A+63.D0)*A+69.D0)*A-168.D0)*U-33600.D0CORA-984
     3-(A+2.D0)*(2.D0*A+3.D0)*((19.D0*A-6.D0)*A-210.D0)*U**2))*F/(2.D0*ACORA-985
     4+3.D0)-(20.D0*(((((4.D0*A+132.D0)*A+35.D0)*A+1509.D0)*A+1812.D0)*ACORA-986
     5+420.D0+(A-2.D0)*(4.D0*A-6.D0)*(2.D0*A-5.D0)*Q)*Q-(A+1.D0)*(4.D0*(CORA-987
     6((((12.D0*A+80.D0)*A+625.D0)*A+1710.D0)*A+1158.D0)*A+630.D0)*Q*U+ACORA-988
     7**2*(2.D0*A+1.D0)*(16800.D0-(20.D0*(((2.D0*A+63.D0)*A+174.D0)*A+42CORA-989
     8.D0)-(3.D0*A+6.D0)*(2.D0*A+3.D0)*(2.D0*A+5.D0)*(2.D0*A+7.D0)*U)*U)CORA-990
     9)+(10.D0*A+20.D0)*((((((4.D0*A+132.D0)*A+371.D0)*A+1173.D0)*A+720.CORA-991
     AD0)*Q-A*(A+1.D0)*(2.D0*A+1.D0)*(840.D0-(2.D0*A+3.D0)*((A+30.D0)*A+CORA-992
     B84.D0)*U))*2.D0*F+(6.D0*A+15.D0)*(A+1.D0)*(((((8.D0*A+60.D0)*A+220CORA-993
     C.D0)*A-30.D0)*A+252.D0)*Q+A*(2.D0*A+1.D0)*(280.D0*A+(2.D0*A+3.D0)*CORA-994
     D(A*(A+2.D0)*(2.D0*A+7.D0)*U-14.D0*(A+3.D0)*((2.D0*A+12.D0)*F+A*(A+CORA-995
     E2.D0)*(2.D0*A+7.D0)/W)/W)))/W)/W)/W)/W                            CORA-996
      GO TO 45                                                          CORA-997
C MULTIPOLE LQ=6 L1=L2-4.                                               CORA-998
   44 BD=181440.D0*DEN*A**2*(A+2.D0)*(A+3.D0)*(A+4.D0)*(2.D0*A+5.D0)*(2.CORA-999
     1D0*A+7.D0)*S*DSQRT(S*Q*((A+1.D0)**2+F**2)*((A+2.D0)**2+F**2)*((A+3CORA1000
     2.D0)**2+F**2))                                                    CORA1001
      B1=-Q*((4320.D0*((7.D0*Q+(A+78.D0)*A+103.D0)*Q+((10.D0*A+255.D0)*ACORA1002
     1+638.D0)*A+420.D0)*Q+2880.D0*A*(A+2.D0)*(A+3.D0)*(A+4.D0)*(2.D0*A+CORA1003
     25.D0)*(2.D0*A+7.D0))/(A-1.D0)-72.D0*((5.D0*((14.D0*A-13.D0)*A+30.DCORA1004
     30)*Q+(((38.D0*A+833.D0)*A+243.D0)*A+678.D0)*A+2520.D0)*Q+A*(((((14CORA1005
     48.D0*A+2304.D0)*A+13879.D0)*A+31194.D0)*A+31336.D0)*A+13440.D0))*UCORA1006
     5+(A+2.D0)*(4.D0*A+10.D0)*(A*((((814.D0*A+8467.D0)*A+23265.D0)*A+30CORA1007
     6834.D0)*A+16632.D0)+(((338.D0*A+1757.D0)*A+8199.D0)*A+8190.D0)*Q-2CORA1008
     74.D0*A*(A+1.D0)*(A+3.D0)*(2.D0*A+3.D0)*(2.D0*A+5.D0)*(2.D0*A+7.D0)CORA1009
     8*U)*U**2)                                                         CORA1010
      B2=-6.D0*Q*F*(2160.D0*((4.D0*A+84.D0)*A+103.D0+7.D0*Q)*Q/(2.D0*A+3CORA1011
     1.D0)+1440.D0*(((2.D0*A+33.D0)*A+226.D0)*A+210.D0)-(A+2.D0)*(12.D0*CORA1012
     2(((62.D0*A+623.D0)*A+126.D0)*A-1260.D0)+180.D0*(4.D0*A-5.D0)*Q+(A+CORA1013
     33.D0)*(2.D0*A+5.D0)*((41.D0*A+591.D0)*A+910.D0)*U)*U)             CORA1014
      B3=-3.D0*F*(((30240.D0*Q+4320.D0*((4.D0*A+84.D0)*A+103.D0)-360.D0*CORA1015
     1(4.D0*A**2+5.D0)*(5.D0*A-6.D0)*U)*Q-24.D0*(((((284.D0*A+3480.D0)*ACORA1016
     2+25.D0)*A-690.D0)*A+396.D0)*A-7560.D0)*U)*Q/(2.D0*A+3.D0)+2880.D0*CORA1017
     3(((2.D0*A+33.D0)*A+226.D0)*A+210.D0)*Q+(2.D0*A+4.D0)*((((98.D0*A-2CORA1018
     4563.D0)*A-7346.D0)*A-19775.D0)*A-13650.D0)*Q*U**2+A*(A+1.D0)*(2.D0CORA1019
     5*A+1.D0)*(480.D0*(((2.D0*A+33.D0)*A+226.D0)*A+840.D0)+(A+2.D0)*((ACORA1020
     6+3.D0)*(2.D0*A+5.D0)*((233.D0*A+1798.D0)*A+2520.D0)*U-4.D0*(((262.CORA1021
     7D0*A+3523.D0)*A+16926.D0)*A+17640.D0))*U)*U)                      CORA1022
      B4=-(4320.D0*(2.D0*A-3.D0)*((7.D0*Q+(A+90.D0)*A+103.D0)*Q+7.D0*(((CORA1023
     12.D0*A+57.D0)*A+118.D0)*A+60.D0))*Q/((A+1.D0)*(2.D0*A+3.D0))-72.D0CORA1024
     2*(2.D0*A-3.D0)*((35.D0*A-50.D0)*Q+((19.D0*A+425.D0)*A-66.D0)*A-840CORA1025
     3.D0)*Q*U+(A+2.D0)*(4.D0*A+10.D0)*(((338.D0*A+1567.D0)*A+8439.D0)*ACORA1026
     4+8190.D0)*Q*U**2+A*(4.D0*A+2.D0)*(1440.D0*(((2.D0*A+33.D0)*A+226.DCORA1027
     50)*A+840.D0)-36.D0*((((74.D0*A+1109.D0)*A+6154.D0)*A+9304.D0)*A+33CORA1028
     660.D0)*U+(A+2.D0)*(2.D0*A+5.D0)*(((407.D0*A+3765.D0)*A+6714.D0)*A+CORA1029
     74536.D0-12.D0*(A+3.D0)*(2.D0*A+3.D0)*(2.D0*A+5.D0)*(2.D0*A+7.D0)*UCORA1030
     8)*U**2))                                                          CORA1031
      IF (LT) GO TO 45                                                  CORA1032
      C(1)=B1+6.D0*Q*(A+2.D0)*((240.D0*(9.D0*Q+(8.D0*A+98.D0)*A+141.D0)*CORA1033
     1Q+4.D0*(((((524.D0*A+7460.D0)*A+40945.D0)*A+88825.D0)*A+77406.D0)*CORA1034
     2A+17640.D0)*U-60.D0*(((6.D0*A-7.D0)*A+23.D0)*A-30.D0)*Q*U-(A+3.D0)CORA1035
     3*(2.D0*A+5.D0)*(480.D0*(A+4.D0)*(2.D0*A+7.D0)+(A+2.D0)*(2.D0*A+3.DCORA1036
     40)*((233.D0*A+1373.D0)*A+420.D0)*U**2))*F+4.D0*(2.D0*A+5.D0)*((A+2CORA1037
     5.D0)*(5.D0*(((116.D0*A+1028.D0)*A+3015.D0)*A+1953.D0)*Q-A*((((292.CORA1038
     6D0*A+2908.D0)*A+9639.D0)*A+5517.D0)*A-756.D0))*U-60.D0*((((14.D0*ACORA1039
     7+153.D0)*A+467.D0)*A+483.D0)+(2.D0*A-3.D0)*Q)*Q+(A+3.D0)*(12.D0*(2CORA1040
     8.D0*A+7.D0)*A*(5.D0*(A+4.D0)*(2.D0*A-3.D0)-(A+1.D0)*(A+2.D0)*(2.D0CORA1041
     9*A+3.D0)*(2.D0*A+5.D0)*U**2)+5.D0*((12.D0*(A+4.D0)*(A-2.D0)*(2.D0*CORA1042
     AA-3.D0)*(2.D0*A+7.D0)-12.D0*((22.D0*A+175.D0)*A+318.D0)*Q)*F+(A+2.CORA1043
     BD0)*((2.D0*A+3.D0)*(2.D0*A+5.D0)*((85.D0*A+409.D0)*A+84.D0)*U*F-6.CORA1044
     CD0*(2.D0*A+7.D0)*(15.D0*((20.D0*A+92.D0)*A+51.D0)*Q-(A+4.D0)*A*((3CORA1045
     D32.D0*A+368.D0)*A+141.D0)-(2.D0*A+3.D0)*(2.D0*A+5.D0)*(4.D0*A*(A+1CORA1046
     E.D0)*(A+3.D0)*U-7.D0*(A+4.D0)*((19.D0*A+3.D0)*F+8.D0*A*(A+1.D0)*(ACORA1047
     F+3.D0)/W)/W))/W))/W))/W)/W                                        CORA1048
      XX=(3.D0*A+6.D0)*(5.D0*F*(120.D0*((((8.D0*A+116.D0)*A+682.D0)*A+14CORA1049
     105.D0)*A+966.D0-(A-2.D0)*(2.D0*A-3.D0)*Q)*Q-(2.D0*A+3.D0)*(4.D0*((CORA1050
     2((140.D0*A+1808.D0)*A+7627.D0)*A+13909.D0)*A+6510.D0)*Q-A*(A+1.D0)CORA1051
     3*(2.D0*A+1.D0)*(12.D0*(((2.D0*A+33.D0)*A+226.D0)*A+840.D0)+(A+2.D0CORA1052
     4)*(A+3.D0)*(10.D0*A+25.D0)*(17.D0*A+84.D0)*U))*U)+(4.D0*A+10.D0)*(CORA1053
     512.D0*(A+2.D0)*(30.D0*(((8.D0*A+62.D0)*A+159.D0)*Q-((((4.D0*A+40.DCORA1054
     60)*A+151.D0)*A+94.D0)*A+84.D0))*Q-(A+1.D0)*(2.D0*A+3.D0)*(5.D0*(((CORA1055
     768.D0*A+564.D0)*A+1285.D0)*A+210.D0)*Q-A**2*(2.D0*A+1.D0)*(((166.DCORA1056
     80*A+1399.D0)*A+3108.D0)+(A+3.D0)*(4.D0*A+10.D0)*(2.D0*A+7.D0)*U))*CORA1057
     9U)+10.D0*(A+3.D0)*(42.D0*(2.D0*A+3.D0)*(((70.D0*A+561.D0)*A+1241.DCORA1058
     A0)*A+510.D0)*Q*F-(2.D0*A+3.D0)*(A+1.D0)*(7.D0*A*(A+2.D0)*(2.D0*A+1CORA1059
     B.D0)*(2.D0*A+5.D0)*(19.D0*A+84.D0)*U*F+(2.D0*A+7.D0)*(1008.D0*A*(ACORA1060
     C+4.D0)*(2.D0*A+1.D0)*F-(6.D0*A+12.D0)*(4.D0*((194.D0*A+797.D0)*A+1CORA1061
     D05.D0)*Q-A*(2.D0*A+1.D0)*(432.D0*A*(A+4.D0)+(2.D0*A+5.D0)*(4.D0*A*CORA1062
     E(A+3.D0)*U-7.D0*(A+4.D0)*(27.D0*F+8.D0*A*(A+3.D0)/W)/W)))/W)))/W)/CORA1063
     FW)/W                                                              CORA1064
      C(2)=B3+2.D0*(1440.D0*(((((8.D0*A+132.D0)*A+906.D0)*A+2133.D0)*A+2CORA1065
     1116.D0)*A+840.D0+(2.D0*A-3.D0)*((A+36.D0)*A+47.D0+3.D0*Q)*Q)*Q-24.CORA1066
     2D0*((((((504.D0*A+7644.D0)*A+39370.D0)*A+119625.D0)*A+174926.D0)*ACORA1067
     3+106596.D0)*A+17640.D0+(10.D0*A-15.D0)*(((8.D0*A-6.D0)*A+13.D0)*A-CORA1068
     430.D0)*Q)*Q*U+(A+1.D0)*(8.D0*A+12.D0)*((A+2.D0)*(2.D0*A+5.D0)*(((3CORA1069
     562.D0*A+3013.D0)*A+8436.D0)*A+1890.D0)*Q*U+A**2*(2.D0*A+1.D0)*(60.CORA1070
     6D0*(((2.D0*A+33.D0)*A+226.D0)*A+840.D0)-(A+2.D0)*(2.D0*A+5.D0)*((7CORA1071
     73.D0*A+822.D0)*A+3024.D0+(A+3.D0)*(12.D0*A+30.D0)*(2.D0*A+7.D0)*U)CORA1072
     8*U))*U+XX)/W                                                      CORA1073
      C(3)=B2+24.D0*Q*(2.D0*(60.D0*(((((8.D0*A+132.D0)*A+906.D0)*A+2133.CORA1074
     1D0)*A+2116.D0)*A+840.D0+(2.D0*A-3.D0)*(3.D0*Q+(A+36.D0)*A+47.D0)*QCORA1075
     2)-(A+2.D0)*(2.D0*A+5.D0)*((((112.D0*A+1218.D0)*A+3459.D0)*A+4242.DCORA1076
     30)*A+1764.D0+5.D0*(2.D0*A-3.D0)**2*Q-3.D0*(A+1.D0)*(A+3.D0)*(2.D0*CORA1077
     4A+3.D0)*(2.D0*A+5.D0)*(2.D0*A+7.D0)*U)*U)+5.D0*(A+2.D0)*(30.D0*(((CORA1078
     5((8.D0*A+116.D0)*A+682.D0)*A+1405.D0)*A+966.D0)-(A-2.D0)*(2.D0*A-3CORA1079
     6.D0)*Q)*F-(A+2.D0)*(2.D0*A+5.D0)*((A+3.D0)*(2.D0*A+3.D0)*(37.D0*A+CORA1080
     7217.D0)*U*F+6.D0*(6.D0*((((4.D0*A+40.D0)*A+151.D0)*A+94.D0)*A+84.DCORA1081
     80)-2.D0*((24.D0*A+186.D0)*A+477.D0)*Q+(2.D0*A+3.D0)*(2.D0*A+5.D0)*CORA1082
     9(2.D0*(A+1.D0)*(A+3.D0)*(2.D0*A+7.D0)*U-7.D0*(A+3.D0)*((11.D0*A+51CORA1083
     A.D0)*F+4.D0*(2.D0*A+7.D0)*(A+1.D0)*(A+3.D0)/W)/W))/W))/W)/W       CORA1084
      C(4)=B4+6.D0*((240.D0*(A-2.D0)*(2.D0*A-3.D0)*(9.D0*Q+(8.D0*A+118.DCORA1085
     10)*A+141.D0)*Q/(2.D0*A+3.D0)-60.D0*(A**2-4.D0)*(2.D0*A-3.D0)*(3.D0CORA1086
     2*A-5.D0)*Q*U-(A+1.D0)*(2.D0*A+1.D0)*(480.D0*(((2.D0*A+33.D0)*A+226CORA1087
     3.D0)*A+840.D0)+(A+2.D0)*((A+3.D0)*(2.D0*A+5.D0)*((233.D0*A+1798.D0CORA1088
     4)*A+2520.D0)*U**2-4.D0*(((262.D0*A+3523.D0)*A+16926.D0)*A+17640.D0CORA1089
     5)*U)))*F-4.D0*(60.D0*(7.D0*(((((4.D0*A+60.D0)*A+291.D0)*A+969.D0)*CORA1090
     6A+1409.D0)*A+690.D0)+(A-2.D0)*(2.D0*A-3.D0)*(2.D0*A-5.D0)*Q)*Q-(A+CORA1091
     71.D0)*(2.D0*A+3.D0)*(60.D0*A*(2.D0*A+1.D0)*(((2.D0*A+33.D0)*A+226.CORA1092
     8D0)*A+840.D0)+(A+2.D0)*(2.D0*A+5.D0)*(5.D0*((58.D0*A+473.D0)*A+130CORA1093
     92.D0)*Q-A*(2.D0*A+1.D0)*((73.D0*A+822.D0)*A+3024.D0+6.D0*(A+3.D0)*CORA1094
     A(2.D0*A+5.D0)*(2.D0*A+7.D0)*U))*U)+5.D0*(A+2.D0)*(2.D0*A+3.D0)*(12CORA1095
     B.D0*(((22.D0*A+273.D0)*A+1001.D0)*A+1590.D0)*Q*F-(A+1.D0)*(12.D0*(CORA1096
     C2.D0*A+1.D0)*(((2.D0*A+33.D0)*A+226.D0)*A+840.D0)*F+(A+2.D0)*(2.D0CORA1097
     D*A+5.D0)*(5.D0*(A+3.D0)*(2.D0*A+1.D0)*(17.D0*A+84.D0)*U*F-6.D0*(15CORA1098
     E.D0*((20.D0*A+164.D0)*A+357.D0)*Q-(2.D0*A+1.D0)*(A*((166.D0*A+1399CORA1099
     F.D0)*A+3108.D0)+(A+3.D0)*(2.D0*A+5.D0)*(A*(4.D0*A+14.D0)*U-7.D0*((CORA1100
     G19.D0*A+84.D0)*F+4.D0*A*(2.D0*A+7.D0)*(A+3.D0)/W)/W)))/W)))/W)/W)/CORA1101
     HW                                                                 CORA1102
   45 A1=-2.D0*A**2*Y*B1*Z                                              CORA1103
      A2=Y*((2.D0*A+1.D0)*((2.D0*(A**2-A)**2*(1.D0+S**2)+Y**2*(A**2-A-.5CORA1104
     1D0))*B1-Y*(P*B2+Q*B3))+(4.D0*A-2.D0)*(A+1.D0)**2*P*Q*B4)*Z        CORA1105
      A3=-Y*((2.D0*A-1.D0)*((2.D0*(A**2+A)**2*(1.D0+S**2)+Y**2*(A**2+A-.CORA1106
     15D0))*B4-Y*(B2+B3))+(4.D0*A+2.D0)*(A-1.D0)**2*B1)*Z               CORA1107
      A4=2.D0*A**2*Y*B4*Z                                               CORA1108
      IF (X.EQ.0.D0) GO TO 48                                           CORA1109
      A1=A1+X*2.D0*P*((A*X-2.D0*F)*B2-2.D0*Q*(A+1.D0)*B4)               CORA1110
      A2=A2+X*P*((2.D0*Q*((2.D0*A+3.D0)*(A+1.D0)*F+((2.D0*A+3.D0)*A-4.D0CORA1111
     1)*E*S)-2.D0*P*(2.D0*A-3.D0)*(2.D0*A+1.D0)*F+(2.D0*A+1.D0)*((((2.D0CORA1112
     2*A-1.D0)*A-8.D0)*A+5.D0)*F-S*((A+1.D0)*((6.D0*A-9.D0)*A+1.D0)*E-S*CORA1113
     3((6.D0*A-7.D0)*(A**2-1.D0)*F-S*((((2.D0*A-5.D0)*A+2.D0)*A-1.D0)*E)CORA1114
     4))))*B2+2.D0*Q*(A+1.D0)*((2.D0*A+3.D0)*(Q+P)+(2.D0*A-1.D0)*(4.D0*ECORA1115
     5*F*S+(1.D0+S**2)*((2.D0*A-4.D0)*A-3.D0)))*B4)                     CORA1116
      A3=A3+X*((2.D0*P*((2.D0*A-3.D0)*(A-1.D0)*E*S+F*((2.D0*A-3.D0)*A-4.CORA1117
     1D0))-2.D0*Q*E*S*(2.D0*A+3.D0)*(2.D0*A-1.D0)-(2.D0*A-1.D0)*((((2.D0CORA1118
     2*A+5.D0)*A+2.D0)*A+1.D0)*F-S*((6.D0*A+7.D0)*(A**2-1.D0)*E-S*((A-1.CORA1119
     3D0)*((6.D0*A+9.D0)*A+1.D0)*F-S*((((2.D0*A+1.D0)*A-8.D0)*A-5.D0)*E)CORA1120
     4))))*B2+(-4.D0*P*Q*((2.D0*A+1.D0)*A+2.D0)-(2.D0*A-1.D0)*(2.D0*S*E*CORA1121
     5F*(2.D0*A+3.D0)*(Q+P+S**2+1.D0)-A**2*(((A+3.D0)*A+2.5D0)*(1.D0+S**CORA1122
     64)-S**2*((2.D0*A-2.D0)*A-7.D0))+A*(A-1.D0)*(2.D0*A+3.D0)*(1.D0+S**CORA1123
     72)*(P+Q)-(S**2-1.D0)*((2.D0*A+4.D0)*A+.5D0)*(Q-P)))*B4)           CORA1124
      A4=A4+X*2.D0*(((A-1.D0)*(P+Q)+Y**2)*B4+(2.D0*E*S-A*X)*B2)         CORA1125
   46 A4=A4+X*2.D0*((2.D0*F+A*X)*B3-2.D0*(A-1.D0)*B1)                   CORA1126
   47 A1=A1+X*2.D0*(((A+1.D0)*(Q+P)-Y**2)*B1-Q*(A*X+2.D0*E*S)*B3)       CORA1127
      A2=A2+X*((-4.D0*P*Q*(2.D0*A**2-A+2.D0)-(2.D0*A+1.D0)*(2.D0*S*E*F*(CORA1128
     12.D0*A-3.D0)*(Q+P+S**2+1.D0)+A**2*(((A-3.D0)*A+2.5D0)*(1.D0+S**4)-CORA1129
     2S**2*((2.D0*A+2.D0)*A-7.D0))+A*(A+1.D0)*(2.D0*A-3.D0)*(1.D0+S**2)*CORA1130
     3(P+Q)+(S**2-1.D0)*((2.D0*A-4.D0)*A+.5D0)*(Q-P)))*B1+Q*(2.D0*P*((2.CORA1131
     4D0*A+3.D0)*(A+1.D0)*E*S+F*((2.D0*A+3.D0)*A-4.D0))-2.D0*Q*E*S*(2.D0CORA1132
     5*A-3.D0)*(2.D0*A+1.D0)-(2.D0*A+1.D0)*((((2.D0*A-5.D0)*A+2.D0)*A-1.CORA1133
     6D0)*F-S*((6.D0*A-7.D0)*(A**2-1.D0)*E-S*((A+1.D0)*((6.D0*A-9.D0)*A+CORA1134
     71.D0)*F-S*((((2.D0*A-1.D0)*A-8.D0)*A+5.D0)*E)))))*B3)             CORA1135
      A3=A3+X*(2.D0*(A-1.D0)*((2.D0*A-3.D0)*(Q+P)+(2.D0*A+1.D0)*(4.D0*E*CORA1136
     1F*S+(1.D0+S**2)*((2.D0*A+4.D0)*A-3.D0)))*B1+(2.D0*Q*((2.D0*A-3.D0)CORA1137
     2*(A-1.D0)*F+((2.D0*A-3.D0)*A-4.D0)*E*S)-2.D0*P*(2.D0*A+3.D0)*(2.D0CORA1138
     3*A-1.D0)*F+(2.D0*A-1.D0)*((((2.D0*A+1.D0)*A-8.D0)*A-5.D0)*F-S*((A-CORA1139
     41.D0)*((6.D0*A+9.D0)*A+1.D0)*E-S*((6.D0*A+7.D0)*(A**2-1.D0)*F-S*((CORA1140
     5((2.D0*A+5.D0)*A+2.D0)*A+1.D0)*E)))))*B3)                         CORA1141
   48 B(1)=A1*DSQRT(((A-1.D0)**2+E**2)*((A-1.D0)**2+F**2))*(2.D0*A-3.D0)CORA1142
     1/(BD*A*(2.D0*A-1.D0))                                             CORA1143
      B(2)=A2/(A*S*BD*(2.D0*A+1.D0))                                    CORA1144
      B(3)=A3/(A*S*BD*(2.D0*A-1.D0))*DSQRT(P*Q)                         CORA1145
      B(4)=A4*DSQRT(P*Q*((A+1.D0)**2+E**2)*((A+1.D0)**2+F**2))*(2.D0*A+3CORA1146
     1.D0)/(BD*A*(2.D0*A+1.D0))                                         CORA1147
      IF (LT) RETURN                                                    CORA1148
      Z=-2.D0*BD*V*W*DSQRT(S)/DEN                                       CORA1149
      BD=A**2*BD*DSQRT(V*W)**3                                          CORA1150
      A1=A1*(A-1.D0)**2                                                 CORA1151
      A4=A4*(A+1.D0)**2                                                 CORA1152
      C(1)=C(1)/Z-(A+Y*W*(.5D0*A-1.D0)/(A-1.D0)**2)*(2.D0*A-1.D0)*A1/BD CORA1153
      C(2)=(C(2)/Z+W*(A1-P*A4)/BD)*DSQRT(Q)                             CORA1154
      C(3)=(C(3)/Z+W*(A1-Q*A4)/BD)*DSQRT(P)                             CORA1155
      C(4)=(C(4)/Z+(A+Y*W*(.5D0*A+1.D0)/(A+1.D0)**2)*(2.D0*A+1.D0)*A4/BDCORA1156
     1)*DSQRT(P*Q)                                                      CORA1157
      GO TO 70                                                          CORA1158
C FOR NEUTRONS ( COULOMB PARAMETER 0.).                                 CORA1159
   49 X=W**2                                                            CORA1160
      GO TO ( 50 , 51 , 54 , 57 , 61 , 65 ) , LQ                        CORA1161
   50 IF ((LL.EQ.1).OR.(U.EQ.0)) GO TO 73                               CORA1162
C MULTIPOLE LQ=1 L1=L2-1 (USELESS).                                     CORA1163
      B(2)=(U+2.D0-2.D0*A)/(DSQRT(S)*U)                                 CORA1164
      B(3)=DSQRT(S)*(2.D0*A+2.D0)/U                                     CORA1165
      IF (LT) RETURN                                                    CORA1166
      C(1)=DSQRT(S)/(U*X)                                               CORA1167
      C(4)=C(1)/S                                                       CORA1168
      RETURN                                                            CORA1169
   51 GO TO ( 52 , 73 , 53 ) , LL                                       CORA1170
C MULTIPOLE LQ=2 L1=L2.                                                 CORA1171
   52 B(2)=1.D0                                                         CORA1172
      RETURN                                                            CORA1173
C MULTIPOLE LQ=2 L1=L2-2.                                               CORA1174
   53 B(2)=(2.D0*A-1.D0)/(3.D0*S)                                       CORA1175
      B(3)=-(A+1.D0)/1.5D0                                              CORA1176
      IF (.NOT.LT) C(2)=-B(2)/(X*V)                                     CORA1177
      RETURN                                                            CORA1178
   54 GO TO ( 73 , 55 , 73 , 56 ) , LL                                  CORA1179
C MULTIPOLE LQ=3 L1=L2-1.                                               CORA1180
   55 B(2)=1.D0/(3.D0*DSQRT(S))                                         CORA1181
      B(3)=-S*B(2)                                                      CORA1182
      IF (.NOT.LT) C(2)=-B(2)/(X*V)                                     CORA1183
      GO TO 70                                                          CORA1184
C MULTIPOLE LQ=3 L1=L2-3.                                               CORA1185
   56 Y=15.D0*S*DSQRT(S)                                                CORA1186
      B(2)=(1.D0-2.D0*U*(A-1.D0))/Y                                     CORA1187
      B(3)=((2.D0*A+1.D0)*U-1.D0)*S/Y                                   CORA1188
      IF (LT) RETURN                                                    CORA1189
      C(1)=(6.D0*A-3.D0)/(Y*X**2*S)                                     CORA1190
      C(2)=(2.D0+U*(2.D0*A+1.D0)-(12.D0*A**2-3.D0)/X)/(Y*X*V)           CORA1191
      C(3)=-3.D0/(Y*X*W)                                                CORA1192
      C(4)=(6.D0*A+3.D0)/(Y*X**2)                                       CORA1193
      GO TO 70                                                          CORA1194
   57 GO TO ( 58 , 73 , 59 , 73 , 60 ) , LL                             CORA1195
C MULTIPOLE LQ=4 L1=L2.                                                 CORA1196
   58 B(2)=1.D0/(3.D0*A+6.D0)                                           CORA1197
      B(3)=-(U+2.D0)*B(2)/(2.D0*S)                                      CORA1198
      IF (LT) RETURN                                                    CORA1199
      C(2)=-B(2)/(2.D0*X*V)                                             CORA1200
      C(3)=C(2)/S                                                       CORA1201
      C(4)=3.D0*C(3)/W                                                  CORA1202
      GO TO 70                                                          CORA1203
C MULTIPOLE LQ=4 L1=L2-2.                                               CORA1204
   59 Y=15.D0*S*(A+2.D0)                                                CORA1205
      B(2)=-(U*(A-1.D0)-3.D0)/Y                                         CORA1206
      B(3)=((A+.5D0)*U-3.D0)*S/Y                                        CORA1207
      IF (LT) RETURN                                                    CORA1208
      C(1)=(3.D0*A+6.D0)/(Y*X**2*S)                                     CORA1209
      C(2)=-(3.D0-(2.D0*A+1.D0)*(U-(6.D0*A+12.D0)/X))/(2.D0*Y*X*V)      CORA1210
      C(3)=-1.5D0/(Y*X*W)                                               CORA1211
      C(4)=-(2.D0*A+1.D0)*C(3)/W                                        CORA1212
      GO TO 70                                                          CORA1213
C MULTIPOLE LQ=4 L1=L2-4.                                               CORA1214
   60 Y=210.D0*(A+2.D0)*S                                               CORA1215
      B(2)=(6.D0-2.D0*(A-1.D0)*U*(2.D0-(2.D0*A+3.D0)*U))/(Y*S)          CORA1216
      B(3)=-(6.D0-U*(4.D0*A-1.D0-U*(2.D0*A+3.D0)*(2.D0*A+1.D0)))/Y      CORA1217
      IF (LT) RETURN                                                    CORA1218
      C1=U-(10.D0*A-5.D0)/X                                             CORA1219
      C(1)=-(6.D0*A+12.D0)*(4.D0+(2.D0*A+3.D0)*C1)/(Y*X**2*S**2)        CORA1220
      C(2)=-(3.D0+(2.D0*A+10.D0)*U-(A+2.D0)*(108.D0*A-6.D0)/X+(2.D0*A+3.CORA1221
     1D0)*(2.D0*A+1.D0)*(U**2-(6.D0*A+12.D0)*C1/X))/(Y*X*V*S)           CORA1222
      C(3)=-(3.D0-(6.D0*A+9.D0)*(U-(10.D0*A+20.D0)/X))/(Y*X*V)          CORA1223
      C(4)=-(2.D0*A+1.D0)*C(3)/W-(30.D0*A+60.D0)/(Y*X**2*S)             CORA1224
      GO TO 70                                                          CORA1225
   61 GO TO ( 73 , 62 , 73 , 63 , 73 , 64 ) , LL                        CORA1226
C MULTIPOLE LQ=5 L1=L2-1.                                               CORA1227
   62 Y=60.D0*(A+2.D0)*(A+3.D0)*DSQRT(S)                                CORA1228
      B(2)=(16.D0-U*(2.D0*A-2.D0))/Y                                    CORA1229
      B(3)=-(16.D0-U*(2.D0*A-10.D0+U*(2.D0*A+1.D0)))/(Y*S)              CORA1230
      IF (LT) RETURN                                                    CORA1231
      C(1)=6.D0*(A+2.D0)/(S*Y*X**2)                                     CORA1232
      C(2)=-(8.D0-(2.D0*A+1.D0)*(U-6.D0*(A+2.D0)/X))/(Y*X*V)            CORA1233
      C(3)=-(8.D0+3.D0*U-30.D0*(A+2.D0)/X)/(Y*X*V*S)                    CORA1234
      C(4)=(6.D0*A-12.D0+(15.D0*U-(60.D0*A+30.D0)*(A+2.D0)/X)/S**2)/(Y*XCORA1235
     1**2)                                                              CORA1236
      GO TO 70                                                          CORA1237
C MULTIPOLE LQ=5 L1=L2-3.                                               CORA1238
   63 Y=420.D0*(A+2.D0)*(A+3.D0)*DSQRT(S)                               CORA1239
      B(2)=(48.D0-U*(2.D0*A-2.D0)*(9.D0-U*(2.D0*A+3.D0)))/(Y*S)         CORA1240
      B(3)=-(48.D0-U*(18.D0*A+6.D0-U*(2.D0*A+1.D0)*(2.D0*A+3.D0)))/Y    CORA1241
      IF (LT) RETURN                                                    CORA1242
      C1=U-(10.D0*A+30.D0)/X                                            CORA1243
      C(1)=(6.D0*A+12.D0)*(3.D0-(2.D0*A+3.D0)*C1)/(Y*X**2*S**2)         CORA1244
      C(2)=((24.D0*A+162.D0)*(A+2.D0)/X-24.D0+(12.D0*A-3.D0)*U-(2.D0*A+1CORA1245
     1.D0)*(2.D0*A+3.D0)*(U**2-(6.D0*A+12.D0)*C1/X))/(Y*X*S*V)          CORA1246
      C(3)=-(24.D0-(6.D0*A+9.D0)*(U-(10.D0*A+20.D0)/X))/(Y*X*V)         CORA1247
      C(4)=-(2.D0*A+1.D0)*C(3)/W-(30.D0*A+60.D0)/(Y*X**2*S)             CORA1248
      GO TO 70                                                          CORA1249
C MULTIPOLE LQ=5 L1=L2-5.                                               CORA1250
   64 Y=3780.D0*(A+2.D0)*(A+3.D0)*S*DSQRT(S)                            CORA1251
      B(2)=(48.D0-U*(2.D0*A-2.D0)*(15.D0-U*(6.D0*A+6.D0-U*(2.D0*A+3.D0)*CORA1252
     1(2.D0*A+5.D0))))/(Y*S)                                            CORA1253
      B(3)=-(48.D0-U*(30.D0*A-6.D0-U*(12.D0*A**2+12.D0*A-9.D0-U*(8.D0*A*CORA1254
     1*3+36.D0*A**2+46.D0*A+15.D0))))/Y                                 CORA1255
      IF (LT) RETURN                                                    CORA1256
      C1=U**2-(10.D0*A+30.D0)*(U-(14.D0*A-7.D0)/X)/X                    CORA1257
      C(1)=(6.D0*A+12.D0)*(3.D0+(6.D0*A+24.D0)*U-(A+3.D0)*(260.D0*A+110.CORA1258
     1D0)/X+(2.D0*A+5.D0)*(2.D0*A+3.D0)*C1)/(Y*X**2*S**2)               CORA1259
      C(2)=-(24.D0-(18.D0*A-9.D0)*U+(6.D0*A+12.D0)*(66.D0*A+183.D0)/X-(2CORA1260
     1.D0*A+3.D0)*(18.D0*U**2-(12.D0*A+24.D0)*((8.D0*A+29.D0)*U-(A+3.D0)CORA1261
     2*(200.D0*A-40.D0)/X)/X+(2.D0*A+5.D0)*(2.D0*A+1.D0)*(U**3-(6.D0*A+1CORA1262
     32.D0)*C1/X)))/(Y*X*V*S)                                           CORA1263
      C1=U**2-(10.D0*A+20.D0)*(U-(14.D0*A+42.D0)/X)/X                   CORA1264
      C(3)=-(24.D0-(12.D0*A+3.D0)*U+(6.D0*A+9.D0)*(2.D0*A+5.D0)*C1-(A+2.CORA1265
     1D0)*(360.D0*A+1170.D0)/X)/(Y*X*V)                                 CORA1266
      C(4)=(18.D0*A-36.D0+(6.D0*A+9.D0)*((6.D0*A+33.D0)*U+(2.D0*A+1.D0)*CORA1267
     1(2.D0*A+5.D0)*C1-(A+2.D0)*(260.D0*A+830.D0)/X))/(Y*X**2*S)        CORA1268
      GO TO 70                                                          CORA1269
   65 GO TO ( 66 , 73 , 67 , 73 , 68 , 73 , 69 ) , LL                   CORA1270
C MULTIPOLE LQ=6 L1=L2.                                                 CORA1271
   66 Y=120.D0*(A+2.D0)*(A+3.D0)*(A+4.D0)*S**2                          CORA1272
      B(2)=(32.D0+U*(32.D0-U*(2.D0*A-2.D0)))/Y                          CORA1273
      B(3)=-(32.D0+U*(48.D0-U*(4.D0*A-14.D0+U*(2.D0*A+1.D0))))/(Y*S)    CORA1274
      IF (LT) RETURN                                                    CORA1275
      C(1)=(A+2.D0)*(12.D0+6.D0*U-(60.D0*A+180.D0)/X)/(Y*X**2*S)        CORA1276
      C(2)=-(16.D0-U*(2.D0*A-12.D0+U*(2.D0*A+1.D0))+(6.D0*A+12.D0)*(2.D0CORA1277
     1*A-4.D0+(2.D0*A+1.D0)*(U-(10.D0*A+30.D0)/X))/X)/(Y*X*V)           CORA1278
      C(3)=-(16.D0+U*(2.D0*A+20.D0+3.D0*U)-(6.D0*A+12.D0)*(4.D0-2.D0*A+5CORA1279
     1.D0*U+(20.D0*A+10.D0)*(A+3.D0)/X)/X)/(Y*X*S*V)                    CORA1280
      C(4)=((12.D0*A-24.D0)*S**2+(6.D0*A+3.D0)*(U**2-(10.D0*A+20.D0)*(2.CORA1281
     1D0+U+(4.D0*A+2.D0)*(A+3.D0)/X)/X))/(Y*X**2*S**2)                  CORA1282
      RETURN                                                            CORA1283
C MULTIPOLE LQ=6 L1=L2-2.                                               CORA1284
   67 Y=840.D0*(A+2.D0)*(A+3.D0)*(A+4.D0)*S                             CORA1285
      B(2)=(160.D0-U*(2.D0*A-2.D0)*(16.D0-U*(2.D0*A+3.D0)))/Y           CORA1286
      B(3)=-(160.D0-U*(32.D0*A-112.D0-U*((4.D0*A-24.D0)*A-10.D0+U*(2.D0*CORA1287
     1A+3.D0)*(2.D0*A+1.D0))))/(Y*S)                                    CORA1288
      IF (LT) RETURN                                                    CORA1289
      C1=U-(10.D0*A+30.D0)/X                                            CORA1290
      C(1)=(60.D0-(12.D0*A+18.D0)*C1)*(A+2.D0)/(Y*X**2*S)               CORA1291
      C(2)=-(80.D0-(26.D0*A+4.D0)*U+(60.D0*A**2-240.D0)/X+(2.D0*A+1.D0)*CORA1292
     1(2.D0*A+3.D0)*(U**2-(6.D0*A+12.D0)*C1/X))/(Y*X*V)                 CORA1293
      C(3)=-(80.D0-U*(6.D0*A-36.D0)+60.D0*(A**2-4.D0)/X-(6.D0*A+9.D0)*(UCORA1294
     1**2-(10.D0*A+20.D0)*(U-(14.D0*A+42.D0)/X)/X))/(Y*X*V*S)           CORA1295
      C(4)=-(105.D0*(1.D0-(10.D0*A+20.D0)/X)-(60.D0*A-15.D0)*S**2+(6.D0*CORA1296
     1A+3.D0)*(2.D0*A+3.D0)*(U*S**2-(10.D0*A+20.D0)*(1.D0+U-(14.D0*A+42.CORA1297
     2D0)/X)/X))/(Y*X**2*S**2)                                          CORA1298
      GO TO 70                                                          CORA1299
C MULTIPOLE LQ=6 L1=L2-4.                                               CORA1300
   68 Y=7560.D0*(A+2.D0)*(A+3.D0)*(A+4.D0)*S                            CORA1301
      B(2)=(480.D0-U*(2.D0*A-2.D0)*(96.D0-U*(24.D0*A+33.D0-U*(2.D0*A+3.DCORA1302
     10)*(2.D0*A+5.D0))))/(Y*S)                                         CORA1303
      B(3)=-(480.D0-U*(192.D0*A+48.D0-U*(2.D0*A+3.D0)*(24.D0*A+6.D0-U*(2CORA1304
     1.D0*A+1.D0)*(2.D0*A+5.D0))))/Y                                    CORA1305
      IF (LT) RETURN                                                    CORA1306
      C1=(U**2-(10.D0*A+30.D0)*(U-(14.D0*A+56.D0)/X)/X)                 CORA1307
      C(1)=(6.D0*A+12.D0)*(30.D0-(12.D0*A+3.D0)*U-(A+3.D0)*(80.D0*A+470.CORA1308
     1D0)/X+(2.D0*A+3.D0)*(2.D0*A+5.D0)*C1)/(Y*X**2*S**2)               CORA1309
      C(2)=-(240.D0-U*(126.D0*A-36.D0)+180.D0*(A**2-4.D0)/X+(2.D0*A+3.D0CORA1310
     1)*((18.D0*A-9.D0)*U**2-(6.D0*A+12.D0)*((2.D0*A-49.D0)*U+(10.D0*A+3CORA1311
     20.D0)*(22.D0*A+109.D0)/X)/X-(2.D0*A+1.D0)*(2.D0*A+5.D0)*(U**3-(6.DCORA1312
     30*A+12.D0)*C1/X)))/(Y*X*W*S**2)                                   CORA1313
      C1=(2.D0*A+5.D0)*(U**2-(10.D0*A+20.D0)*(U-(14.D0*A+42.D0)/X)/X)   CORA1314
      C(3)=-(240.D0-(66.D0*A+84.D0)*U+180.D0*(A**2-4.D0)/X+(6.D0*A+9.D0)CORA1315
     1*C1)/(Y*X*V)                                                      CORA1316
      C(4)=((A-2.D0)*(180.D0-(72.D0*A+108.D0)*U)+(6.D0*A+9.D0)*((2.D0*A+CORA1317
     11.D0)*C1-(A+2.D0)*(80.D0*A+740.D0)/X))/(Y*X**2*S)                 CORA1318
      GO TO 70                                                          CORA1319
C MULTIPOLE LQ=6 L1=L2-6.                                               CORA1320
   69 Y=83160.D0*(A+2.D0)*(A+3.D0)*(A+4.D0)*S**2                        CORA1321
      B(2)=(480.D0-U*(2.D0*A-2.D0)*(144.D0-U*(54.D0*A+57.D0-U*(2.D0*A+5.CORA1322
     1D0)*(8.D0*A+6.D0-U*(2.D0*A+7.D0)*(2.D0*A+3.D0)))))/(Y*S)          CORA1323
      B(3)=-(480.D0-U*(288.D0*A-48.D0-U*(108.D0*A**2+120.D0*A-78.D0-U*(2CORA1324
     1.D0*A+3.D0)*(2.D0*A+5.D0)*(8.D0*A-5.D0-U*(2.D0*A+1.D0)*(2.D0*A+7.DCORA1325
     20)))))/Y                                                          CORA1326
      IF (LT) RETURN                                                    CORA1327
      C1=U**3-(10.D0*A+30.D0)*(U**2-(14.D0*A+56.D0)*(U-(18.D0*A-9.D0)/X)CORA1328
     1/X)/X                                                             CORA1329
      C(1)=(6.D0*A+12.D0)*(30.D0-(18.D0*A-3.D0)*U+(10.D0*A+30.D0)*(118.DCORA1330
     10*A+457.D0)/X-(2.D0*A+5.D0)*((4.D0*A+36.D0)*U**2-(10.D0*A+30.D0)*(CORA1331
     2(24.D0*A+106.D0)*U-(14.D0*A+56.D0)*(52.D0*A+6.D0)/X)/X+(2.D0*A+7.DCORA1332
     30)*(2.D0*A+3.D0)*C1))/(Y*X**2*S**2)                               CORA1333
      C(2)=-(240.D0-(174.D0*A-84.D0)*U+((72.D0*A+36.D0)*A-63.D0)*U**2+(1CORA1334
     180.D0*(A**2-4.D0)-(6.D0*A+12.D0)*(((136.D0*A+772.D0)*A+1377.D0)*U-CORA1335
     2(10.D0*A+30.D0)*((712.D0*A+3370.D0)*A+2193.D0)/X))/X-(2.D0*A+5.D0)CORA1336
     3*(2.D0*A+3.D0)*((2.D0*A-26.D0)*U**3+(6.D0*A+12.D0)*((14.D0*A+82.D0CORA1337
     4)*U**2-(10.D0*A+30.D0)*((38.D0*A+166.D0)*U-(14.D0*A+56.D0)*(70.D0*CORA1338
     5A-19.D0)/X)/X)/X-(2.D0*A+7.D0)*(2.D0*A+1.D0)*(U**4-6.D0*(A+2.D0)*CCORA1339
     61/X)))/(Y*X*V*S)                                                  CORA1340
      C1=U**3-(10.D0*A+20.D0)*(U**2-(14.D0*A+42.D0)*(U-(18.D0*A+72.D0)/XCORA1341
     1)/X)/X                                                            CORA1342
      C(3)=-(240.D0-(114.D0*A+36.D0)*U+180*(A**2-4.D0)/X+(6.D0*A+15.D0)*CORA1343
     1((6.D0*A-1.D0)*U**2+(10.D0*A+20.D0)*((10.D0*A+57.D0)*U-(14.D0*A+42CORA1344
     2.D0)*(34.D0*A+141.D0)/X)/X-(2.D0*A+7.D0)*(2.D0*A+3.D0)*C1))/(Y*X*VCORA1345
     3)                                                                 CORA1346
      C(4)=(180.D0*A-360.D0-(108.D0*A+72.D0)*(A-2.D0)*U+(60.D0*A+120.D0)CORA1347
     1*((118.D0*A+796.D0)*A+1401.D0)/X-(2.D0*A+5.D0)*(6.D0*A+9.D0)*((4.DCORA1348
     20*A+47.D0)*U**2-(10.D0*A+20.D0)*((24.D0*A+117.D0)*U-(14.D0*A+42.D0CORA1349
     3)*(52.D0*A+215.D0)/X)/X+(2.D0*A+7.D0)*(2.D0*A+1.D0)*C1))/(Y*X**2*SCORA1350
     4)                                                                 CORA1351
   70 IF (LL.EQ.L2-L1+1) RETURN                                         CORA1352
      C1=C(3)                                                           CORA1353
      C(3)=C(2)                                                         CORA1354
      C(2)=C1                                                           CORA1355
      RETURN                                                            CORA1356
   71 WRITE (MW,1000) LQ,LL                                             CORA1357
      RETURN                                                            CORA1358
   72 WRITE (MW,1001) LQ,LL                                             CORA1359
      RETURN                                                            CORA1360
   73 WRITE (MW,1002) LQ,LL                                             CORA1361
      RETURN                                                            CORA1362
 1000 FORMAT (' IN CORA: LQ =',I3,' TOO LARGE OR L2-L1+1=',I3,' OUT OF BCORA1363
     1OUNDS.')                                                          CORA1364
 1001 FORMAT (' IN CORA: COEFFICIENTS SUPPRESSED FOR LQ =',I3,' AND L1-LCORA1365
     12+1=',I3)                                                         CORA1366
 1002 FORMAT (' IN CORA: COEFFICIENTS NOT GIVEN BY THE METHOD FOR LQ =',CORA1367
     1I3,' AND L2-L1+1=',I3)                                            CORA1368
      END                                                               CORA1369
