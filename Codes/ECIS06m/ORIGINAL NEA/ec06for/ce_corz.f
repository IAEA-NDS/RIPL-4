C 27/06/06                                                      ECIS06  CORZ-000
      SUBROUTINE CORZ(XI,XF,T,FGI,FGF,LM1,LM2,LM,W)                     CORZ-001
C    COMPUTATION OF THE INTEGRALS FROM 1 TO INFINITY OF PRODUCTS OF     CORZ-002
C REGULAR OR IRREGULAR BESSEL FUNCTIONS OF SAME L-VALUE AND ARGUMENT    CORZ-003
C XI*R AND XF*R RESPECTIVELY, DIVIDED BY R**3 FOR L=1 TO LM, AND OF     CORZ-004
C INTEGRALS FROM 0 TO INFINITY OF REGULAR FUNCTIONS OF THESE ARGUMENTS  CORZ-005
C DIVIDED BY R**3.                                                      CORZ-006
C    THE LIMIT OF STABILITY OF THE UPWARDS RECURRENCE HAS               CORZ-007
C     BEEN FOUND TO BE   LM2*(1-XI/XF) < 3 .                            CORZ-008
C INPUT:     XI,XF:   WAVE NUMBERS MULTIPLIED BY THE MATCHING RADIUS.   CORZ-009
C            FGI,FGF: REGULAR AND IRREGULAR BESSEL FUNCTIONS.           CORZ-010
C            LM1:     MAXIMUM NUMBER OF INTEGRALS FROM 1 TO INFINITY.   CORZ-011
C            LM2:     NUMBER OF INTEGRALS FROM 0 TO INFINITY.           CORZ-012
C            LM:      ACTUAL NUMBER OF INTEGRALS FROM 1 TO INFINITY.    CORZ-013
C OUTPUT:    T:       INTEGRALS FROM 0 TO INFINITY.                     CORZ-014
C            W:       INTEGRALS FROM RS TO INFINITY.                    CORZ-015
C***********************************************************************CORZ-016
      IMPLICIT REAL*8 (A-H,O-Z)                                         CORZ-017
      DIMENSION T(*),FGI(LM1,*),FGF(LM1,*),W(LM1,*),B(4),C(2),D(2)      CORZ-018
C INTEGRALS FROM 0 TO INFINITY.                                         CORZ-019
C THE RECURRENCE RELATION STARTS FROM L=0 WITH L*(0,3;L)=1/3 FOR L=0.   CORZ-020
      T(1)=1.D0/3.D0                                                    CORZ-021
      IF (XI.NE.XF) GO TO 2                                             CORZ-022
C WHEN KI=KF   (0,3;L)=1/(2*L*(L+1)).                                   CORZ-023
      A1=0.D0                                                           CORZ-024
      DO 1 I=2,LM2                                                      CORZ-025
      A1=A1+1.D0                                                        CORZ-026
    1 T(I)=.5D0/(A1*(A1+1.D0))                                          CORZ-027
      GO TO 11                                                          CORZ-028
    2 X=XI/XF                                                           CORZ-029
      IF (X.GT.1.D0) X=1.D0/X                                           CORZ-030
C (0,3;L) = X**L GT(1/2) GT(L+1)/2/GT(L+3/2) 2F1(L,-1/2;L+3/2;X**2).    CORZ-031
      IF (X.GE..9D0) GO TO 4                                            CORZ-032
C DIRECT EXPANSION OF THE 2F1 FUNCTION.                                 CORZ-033
      A1=X**2                                                           CORZ-034
      T(2)=X/3.D0                                                       CORZ-035
      A2=T(2)                                                           CORZ-036
      A3=0.D0                                                           CORZ-037
      DO 3 J=1,2000                                                     CORZ-038
      A3=A3+1.D0                                                        CORZ-039
      A2=A2*A1*(A3-1.5D0)/(A3+1.5D0)                                    CORZ-040
      T(2)=T(2)+A2                                                      CORZ-041
      IF (DABS(A2).LT.1.D-15*DABS(T(2))) GO TO 6                        CORZ-042
    3 CONTINUE                                                          CORZ-043
      GO TO 6                                                           CORZ-044
C ANALYTIC CONTINUATION OF THE 2F1 FUNCTION - ERDELYI,.... 15.3.11.     CORZ-045
    4 A1=1.D0-X**2                                                      CORZ-046
      A2=X*A1**2/16.D0                                                  CORZ-047
      A3=DLOG(A1/4.D0)+2.D0                                             CORZ-048
      T(2)=X*(2.D0+A1)/8.D0+A2*A3                                       CORZ-049
      A4=0.D0                                                           CORZ-050
      DO 5 J=1,2000                                                     CORZ-051
      A4=A4+1.D0                                                        CORZ-052
      A3=A3-.5D0/((A4+.5D0)*A4)                                         CORZ-053
      A2=A2*A1*(A4+.5D0)/A4                                             CORZ-054
      A5=A2*A3                                                          CORZ-055
      T(2)=T(2)+A5                                                      CORZ-056
      IF (DABS(A5).LT.1.D-15*DABS(T(2))) GO TO 6                        CORZ-057
    5 CONTINUE                                                          CORZ-058
C RECURRENCE RELATION FOR  (0,3;L) :                                    CORZ-059
C  2*(L-1)*(0,3;L-1)-(2*L+1)*(X+1/X)*(0,3;L)+2*(L+2)*(0,3;L+1)=0.       CORZ-060
C UPWARDS RECURRENCE.                                                   CORZ-061
    6 A1=X+1.D0/X                                                       CORZ-062
      XX=LM2*(1-X)                                                      CORZ-063
      IF (XX.GT.3.D0) GO TO 8                                           CORZ-064
      T(3)=.5D0*(A1*T(2)-T(1))                                          CORZ-065
      A2=0.D0                                                           CORZ-066
      DO 7 I=4,LM2                                                      CORZ-067
      A2=A2+1.D0                                                        CORZ-068
    7 T(I)=((A2+1.5D0)*A1*T(I-1)-A2*T(I-2))/(A2+3.D0)                   CORZ-069
      GO TO 11                                                          CORZ-070
C DOWNWARDS RECURRENCE.                                                 CORZ-071
    8 LK=4*LM2                                                          CORZ-072
      A2=DFLOAT(LK-1)                                                   CORZ-073
      A3=A2/(A2+2.D0)                                                   CORZ-074
      DO 9 I=3,LK                                                       CORZ-075
      J=LK-I+3                                                          CORZ-076
      A2=A2-1.D0                                                        CORZ-077
      A3=A2/((A2+1.5D0)*A1-(A2+3.D0)*A3)                                CORZ-078
      IF (J.LE.LM2) T(J)=A3                                             CORZ-079
    9 CONTINUE                                                          CORZ-080
      DO 10 I=3,LM2                                                     CORZ-081
   10 T(I)=T(I)*T(I-1)                                                  CORZ-082
   11 B(1)=XI+XF                                                        CORZ-083
      B(2)=XI-XF                                                        CORZ-084
C COMPUTATION OF THE INTEGRALS FROM 1 TO INFINITY OF EXP(I*B*R)/R**3 DR.CORZ-085
      DO 17 I=1,2                                                       CORZ-086
      IF (B(I).LT.5.D0) GO TO 13                                        CORZ-087
C USE OF PADE APPROXIMANT.                                              CORZ-088
      A1=0.D0                                                           CORZ-089
      A2=0.D0                                                           CORZ-090
      A3=21.D0                                                          CORZ-091
      DO 12 J=1,20                                                      CORZ-092
      A3=A3-1.D0                                                        CORZ-093
      A2=A2-B(I)                                                        CORZ-094
      B1=A1**2+A2**2                                                    CORZ-095
      B2=A1*A3/B1+1.D0                                                  CORZ-096
      B3=-A2*A3/B1                                                      CORZ-097
      B1=B2**2+B3**2                                                    CORZ-098
      A1=B2*(A3+2.D0)/B1                                                CORZ-099
   12 A2=-B3*(A3+2.D0)/B1                                               CORZ-100
      A2=A2-B(I)                                                        CORZ-101
      B1=A1**2+A2**2                                                    CORZ-102
      C(I)=(DCOS(B(I))*A1+DSIN(B(I))*A2)/B1                             CORZ-103
      D(I)=(DSIN(B(I))*A1-DCOS(B(I))*A2)/B1                             CORZ-104
      GO TO 17                                                          CORZ-105
   13 IF (B(I).EQ.0.D0) GO TO 16                                        CORZ-106
C USE OF THE TAYLOR EXPANSION.                                          CORZ-107
      A1=DABS(B(I))                                                     CORZ-108
      C(I)=.5D0*(1.D0-A1**2*(.922784335098467D0-DLOG(A1)))              CORZ-109
      D(I)=A1*(1.D0-.78539816339744831D0*A1)                            CORZ-110
      A2=A1**2/2.D0                                                     CORZ-111
      A3=0.D0                                                           CORZ-112
      DO 14 J=2,20                                                      CORZ-113
      A3=A3+1.D0                                                        CORZ-114
      A2=A2*A1/(2.D0*A3+1.D0)                                           CORZ-115
      D(I)=D(I)+A2/(2.D0*A3-1.D0)                                       CORZ-116
      A2=-A2*A1/(2.D0*A3+2.D0)                                          CORZ-117
      C(I)=C(I)+A2/(2.D0*A3)                                            CORZ-118
      IF (DABS(A2).LT.1.D-16) GO TO 15                                  CORZ-119
   14 CONTINUE                                                          CORZ-120
   15 IF (A1.NE.B(I)) D(I)=-D(I)                                        CORZ-121
      GO TO 17                                                          CORZ-122
   16 C(I)=.5D0                                                         CORZ-123
      D(I)=0.D0                                                         CORZ-124
   17 CONTINUE                                                          CORZ-125
C STORAGE OF THE TWO FIRST INTEGRALS.                                   CORZ-126
      A1=2.D0*XI*XF                                                     CORZ-127
      DO 18 I=1,2                                                       CORZ-128
      W(I,1)=(C(2)-C(1))/A1                                             CORZ-129
      W(I,2)=(D(1)-D(2))/A1                                             CORZ-130
      W(I,3)=(D(1)+D(2))/A1                                             CORZ-131
      W(I,4)=(C(2)+C(1))/A1                                             CORZ-132
      IF (I.EQ.2) GO TO 18                                              CORZ-133
      C(1)=C(1)*B(2)**2+DCOS(B(1))+B(1)*DSIN(B(1))                      CORZ-134
      D(1)=D(1)*B(2)**2+DSIN(B(1))-B(1)*DCOS(B(1))                      CORZ-135
      C(2)=C(2)*B(1)**2+DCOS(B(2))+B(2)*DSIN(B(2))                      CORZ-136
      D(2)=D(2)*B(1)**2+DSIN(B(2))-B(2)*DCOS(B(2))                      CORZ-137
      A1=2.D0*A1**2                                                     CORZ-138
   18 CONTINUE                                                          CORZ-139
C UPWARDS RECURRENCE.                                                   CORZ-140
      A1=XI/XF+XF/XI                                                    CORZ-141
      A2=2.D0*XI*XF                                                     CORZ-142
      A3=2.D0                                                           CORZ-143
      DO 20 J=3,LM                                                      CORZ-144
      A3=A3+1.D0                                                        CORZ-145
      B(2)=(FGI(J-2,3)*FGF(J-2,1)-FGI(J,3)*FGF(J,1))/A2                 CORZ-146
      B(3)=(FGI(J-2,1)*FGF(J-2,3)-FGI(J,1)*FGF(J,3))/A2                 CORZ-147
      B(4)=(FGI(J-2,3)*FGF(J-2,3)-FGI(J,3)*FGF(J,3))/A2                 CORZ-148
      DO 19 I=2,4                                                       CORZ-149
   19 W(J,I)=((A3-1.5D0)*A1*W(J-1,I)-(A3-3.D0)*W(J-2,I)-B(I))/A3        CORZ-150
   20 CONTINUE                                                          CORZ-151
      LN=MAX0(LM,IDINT(DMAX1(XI,XF)))+50                                CORZ-152
      B2=0.D0                                                           CORZ-153
      B3=0.D0                                                           CORZ-154
      C2=0.D0                                                           CORZ-155
      C3=1.D-15                                                         CORZ-156
      D2=0.D0                                                           CORZ-157
      D3=1.D-15                                                         CORZ-158
      A3=LN                                                             CORZ-159
      DO 24 I=2,LN                                                      CORZ-160
      J=LN-I+2                                                          CORZ-161
      A3=A3-1.D0                                                        CORZ-162
      B1=B2                                                             CORZ-163
      B2=B3                                                             CORZ-164
      C1=C2                                                             CORZ-165
      C2=C3                                                             CORZ-166
      D1=D2                                                             CORZ-167
      D2=D3                                                             CORZ-168
      C3=(2.D0*A3+3.D0)*C2/XI-C1                                        CORZ-169
      D3=(2.D0*A3+3.D0)*D2/XF-D1                                        CORZ-170
      B3=((A3+1.5D0)*A1*B2-(A3+3.D0)*B1-(C3*D3-C1*D1)/A2)/A3            CORZ-171
      IF (J.LE.LM) W(J,1)=B3                                            CORZ-172
      IF (D3.LT.1.D15) GO TO 21                                         CORZ-173
      D2=D2*1.D-30                                                      CORZ-174
      D3=D3*1.D-30                                                      CORZ-175
      GO TO 22                                                          CORZ-176
   21 IF (C3.LT.1.D15) GO TO 24                                         CORZ-177
      C2=C2*1.D-30                                                      CORZ-178
      C3=C3*1.D-30                                                      CORZ-179
   22 B2=B2*1.D-30                                                      CORZ-180
      B3=B3*1.D-30                                                      CORZ-181
      IF (J.GT.LM) GO TO 24                                             CORZ-182
      DO 23 K=J,LM                                                      CORZ-183
   23 W(K,1)=W(K,1)*1.D-30                                              CORZ-184
   24 CONTINUE                                                          CORZ-185
      A1=FGI(2,1)*FGF(2,1)/(C3*D3)                                      CORZ-186
      DO 25 I=2,LM                                                      CORZ-187
   25 W(I,1)=W(I,1)*A1+T(I)                                             CORZ-188
      RETURN                                                            CORZ-189
      END                                                               CORZ-190
