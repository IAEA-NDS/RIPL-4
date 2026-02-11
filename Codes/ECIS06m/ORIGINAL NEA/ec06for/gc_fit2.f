C 07/03/07                                                      ECIS06  FIT2-000
      SUBROUTINE FIT2(M,N,A,D,IP,IR,LLO)                                FIT2-001
C IF LLO=.TRUE. EX SUBROUTINE LILESQ WRITTEN BY SCHWEIMER..             FIT2-002
C  LINEAR LEAST SQUARES PROBLEM ||B-A*D||=MIN(D)                        FIT2-003
C  SOLVED BY HOUSEHOLDER TRANSFORMATIONS.                               FIT2-004
C INPUT:     M:       NUMBER OF ROWS OF A AND B.                        FIT2-005
C            N:       NUMBER OF COLUMNS OF A AND ROWS OF D.             FIT2-006
C            A:       M*N MATRIX FOLLOWED BY THE VECTOR B OF M          FIT2-007
C                     COMPONENTS (DESTROYED).                           FIT2-008
C OUTPUT:    A:       THE UPPER PART CONTAINS THE TRANSFORMED INPUT A   FIT2-009
C                     A(2,1) CONTAINS THE SQUARE OF THE EUCLIDEAN NORM. FIT2-010
C            D:       VECTOR OF VARIABLES, THE REDUNDANT VARIABLES ARE  FIT2-011
C                     SET TO ZERO. THE ||D||=MIN IS NOT USED BECAUSE THEFIT2-012
C                     COMPONENTS OF D ARE ASSUMED NOT COMMENSURABLE.    FIT2-013
C            IP:      PERMUTATION VECTOR OF N COMPONENTS, IT CONTAINS   FIT2-014
C                     THE COLUMN LABELS OF MATRIX A ORDERED ACCORDING   FIT2-015
C                     THEIR IMPORTANCE IN REDUCING THE EUCLIDEAN NORM.  FIT2-016
C            IR:      ERROR CODE: IR=0 NO ERROR, IR=-1 ALL COMPONENTS   FIT2-017
C                     OF D ARE ZERO AND MAY BE REDUNDANT, IR>0 THE      FIT2-018
C                     FIRST IR COMPONENTS OF IP CONTAIN THE LABELS OF   FIT2-019
C                     THE NONZERO COMPONENTS OF D, THE REMAINING        FIT2-020
C                     COMPONENTS OF D ARE ZERO AND MAY BE REDUNDANT.    FIT2-021
C  NOTE: ALL ARITHMETIC OPERATIONS ARE PERFORMED IN DOUBLE PRECISION,   FIT2-022
C  AN ITERATIVE IMPROVEMENT IS IMPOSSIBLE WITHOUT SAVING A AND B.       FIT2-023
C  THE ROUND OFF ERROR OF ||B-A*D||**2 IS APPROXIMATELY GIVEN BY        FIT2-024
C  ||B(INITIAL)||**2 - ||B(TRANSFORMED)||**2                            FIT2-025
C                                                                       FIT2-026
C IF LLO=.FALSE. EX SUBROUTINE INVATA WRITTEN BY SCHWEIMER.             FIT2-027
C  INVERSION OF THE PRODUCT MATRIX A(TRANSPOSED)*A                      FIT2-028
C  THE MATRIX A IS REDUCED TO AN UPPER TRIANGULAR MATRIX R BY           FIT2-029
C  HOUSEHOLDER TRANSFORMATIONS. THE REMAINING COMPUTATION IS STRAIGHT   FIT2-030
C  FORWARD.                                                             FIT2-031
C INPUT:     M:       NUMBER OF COLUMNS OF MATRIX A.                    FIT2-032
C            N:       NUMBER OF ROWS OF MATRIX A, M >= N > 0.           FIT2-033
C            A:       INPUT MATRIX (DESTROYED).                         FIT2-034
C OUTPUT:    A:       TRIANGULAR MATRIX R, R=A(I,J) I<=J=1,N.           FIT2-035
C            D:       VECTOR OF LENGTH (N*(N+1))/2, IT CONTAINS THE     FIT2-036
C                     UPPER TRIANGULAR PART OF THE INVERSE OF A(T)*A.   FIT2-037
C            IP:      PERMUTATION VECTOR OF LENGTH N, ITS FIRST IR      FIT2-038
C                     COMPONENTS CONTAIN THE LABELS OF THE USEFUL       FIT2-039
C                     COLUMNS OF A, THE LAST COMPONENTS CONTAIN THE     FIT2-040
C                     LABELS OF THE COLUMNS WHICH ARE LINEAR            FIT2-041
C                     COMBINATIONS OF THE FIRST.                        FIT2-042
C            IR:      ERROR CODE: IR=0 NO ERROR (RANK OF MATRIX A IS N),FIT2-043
C                     IR=-1 RANK OF MATRIX A IS ZERO, IR>0 RANK OF      FIT2-044
C                     MATRIX A IS IR, THE INVERSE OF A(T)*A IS COMPUTED FIT2-045
C                     CONSIDERING THE IR COLUMNS OF A INDICATED BY THE  FIT2-046
C                     FIRST IR COMPONENTS OF IP.                        FIT2-047
C***********************************************************************FIT2-048
      IMPLICIT REAL*8 (A-H,O-Z)                                         FIT2-049
      LOGICAL LLO                                                       FIT2-050
      DIMENSION A(M,*),D(*),IP(N)                                       FIT2-051
      IR=0                                                              FIT2-052
      N1=N                                                              FIT2-053
      IF (LLO) N1=N+1                                                   FIT2-054
      DO 1 J=1,N                                                        FIT2-055
    1 IP(J)=J                                                           FIT2-056
C ROTATION LOOP.                                                        FIT2-057
      DO 13 K=1,N                                                       FIT2-058
C PIVOT ELEMENT    COLUMN J WHICH GENERATES THE LARGEST NEW A(*,M) AND  FIT2-059
C AND LINE I OF THE LARGEST ELEMENT OF COLUMN J IF LLO=.TRUE.           FIT2-060
C COLUMN AND LINE OF LARGEST ELEMENT IF LLO=.FALSE.                     FIT2-061
      U=0.D0                                                            FIT2-062
      DO 7 J=K,N                                                        FIT2-063
      C=0.D0                                                            FIT2-064
      DO 2 I=K,M                                                        FIT2-065
      IF (DABS(A(I,J)).LE.DABS(C)) GO TO 2                              FIT2-066
      L2=I                                                              FIT2-067
      C=A(I,J)                                                          FIT2-068
    2 CONTINUE                                                          FIT2-069
      IF (C.EQ.0.D0) GO TO 7                                            FIT2-070
      S=0.D0                                                            FIT2-071
      IF (LLO) GO TO 4                                                  FIT2-072
      IF (DABS(C).LT.U) GO TO 7                                         FIT2-073
      U=DABS(C)                                                         FIT2-074
      DO 3 I=K,M                                                        FIT2-075
      V=A(I,J)/C                                                        FIT2-076
    3 S=S+V*V                                                           FIT2-077
      GO TO 6                                                           FIT2-078
    4 T=0.D0                                                            FIT2-079
      DO 5 I=K,M                                                        FIT2-080
      V=A(I,J)/C                                                        FIT2-081
      S=S+V*V                                                           FIT2-082
    5 T=T+V*A(I,N1)                                                     FIT2-083
      IF (U.GE.T*(T/S)) GO TO 7                                         FIT2-084
      U=T*(T/S)                                                         FIT2-085
    6 SIG=C*DSQRT(S)                                                    FIT2-086
      L=J                                                               FIT2-087
      L1=L2                                                             FIT2-088
    7 CONTINUE                                                          FIT2-089
      IF (U.EQ.0.D0) GO TO 14                                           FIT2-090
C PERMUTE COLUMNS OF A(K).                                              FIT2-091
      I=IP(L)                                                           FIT2-092
      IP(L)=IP(K)                                                       FIT2-093
      IP(K)=I                                                           FIT2-094
      DO 8 I=1,M                                                        FIT2-095
      C=A(I,L)                                                          FIT2-096
      A(I,L)=A(I,K)                                                     FIT2-097
    8 A(I,K)=C                                                          FIT2-098
C PERMUTE LINES OF A(K).                                                FIT2-099
      DO 9 J=K,N1                                                       FIT2-100
      C=A(K,J)                                                          FIT2-101
      A(K,J)=A(L1,J)                                                    FIT2-102
    9 A(L1,J)=C                                                         FIT2-103
C ROTATION OF THE LOWER COLUMN FRAGMENT OF A(K) AND B(K).               FIT2-104
      U=SIG+A(K,K)                                                      FIT2-105
      V=A(K,K)/SIG                                                      FIT2-106
      A(K,K)=-SIG                                                       FIT2-107
      L=K+1                                                             FIT2-108
      IF (L.GT.M) A(K,L)=-A(K,L)                                        FIT2-109
      IF (L.GT.N1.OR.L.GT.M) GO TO 13                                   FIT2-110
      DO 12 J=L,N1                                                      FIT2-111
      S=V*A(K,J)                                                        FIT2-112
      DO 10 I=L,M                                                       FIT2-113
      T=A(I,K)/SIG                                                      FIT2-114
   10 S=S+T*A(I,J)                                                      FIT2-115
      T=(A(K,J)+S)/U                                                    FIT2-116
      A(K,J)=-S                                                         FIT2-117
      DO 11 I=L,M                                                       FIT2-118
   11 A(I,J)=A(I,J)-T*A(I,K)                                            FIT2-119
   12 CONTINUE                                                          FIT2-120
   13 CONTINUE                                                          FIT2-121
C END OF ROTATION LOOP.                                                 FIT2-122
      K=N                                                               FIT2-123
      GO TO 15                                                          FIT2-124
   14 K=K-1                                                             FIT2-125
      IR=K                                                              FIT2-126
   15 IF (LLO) GO TO 22                                                 FIT2-127
      IF (K.EQ.0) GO TO 29                                              FIT2-128
C INVERSE OF THE TRIANGULAR MATRIX R STORED IN D.                       FIT2-129
      DO 18 J=1,K                                                       FIT2-130
      D(J)=A(J,J)                                                       FIT2-131
      A(J,J)=1.D0/D(J)                                                  FIT2-132
      IF (J.EQ.1) GO TO 18                                              FIT2-133
      I=J                                                               FIT2-134
      DO 17 L=2,J                                                       FIT2-135
      I1=I                                                              FIT2-136
      I=I-1                                                             FIT2-137
      S=0.D0                                                            FIT2-138
      DO 16 L1=I1,J                                                     FIT2-139
   16 S=S+A(I,L1)*A(J,L1)                                               FIT2-140
   17 A(J,I)=-S/D(I)                                                    FIT2-141
   18 CONTINUE                                                          FIT2-142
C INVERSE OF THE PRODUCT MATRIX.                                        FIT2-143
      DO 21 J=1,K                                                       FIT2-144
      DO 20 I=1,J                                                       FIT2-145
      L=MAX0(IP(J),IP(I))                                               FIT2-146
      IJ=IP(I)+IP(J)+(L*(L-3))/2                                        FIT2-147
      S=0.D0                                                            FIT2-148
      DO 19 L1=J,K                                                      FIT2-149
   19 S=S+A(L1,I)*A(L1,J)                                               FIT2-150
   20 D(IJ)=S                                                           FIT2-151
   21 CONTINUE                                                          FIT2-152
      GO TO 30                                                          FIT2-153
C SQUARE OF THE EUCLIDEAN NORM.                                         FIT2-154
   22 S=0.D0                                                            FIT2-155
      L=K+1                                                             FIT2-156
      IF (K.EQ.M) GO TO 24                                              FIT2-157
      DO 23 I=L,M                                                       FIT2-158
   23 S=S+A(I,N1)*A(I,N1)                                               FIT2-159
   24 A(2,1)=S                                                          FIT2-160
      IF (K.EQ.N) GO TO 26                                              FIT2-161
C COMPONENTS OF D WHICH DO NOT REDUCE THE EUCLIDEAN NORM.               FIT2-162
      DO 25 J=L,N                                                       FIT2-163
      IJ=IP(J)                                                          FIT2-164
   25 D(IJ)=0.D0                                                        FIT2-165
      IF (K.EQ.0) GO TO 29                                              FIT2-166
C COMPUTATION OF D.                                                     FIT2-167
   26 IJ=IP(K)                                                          FIT2-168
      D(IJ)=A(K,N1)/A(K,K)                                              FIT2-169
      IF (K.EQ.1) GO TO 30                                              FIT2-170
      DO 28 J=2,K                                                       FIT2-171
      L=K+2-J                                                           FIT2-172
      S=A(L-1,N1)                                                       FIT2-173
      DO 27 I=L,K                                                       FIT2-174
      IJ=IP(I)                                                          FIT2-175
   27 S=S-A(L-1,I)*D(IJ)                                                FIT2-176
      IJ=IP(L-1)                                                        FIT2-177
   28 D(IJ)=S/A(L-1,L-1)                                                FIT2-178
      GO TO 30                                                          FIT2-179
C ERROR CODE.                                                           FIT2-180
   29 IR=-1                                                             FIT2-181
   30 RETURN                                                            FIT2-182
      END                                                               FIT2-183
