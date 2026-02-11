C 09/03/07                                                      ECIS06  LINS-000
      SUBROUTINE LINS(A,IA,B,IB,X,IX,Y,IY,N,M,K,IER)                    LINS-001
C LINS: JLSB2 VERSION 01 18/12/68 MATH L003 OF SACLAY                   LINS-002
C SOLUTION OF DOUBLE PRECISION COMPLEX LINEAR SYSTEMS WITH REAL AND     LINS-003
C IMAGINARY PARTS IN DIFFERENT AREAS.                                   LINS-004
C INPUT:     A:       REAL COMPONENT OF MATRIX.                         LINS-005
C            IA:      FIRST DIMENSION OF AREA A.                        LINS-006
C            B:       IMAGINARY COMPONENT OF MATRIX.                    LINS-007
C            IB:      FIRST DIMENSION OF AREA B.                        LINS-008
C            X:       REAL SECOND MEMBERS.                              LINS-009
C            IX:      FIRST DIMENSION OF AREA X.                        LINS-010
C            Y:       IMAGINARY SECOND MEMBERS.                         LINS-011
C            IY:      FIRST DIMENSION OF AREA Y.                        LINS-012
C            N:       ORDER OF THE SYSTEM.                              LINS-013
C            M:       NUMBER OF SECOND MEMBERS.                         LINS-014
C OUTPUT:    X,Y:     SECOND MEMBERS ARE REPLACED BY SOLUTIONS.         LINS-015
C            IER:     RETURNS DIAGNOSTIC: 0 NON SINGULAR MATRIX,        LINS-016
C                                         1 SINGULAR MATRIX,            LINS-017
C                                         2 QUASI SINGULAR MATRIX.      LINS-018
C WORKING AREAS:                                                        LINS-019
C            K:       WORKING FIELD OF DIMENSION >/= N.                 LINS-020
C***********************************************************************LINS-021
      IMPLICIT REAL*8 (A-H,O-Z)                                         LINS-022
      DIMENSION A(IA,*),B(IB,*),X(IX,*),Y(IY,*),K(*)                    LINS-023
      COMMON /INOUT/ MR,MW,MS                                           LINS-024
      IER=0                                                             LINS-025
      DO 1 I=1,N                                                        LINS-026
    1 K(I)=I                                                            LINS-027
      DO 14 I=1,N                                                       LINS-028
      AMAX=DABS(A(I,I))+DABS(B(I,I))                                    LINS-029
      JMAX=I                                                            LINS-030
      I1=I+1                                                            LINS-031
      IF (I.EQ.N) GO TO 3                                               LINS-032
      DO 2 J=I1,N                                                       LINS-033
      AMAY=DABS(A(I,J))+DABS(B(I,J))                                    LINS-034
      IF (AMAX.GE.AMAY) GO TO 2                                         LINS-035
      AMAX=AMAY                                                         LINS-036
      JMAX=J                                                            LINS-037
    2 CONTINUE                                                          LINS-038
    3 IF (AMAX.EQ.0.D0) GO TO 21                                        LINS-039
      IF (JMAX.EQ.I) GO TO 5                                            LINS-040
      DO 4 I2=1,N                                                       LINS-041
      AUX=B(I2,I)                                                       LINS-042
      B(I2,I)=B(I2,JMAX)                                                LINS-043
      B(I2,JMAX)=AUX                                                    LINS-044
      AUX=A(I2,I)                                                       LINS-045
      A(I2,I)=A(I2,JMAX)                                                LINS-046
    4 A(I2,JMAX)=AUX                                                    LINS-047
      NAB=K(JMAX)                                                       LINS-048
      K(JMAX)=K(I)                                                      LINS-049
      K(I)=NAB                                                          LINS-050
    5 IF (I.EQ.1) GO TO 7                                               LINS-051
      SR=0.D0                                                           LINS-052
      SI=0.D0                                                           LINS-053
      T=0.D0                                                            LINS-054
      I4=I-1                                                            LINS-055
      DO 6 IT=1,I4                                                      LINS-056
      PR=A(IT,I)*A(I,IT)-B(IT,I)*B(I,IT)                                LINS-057
      PI=A(IT,I)*B(I,IT)+A(I,IT)*B(IT,I)                                LINS-058
      SR=SR+PR                                                          LINS-059
      SI=SI+PI                                                          LINS-060
    6 T=T+DABS(PR)+DABS(PI)                                             LINS-061
      ERA=1.D-16*(T+DABS(A(I,I)-SR)+DABS(B(I,I)-SI))                    LINS-062
      IF (AMAX.GT.ERA) GO TO 7                                          LINS-063
      IER=2                                                             LINS-064
      WRITE (MW,1000)                                                   LINS-065
      GO TO 22                                                          LINS-066
    7 IF (I.EQ.N) GO TO 9                                               LINS-067
      DO 8 J=I1,N                                                       LINS-068
      AA=A(I,J)                                                         LINS-069
      BB=B(I,J)                                                         LINS-070
      AI=A(I,I)                                                         LINS-071
      BI=B(I,I)                                                         LINS-072
      D=AI*AI+BI*BI                                                     LINS-073
      A(I,J)=(AA*AI+BB*BI)/D                                            LINS-074
    8 B(I,J)=(BB*AI-AA*BI)/D                                            LINS-075
    9 DO 10 J=1,M                                                       LINS-076
      AA=X(I,J)                                                         LINS-077
      BB=Y(I,J)                                                         LINS-078
      AI=A(I,I)                                                         LINS-079
      BI=B(I,I)                                                         LINS-080
      D=AI*AI+BI*BI                                                     LINS-081
      X(I,J)=(AA*AI+BB*BI)/D                                            LINS-082
   10 Y(I,J)=(BB*AI-AA*BI)/D                                            LINS-083
      IF (I.EQ.N) GO TO 14                                              LINS-084
      DO 13 I3 =I1,N                                                    LINS-085
      DO 11 J3 =I1,N                                                    LINS-086
      B(I3,J3)=B(I3,J3)-A(I3,I)*B(I,J3)-B(I3,I)*A(I,J3)                 LINS-087
   11 A(I3,J3)=A(I3,J3)-A(I3,I)*A(I,J3)+B(I3,I)*B(I,J3)                 LINS-088
      DO 12 J3 =1,M                                                     LINS-089
      Y(I3,J3)=Y(I3,J3)-A(I3,I)*Y(I,J3)-B(I3,I)*X(I,J3)                 LINS-090
   12 X(I3,J3)=X(I3,J3)-A(I3,I)*X(I,J3)+B(I3,I)*Y(I,J3)                 LINS-091
   13 CONTINUE                                                          LINS-092
   14 CONTINUE                                                          LINS-093
      IF (N.EQ.1) GO TO 22                                              LINS-094
      DO 17 KC=1,M                                                      LINS-095
      J=N                                                               LINS-096
   15 I=J-1                                                             LINS-097
   16 X(I,KC)=X(I,KC)-X(J,KC)*A(I,J)+Y(J,KC)*B(I,J)                     LINS-098
      Y(I,KC)=Y(I,KC)-X(J,KC)*B(I,J)-Y(J,KC)*A(I,J)                     LINS-099
      I=I-1                                                             LINS-100
      IF (I.NE.0) GO TO 16                                              LINS-101
      J=J-1                                                             LINS-102
      IF (J.NE.1) GO TO 15                                              LINS-103
   17 CONTINUE                                                          LINS-104
      DO 20 I=1,N                                                       LINS-105
   18 J=K(I)                                                            LINS-106
      IF (J.LE.I) GO TO 20                                              LINS-107
      K(I)=K(J)                                                         LINS-108
      K(J)=J                                                            LINS-109
      DO 19 MP=1,M                                                      LINS-110
      AUX=Y(J,MP)                                                       LINS-111
      Y(J,MP)=Y(I,MP)                                                   LINS-112
      Y(I,MP)=AUX                                                       LINS-113
      AUX=X(J,MP)                                                       LINS-114
      X(J,MP)=X(I,MP)                                                   LINS-115
   19 X(I,MP)=AUX                                                       LINS-116
      GO TO 18                                                          LINS-117
   20 CONTINUE                                                          LINS-118
      GO TO 22                                                          LINS-119
   21 IER=1                                                             LINS-120
      WRITE (MW,1001)                                                   LINS-121
   22 RETURN                                                            LINS-122
 1000 FORMAT (' *****  LINS  ***** QUASI SINGULAR MATRIX.')             LINS-123
 1001 FORMAT (' *****  LINS  *****       SINGULAR MATRIX.')             LINS-124
      END                                                               LINS-125
