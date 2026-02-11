C 09/02/06                                                      ECIS06  EMRO-000
      SUBROUTINE EMRO(IA,I1,I2,X2,X3,RES,NX)                            EMRO-001
C COMPUTATION OF ROTATION MATRIX ELEMENTS.                              EMRO-002
C INPUT:     IA:      TWICE THE FIRST J VALUE.                          EMRO-003
C            I1:      TWICE M1.                                         EMRO-004
C            I2:      TWICE M2.                                         EMRO-005
C            X2:      COSINES OF THETA/2.                               EMRO-006
C            X3:      SINUS OF THETA/2.                                 EMRO-007
C            NX:      NUMBER OF VALUES REQUESTED.                       EMRO-008
C OUTPUT:    RES:     MATRIX ELEMENTS:                                  EMRO-009
C                      J          J+1                     J+NX-1        EMRO-010
C                     R (THETA), R (THETA) ............. R (THETA)      EMRO-011
C                      M1,M2      M1,M2                   M1,M2         EMRO-012
C                     USING THE RECURRENCE RELATION OBTAINED FROM       EMRO-013
C                                   J                  J'               EMRO-014
C                     COS(THETA) * R (THETA)  = SUM   R (THETA)         EMRO-015
C                                   M1,M2       ON J'  M1,M2            EMRO-016
C                                                                       EMRO-017
C***********************************************************************EMRO-018
      IMPLICIT REAL*8 (O-Z)                                             EMRO-019
      DIMENSION RES(*)                                                  EMRO-020
      NR=1                                                              EMRO-021
      JA=IA                                                             EMRO-022
      X1=X2*X2-X3*X3                                                    EMRO-023
      M1=IABS(I1+I2)/2                                                  EMRO-024
      M2=IABS(I1-I2)/2                                                  EMRO-025
      MJ=M1+M2                                                          EMRO-026
      NJ=MJ                                                             EMRO-027
    1 IF (MJ.LE.JA) GO TO 2                                             EMRO-028
      RES(NR)=0.D0                                                      EMRO-029
      NR=NR+1                                                           EMRO-030
      IF (NR.GT.NX) RETURN                                              EMRO-031
      JA=JA+2                                                           EMRO-032
      GO TO 1                                                           EMRO-033
    2 S1=0.D0                                                           EMRO-034
      S2=0.D0                                                           EMRO-035
      S3=1.D0                                                           EMRO-036
      IF (MJ.EQ.0) GO TO 11                                             EMRO-037
C THE MAGNETIC QUANTUM NUMBERS ARE NOT BOTH ZEROS.                      EMRO-038
      IF (M1.EQ.M2) GO TO 8                                             EMRO-039
C NO MAGNETIC QUANTUM NUMBER IS ZERO.                                   EMRO-040
      IF (M2.EQ.0) GO TO 4                                              EMRO-041
      DO 3 I=1,M2                                                       EMRO-042
    3 S3=S3*X3*DSQRT(DFLOAT(M1+I)/DFLOAT(I))                            EMRO-043
    4 IF (M1.EQ.0) GO TO 6                                              EMRO-044
      DO 5 M=1,M1                                                       EMRO-045
    5 S3=S3*X2                                                          EMRO-046
    6 IF ((I1.GT.I2).AND.(2*(M2/2).NE.M2)) S3=-S3                       EMRO-047
      KX=1                                                              EMRO-048
      X4=0.D0                                                           EMRO-049
      Z1=DFLOAT(M1)                                                     EMRO-050
      Z2=DFLOAT(M2)                                                     EMRO-051
      SL=(Z1-Z2)/DFLOAT(NJ+2)                                           EMRO-052
    7 IF (JA.LE.NJ) GO TO 13                                            EMRO-053
C RECURRENCE.                                                           EMRO-054
      NJ=NJ+2                                                           EMRO-055
      S2=S3                                                             EMRO-056
      Z1=Z1+1.D0                                                        EMRO-057
      Z2=Z2+1.D0                                                        EMRO-058
      MJ=MJ+1                                                           EMRO-059
      X4=X4+1.D0                                                        EMRO-060
      Y2=DSQRT(Z1*Z2*X4*DFLOAT(MJ))                                     EMRO-061
      SJ=DFLOAT(NJ)                                                     EMRO-062
      S3=(SJ-1.D0)*SJ*((X1-SL)*S2-S1)/(2.D0*Y2)                         EMRO-063
      S1=2.D0*Y2*S2/((SJ+1.D0)*SJ)                                      EMRO-064
      SL=SL*(SJ-2.D0)/(SJ+2.D0)                                         EMRO-065
      GO TO 7                                                           EMRO-066
C A MAGNETIC QUANTUM NUMBER IS ZERO.                                    EMRO-067
    8 DO 9 I=1,M1                                                       EMRO-068
    9 S3=S3*X3*X2*DSQRT(DFLOAT(M1+I)/DFLOAT(I))                         EMRO-069
      KX=2                                                              EMRO-070
      IF (I1.GT.I2.AND.(2*(M2/2).NE.M2)) S3=-S3                         EMRO-071
      X4=0.D0                                                           EMRO-072
   10 IF (JA.LE.NJ) GO TO 13                                            EMRO-073
C RECURRENCE.                                                           EMRO-074
      NJ=NJ+2                                                           EMRO-075
      S2=S3                                                             EMRO-076
      MJ=MJ+1                                                           EMRO-077
      X4=X4+1.D0                                                        EMRO-078
      Y2=DSQRT(X4*DFLOAT(MJ))                                           EMRO-079
      S3=(DFLOAT(NJ-1)*X1*S2-S1)/Y2                                     EMRO-080
      S1=Y2*S2                                                          EMRO-081
      GO TO 10                                                          EMRO-082
C THE MAGNETIC QUANTUM NUMBERS ARE BOTH ZEROS.                          EMRO-083
   11 KX=3                                                              EMRO-084
   12 IF (JA.LE.NJ) GO TO 13                                            EMRO-085
C RECURRENCE.                                                           EMRO-086
      NJ=NJ+2                                                           EMRO-087
      SJ=DFLOAT(NJ/2)                                                   EMRO-088
      S2=S3                                                             EMRO-089
      S3=((2.D0*SJ-1.D0)*X1*S2-S1)/SJ                                   EMRO-090
      S1=S2*SJ                                                          EMRO-091
      GO TO 12                                                          EMRO-092
C STORAGE.                                                              EMRO-093
   13 RES(NR)=S3                                                        EMRO-094
      NR=NR+1                                                           EMRO-095
      IF (NR.GT.NX) RETURN                                              EMRO-096
      JA=JA+2                                                           EMRO-097
      GO TO ( 7 , 10 , 12 ),KX                                          EMRO-098
      END                                                               EMRO-099
