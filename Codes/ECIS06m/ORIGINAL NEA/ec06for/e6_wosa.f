C 15/11/06                                                      ECIS06  WOSA-000
      SUBROUTINE WOSA(V,A,P,CC,R,N,LX)                                  WOSA-001
C COMPUTATION OF WOODS-SAXON FORM-FACTORS AT SOME POWER AND THEIR       WOSA-002
C DERIVATIVES.                                                          WOSA-003
C LIMITED TO THE SIXTH DERIVATIVE (THE FOURTH ONE ONLY IS USED IN ECIS).WOSA-004
C INPUT:     A:       DIFFUSENESS.                                      WOSA-005
C            P:       THE WOOD-SAXON PARAMETERS ARE AT THE POWER (1+P). WOSA-006
C            CC:      VALUE OF EXP(-(R-RR/A) WHERE R IS THE RADIUS.     WOSA-007
C            R:       VALUE OF RADIUS OF THE FORM FACTOR.               WOSA-008
C            N:       NUMBER OF VALUES (ONE PLUS HIGHEST DERIVATIVE).   WOSA-009
C            LX:      LOGICAL FALSE FOR SYMMETRISED FORM-FACTORS.       WOSA-010
C OUTPUT:    V:       FORM-FACTORS AND ITS DERIVATIVES.                 WOSA-011
C***********************************************************************WOSA-012
      IMPLICIT REAL*8 (A-H,O-Z)                                         WOSA-013
      LOGICAL LX                                                        WOSA-014
      DIMENSION V(7),VV(7,2)                                            WOSA-015
      COMMON /INOUT/ MR,MW,MS                                           WOSA-016
      IF (N.GT.7) GO TO 6                                               WOSA-017
      C=CC                                                              WOSA-018
      K=1                                                               WOSA-019
C FORM-FACTOR AND ITS DERIVATIVES.                                      WOSA-020
    1 B=1.D0/(1.D0+C)                                                   WOSA-021
      V(1)=B**(P+1.D0)                                                  WOSA-022
      IF (N.EQ.1) GO TO 2                                               WOSA-023
      V(2)=B*(P+1.D0)*C*V(1)/A                                          WOSA-024
      IF (N.EQ.2) GO TO 2                                               WOSA-025
      V(3)=B*(P*C+C-1.D0)*V(2)/A                                        WOSA-026
      IF (N.EQ.3) GO TO 2                                               WOSA-027
      V(4)=B*((P*C+C-2.D0)*V(3)-V(2)/A)/A                               WOSA-028
      IF (N.EQ.4) GO TO 2                                               WOSA-029
      V(5)=B*((P*C+2.D0*C-2.D0)*V(4)-(P*C+C+1.D0)*V(3)/A)/A             WOSA-030
      IF (N.EQ.5) GO TO 2                                               WOSA-031
      V(6)=B*((P*C+2.D0*C-3.D0)*V(5)-((P*C+C+3.D0)*V(4)+V(3)/A)/A)/A    WOSA-032
      IF (N.EQ.6) GO TO 2                                               WOSA-033
      V(7)=B*((P*C+3.D0*C-3.D0)*V(6)-((2.D0*P*C+3.D0*C+3.D0)*V(5)-(P*C+CWOSA-034
     1-1.D0)*V(4)/A)/A)/A                                               WOSA-035
    2 IF (LX) RETURN                                                    WOSA-036
C STORAGE OF INTERMEDIATE RESULTS FOR SYMMETRISATION.                   WOSA-037
      DO 3 I=1,N                                                        WOSA-038
    3 VV(I,K)=V(I)                                                      WOSA-039
      C=DEXP(-2.D0*DABS(R)/A)/C                                         WOSA-040
      K=K+1                                                             WOSA-041
      IF (K.EQ.2) GO TO 1                                               WOSA-042
C SYMMETRISATION.                                                       WOSA-043
      DO 5 L=1,7                                                        WOSA-044
      V(L)=VV(L,1)*VV(1,2)                                              WOSA-045
      IF (L.EQ.1) GO TO 5                                               WOSA-046
      A1=1.D0                                                           WOSA-047
      DO 4 J=2,L                                                        WOSA-048
      A1=-A1*DFLOAT(L-J+1)/DFLOAT(J-1)                                  WOSA-049
    4 V(L)=V(L)+A1*VV(L-J+1,1)*VV(J,2)                                  WOSA-050
    5 CONTINUE                                                          WOSA-051
      RETURN                                                            WOSA-052
    6 WRITE (MW,1000) N                                                 WOSA-053
      STOP                                                              WOSA-054
 1000 FORMAT (5X,I5,' VALUES ASKED IN SUBROUTINE WOSA TOO LARGE (LIMITATWOSA-055
     1ION TO 7.'///' IN WOSA  ...  STOP  ...')                          WOSA-056
      END                                                               WOSA-057
