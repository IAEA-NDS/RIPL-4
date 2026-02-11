C 13/04/06                                                      ECIS06  DJ6J-000
      FUNCTION DJ6J(J1,J2,J3,J4,J5,J6,FAC,NFA)                          DJ6J-001
C COMPUTATION OF SIX-J COEFFICIENTS.                                    DJ6J-002
C INPUT:     J1..J6:  DOUBLE VALUE OF ANGULAR MOMENTA.                  DJ6J-003
C            FAC:     TABLE OF LOGARITHM OF FACTORIALS.                 DJ6J-004
C            NFA:     LENGTH OF FAC.                                    DJ6J-005
C OUTPUT:                                                               DJ6J-006
C                                  ( J1  J2  J3 )                       DJ6J-007
C            6-J COEFFICIENTS      )            (                       DJ6J-008
C                                  ( J4  J5  J6 )                       DJ6J-009
C                                                                       DJ6J-010
C***********************************************************************DJ6J-011
      IMPLICIT REAL*8(A-F)                                              DJ6J-012
      DIMENSION IX(6),IA(3,4),IY(4),FAC(*)                              DJ6J-013
      COMMON /INOUT/ MR,MW,MS                                           DJ6J-014
      DATA IA /1,2,3,1,5,6,4,2,6,4,5,3/                                 DJ6J-015
      DJ6J=0.D0                                                         DJ6J-016
      IX(1)=J1                                                          DJ6J-017
      IX(2)=J2                                                          DJ6J-018
      IX(3)=J3                                                          DJ6J-019
      IX(4)=J4                                                          DJ6J-020
      IX(5)=J5                                                          DJ6J-021
      IX(6)=J6                                                          DJ6J-022
C THE QUANTUM NUMBERS MULTIPLIED BY 2 ARE IN THE TABLE IX.              DJ6J-023
C SEARCH FOR A ZERO QUANTUM NUMBER.                                     DJ6J-024
      DO 1 I=1,6                                                        DJ6J-025
      IF (IX(I).LT.0) GO TO 14                                          DJ6J-026
      IF (IX(I).EQ.0) GO TO ( 5 , 6 , 7 , 8 , 9, 10 ),I                 DJ6J-027
    1 CONTINUE                                                          DJ6J-028
C GENERAL CASE.                                                         DJ6J-029
C CHECK OF THE TRIANGULAR RELATIONS AND COMPUTATION OF DELTA.           DJ6J-030
      DO 2 K=1,4                                                        DJ6J-031
      IZ1=IA(1,K)                                                       DJ6J-032
      IZ2=IA(2,K)                                                       DJ6J-033
      IZ3=IA(3,K)                                                       DJ6J-034
      IX1=IX(IZ1)                                                       DJ6J-035
      IX2=IX(IZ2)                                                       DJ6J-036
      IX3=IX(IZ3)                                                       DJ6J-037
      N=IX1+IX2+IX3+2                                                   DJ6J-038
      I1=N/2                                                            DJ6J-039
      IF (2*I1.NE.N) GO TO 12                                           DJ6J-040
      IF (I1.GT.NFA) GO TO 13                                           DJ6J-041
      N1=I1-IX3                                                         DJ6J-042
      N2=I1-IX2                                                         DJ6J-043
      N3=I1-IX1                                                         DJ6J-044
      IF ((I1.LE.0).OR.(N1.LE.0).OR.(N2.LE.0).OR.(N3.LE.0)) GO TO 15    DJ6J-045
      IY(K)=I1+1                                                        DJ6J-046
    2 DJ6J=DJ6J+FAC(N1)-FAC(I1+1)+FAC(N2)+FAC(N3)                       DJ6J-047
      N1=(IX(1)+IX(2)+IX(4)+IX(5))/2                                    DJ6J-048
      N2=(IX(1)+IX(3)+IX(4)+IX(6))/2                                    DJ6J-049
      N3=(IX(2)+IX(3)+IX(5)+IX(6))/2                                    DJ6J-050
C K1 AND K2 ARE THE LIMITS OF THE SUM.                                  DJ6J-051
C K1, L1, L2, L3, L4, M1, M2, M3 ARE THE FACTORIALS OF THE FIRST TERM.  DJ6J-052
      K1=MAX0(IY(1),IY(2),IY(3),IY(4))                                  DJ6J-053
      K2=MIN0(N1,N2,N3)+2                                               DJ6J-054
      L1=K1-IY(1)+1                                                     DJ6J-055
      L2=K1-IY(2)+1                                                     DJ6J-056
      L3=K1-IY(3)+1                                                     DJ6J-057
      L4=K1-IY(4)+1                                                     DJ6J-058
      M1=N1-K1+3                                                        DJ6J-059
      M2=N2-K1+3                                                        DJ6J-060
      M3=N3-K1+3                                                        DJ6J-061
      DJ6J=DEXP(.5D0*DJ6J+FAC(K1)-FAC(L1)-FAC(L2)-FAC(L3)-FAC(L4)-FAC(M1DJ6J-062
     1)-FAC(M2)-FAC(M3))                                                DJ6J-063
      IF (2*(K1/2).NE.K1) DJ6J=-DJ6J                                    DJ6J-064
      IF (K2.EQ.K1) GO TO 4                                             DJ6J-065
      A2=DJ6J                                                           DJ6J-066
      K=K2-K1                                                           DJ6J-067
    3 A1=DFLOAT((M1-K)*(M2-K)*(M3-K))                                   DJ6J-068
      K=K-1                                                             DJ6J-069
      A3=DFLOAT((L1+K)*(L2+K)*(L3+K)*(L4+K))                            DJ6J-070
      DJ6J=A2-DFLOAT(K1+K)*A1*DJ6J/A3                                   DJ6J-071
      IF (K.GT.0) GO TO 3                                               DJ6J-072
    4 RETURN                                                            DJ6J-073
C ONE QUANTUM NUMBER IS ZERO; VERIFICATION OF THE TRIANGULAR RELATION.  DJ6J-074
    5 IF (IX(2).NE.IX(3).OR.IX(5).NE.IX(6)) GO TO 15                    DJ6J-075
      IX(6)=IX(4)                                                       DJ6J-076
      IX(1)=IX(2)                                                       DJ6J-077
      IX(4)=IX(5)                                                       DJ6J-078
      GO TO 11                                                          DJ6J-079
    6 IF (IX(1).NE.IX(3).OR.IX(4).NE.IX(6)) GO TO 15                    DJ6J-080
      IX(6)=IX(5)                                                       DJ6J-081
      GO TO 11                                                          DJ6J-082
    7 IF (IX(1).NE.IX(2).OR.IX(4).NE.IX(5)) GO TO 15                    DJ6J-083
      IF (IX(4).NE.IX(5)) GO TO 15                                      DJ6J-084
      GO TO 11                                                          DJ6J-085
    8 IF (IX(2).NE.IX(6).OR.IX(3).NE.IX(5)) GO TO 15                    DJ6J-086
      IX(6)=IX(1)                                                       DJ6J-087
      IX(1)=IX(5)                                                       DJ6J-088
      IX(4)=IX(2)                                                       DJ6J-089
      GO TO 11                                                          DJ6J-090
    9 IF (IX(1).NE.IX(6).OR.IX(3).NE.IX(4)) GO TO 15                    DJ6J-091
      IX(6)=IX(2)                                                       DJ6J-092
      GO TO 11                                                          DJ6J-093
   10 IF (IX(1).NE.IX(5).OR.IX(2).NE.IX(4)) GO TO 15                    DJ6J-094
      IX(6)=IX(3)                                                       DJ6J-095
C VALUE OF 6-J COEFFICIENT WITH ONE QUANTUM NUMBER ZERO.                DJ6J-096
   11 IF (MIN0(IX(1),IX(4),IX(6)).LT.0) GO TO 14                        DJ6J-097
      IF (IX(6).GT.IX(1)+IX(4).OR.IX(6).LT.IABS(IX(1)-IX(4))) GO TO 15  DJ6J-098
      IF (MAX0(IX(1),IX(4)).GT.NFA) GO TO 13                            DJ6J-099
      K=IX(1)+IX(4)+IX(6)                                               DJ6J-100
      N=K/2                                                             DJ6J-101
      IF (2*N.NE.K) GO TO 12                                            DJ6J-102
      DJ6J=DFLOAT(1-2*MOD(N,2))/DSQRT(DFLOAT((IX(1)+1)*(IX(4)+1)))      DJ6J-103
      RETURN                                                            DJ6J-104
   12 WRITE (MW,1000)                                                   DJ6J-105
      GO TO 15                                                          DJ6J-106
   13 WRITE (MW,1001)                                                   DJ6J-107
      GO TO 15                                                          DJ6J-108
   14 WRITE (MW,1002)                                                   DJ6J-109
   15 DJ6J=0.D0                                                         DJ6J-110
      RETURN                                                            DJ6J-111
 1000 FORMAT (' VIOLATION OF THE INTEGER/HALF-INTEGER RULE BETWEEN QUANTDJ6J-112
     1UM NUMBERS IN DJ6J.')                                             DJ6J-113
 1001 FORMAT (' FACTORIAL TOO LARGE FOR DJ6J.')                         DJ6J-114
 1002 FORMAT (' NEGATIVE ANGULAR MOMENTUM IN DJ6J.')                    DJ6J-115
      END                                                               DJ6J-116
