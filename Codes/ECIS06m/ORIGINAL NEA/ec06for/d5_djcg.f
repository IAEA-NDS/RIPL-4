C 13/04/06                                                      ECIS06  DJCG-000
      FUNCTION DJCG(J1,J2,J3,M1,M2,FAC,NFA)                             DJCG-001
C CLEBSCH-GORDAN COEFFICIENTS ( J1  J2  M1  M2 | J3  M1+M2 ).           DJCG-002
C INPUT:     J1,..M2: INTEGER DOUBLED VALUES.                           DJCG-003
C            FAC:     TABLE OF LOGARITHM OF FACTORIALS.                 DJCG-004
C            NFA:     LENGTH OF FAC.                                    DJCG-005
C***********************************************************************DJCG-006
      IMPLICIT REAL*8(A-F)                                              DJCG-007
      DIMENSION FAC(*)                                                  DJCG-008
      COMMON /INOUT/ MR,MW,MS                                           DJCG-009
      DJCG=0.D0                                                         DJCG-010
      M3=M1+M2                                                          DJCG-011
      IF (J1+J2+J3.GT.2*NFA) GO TO 18                                   DJCG-012
      IF ((J1.LT.0).OR.(J2.LT.0).OR.(J3.LT.0)) GO TO 16                 DJCG-013
      IY1=J1+M1+2                                                       DJCG-014
      IX1=IY1/2                                                         DJCG-015
      IY2=J2+M2+2                                                       DJCG-016
      IX2=IY2/2                                                         DJCG-017
      IY3=J3-M3+2                                                       DJCG-018
      IX3=IY3/2                                                         DJCG-019
      IF ((IX1.LE.0).OR.(IX2.LE.0).OR.(IX3.LE.0)) RETURN                DJCG-020
      IF ((2*IX1.NE.IY1).OR.(2*IX2.NE.IY2).OR.(2*IX3.NE.IY3)) GO TO 17  DJCG-021
      IY1=IX1-M1                                                        DJCG-022
      IY2=IX2-M2                                                        DJCG-023
      IY3=IX3+M3                                                        DJCG-024
      IF ((IY1.LE.0).OR.(IY2.LE.0).OR.(IY3.LE.0)) RETURN                DJCG-025
C AT THIS POINT IX1,IX2,IX3 ARE TWICE J+M AND IY1,IY2,IY3 TWICE J-M     DJCG-026
C SEARCH FOR A ZERO ARGUMENT.                                           DJCG-027
      IF (J3.EQ.0) GO TO 4                                              DJCG-028
      IF ((J1.EQ.0).OR.(J2.EQ.0)) GO TO 5                               DJCG-029
      IF (M3.EQ.0) GO TO 6                                              DJCG-030
      IF (M2.EQ.0) GO TO 7                                              DJCG-031
      IF (M1.EQ.0) GO TO 8                                              DJCG-032
C GENERAL CASE.                                                         DJCG-033
    1 NX=IX1+IX2+IX3                                                    DJCG-034
      IZ1=NX-IX1-IY1                                                    DJCG-035
      IZ2=NX-IX2-IY2                                                    DJCG-036
      IZ3=NX-IX3-IY3                                                    DJCG-037
      IF ((IZ1.LE.0).OR.(IZ2.LE.0).OR.(IZ3.LE.0)) RETURN                DJCG-038
      NXY=NX-1                                                          DJCG-039
      I1=IX2-IY3                                                        DJCG-040
      I2=IY1-IX3                                                        DJCG-041
C K1 AND K2 ARE THE LIMITS OF THE SUM.                                  DJCG-042
C M1,M2,M3,K1,K3,K4 HAVE THEIR FACTORIALS IN THE DENOMINATOR            DJCG-043
C NXY=J1+J2+J3+1      IZ1,IZ2,IZ3  ARE J1+J2-J3, J1-J2+J3 AND -J1+J2+J3.DJCG-044
      K1=MAX0(I1,I2,0)+1                                                DJCG-045
      K2=MIN0(IY1,IX2,IZ3)                                              DJCG-046
      K3=K1-I1                                                          DJCG-047
      K4=K1-I2                                                          DJCG-048
      N1=IY1-K1+1                                                       DJCG-049
      N2=IX2-K1+1                                                       DJCG-050
      N3=IZ3-K1+1                                                       DJCG-051
      DJCG=DEXP(0.5D0*(FAC(IX3+IY3)-FAC(IX3+IY3-1)-FAC(NXY)+FAC(IZ1)+FACDJCG-052
     1(IZ2)+FAC(IZ3)+FAC(IX1)+FAC(IX2)+FAC(IX3)+FAC(IY1)+FAC(IY2)+FAC(IYDJCG-053
     23))-FAC(N1)-FAC(N2)-FAC(N3)-FAC(K1)-FAC(K3)-FAC(K4))              DJCG-054
      IF (2*(K1/2).EQ.K1) DJCG=-DJCG                                    DJCG-055
      IF (K1.EQ.K2) GO TO 3                                             DJCG-056
      A4=DJCG                                                           DJCG-057
      K=K2-K1                                                           DJCG-058
      K3=K2-I1                                                          DJCG-059
      K4=K2-I2                                                          DJCG-060
      N1=IY1-K2                                                         DJCG-061
      N2=IX2-K2                                                         DJCG-062
      N3=IZ3-K2                                                         DJCG-063
C K2,K3,K4,N1,N2,N3 ARE THE ARGUMENTS OF THE FACTORIALS IN THE LAST TERMDJCG-064
      DO 2 I=1,K                                                        DJCG-065
      A1=DFLOAT((K2-I)*(K3-I)*(K4-I))                                   DJCG-066
      A2=DFLOAT((N1+I)*(N2+I)*(N3+I))                                   DJCG-067
    2 DJCG=A4-DJCG*A2/A1                                                DJCG-068
    3 RETURN                                                            DJCG-069
C J1,J2 OR J3  IS ZERO.                                                 DJCG-070
    4 IF (J1.NE.J2) RETURN                                              DJCG-071
      A1=DFLOAT(J1+1)                                                   DJCG-072
      DJCG=1.D0/DSQRT(A1)                                               DJCG-073
      IF (MOD(IY1,2).EQ.0) DJCG=-DJCG                                   DJCG-074
      RETURN                                                            DJCG-075
    5 IF (J1+J2.NE.J3) RETURN                                           DJCG-076
      DJCG=1.D0                                                         DJCG-077
      RETURN                                                            DJCG-078
C M1,M2 OR M3 IS ZERO; IF THE OTHERS M ARE LARGER THAN 1/2,GENERAL CASE.DJCG-079
    6 IF (IABS(M1)-1) 9 , 10 , 1                                        DJCG-080
    7 IF (IABS(M1).GT.1) GO TO 1                                        DJCG-081
      GO TO 11                                                          DJCG-082
    8 IF (IABS(M3).GT.1) GO TO 1                                        DJCG-083
      GO TO 12                                                          DJCG-084
C ALL THE M ARE ZEROS.                                                  DJCG-085
    9 N5=IX1+IX2+IX3-1                                                  DJCG-086
      IF (2*(N5/2).NE.N5) RETURN                                        DJCG-087
      N2=IX1+IX2-IX3                                                    DJCG-088
      N3=IX2+IX3-IX1                                                    DJCG-089
      N4=IX3+IX1-IX2                                                    DJCG-090
      A1=FAC(2*IX3)-FAC(2*IX3-1)                                        DJCG-091
      KC=1                                                              DJCG-092
      GO TO 14                                                          DJCG-093
C ONE M IS ZERO AND THE OTHERS +-1/2  FORMULA OF DCGS.                  DJCG-094
   10 IQ=IX3-1                                                          DJCG-095
      IZ1=IX1+IX2-2                                                     DJCG-096
      IZ2=IX1-IY2                                                       DJCG-097
      LZ1=J1                                                            DJCG-098
      LZ2=J2                                                            DJCG-099
      LW2=M2                                                            DJCG-100
      LW1=IY1                                                           DJCG-101
      GO TO 13                                                          DJCG-102
   11 IQ=IX2-1                                                          DJCG-103
      IZ1=IX1+IX3-2                                                     DJCG-104
      IZ2=IX3-IY1                                                       DJCG-105
      LZ1=J1                                                            DJCG-106
      LZ2=J3                                                            DJCG-107
      LW2=M1                                                            DJCG-108
      LW1=IY3                                                           DJCG-109
      GO TO 13                                                          DJCG-110
   12 IQ=IX1-1                                                          DJCG-111
      IZ1=IX2+IX3-2                                                     DJCG-112
      IZ2=IX2-IY3                                                       DJCG-113
      LZ1=J2                                                            DJCG-114
      LZ2=J3                                                            DJCG-115
      LW2=-M3                                                           DJCG-116
      LW1=IY2                                                           DJCG-117
   13 N2=IZ1-IQ+1                                                       DJCG-118
      N3=IQ+IZ2+1                                                       DJCG-119
      N4=IQ-IZ2+1                                                       DJCG-120
      IF (N2.LT.1.OR.N3.LT.1.OR.N4.LT.1) RETURN                         DJCG-121
      N5=IZ1+IQ+2                                                       DJCG-122
      A1=FAC(J3+2)-FAC(J3+1)-FAC(LZ1+2)+FAC(LZ1+1)-FAC(LZ2+2)+FAC(LZ2+1)DJCG-123
      KC=2                                                              DJCG-124
C SIMPLE FORMULA.                                                       DJCG-125
   14 IF (N5-1.GT.NFA) GO TO 18                                         DJCG-126
      L1=(N5+1)/2                                                       DJCG-127
      L2=(N2+1)/2                                                       DJCG-128
      L3=(N3+1)/2                                                       DJCG-129
      L4=(N4+1)/2                                                       DJCG-130
      DJCG=DEXP(0.5D0*(A1+FAC(N2)+FAC(N3)+FAC(N4)-FAC(N5))+FAC(L1)-FAC(LDJCG-131
     12)-FAC(L3)-FAC(L4))                                               DJCG-132
      IF (KC.EQ.2) GO TO 15                                             DJCG-133
      IF (MOD(L1+IX1-IX2,2).EQ.0) DJCG=-DJCG                            DJCG-134
      RETURN                                                            DJCG-135
   15 IF (LW2.GT.0) L4=L4+N5+1                                          DJCG-136
      DJCG=2.D0*DJCG                                                    DJCG-137
      IF (MOD(L4+LW1+IX1-IY2,2).NE.0) DJCG=-DJCG                        DJCG-138
      RETURN                                                            DJCG-139
   16 WRITE (MW,1000)                                                   DJCG-140
      RETURN                                                            DJCG-141
   17 WRITE (MW,1001)                                                   DJCG-142
      RETURN                                                            DJCG-143
   18 WRITE (MW,1002)                                                   DJCG-144
      RETURN                                                            DJCG-145
 1000 FORMAT (' NEGATIVE ANGULAR MOMENTUM IN DJCG.')                    DJCG-146
 1001 FORMAT (' INTEGER/HALF-INTEGER RULE BETWEEN QUANTUM NUMBERS TRANSGDJCG-147
     1RESSED IN DJCG.')                                                 DJCG-148
 1002 FORMAT (' FACTORIAL TOO LARGE IN DJCG.')                          DJCG-149
      END                                                               DJCG-150
