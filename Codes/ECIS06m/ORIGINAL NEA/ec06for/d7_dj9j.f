C 27/06/06                                                      ECIS06  DJ9J-000
      FUNCTION DJ9J(J1,J2,J3,J4,J5,J6,J7,J8,J9,AA,ITO)                  DJ9J-001
C COMPUTATION OF 9-J COEFFICIENTS.                                      DJ9J-002
C INPUT:     J1,..J9: INTEGER DOUBLED VALUES.                           DJ9J-003
C OUTPUT:                             | J1  J2  J3 |                    DJ9J-004
C                                     |            |                    DJ9J-005
C   9J(J1,J2,J3,J4,J5,J6,J7,J8,J9) =  | J4  J5  J6 |                    DJ9J-006
C                                     |            |                    DJ9J-007
C                                     | J7  J8  J9 |                    DJ9J-008
C                                                                       DJ9J-009
C WORKING AREA:                                                         DJ9J-010
C            AA:      FOR COMPUTATION OF UNNORMALISED 6-J COEFFICIENTS  DJ9J-011
C                     BY RECURRENCE.                                    DJ9J-012
C            ITO:     LENGTH OF AA.                                     DJ9J-013
C***********************************************************************DJ9J-014
      IMPLICIT REAL*8 (A-H,O-Z)                                         DJ9J-015
      DIMENSION AA(*),J(5),IX(9),JX(3,3),JS(3),JF(3),IA(5,3),IB(3)      DJ9J-016
      EQUIVALENCE (IX(1),JX(1,1)),(IB(1),IB1),(IB(2),IB2),(IB(3),IB3)   DJ9J-017
      COMMON /INOUT/ MR,MW,MS                                           DJ9J-018
      DATA IA /1,2,3,6,9,6,4,5,8,2,8,9,7,1,4/                           DJ9J-019
      DJ9J=0.D0                                                         DJ9J-020
      IX(1)=J1                                                          DJ9J-021
      IX(2)=J2                                                          DJ9J-022
      IX(3)=J3                                                          DJ9J-023
      IX(4)=J4                                                          DJ9J-024
      IX(5)=J5                                                          DJ9J-025
      IX(6)=J6                                                          DJ9J-026
      IX(7)=J7                                                          DJ9J-027
      IX(8)=J8                                                          DJ9J-028
      IX(9)=J9                                                          DJ9J-029
      DO 2 I=1,3                                                        DJ9J-030
      DO 1 K=1,3                                                        DJ9J-031
      IF (JX(K,I).LT.0) GO TO 13                                        DJ9J-032
    1 CONTINUE                                                          DJ9J-033
C CHECK OF TRIANGULAR RELATIONS.                                        DJ9J-034
      IF (MOD(JX(1,I)+JX(2,I)+JX(3,I),2).NE.0.OR.MOD(JX(I,1)+JX(I,2)+JX(DJ9J-035
     1I,3),2).NE.0) GO TO 11                                            DJ9J-036
      IF (IABS(JX(1,I)-JX(2,I)).GT.JX(3,I).OR.JX(1,I)+JX(2,I).LT.JX(3,I)DJ9J-037
     1) RETURN                                                          DJ9J-038
      IF (IABS(JX(I,1)-JX(I,2)).GT.JX(I,3).OR.JX(I,1)+JX(I,2).LT.JX(I,3)DJ9J-039
     1) RETURN                                                          DJ9J-040
    2 CONTINUE                                                          DJ9J-041
C SEARCH OF THE CONFIGURATION FOR WHICH THE SUM ON PRODUCTS OF 6-J      DJ9J-042
C COEFFICIENTS IS THE SMALLEST ONE.                                     DJ9J-043
      K1=MIN0(IX(3),IX(5),IX(7))                                        DJ9J-044
      K2=MIN0(IX(2),IX(4),IX(9))                                        DJ9J-045
      K3=MIN0(IX(1),IX(6),IX(8))                                        DJ9J-046
      KT=MAX0(K1,K2,K3)                                                 DJ9J-047
      IF (KT.EQ.K1) GO TO 6                                             DJ9J-048
      IF (K2.GT.K3) GO TO 4                                             DJ9J-049
      DO 3 I=1,3                                                        DJ9J-050
      II=IX(I)                                                          DJ9J-051
      IX(I)=IX(I+3)                                                     DJ9J-052
      IX(I+3)=IX(I+6)                                                   DJ9J-053
    3 IX(I+6)=II                                                        DJ9J-054
      GO TO 6                                                           DJ9J-055
    4 DO 5 I=1,3                                                        DJ9J-056
      II=IX(I+6)                                                        DJ9J-057
      IX(I+6)=IX(I+3)                                                   DJ9J-058
      IX(I+3)=IX(I)                                                     DJ9J-059
    5 IX(I)=II                                                          DJ9J-060
    6 DO 7 K=1,3                                                        DJ9J-061
      IA1=IA(1,K)                                                       DJ9J-062
      IA2=IA(2,K)                                                       DJ9J-063
      IA4=IA(4,K)                                                       DJ9J-064
      IA5=IA(5,K)                                                       DJ9J-065
      JS(K)=MAX0(IABS(IX(IA1)-IX(IA5)),IABS(IX(IA4)-IX(IA2)))+1         DJ9J-066
    7 JF(K)=MIN0(IX(IA1)+IX(IA5),IX(IA2)+IX(IA4))+1                     DJ9J-067
      KF=MIN0(JF(1),JF(2),JF(3))                                        DJ9J-068
      KT=2+(KF-MAX0(JS(1),JS(2),JS(3)))/2                               DJ9J-069
C LOOP ON THE THREE 6-J COEFFICIENTS.                                   DJ9J-070
      ITX=0                                                             DJ9J-071
      AZ=1.D0                                                           DJ9J-072
      DO 9 K=1,3                                                        DJ9J-073
      IJ=ITX                                                            DJ9J-074
      IB(K)=ITX+(JF(K)-KF)/2                                            DJ9J-075
      DO 8 I=1,5                                                        DJ9J-076
      IA1=IA(I,K)                                                       DJ9J-077
    8 J(I)=IX(IA1)                                                      DJ9J-078
      JT=2+(JF(K)-JS(K))/2                                              DJ9J-079
      ITX=ITX+JT                                                        DJ9J-080
      IF (ITX.GT.ITO) GO TO 12                                          DJ9J-081
      AT=DFLOAT(JF(K))                                                  DJ9J-082
      CALL DX6J(AA(1+IJ),AT,J,JT)                                       DJ9J-083
    9 AZ=AZ*AT*DFLOAT(J(3)+1)                                           DJ9J-084
C SUMMATION ON PRODUCTS OF 6-J.                                         DJ9J-085
      AF=DFLOAT(KF)                                                     DJ9J-086
      DO 10 I=2,KT                                                      DJ9J-087
      DJ9J=DJ9J+AF*AA(I+IB1)*AA(I+IB2)*AA(I+IB3)                        DJ9J-088
   10 AF=AF-2.D0                                                        DJ9J-089
C NORMALISATION.                                                        DJ9J-090
      DJ9J=DJ9J/DSQRT(AZ)                                               DJ9J-091
      RETURN                                                            DJ9J-092
   11 WRITE (MW,1000)                                                   DJ9J-093
      RETURN                                                            DJ9J-094
   12 WRITE (MW,1001)                                                   DJ9J-095
      RETURN                                                            DJ9J-096
   13 WRITE (MW,1002)                                                   DJ9J-097
      RETURN                                                            DJ9J-098
 1000 FORMAT (' INTEGER/HALF-INTEGER RULE BETWEEN QUANTUM NUMBERS TRANSGDJ9J-099
     1RESSED IN DJ9J.')                                                 DJ9J-100
 1001 FORMAT (' TOO MANY 6-J NEEDED IN DJ9J.')                          DJ9J-101
 1002 FORMAT (' NEGATIVE ANGULAR MOMENTUM IN DJ9J.')                    DJ9J-102
      END                                                               DJ9J-103
