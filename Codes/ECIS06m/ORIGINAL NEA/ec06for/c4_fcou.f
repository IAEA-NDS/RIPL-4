C 19/11/05                                                      ECIS06  FCOU-000
      SUBROUTINE FCOU(L,ETA,RO,F,FP,G,GP,IEXP,SIGMA)                    FCOU-001
C COULOMB FUNCTIONS FOR RHO > 0 AND -500 < ETA < 500  FROM - DFCOUL     FCOU-002
C FOR THE MATHEMATICAL DESCRIPTION SEE CHR. BARDIN ET. AL. CEA-N-906    FCOU-003
C OR BARDIN ET AL. COMP. PHYSICS COMM. VOL 3 (1972) PAGES 73-87         FCOU-004
C NUMERICAL PRECISION: AT LEAST 8 SIGNIFICANT FIGURES, EXCEPT FOR       FCOU-005
C THE REGION -8 < ETA < -6 AND RO < 125/6                               FCOU-006
C THE COMPUTATION OF THE PHASE SHIFTS HAS BEEN REMOVED EXCEPT FOR L=0   FCOU-007
C INPUT:     L:       FINAL ANGULAR MOMENTUM REQUIRED.                  FCOU-008
C            ETA:     COULOMB PARAMETER.                                FCOU-009
C            RO:      RADIUS.                                           FCOU-010
C OUTPUT:    F:       REGULAR FUNCTION.                                 FCOU-011
C            FP:      DERIVATIVE OF THE REGULAR FUNCTION.               FCOU-012
C            G:       IRREGULAR FUNCTION.                               FCOU-013
C            GP:      DERIVATIVE OF THE IRREGULAR FUNCTION.             FCOU-014
C            IEXP:    POWERS OF 10 BY WHICH TO MULTIPLY THE IRREGULAR   FCOU-015
C                     FUNCTIONS AND DIVIDE THE REGULAR , IF THIS POWER  FCOU-016
C                     EXCEEDS 15.                                       FCOU-017
C            SIGMA:   COULOMB PHASE-SHIFT FOR L=0.                      FCOU-018
C                                                                       FCOU-019
C SUBROUTINES FCOU, FCZ0, YFRI, YFCL, YFAS, YFIR AND FUNCTIONS SIGM AND FCOU-020
C PSI WILL BE KEPT THE SAME FOR DWBA AND ECIS.                          FCOU-021
C***********************************************************************FCOU-022
      IMPLICIT REAL*8 (A-H,O-Z)                                         FCOU-023
      DIMENSION F(*),FP(*),G(*),GP(*),IEXP(*),SIGMA(*)                  FCOU-024
      COMMON /INOUT/ MR,MW,MS                                           FCOU-025
      IF (RO.GT.0.D0.AND.DABS(ETA).LE.500.D0) GO TO 1                   FCOU-026
      WRITE (MW,1000) ETA,RO                                            FCOU-027
      STOP                                                              FCOU-028
C COMPUTATION OF THE COULOMB FUNCTIONS FOR L=0.                         FCOU-029
    1 CALL FCZ0(ETA,RO,F1,FP1,G1,GP1,IEXP1,SIGMA1)                      FCOU-030
      F(1)=F1                                                           FCOU-031
      FP(1)=FP1                                                         FCOU-032
      G(1)=G1                                                           FCOU-033
      GP(1)=GP1                                                         FCOU-034
      IEXP(1)=IEXP1                                                     FCOU-035
      SIGMA(1)=SIGMA1                                                   FCOU-036
      IF (L.LE.0) RETURN                                                FCOU-037
      LINF=0                                                            FCOU-038
      LIN=1                                                             FCOU-039
      IND=0                                                             FCOU-040
      L1=L+1                                                            FCOU-041
      ETAC=ETA*ETA                                                      FCOU-042
      IF ((ETA.GT.0.D0.AND.RO.LT.ETA+ETA.OR.RO.LT.ETA+DSQRT(ETAC+1.D0)))FCOU-043
     1 GO TO 7                                                          FCOU-044
      IF (RO.GE.ETA+DSQRT(ETAC+DFLOAT(L*(L+1)))) GO TO 5                FCOU-045
    2 ROINF=ETA+DSQRT(ETAC+DFLOAT(LINF*(LINF+1)))                       FCOU-046
      IF (RO.LT.ROINF) GO TO 3                                          FCOU-047
      IF (LINF.GE.L) GO TO 4                                            FCOU-048
      LINF=LINF+1                                                       FCOU-049
      GO TO 2                                                           FCOU-050
    3 IND=1                                                             FCOU-051
    4 LIN=LINF+1                                                        FCOU-052
    5 XM=1.D0                                                           FCOU-053
      IF (IND.EQ.0) LIN=L1                                              FCOU-054
C UPWARD RECURSION FOR THE REGULAR AND IRREGULAR FUNCTIONS.             FCOU-055
      DO 6 J=2,LIN                                                      FCOU-056
      ZIG=(DSQRT(ETAC+XM*XM))/XM                                        FCOU-057
      ZAG=ETA/XM+XM/RO                                                  FCOU-058
      F(J)=(ZAG*F(J-1)-FP(J-1))/ZIG                                     FCOU-059
      FP(J)=ZIG*F(J-1)-ZAG*F(J)                                         FCOU-060
      G(J)=(ZAG*G(J-1)-GP(J-1))/ZIG                                     FCOU-061
      GP(J)=ZIG*G(J-1)-ZAG*G(J)                                         FCOU-062
      IEXP(J)=IEXP(1)                                                   FCOU-063
    6 XM=XM+1.D0                                                        FCOU-064
      IF (IND.EQ.0) RETURN                                              FCOU-065
C DESCENDING RECURSION FOR THE REGULAR FUNCTION.                        FCOU-066
    7 FTEST=F(LIN)                                                      FCOU-067
      FPTEST=FP(LIN)                                                    FCOU-068
      LMAX=LINF+25+IDINT(5.D0*DABS(ETA))                                FCOU-069
      IF (LMAX.LT.L) LMAX=L                                             FCOU-070
      FI=1.D0                                                           FCOU-071
      FPI=1.D0                                                          FCOU-072
C ANGULAR MOMENTUM GREATER THAN THE MAXIMUM REQUIRED.                   FCOU-073
    8 XM=DFLOAT(LMAX+1)                                                 FCOU-074
      ZIG=(DSQRT(ETAC+XM*XM))/XM                                        FCOU-075
      ZAG=ETA/XM+XM/RO                                                  FCOU-076
      FL=(ZAG*FI+FPI)/ZIG                                               FCOU-077
      FPL=ZAG*FL-ZIG*FI                                                 FCOU-078
      IF (DABS(FL).LT.1.D15.AND.DABS(FPL).LT.1.D15) GO TO 9             FCOU-079
      FL=FL*1.D-15                                                      FCOU-080
      FPL=FPL*1.D-15                                                    FCOU-081
    9 FI=FL                                                             FCOU-082
      FPI=FPL                                                           FCOU-083
      IF (LMAX.LE.L) GO TO 11                                           FCOU-084
   10 LMAX=LMAX-1                                                       FCOU-085
      GO TO 8                                                           FCOU-086
C ANGULAR MOMENTUM SMALLER THAN THAT REQUIRED.                          FCOU-087
   11 F(LMAX+1)=FL                                                      FCOU-088
      FP(LMAX+1)=FPL                                                    FCOU-089
      IF (LMAX.GT.LINF) GO TO 10                                        FCOU-090
      FACT=FTEST/F(LIN)                                                 FCOU-091
      FACTP=FPTEST/FP(LIN)                                              FCOU-092
      INDICE=IEXP(1)/15                                                 FCOU-093
      XM=DFLOAT(LINF)                                                   FCOU-094
C NORMALISATION OF THE RESULTS OF THE DESCENDING RECURSION AND          FCOU-095
C UPWARDS RECURSION FOR THE IRREGULAR FUNCTION.                         FCOU-096
      DO 13 J=LIN,L1                                                    FCOU-097
      F(J)=F(J)*FACT                                                    FCOU-098
      FP(J)=FP(J)*FACTP                                                 FCOU-099
      IF (J.EQ.1) GO TO 13                                              FCOU-100
      ZIG=(DSQRT(ETAC+XM*XM))/XM                                        FCOU-101
      ZAG=ETA/XM+XM/RO                                                  FCOU-102
      G(J)=(ZAG*G(J-1)-GP(J-1))/ZIG                                     FCOU-103
      GP(J)=ZIG*G(J-1)-ZAG*G(J)                                         FCOU-104
      IF (DABS(G(J)).LT.1.D15.AND.DABS(GP(J)).LT.1.D15) GO TO 12        FCOU-105
      G(J)=G(J)/1.D15                                                   FCOU-106
      GP(J)=GP(J)/1.D15                                                 FCOU-107
      INDICE=INDICE+1                                                   FCOU-108
   12 IEXP(J)=INDICE*15                                                 FCOU-109
      A=DLOG10(DABS(FP(J)))+DLOG10(DABS(G(J)))                          FCOU-110
      B=0.D0                                                            FCOU-111
      IF (A.GE.0.D0) B=1.D0                                             FCOU-112
      I1=IDINT(B+A)                                                     FCOU-113
      I2=IDINT(B+DLOG10(DABS(GP(J)))+DLOG10(DABS(F(J))))                FCOU-114
      F(J)=F(J)*1.D1**(-I2)                                             FCOU-115
      FP(J)=FP(J)*1.D1**(-I1)                                           FCOU-116
   13 XM=XM+1.D0                                                        FCOU-117
      RETURN                                                            FCOU-118
 1000 FORMAT ('  FCOU   ***  ETA =',1P,D13.5,',  RHO =',D13.5,'   ARGUMEFCOU-119
     1NT OUT OFF RANGE.')                                               FCOU-120
      END                                                               FCOU-121
