C 08/04/06                                                      ECIS06  FIT1-000
      SUBROUTINE FIT1(KE,I,W)                                           FIT1-001
C MINIMISATION OF A FUNCTION F(X) OF ONE VARIABLE X.                    FIT1-002
C INPUT:     AM:      TABLES OF MULTIPOLES.                             FIT1-003
C CALLING SEQUENCE                                                      FIT1-004
C    KE=0                                                               FIT1-005
C    I(2)=MAXIMUM NUMBER OF FUNCTION EVALUATIONS                        FIT1-006
C    W(1)=START VALUE OF X                                              FIT1-007
C    W(3)=FIRST STEP SIZE                                               FIT1-008
C    W(4)=ABSOLUTE SEARCH ACCURACY                                      FIT1-009
C    W(5)=RELATIVE SEARCH ACCURACY                                      FIT1-010
C  1 W(2)=FUNCTION VALUE F(X) AT X=W(1)                                 FIT1-011
C    CALL FIT1(KE,I,W)                                                  FIT1-012
C    IF(KE.EQ.1) GO TO 1                                                FIT1-013
C    XMIN=W(1)                                                          FIT1-014
C    FMIN=W(2)                                                          FIT1-015
C KE = ERROR CODE: KE=0 NO ERRORS, KE=                                  FIT1-016
C  2 MAXIMUM NUMBER OF FUNCTION EVALUATIONS                             FIT1-017
C  3 ROUNDING ERRORS, PROBABLY BECAUSE BOTH W(4) AND W(5) ARE TOO SMALL FIT1-018
C THE WORKING FIELDS I AND W HAVE THE LENGTH 3 AND 11 RESPECTIVELY      FIT1-019
C THEY CONTAIN ALL INFORMATION FOR THE CONTINUATION OF THE SEARCH       FIT1-020
C THEREFORE A SEARCH WITHIN ANOTHER SEARCH CAN BE DONE JUST CHANGING    FIT1-021
C THE WORKING FIELDS                                                    FIT1-022
C IF 2 FUNCTION VALUES F1 AND F2 ARE KNOWN FOR X = X1 AND X2 RESPECTIVE FIT1-023
C LY WITH X1 NE X2 ENTER THE CALLING SEQUENCE AFTER DEFINING :          FIT1-024
C KE = 1; I(1) = 3; W(6) = X1; W(7) = X2; W(9) = F1; W(10) = F2 AND     FIT1-025
C W(1) = USERS CHOICE                                                   FIT1-026
C WORKING FIELD VARIABLES:                                              FIT1-027
C I(1): CURRENT NUMBER OF FUNCTION EVALUATIONS                          FIT1-028
C I(2): MAXIMUM NUMBER OF FUNCTION EVALUATIONS                          FIT1-029
C I(3): MINIMUM POINTER, THE MINIMUM FUNCTION VALUE IS AT W(7+I(3))     FIT1-030
C W(1): CURRENT VALUE OF X                                              FIT1-031
C W(2): USER SUPPLIED FUNCTION VALUE                                    FIT1-032
C W(3): FIRST STEP SIZE                                                 FIT1-033
C W(4 AND 5): SEARCH ACCURACIES                                         FIT1-034
C W(6, 7 AND 8): X1, X2 AND X3 WITH X1 < X2 < X3                        FIT1-035
C W(9, 10 AND 11): FUNCTION VALUES AT X1, X2 AND X3 RESPECTIVELY        FIT1-036
C***********************************************************************FIT1-037
      IMPLICIT REAL*8 (A-H,O-Z)                                         FIT1-038
      DIMENSION I(3),W(11)                                              FIT1-039
      IF (KE.EQ.1) GO TO 2                                              FIT1-040
      KE=1                                                              FIT1-041
      I(1)=1                                                            FIT1-042
      I(3)=-1                                                           FIT1-043
      W(6)=W(1)                                                         FIT1-044
      W(9)=W(2)                                                         FIT1-045
    1 W(1)=W(1)+W(3)                                                    FIT1-046
      GO TO 11                                                          FIT1-047
    2 IF (I(1).GT.2) GO TO 3                                            FIT1-048
      I(3)=0                                                            FIT1-049
      W(7)=W(1)                                                         FIT1-050
      W(10)=W(2)                                                        FIT1-051
      IF (W(2).LE.W(9)) GO TO 1                                         FIT1-052
      I(3)=-1                                                           FIT1-053
      W(1)=W(6)-W(3)                                                    FIT1-054
      GO TO 11                                                          FIT1-055
    3 IF (I(1).GT.3) GO TO 5                                            FIT1-056
      W(8)=W(1)                                                         FIT1-057
      W(11)=W(2)                                                        FIT1-058
C ORDERING OF THE 3 FIRST VALUES OF X: W(6) < W(7) < W(8).              FIT1-059
      DO 4 J=1,3                                                        FIT1-060
      K=7-MOD(J,2)                                                      FIT1-061
      IF (W(K).LE.W(K+1)) GO TO 4                                       FIT1-062
      W(1)=W(K)                                                         FIT1-063
      W(K)=W(K+1)                                                       FIT1-064
      W(K+1)=W(1)                                                       FIT1-065
      K=K+3                                                             FIT1-066
      W(1)=W(K)                                                         FIT1-067
      W(K)=W(K+1)                                                       FIT1-068
      W(K+1)=W(1)                                                       FIT1-069
    4 CONTINUE                                                          FIT1-070
      I(3)=0                                                            FIT1-071
      IF (W(9).LT.W(10).AND.W(9).LT.W(11)) I(3)=-1                      FIT1-072
      IF (W(11).LT.W(10).AND.W(11).LT.W(9)) I(3)=1                      FIT1-073
      GO TO 9                                                           FIT1-074
C SORT IN THE NEW VALUES OF X AND F.                                    FIT1-075
    5 IF (I(3).EQ.0) GO TO 6                                            FIT1-076
      J=I(3)                                                            FIT1-077
      W(7-J)=W(7)                                                       FIT1-078
      W(10-J)=W(10)                                                     FIT1-079
      IF ((W(7+J)-W(1))*(W(1)-W(7)).GT.0.D0) GO TO 7                    FIT1-080
      W(7)=W(7+J)                                                       FIT1-081
      W(10)=W(10+J)                                                     FIT1-082
      W(7+J)=W(1)                                                       FIT1-083
      W(10+J)=W(2)                                                      FIT1-084
      IF (W(2).GE.W(10)) I(3)=0                                         FIT1-085
      GO TO 9                                                           FIT1-086
    6 J=-1                                                              FIT1-087
      IF (W(1).LT.W(7)) J=1                                             FIT1-088
      IF (W(2).GT.W(10)) GO TO 8                                        FIT1-089
      W(7+J)=W(7)                                                       FIT1-090
      W(10+J)=W(10)                                                     FIT1-091
    7 W(7)=W(1)                                                         FIT1-092
      W(10)=W(2)                                                        FIT1-093
      I31=I(3)+10                                                       FIT1-094
      IF (W(2).LE.W(I31)) I(3)=0                                        FIT1-095
      GO TO 9                                                           FIT1-096
    8 W(7-J)=W(1)                                                       FIT1-097
      W(10-J)=W(2)                                                      FIT1-098
    9 J=7+I(3)                                                          FIT1-099
C ERROR TESTS.                                                          FIT1-100
      IF (W(6).EQ.W(7).OR.W(7).EQ.W(8).OR.(W(9).EQ.W(10).AND.W(10).EQ.W(FIT1-101
     111))) GO TO 14                                                    FIT1-102
      IF (I(1).GE.I(2)) GO TO 15                                        FIT1-103
      IF (I(3).EQ.0) GO TO 10                                           FIT1-104
C STEP SIZE LIMITATION.                                                 FIT1-105
      W(1)=W(J)+(2*I(3))*DABS(W(6)-W(8))                                FIT1-106
      GO TO 11                                                          FIT1-107
C PREDICTION OF THE POSITION OF THE MINIMUM (PARABOLIC APPROXIMATION).  FIT1-108
   10 W(1)=((W(9)-W(10))/(W(6)-W(7))-(W(10)-W(11))/(W(7)-W(8)))/(W(6)-W(FIT1-109
     18))                                                               FIT1-110
      W(1)=.5D0*(W(6)+W(8)+(W(11)-W(9))/(W(1)*(W(6)-W(8))))             FIT1-111
C TEST OF CONVERGENCE.                                                  FIT1-112
      W(2)=DABS(W(1)-W(J))                                              FIT1-113
      IF (W(2).LT.DABS(W(4)).OR.W(2).LT.DABS(W(5)*W(J))) GO TO 12       FIT1-114
   11 I(1)=I(1)+1                                                       FIT1-115
      RETURN                                                            FIT1-116
   12 KE=0                                                              FIT1-117
   13 I37=I(3)+7                                                        FIT1-118
      W(1)=W(I37)                                                       FIT1-119
      I31=I(3)+10                                                       FIT1-120
      W(2)=W(I31)                                                       FIT1-121
      RETURN                                                            FIT1-122
   14 KE=KE+1                                                           FIT1-123
   15 KE=KE+1                                                           FIT1-124
      GO TO 13                                                          FIT1-125
      END                                                               FIT1-126
