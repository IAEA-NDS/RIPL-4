C 02/04/06                                                      ECIS06  COCN-000
      SUBROUTINE COCN(LA,LB,JA,JB,IW,IV,IJ,IP,RB,RC,NT)                 COCN-001
C COMPUTATION OF COEFFICIENT FOR ANGULAR DISTRIBUTION OF COMPOUND       COCN-002
C NUCLEUS CROSS-SECTIONS.                                               COCN-003
C INPUT:     LA,LB,JA,JB,IW,IV,IJ: INTEGER DOUBLE VALUES.               COCN-004
C            IP:      NUMBER OF COEFFICIENTS REQUESTED.                 COCN-005
C            NT:      SIZE OF THE WORKING SPACE RC.                     COCN-006
C OUTPUT:    RB:      IN RB(I), VALUE OF:                               COCN-007
C                     ( LL  LA  LB )  ( JA  JB  LL )  ( LA  LB  LL )    COCN-008
C                     (            )  )            (  )            (    COCN-009
C                     (  0   0   0 )  ( IJ  IJ  IW )  ( JB  JA  IV )    COCN-010
C                     *(2 IJ+1)*SQRT((2 LA+1)(2 LB+1)(2 JA+1)(2 JB+1))  COCN-011
C                     FOR EVEN VALUES OF LL, FROM LL=0 WITH POSITIVE    COCN-012
C                     OR NULL VALUE FOR LL=0.                           COCN-013
C WORKING AREA:                                                         COCN-014
C             RC:     FOR ALL THE NON-ZERO 6-J COEFFICIENTS WITH EVEN   COCN-015
C                     AND ODD VALUES OF LL, STARTING WITH THE LARGEST   COCN-016
C                     VALUE OF LL IN RC(2).                             COCN-017
C***********************************************************************COCN-018
      IMPLICIT REAL*8 (A-H,O-Z)                                         COCN-019
      COMMON /INOUT/ MR,MW,MS                                           COCN-020
      DIMENSION RB(*),RC(*),J(5)                                        COCN-021
      LD=IABS(LA-LB)                                                    COCN-022
      LP=(LA+LB)/4+1                                                    COCN-023
      LM=LD/4+1                                                         COCN-024
      AD=DFLOAT(LD**2)                                                  COCN-025
      AP=DFLOAT(LA+LB+2)**2                                             COCN-026
      JT=MAX0(IP,LP)                                                    COCN-027
      IF (JT.GT.NT) GO TO 7                                             COCN-028
      TI=0.D0                                                           COCN-029
C RECURRENCE COMPUTATION OF 3-J COEFFICIENTS FOR EVEN VALUES OF LL.     COCN-030
      DO 2 I=1,JT                                                       COCN-031
      AL=DFLOAT(4*I-4)                                                  COCN-032
      RC(I)=0.D0                                                        COCN-033
      IF (I.LT.LM.OR.I.GT.LP) GO TO 2                                   COCN-034
      IF (I.GT.LM) GO TO 1                                              COCN-035
      RC(I)=1.D0                                                        COCN-036
      GO TO 2                                                           COCN-037
    1 RC(I)=-RC(I-1)*DSQRT(((AL-2.D0)**2-AD)*(AP-(AL-2.D0)**2)/((AL**2-ACOCN-038
     1D)*(AP-AL**2)))                                                   COCN-039
    2 TI=TI+RC(I)**2*(AL+1.D0)                                          COCN-040
      TI=DSQRT(TI)                                                      COCN-041
      DO 3 I=1,IP                                                       COCN-042
    3 RB(I)=RC(I)                                                       COCN-043
C QUANTUM NUMBERS OF THE FIRST 6-J COEFFICIENT.                         COCN-044
      J(1)=LA                                                           COCN-045
      J(2)=JA                                                           COCN-046
      J(3)=IV                                                           COCN-047
      J(4)=JB                                                           COCN-048
      J(5)=LB                                                           COCN-049
      DO 5 K=1,2                                                        COCN-050
      JI=MAX0(IABS(J(1)-J(5)),IABS(J(2)-J(4)))                          COCN-051
      JF=MIN0(J(1)+J(5),J(2)+J(4))                                      COCN-052
      JT=(JF-JI)/2+2                                                    COCN-053
      IF (JT.GT.NT) GO TO 7                                             COCN-054
      AT=DFLOAT(JF+1)                                                   COCN-055
      CALL DX6J(RC,AT,J,JT)                                             COCN-056
      TI=TI*DSQRT(AT*DFLOAT(J(3)+1))                                    COCN-057
      DO 4 I=1,IP                                                       COCN-058
      LL=4*I-4                                                          COCN-059
      LK=1                                                              COCN-060
      IF (LL.GE.JI.AND.LL.LE.JF) LK=2+(JF-LL)/2                         COCN-061
      RB(I)=RC(LK)*RB(I)                                                COCN-062
    4 CONTINUE                                                          COCN-063
C QUANTUM NUMBERS OF THE SECOND 6-J COEFFICIENT.                        COCN-064
      J(1)=IJ                                                           COCN-065
      J(3)=IW                                                           COCN-066
    5 J(5)=IJ                                                           COCN-067
C NORMALISATION.                                                        COCN-068
      TI=DSQRT(DFLOAT(LA+1)*DFLOAT(LB+1)*DFLOAT(JA+1)*DFLOAT(JB+1))*DFLOCOCN-069
     1AT(IJ+1)/TI                                                       COCN-070
      IF (RB(1).LT.0.D0) TI=-TI                                         COCN-071
      DO 6 I=1,IP                                                       COCN-072
    6 RB(I)=RB(I)*TI                                                    COCN-073
      RETURN                                                            COCN-074
    7 WRITE (MW,1000) NT,JT                                             COCN-075
      STOP                                                              COCN-076
 1000 FORMAT (' WORKING SPACE TOO SMALL IN COCN:',I5,' AVAILABLE',I6,' RCOCN-077
     1EQUESTED.'///' IN COCN  ...  STOP  ...')                          COCN-078
      END                                                               COCN-079
