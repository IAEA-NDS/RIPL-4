C 27/06/06                                                      ECIS06  COPO-000
      SUBROUTINE COPO(W,V,Q,ISM,H,L,VAL,VAC,CCZ,ZT,LT,LZ)               COPO-001
C COMPUTES THE COULOMB POTENTIAL OF THE CHARGE DISTRIBUTION STORED IN V.COPO-002
C INPUT:     V:       CHARGE DISTRIBUTION UNNORMALISED.                 COPO-003
C            H:       STEP SIZE.                                        COPO-004
C            ISM:     NUMBER OF POINTS.                                 COPO-005
C            L:       ANGULAR MOMENTUM, RETURN 0 IF L IS NEGATIVE.      COPO-006
C            VAL:     PRODUCT OF CHARGES IS IN VAL(1).                  COPO-007
C            VAC:     PARAMETER OF CHARGE DISTRIBUTION.                 COPO-008
C            LT:      LOGICAL  .TRUE.  FORM FACTOR ALREADY NORMALISED   COPO-009
C                    .FALSE.  NORMALISATION STORED FROM THE CENTRAL.    COPO-010
C            LZ:      LOGICAL  .TRUE. TO USE THE NORMALISATION STORED   COPO-011
C                     EVEN IF L=0.                                      COPO-012
C OUTPUT:    W:       FORM FACTOR WHICH CAN BE AT THE SAME PLACE AS V.  COPO-013
C            ZT:      NORMALISATION FACTOR.                             COPO-014
C WORKING AREA:                                                         COPO-015
C            Q:       WHICH CAN BE AT THE SAME PLACE AS W.              COPO-016
C***********************************************************************COPO-017
      IMPLICIT REAL*8 (A-H,O-Z)                                         COPO-018
      DIMENSION W(*),V(*),Q(ISM,5),VAL(*)                               COPO-019
      LOGICAL LT,LZ                                                     COPO-020
      IF (L.GE.0) GO TO 2                                               COPO-021
      DO 1 IS=1,ISM                                                     COPO-022
    1 V(IS)=0.D0                                                        COPO-023
      RETURN                                                            COPO-024
    2 RR=0.D0                                                           COPO-025
      DO 3 IS=1,ISM                                                     COPO-026
      RR=RR+H                                                           COPO-027
      V(IS)=V(IS)*(1.D0+VAC*RR*RR)                                      COPO-028
      Q(IS,4)=RR**L                                                     COPO-029
      Q(IS,3)=Q(IS,4)/RR                                                COPO-030
      Q(IS,5)=Q(IS,4)*RR                                                COPO-031
    3 Q(IS,1)=Q(IS,5)*RR                                                COPO-032
      V(ISM)=0.D0                                                       COPO-033
      V(ISM-1)=0.D0                                                     COPO-034
      Q(1,1)=Q(1,1)*V(1)                                                COPO-035
      Q(ISM,2)=0.D0                                                     COPO-036
      DO 4 IS=2,ISM                                                     COPO-037
      JS=ISM+1-IS                                                       COPO-038
      Q(IS,1)=Q(IS-1,1)+Q(IS,1)*V(IS)                                   COPO-039
    4 Q(JS,2)=Q(JS+1,2)+V(JS+1)/Q(JS+1,3)                               COPO-040
      C=DFLOAT(2*L+1)                                                   COPO-041
      ZZ=1.D0/C                                                         COPO-042
      IF (LT) GO TO 6                                                   COPO-043
      IF (L.NE.0.OR.LZ) GO TO 5                                         COPO-044
      ZT=CCZ*VAL(1)/Q(ISM,1)                                            COPO-045
    5 ZZ=ZT*ZZ                                                          COPO-046
    6 C=C*H/12.D0                                                       COPO-047
      DO 7 IS=1,ISM                                                     COPO-048
    7 W(IS)=(Q(IS,1)/Q(IS,5)+Q(IS,2)*Q(IS,4)+C*V(IS))*ZZ                COPO-049
      RETURN                                                            COPO-050
      END                                                               COPO-051
