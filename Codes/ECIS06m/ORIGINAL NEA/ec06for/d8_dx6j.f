C 27/06/06                                                      ECIS06  DX6J-000
      SUBROUTINE DX6J(AA,AT,J,JT)                                       DX6J-001
C RECURRENCE COMPUTATION OF UNNORMALISED 6-J COEFFICIENTS.              DX6J-002
C                       ( J(1)  J(2)  J(3) )                            DX6J-003
C 6-J COEFFICIENTS      )                  ( FOR ALL VALUES OF JJ.      DX6J-004
C                       ( J(4)  J(5)   JJ  )                            DX6J-005
C INPUT:     J:       INTEGER DOUBLED VALUES.                           DX6J-006
C            JT:      NUMBER OF 6J COEFFICIENTS PLUS ONE.               DX6J-007
C OUTPUT:    AA:      UNNORMALISED 6J COEFFICIENTS IN AA(2) TO AA(JT)   DX6J-008
C                     STARTING FROM THE LARGEST VALUE OF JJ; AA(1)=0.   DX6J-009
C            AT:      NORMALISATION: THE VALUES AA MUST BE DIVIDED BY   DX6J-010
C                     (-)**(J(1)+J(2)+J(4)+J(5))*SQRT((J(3)+1)*AT).     DX6J-011
C***********************************************************************DX6J-012
      IMPLICIT REAL*8 (A-H,O-Z)                                         DX6J-013
      DIMENSION AA(*),J(5)                                              DX6J-014
      AA(1)=0.D0                                                        DX6J-015
      AA(2)=1.D0                                                        DX6J-016
      IF (JT.LE.2) RETURN                                               DX6J-017
      AL=AT                                                             DX6J-018
      C2=0.D0                                                           DX6J-019
      BK1=DFLOAT(J(1)-J(5))**2                                          DX6J-020
      BK3=DFLOAT(J(2)-J(4))**2                                          DX6J-021
      BK2=DFLOAT(J(1)+J(5)+2)**2                                        DX6J-022
      BK4=DFLOAT(J(2)+J(4)+2)**2                                        DX6J-023
      D1=DFLOAT(J(1)-J(5))*DFLOAT(J(1)+J(5)+2)*DFLOAT(J(4)-J(2))*DFLOAT(DX6J-024
     1J(2)+J(4)+2)/16.D0                                                DX6J-025
      D2=(BK1+BK2+BK3+BK4-DFLOAT(4*J(3)*(J(3)+2)))/8.D0-1.D0            DX6J-026
      BK=(AL+1.D0)**2                                                   DX6J-027
      DO 2 I=3,JT                                                       DX6J-028
      C1=C2                                                             DX6J-029
      BK=BK-AL*4.D0                                                     DX6J-030
      C2=.03125D0*DSQRT((BK3-BK)*(BK1-BK)*(BK2-BK)*(BK4-BK))            DX6J-031
      D4=.5D0*(AL+.5D0*BK-1.D0)                                         DX6J-032
      AA(I)=-(AL*(D1+(D2-D4)*D4)*AA(I-1)+(AL-1.D0)*C1*AA(I-2))/(C2*(AL+1DX6J-033
     1.D0))                                                             DX6J-034
      AL=AL-2.D0                                                        DX6J-035
      AT=AT+AL*AA(I)*AA(I)                                              DX6J-036
      IF (AT.LT.1.D12) GO TO 2                                          DX6J-037
      AT=AT*1.D-24                                                      DX6J-038
      DO 1 II=2,I                                                       DX6J-039
    1 AA(II)=AA(II)*1.D-12                                              DX6J-040
    2 CONTINUE                                                          DX6J-041
      RETURN                                                            DX6J-042
      END                                                               DX6J-043
