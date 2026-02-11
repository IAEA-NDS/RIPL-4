C 12/01/07                                                      ECIS06  DVWI-000
      FUNCTION DVWI(E,EO,NV,B,CS,CR)                                    DVWI-001
C CONTRIBUTION TO THE REAL POTENTIAL OF A DISPERSIVE IMAGINARY TERM     DVWI-002
C '(E-EO)**NV/((E-EO)**NV+BV**NV)' WITH A DAMPING FACTOR                DVWI-003
C 'EXP[-CS*|E-EO|-CR*(E-EO)]'. WITH 'EY=E-E0' FOR POSITIVE VALUES AND   DVWI-004
C 'EY=-E-E0' FOR NEGATIVE VALUES, THE POLES ARE:                        DVWI-005
C 'EY+B*DEXP( I *PI*(2*J-1)/NV)' FOR J=1 TO NV AND THEIR RESIDUES ARE:  DVWI-006
C '2 B*DEXP( I *PI*(2*J-1)/NV)/NV' EVENTUALLY MULTIPLIED BY THE DAMPING DVWI-007
C FACTOR.                                                               DVWI-008
C INPUT:     E:       ENERGY MINUS THE FERMI ENERGY.                    DVWI-009
C            EO:      THRESHOLD ENERGY MINUS THE FERMI ENERGY.          DVWI-010
C            NV:      POWER IN THE EXPRESSION GIVEN ABOVE.              DVWI-011
C            B:       CONSTANT IN THE EXPRESSION GIVEN ABOVE.           DVWI-012
C            CS:      EXPONENTIAL DECREASE IN |E-EF|.                   DVWI-013
C            CR:      EXPONENTIAL DECREASE IN (E-EF).                   DVWI-014
C***********************************************************************DVWI-015
      IMPLICIT REAL*8 (A-H,O-Z)                                         DVWI-016
      DATA PI /3.1415926535897932D0/                                    DVWI-017
      M=NV/2                                                            DVWI-018
      IF ((CR.NE.0.D0).OR.(CS.NE.0.D0)) GO TO 2                         DVWI-019
C WITHOUT DAMPING. NO INTEGRAL BETWEEN -E-EO AND E-EO.                  DVWI-020
      DVWI=(DLOG(DABS(E+EO)/B)*(E+EO)**NV/((E+EO)**NV+B**NV)-DLOG(DABS(EDVWI-021
     1-EO)/B)*(E-EO)**NV/((E-EO)**NV+B**NV))                            DVWI-022
C LOOP ON POLES IN THE COMPLEX UPPER PLANE.                             DVWI-023
      DO 1 J=1,M                                                        DVWI-024
      BT=PI*DFLOAT(2*J-1)/DFLOAT(NV)                                    DVWI-025
      H=((E-EO)**2+2.D0*B*(E-EO)*DCOS(BT)+B**2)*M/B                     DVWI-026
      G=((E+EO)**2-2.D0*B*(E+EO)*DCOS(BT)+B**2)*M/B                     DVWI-027
    1 DVWI=DVWI+((E-EO)/H+(E+EO)/G)*DSIN(BT)*BT                         DVWI-028
      GO TO 4                                                           DVWI-029
    2 H=DFLOAT(M)*(E-EO)**NV/((E-EO)**NV+B**NV)*DREI((CS+CR)*(EO-E))    DVWI-030
      G=-DFLOAT(M)*(E+EO)**NV/((E+EO)**NV+B**NV)*DREI((CS-CR)*(EO+E))   DVWI-031
C LOOP ON POLES IN THE COMPLEX UPPER PLANE.                             DVWI-032
      DO 3 J=1,M                                                        DVWI-033
      BT=PI*DFLOAT(2*J-1)/DFLOAT(NV)                                    DVWI-034
      AT=B*DCOS(BT)                                                     DVWI-035
      AZ=B*DSIN(BT)                                                     DVWI-036
      CALL DCEI((CS+CR)*AT,(CS+CR)*AZ,BR,BI)                            DVWI-037
      H=H+((BR*AT-BI*AZ)*(E-EO)+BR*B**2)/((E-EO+AT)**2+AZ**2)           DVWI-038
      CALL DCEI((CS-CR)*AT,(CS-CR)*AZ,BR,BI)                            DVWI-039
    3 G=G+((BR*AT-BI*AZ)*(E+EO)-BR*B**2)/((E+EO-AT)**2+AZ**2)           DVWI-040
      DVWI=(DEXP(-(CS+CR)*EO)*H+DEXP(-(CS-CR)*EO)*G)/DFLOAT(M)          DVWI-041
    4 DVWI=DVWI/PI                                                      DVWI-042
      RETURN                                                            DVWI-043
      END                                                               DVWI-044
