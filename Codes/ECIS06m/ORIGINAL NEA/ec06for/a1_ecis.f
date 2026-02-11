C 29/01/07                                                      ECIS06  ECIS-000
C                                                                       ECIS-001
C THE COMMON /INOUT/ IS USED IN ANY SUBROUTINE WITH STANDARD INPUT OR   ECIS-002
C                            OUTPUT AND SAVE OR RESTART A SEARCH.       ECIS-003
C FOR THE COMMON  /DCONS/ SEE CALC.                                     ECIS-004
C                                                                       ECIS-005
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /INOUT/:                     ECIS-006
C  MR:        STANDARD INPUT FILE (USUAL VALUE 5).                      ECIS-007
C  MW:        STANDARD OUTPUT FILE (USUAL VALUE 6).                     ECIS-008
C  MS:        FILE TO SAVE OR RESTART A SEARCH IN REST (USUAL VALUE 8). ECIS-009
C                                                                       ECIS-010
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     ECIS-011
C  CM:        ATOMIC MASS IN MEV.                                       ECIS-012
C  CHB:       PLANCK CONSTANT X SPEED OF LIGHT /(2*PI) IN MEV*FERMI.    ECIS-013
C  CZ:        ELECTRIC CONSTANT.                                        ECIS-014
C   DEFINED:  CM,CHB,CZ                                                 ECIS-015
C                                                                       ECIS-016
C CONSTANTS COMPUTED FROM THE FUNDAMENTAL CONSTANTS, ATOMIC MASS, HBAR*CECIS-017
C AND ALPHA, AS GIVEN IN THE JOURNAL OF PHYSICS G, VOLUME 33, PAGE 97,  ECIS-018
C (JULY 2006) REFERRING FOR THESE VALUES TO THE 2002 CODATA SET WHICH   ECIS-019
C MAY BE FOUND AT HTTP://PHYSICS.NIST.GOV/CONSTANTS.                    ECIS-020
C     CM=931.494043 +/- 0.000080 MEV/C**2                               ECIS-021
C     CHB=197.326968 +/- 0.000017 MEV FM                                ECIS-022
C     CZ=137.03599911 +/- 0.00000046 WITHOUT DIMENSION                  ECIS-023
C                                                                       ECIS-024
C***********************************************************************ECIS-025
      IMPLICIT REAL*8 (A-H,O-Z)                                         ECIS-026
      PARAMETER (IDMX=600000)                                           ECIS-027
      CHARACTER*4 CW(2,IDMX)                                            ECIS-028
      DIMENSION NW(2,IDMX),DW(IDMX)                                     ECIS-029
      EQUIVALENCE (DW,NW,CW,W)                                          ECIS-030
      COMMON DW                                                         ECIS-031
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XE,XM,XN,XZ                   ECIS-032
      COMMON /INOUT/ MR,MW,MS                                           ECIS-033
      CM=931.494043D0                                                   ECIS-034
      CHB=197.326968D0                                                  ECIS-035
      CZ=137.03599911D0                                                 ECIS-036
      MR=5                                                              ECIS-037
      MW=6                                                              ECIS-038
      MS=8                                                              ECIS-039
      CALL CALC(NW,CW,DW,IDMX)                                          ECIS-040
      STOP                                                              ECIS-041
      END                                                               ECIS-042
