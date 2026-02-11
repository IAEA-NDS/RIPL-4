C 20/08/06                                                      ECIS06  GGDR-000
      SUBROUTINE GGDR(IPI,WV,SCN,LO)                                    GGDR-001
C CALCULATIONS RELATED TO GIANT DIPOLE RESONANCE.                       GGDR-002
C INPUT:     IPI(3,*):MULTIPLICITY OF THE TARGET.                       GGDR-003
C            WV(J,*): MASS OF PARTICLE AND TARGET FOR J=1 AND 2,        GGDR-004
C                     CENTRE OF MASS ENERGY IN MEV FOR J=3.             GGDR-005
C            SCN:     DESCRIPTIONS OF LEVEL DENSITIES:                  GGDR-006
C                     1-SA  2-UX   3-TAU  4-SG   5-E0   6-EX   7-NZ     GGDR-007
C            LO(I):   LOGICAL CONTROLS:                                 GGDR-008
C               LO(81) =.TRUE. HAUSER-FESHBACH CORRECTIONS.             GGDR-009
C               LO(82) =.TRUE. OLD SIMPLIFIED COMPOUND NUCLEUS.         GGDR-010
C               LO(86) =.TRUE. GAMMA EMISSION IN COMPOUND NUCLEUS.      GGDR-011
C               LO(112)=.TRUE. SPIN-ORBIT OR COMPOUND NUCLEUS PARAMETERSGGDR-012
C                              ARE CHANGED IN SEARCH.                   GGDR-013
C               LO(115)=.TRUE. FIRST COMPUTATION FOR THIS ENERGY.       GGDR-014
C               LO(116)=.TRUE. NO OUTPUT.                               GGDR-015
C OUTPUT:    IN COMMON /NCOMP/.                                         GGDR-016
C                                                                       GGDR-017
C FOR THE COMMON  /NCOMP/ SEE CALX.                                     GGDR-018
C                                                                       GGDR-019
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /NCOMP/:                     GGDR-020
C  NRD:       NUMBER OF GAMMA TRANSMISSION COEFFICIENTS.                GGDR-021
C  NCONT:     NUMBER OF CONTINUUM FOR COMPOUND NUCLEUS.                 GGDR-022
C  BZ:        HAUSER-FESHBACH AND MOLDAUER'S PARAMETERS DESCRIBED BELOW.GGDR-023
C   BZ(1):    SQUARE ROOT OF ELASTIC ENHANCEMENT.                       GGDR-024
C   BZ(2):    IF LO(82)=.TRUE., SPIN CUT-OFF PARAMETER.                 GGDR-025
C             IF LO(82)=.FALSE., PARTICLE DEGREES OF FREEDOM.           GGDR-026
C   BZ(3):    SQUARE ROOT OF LEVEL DENSITY PARAMETER. IF LO(82)=LO(87)= GGDR-027
C             .FALSE., PARAMETER BZ(3) IN MOLDAUER'S FORMULA OF INPUT   GGDR-028
C             DESCRIPTION.                                              GGDR-029
C   BZ(4):    IF LO(82)=LO(87)=.FALSE., PARAMETER BZ(4) IN SAME FORMULA.GGDR-030
C   BZ(5):    IF LO(82)=LO(87)=.FALSE., PARAMETER BZ(5) IN SAME FORMULA.GGDR-031
C  TG0:       SLOW S-WAVE NEUTRON GAMMA WIDTH/SPACING FOR NORMALISATION.GGDR-032
C  BN:        NEUTRON SEPARATION ENERGY.                                GGDR-033
C  FNUG:      RADIATIVE DEGREES OF FREEDOM.                             GGDR-034
C  EGD:       ENERGY OF THE GIANT DIPOLE RESONANCE.                     GGDR-035
C  GGD:       RESONANCE WIDTH.                                          GGDR-036
C  TG1:       DERIVED DATA FOR GAMMA IN COMPOUND NUCLEUS.               GGDR-037
C  SGSQ:      DERIVED DATA FOR GAMMA IN COMPOUND NUCLEUS.               GGDR-038
C   DEFINED:  BZ,EGD,GGD,SGSQ.                                          GGDR-039
C   USED:     NRD,NCONT,BZ,TGO,BN,FNUG,EGD,GGD,TG1,SGSQ.                GGDR-040
C                                                                       GGDR-041
C***********************************************************************GGDR-042
      IMPLICIT REAL*8 (A-H,O-Z)                                         GGDR-043
      LOGICAL LO(150)                                                   GGDR-044
      DIMENSION IPI(11,*),WV(22,*),SCN(7)                               GGDR-045
      COMMON /NCOMP/ NSP(3),NFISS,NRD,NCONT,NCOJ,NCONS,NIE,NCOLX,NDP,NDQGGDR-046
     1,ACN1,ACN2,AZ(6),BZ(5),TG0,BN,FNUG,EGD,GGD,TG1,SGSQ               GGDR-047
      COMMON /INOUT/ MR,MW,MS                                           GGDR-048
      IF (.NOT.LO(81)) RETURN                                           GGDR-049
      BZ(2)=DABS(BZ(2))                                                 GGDR-050
      IF (LO(82)) BZ(3)=BZ(1)*BZ(1)-1.D0                                GGDR-051
      IF (.NOT.(LO(112).OR.LO(115))) RETURN                             GGDR-052
C  GAMMA RAY INPUT.                                                     GGDR-053
      IF ((.NOT.LO(86)).OR.(NRD.NE.0)) RETURN                           GGDR-054
      XFR=0.5D0                                                         GGDR-055
      NA=IDINT(WV(2,1)+WV(1,1)+.5D0)                                    GGDR-056
      AA=DFLOAT(NA)                                                     GGDR-057
      NZ=IDINT(SCN(7)+.1D0)                                             GGDR-058
      NN=NA-NZ                                                          GGDR-059
      ZZ=DFLOAT(NZ*NN)                                                  GGDR-060
      KGD=0                                                             GGDR-061
      IF (EGD.EQ.0.D0) EGD=163.D0*DSQRT(ZZ)/(AA**1.3333333D0)           GGDR-062
      IF (EGD.LT.0.D0) GO TO 1                                          GGDR-063
      KGD=1                                                             GGDR-064
      IF (GGD.LE.0.D0) GGD=5.D0                                         GGDR-065
    1 IF (.NOT.LO(116)) WRITE (MW,1000) NA,NZ,NN,BN,FNUG,SCN(4)         GGDR-066
      ATG0=DABS(TG0)                                                    GGDR-067
      IF (.NOT.LO(116)) WRITE (MW,1001) ATG0                            GGDR-068
      IF (KGD.EQ.0.AND..NOT.LO(116)) WRITE (MW,1002)                    GGDR-069
      IF (KGD.EQ.1.AND..NOT.LO(116)) WRITE (MW,1003) EGD,GGD,XFR        GGDR-070
      SGSQ=2.D0*SCN(4)**2                                               GGDR-071
      ROJ=0.D0                                                          GGDR-072
      ECM=0.D0                                                          GGDR-073
    2 TG=0.D0                                                           GGDR-074
      ELIM=BN+ECM                                                       GGDR-075
      ELO=0.D0                                                          GGDR-076
    3 EHI=ELO+0.05D0                                                    GGDR-077
      IF (EHI.GE.ELIM) EHI=ELIM                                         GGDR-078
      EPS=(EHI+ELO)/2.D0                                                GGDR-079
      ESQ=EPS*EPS                                                       GGDR-080
      TINT=(EHI-ELO)*EPS*ESQ                                            GGDR-081
      EEE=-EGD*EGD+ESQ                                                  GGDR-082
      IF (KGD.EQ.1) TINT=TINT*EPS/(ESQ*GGD*GGD+EEE*EEE)                 GGDR-083
      EXC=BN+ECM-EPS                                                    GGDR-084
      IF (EXC.GT.SCN(6)) GO TO 4                                        GGDR-085
      EBYT=(EXC-SCN(5))/SCN(3)                                          GGDR-086
      TINT=TINT*DEXP(EBYT)/SCN(3)                                       GGDR-087
      GO TO 5                                                           GGDR-088
    4 EXD=EXC+SCN(2)-SCN(6)                                             GGDR-089
      EBYT=DSQRT(SCN(1)*EXD)                                            GGDR-090
      SCNB=DEXP(2.D0*DSQRT(SCN(1)*SCN(2))-(SCN(6)-SCN(5))/SCN(3))*SCN(3)GGDR-091
     1/(SCN(1)**.5D0*SCN(2)**1.5D0)                                     GGDR-092
      TINT=TINT*DEXP(2.D0*EBYT)/(EBYT*EXD*SCNB)                         GGDR-093
    5 TG=TG+TINT                                                        GGDR-094
      IF (EHI.EQ.ELIM) GO TO 6                                          GGDR-095
      ELO=EHI                                                           GGDR-096
      GO TO 3                                                           GGDR-097
    6 IF (ECM.NE.0.D0) GO TO 9                                          GGDR-098
      TEMP=TG                                                           GGDR-099
      JD=IPI(3,1)-1                                                     GGDR-100
      N1=IABS(JD-1)                                                     GGDR-101
      N2=JD+1                                                           GGDR-102
      DO 8 I=N1,N2,2                                                    GGDR-103
      N3=IABS(I-2)                                                      GGDR-104
      N4=I+2                                                            GGDR-105
      DO 7 J=N3,N4,2                                                    GGDR-106
      DEX=-DFLOAT(J+1)**2/(SGSQ*4.D0)                                   GGDR-107
    7 ROJ=ROJ+DEXP(DEX)*DFLOAT(J+1)/SGSQ                                GGDR-108
    8 CONTINUE                                                          GGDR-109
      CTG=ATG0/(TEMP*ROJ)                                               GGDR-110
      ECM=WV(3,1)                                                       GGDR-111
      GO TO 2                                                           GGDR-112
    9 TG1=TG*CTG                                                        GGDR-113
      RETURN                                                            GGDR-114
 1000 FORMAT (/'    GAMMA RAY CHANNEL PARAMETERS :'//'   A=',I3,'  Z=',IGGDR-115
     13,'  N=',I3,5X,F9.3,' MEV  NEUTRON BINDING',5X,F6.2,' RADIATIVE D.GGDR-116
     2 OF F.',5X,'SIGMA=',F6.3)                                         GGDR-117
 1001 FORMAT ('   NORMALISED TO SLOW S-WAVE NEUTRON GAMMA WIDTHS/SPACINGGGDR-118
     1S =',D12.4)                                                       GGDR-119
 1002 FORMAT ('   E1 STRONG COUPLING MODEL.')                           GGDR-120
 1003 FORMAT ('   E1 GIANT RESONANCE AT ',F7.2,' MEV   WIDTH=',F7.2,' MEGGDR-121
     1V',5X,'EXCHANGE FRACTION=',F5.2)                                  GGDR-122
      END                                                               GGDR-123
