C 28/02/07                                                      ECIS06  DERI-000
      SUBROUTINE DERI(X,Y,H,N,LT)                                       DERI-001
C NUMERICAL DERIVATION OF THE FUNCTION Y KNOWN IN N POINTS WITH STEP H  DERI-002
C IT NEEDS AT LEAST 7 POINTS AND RETURNS IN X THE VALUE OF -D(Y)/DR     DERI-003
C IF LT IS FALSE, THE RESULT IS DIVIDED BY R..                          DERI-004
C INPUT:     Y:       FUNCTION TO BE DERIVED.                           DERI-005
C            H:       STEP SIZE OF Y.                                   DERI-006
C            N:       NUMBER OF POINTS OF Y.                            DERI-007
C            LT:      IF LT=.FALSE., THE RESULT IS DIVIDED BY R.        DERI-008
C OUTPUT:    X:       RESULT.                                           DERI-009
C***********************************************************************DERI-010
      IMPLICIT REAL*8 (A-H,O-Z)                                         DERI-011
      LOGICAL LT                                                        DERI-012
      DIMENSION X(*),Y(*)                                               DERI-013
      COMMON /INOUT/ MR,MW,MS                                           DERI-014
      N3=N-3                                                            DERI-015
      IF (N3.GE.4) GO TO 1                                              DERI-016
      WRITE (MW,1000) N                                                 DERI-017
      STOP                                                              DERI-018
    1 HH=-H*60.D0                                                       DERI-019
      X(1)=(-147.D0*Y(1)+360.D0*Y(2)-450.D0*Y(3)+400.D0*Y(4)-225.D0*Y(5)DERI-020
     1+72.D0*Y(6)-10.D0*Y(7))/HH                                        DERI-021
      X(2)=(-10.D0*Y(1)-77.D0*Y(2)+150.D0*Y(3)-100.D0*Y(4)+50.D0*Y(5)-15DERI-022
     1.D0*Y(6)+2.D0*Y(7))/HH                                            DERI-023
      X(3)=(2.D0*Y(1)-24.D0*Y(2)-35.D0*Y(3)+80.D0*Y(4)-30.D0*Y(5)+8.D0*YDERI-024
     1(6)-Y(7))/HH                                                      DERI-025
      DO 2 I=4,N3                                                       DERI-026
    2 X(I)=(45.D0*(Y(I+1)-Y(I-1))-9.D0*(Y(I+2)-Y(I-2))+Y(I+3)-Y(I-3))/HHDERI-027
      X(N-2)=(Y(N-6)-8.D0*Y(N-5)+30.D0*Y(N-4)-80.D0*Y(N3)+35.D0*Y(N-2)+2DERI-028
     14.D0*Y(N-1)-2.D0*Y(N))/HH                                         DERI-029
      X(N-1)=(-2.D0*Y(N-6)+15.D0*Y(N-5)-50.D0*Y(N-4)+100.D0*Y(N3)-150.D0DERI-030
     1*Y(N-2)+77.D0*Y(N-1)+10.D0*Y(N))/HH                               DERI-031
      X(N)=(10.D0*Y(N-6)-72.D0*Y(N-5)+225.D0*Y(N-4)-400.D0*Y(N3)+450.D0*DERI-032
     1Y(N-2)-360.D0*Y(N-1)+147.D0*Y(N))/HH                              DERI-033
      IF (LT) RETURN                                                    DERI-034
      R=0.D0                                                            DERI-035
      DO 3 I=1,N                                                        DERI-036
      R=R+H                                                             DERI-037
    3 X(I)=X(I)/R                                                       DERI-038
      RETURN                                                            DERI-039
 1000 FORMAT (5X,I5,' POINTS INSUFFICIENT FOR DERIVATION'///' IN DERI  .DERI-040
     1..  STOP  ...')                                                   DERI-041
      END                                                               DERI-042
