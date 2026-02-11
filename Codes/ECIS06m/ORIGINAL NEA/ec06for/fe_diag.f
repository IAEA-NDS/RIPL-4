C 02/04/06                                                      ECIS06  DIAG-000
      SUBROUTINE DIAG(ZR,ZI,XR,XI,N,NC,EPS,AX,IER)                      DIAG-001
C DIAGONALISATION OF A HERMITIAN COMPLEX MATRIX BY AN EXTENSION OF THE  DIAG-002
C JACOBI'S METHOD.                                                      DIAG-003
C INPUT:     ZR,ZI:  REAL AND IMAGINARY PARTS OF THE MATRIX.            DIAG-004
C            XR,XI:  REAL AND IMAGINARY PARTS OF THE UNIT MATRIX.       DIAG-005
C            N:      FIRST DIMENSION OF ZR,ZI,XR AND XI.                DIAG-006
C            NC:     DIMENSION OF THE MATRIX.                           DIAG-007
C            EPS:    VALUE BELOW WHICH MATRIX ELEMENTS ARE SET TO 0.    DIAG-008
C OUTPUT:    ZR,ZI:  THE EIGENVALUES ARE ON THE DIAGONAL OF ZR.         DIAG-009
C                    ALL THE OTHER ELEMENTS ARE 0, IF PROCESS SUCCEEDED.DIAG-010
C            XR,XI:  EIGENVECTORS.                                      DIAG-011
C            AX:     SQUARE OF NORM OF THE LARGEST NON DIAGONAL ELEMENT.DIAG-012
C            IER:    RETURNS 0 OR -1 AFTER 4*NC**2 ROTATIONS.           DIAG-013
C***********************************************************************DIAG-014
      IMPLICIT REAL*8 (A-H,O-Z)                                         DIAG-015
      DIMENSION ZR(N,*),ZI(N,*),XR(N,*),XI(N,*)                         DIAG-016
      IER=0                                                             DIAG-017
      NT=0                                                              DIAG-018
    1 NT=NT+1                                                           DIAG-019
      IF (NT.GT.4*NC*NC) GO TO 6                                        DIAG-020
      AX=0.D0                                                           DIAG-021
      L=1                                                               DIAG-022
      M=2                                                               DIAG-023
C SYMMETRISATION AND SEARCH FOR THE LARGEST NON DIAGONAL ELEMENT.       DIAG-024
      DO 3 I=1,NC                                                       DIAG-025
      DO 2 J=I,NC                                                       DIAG-026
      IF (ZR(J,I).EQ.0.D0) ZR(I,J)=0.D0                                 DIAG-027
      IF (ZI(J,I).EQ.0.D0) ZI(I,J)=0.D0                                 DIAG-028
      IF (ZR(I,J).EQ.0.D0) ZR(J,I)=0.D0                                 DIAG-029
      IF (ZI(I,J).EQ.0.D0) ZI(J,I)=0.D0                                 DIAG-030
      AR=(ZR(I,J)+ZR(J,I))/2.D0                                         DIAG-031
      AI=(ZI(I,J)-ZI(J,I))/2.D0                                         DIAG-032
      ZR(J,I)=AR                                                        DIAG-033
      ZR(I,J)=AR                                                        DIAG-034
      ZI(I,J)=AI                                                        DIAG-035
      ZI(J,I)=-AI                                                       DIAG-036
      IF (I.EQ.J) GO TO 2                                               DIAG-037
      AY=ZR(I,J)**2+ZI(I,J)**2                                          DIAG-038
      IF (AX.GT.AY) GO TO 2                                             DIAG-039
      AX=AY                                                             DIAG-040
      L=I                                                               DIAG-041
      M=J                                                               DIAG-042
    2 CONTINUE                                                          DIAG-043
    3 CONTINUE                                                          DIAG-044
      IF (AX.EQ.0.D0) RETURN                                            DIAG-045
C ELEMENTARY TRANSFORMATION.                                            DIAG-046
      U=DATAN2(-ZI(L,M),ZR(L,M))/2.D0                                   DIAG-047
      V=DATAN2(2.D0*DSQRT(ZR(L,M)**2+ZI(L,M)**2),ZR(M,M)-ZR(L,L))/2.D0  DIAG-048
      UC=DCOS(U)                                                        DIAG-049
      US=DSIN(U)                                                        DIAG-050
      TC=DCOS(V)                                                        DIAG-051
      TS=-DSIN(V)                                                       DIAG-052
      UCC=UC*TC                                                         DIAG-053
      UCS=UC*TS                                                         DIAG-054
      USC=US*TC                                                         DIAG-055
      USS=US*TS                                                         DIAG-056
C TRANSFORMATION OF ROWS.                                               DIAG-057
      DO 4 I=1,NC                                                       DIAG-058
      AR=XR(I,L)*UCC+XI(I,L)*USC+XR(I,M)*UCS-XI(I,M)*USS                DIAG-059
      BR=-XR(I,L)*UCS-XI(I,L)*USS+XR(I,M)*UCC-XI(I,M)*USC               DIAG-060
      AI=XI(I,L)*UCC-XR(I,L)*USC+XI(I,M)*UCS+XR(I,M)*USS                DIAG-061
      BI=-XI(I,L)*UCS+XR(I,L)*USS+XI(I,M)*UCC+XR(I,M)*USC               DIAG-062
      XR(I,L)=AR                                                        DIAG-063
      XR(I,M)=BR                                                        DIAG-064
      XI(I,L)=AI                                                        DIAG-065
      XI(I,M)=BI                                                        DIAG-066
      AR=ZR(I,L)*UCC+ZI(I,L)*USC+ZR(I,M)*UCS-ZI(I,M)*USS                DIAG-067
      BR=-ZR(I,L)*UCS-ZI(I,L)*USS+ZR(I,M)*UCC-ZI(I,M)*USC               DIAG-068
      AI=ZI(I,L)*UCC-ZR(I,L)*USC+ZI(I,M)*UCS+ZR(I,M)*USS                DIAG-069
      BI=-ZI(I,L)*UCS+ZR(I,L)*USS+ZI(I,M)*UCC+ZR(I,M)*USC               DIAG-070
      ZR(I,L)=AR                                                        DIAG-071
      ZR(I,M)=BR                                                        DIAG-072
      ZI(I,L)=AI                                                        DIAG-073
    4 ZI(I,M)=BI                                                        DIAG-074
C TRANSFORMATION OF COLUMNS.                                            DIAG-075
      DO 5 I=1,NC                                                       DIAG-076
      AR=ZR(L,I)*UCC-ZI(L,I)*USC+ZR(M,I)*UCS+ZI(M,I)*USS                DIAG-077
      BR=-ZR(L,I)*UCS+ZI(L,I)*USS+ZR(M,I)*UCC+ZI(M,I)*USC               DIAG-078
      AI=ZI(L,I)*UCC+ZR(L,I)*USC+ZI(M,I)*UCS-ZR(M,I)*USS                DIAG-079
      BI=-ZI(L,I)*UCS-ZR(L,I)*USS+ZI(M,I)*UCC-ZR(M,I)*USC               DIAG-080
      IF (DABS(AR).LT.EPS) AR=0.D0                                      DIAG-081
      IF (DABS(BR).LT.EPS) BR=0.D0                                      DIAG-082
      IF (DABS(AI).LT.EPS) AI=0.D0                                      DIAG-083
      IF (DABS(BI).LT.EPS) BI=0.D0                                      DIAG-084
      ZR(L,I)=AR                                                        DIAG-085
      ZR(M,I)=BR                                                        DIAG-086
      ZI(L,I)=AI                                                        DIAG-087
    5 ZI(M,I)=BI                                                        DIAG-088
      GO TO 1                                                           DIAG-089
    6 IER=-1                                                            DIAG-090
      RETURN                                                            DIAG-091
      END                                                               DIAG-092
