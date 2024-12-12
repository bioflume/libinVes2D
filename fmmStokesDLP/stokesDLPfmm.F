#include <fintrf.h>
C
C stokesDLPfmm.F 
C
C This subroutine is the main gateway to MATLAB.  When a MEX function
C  is executed MATLAB calls the MEXFUNCTION subroutine in the corresponding
C  MEX file.  
C
C
      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER PLHS(*), PRHS(*)

C-----------------------------------------------------------------------
C

      INTEGER NLHS, NRHS
      INTEGER ERROR
      
C
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER MXCREATEDOUBLEMATRIX, MXGETPR

C-----------------------------------------------------------------------
C

      MWSIZE MXGETM, MXGETN
C
C KEEP THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS FOR USE
C IN ALL YOUR FORTRAN MEX FILES.
C---------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER F1,F2,XS,YS,DIR1,DIR2
      MWPOINTER U1,U2

C-----------------------------------------------------------------------
C

      MWSIZE MS, NS, NEL
      REAL *8, ALLOCATABLE :: RF1(:),RF2(:),RXS(:),RYS(:)
      REAL *8, ALLOCATABLE :: RDIR1(:),RDIR2(:)
      REAL *8, ALLOCATABLE :: RU1(:),RU2(:)

#if defined MSWIND
C     For Windows only!
C     This resets the floating point exception to allow divide by zero,
C     overflow and invalid numbers. 
C     
      INTEGER(2) CONTROL
      CALL GETCONTROLFPQQ(CONTROL)
      CONTROL = CONTROL .OR. FPCW$ZERODIVIDE
      CONTROL = CONTROL .OR. FPCW$INVALID
      CONTROL = CONTROL .OR. FPCW$OVERFLOW
      CALL SETCONTROLFPQQ(CONTROL)
#endif

C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
      IF (NRHS .NE. 6) THEN
        CALL MEXERRMSGTXT('FMM requires six input arguments')
      ELSEIF (NLHS .NE. 2) THEN
        CALL MEXERRMSGTXT('FMM Requires two output argument')
      ENDIF
      

c     GET THE NUMBER OF CHARGES/SOURCES.
      MS = MXGETM(PRHS(1))
      NS = MXGETN(PRHS(1));
      MS = MAX(MS,NS)
c     MS is the number of sources

      ALLOCATE(RF1(MS),STAT=ERROR)
      ALLOCATE(RF2(MS),STAT=ERROR)
      ALLOCATE(RXS(MS),STAT=ERROR)
      ALLOCATE(RYS(MS),STAT=ERROR)
      ALLOCATE(RDIR1(MS),STAT=ERROR)
      ALLOCATE(RDIR2(MS),STAT=ERROR)
      ALLOCATE(RU1(MS),STAT=ERROR)
      ALLOCATE(RU2(MS),STAT=ERROR)
C
C CREATE A MATRIX FOR RETURN ARGUMENT
      PLHS(1) = MXCREATEDOUBLEMATRIX(MS,1,0)
      PLHS(2) = MXCREATEDOUBLEMATRIX(MS,1,0)
C
C ASSIGN POINTERS TO THE VARIOUS PARAMETERS
C
      F1 = MXGETPR(PRHS(1))
      F2 = MXGETPR(PRHS(2))
      XS = MXGETPR(PRHS(3))
      YS = MXGETPR(PRHS(4))
      DIR1 = MXGETPR(PRHS(5))
      DIR2 = MXGETPR(PRHS(6))
C
      U1 = MXGETPR(PLHS(1))
      U2 = MXGETPR(PLHS(2))
C
C COPY RIGHT HAND ARGUMENTS TO LOCAL ARRAYS OR VARIABLES
      NEL = MS
      CALL MXCOPYPTRTOREAL8(F1, RF1, NEL)
      CALL MXCOPYPTRTOREAL8(F2, RF2, NEL)
      CALL MXCOPYPTRTOREAL8(XS, RXS, NEL)
      CALL MXCOPYPTRTOREAL8(YS, RYS, NEL)
      CALL MXCOPYPTRTOREAL8(DIR1, RDIR1, NEL)
      CALL MXCOPYPTRTOREAL8(DIR2, RDIR2, NEL)

C
C DO THE ACTUAL COMPUTATIONS IN A SUBROUTINE
C       CREATED ARRAYS.
C
      CALL stokesDLP(RF1,RF2,RXS,RYS,MS,RDIR1,RDIR2,
     $        RU1,RU2)

C
C COPY OUTPUT WHICH IS STORED IN LOCAL ARRAY TO MATRIX OUTPUT
      NEL = MS
      CALL MXCOPYREAL8TOPTR(RU1, U1, NEL)
      CALL MXCOPYREAL8TOPTR(RU2, U2, NEL)
C



      RETURN
      END
