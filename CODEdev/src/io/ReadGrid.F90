!> \file: ReadGrid.F90
!> \author: Sayop Kim
!> \brief: Provides routines to read grid file

MODULE ReadGrid_m
   USE Parameters_m, ONLY: wp

   IMPLICIT NONE

   PUBLIC :: INGRID

CONTAINS

!-----------------------------------------------------------------------------!
  SUBROUTINE INGRID()
!-----------------------------------------------------------------------------!
! Read grid file
! X: X locations of grid points
! Y: Y locations of grid points
! IMAX: maximum number of points in the i-direction
! JMAX: maximum number of points in the j-direction
!-----------------------------------------------------------------------------!
    USE SimulationVars_m, ONLY: XP, IMAX, JMAX, IMIN, JMIN
    USE io_m, ONLY: IOunit, gridFile

    IMPLICIT NONE
    INTEGER :: ios, I, J
    REAL(KIND=wp) :: var

    OPEN(IOunit, FILE = gridFile, FORM = 'FORMATTED', ACTION = 'READ', &
         STATUS = 'OLD', IOSTAT = ios)
    IF(ios /= 0) THEN
      WRITE(*,'(a)') ""
      WRITE(*,'(a)') "Fatal error: Could not open the grid file."
      RETURN
    ELSE
      WRITE(*,'(a)') ""
      WRITE(*,'(a,1x,a)') "Reading grid file:", gridFile
    ENDIF

    DO I = IMIN, IMAX
      DO J = JMIN, JMAX
        READ(IOunit,'(1x,e17.10,5x,e17.10,5x,e17.10)') XP(1,I,J), XP(2,I,J), var
        !READ(IOunit,'(1x,e11.4,5x,e11.4)') XP(1,I,J), XP(2,I,J)
      END DO  ! J-LOOP
      READ(IOunit,*)
    END DO  ! I-LOOP

  END SUBROUTINE INGRID

END MODULE ReadGrid_m
