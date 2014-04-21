!> \file parameters.F90
!> \author Sayop Kim
!> \brief Provides parameters and physical constants for use throughout the
!! code.
MODULE Parameters_m
   INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(8)
   ! Number of total elements in state vector
   INTEGER, PARAMETER :: nelm = 4

   CHARACTER(LEN=10), PARAMETER :: CODE_VER_STRING = "V.001.001"
   REAL(KIND=wp), PARAMETER :: PI = 3.14159265358979323846264338_wp

END MODULE Parameters_m

