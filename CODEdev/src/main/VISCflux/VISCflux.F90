!> \file: SimulationSetup.F90
!> \author: Sayop Kim

MODULE VISCflux_m
   USE Parameters_m, ONLY: wp

   IMPLICIT NONE

CONTAINS
!-----------------------------------------------------------------------------!
  SUBROUTINE VISCflux()
!-----------------------------------------------------------------------------!
    USE SimulationVars_m, ONLY: imin, jmin, imax, jmax, incell, jncell, &
                                DF
    USE GridJacobian_m, ONLY: PIPX, PIPY, PJPX, PJPY, JACOBIAN

    INTEGER i, j, n
    REAL(KIND=wp), DIMENSION(4,incell,jncell) :: FPhalf, GPhalf

    FPhalf = 0.0_wp
    GPhalf = 0.0_wp
    !
    ! Set transformed flux along j-const lines at every (i+1/2) points
    !
    DO j = jmin + 1, jmax - 1
      DO i = imin, imax - 1


        FPhalf(1,i,j) = 0.0_wp
      END DO
    END DO
    
    !
    ! Set transformed flux in along i-const lines at every (j+1/2) points
    !
    DO i = imin + 1, imax - 1
      DO j = jmin, jmax - 1


        GPhalf(1,i,j) = 0.0_wp
      END DO
    END DO

    ! Update DF vector
    DO i = imin + 1, imax - 1
      DO j = jmin + 1, jmax - 1
        DO n = 1, 4
          DF(n,i,j) = DF(n,i,j) - ( &
                      FPhalf(n,i,j) - FPhalf(n,i-1,j) + &
                      GPhalf(n,i,j) - GPhalf(n,i,j-1) )
        END DO
      END DO
    END DO

  END SUBROUTINE VISCflux

  FUNCTION SutherlandLaw(i,j,isu) RESULT(SU)
    USE SimulationVars_m, ONLY: Temp, imin, jmin, imax, jmax
    IMPLICIT NONE
    INTEGER :: i, j, isu
    REAL(KIND=wp) :: SU

  END FUNCTION SutherlandLaw
END MODULE VISCflux_m
