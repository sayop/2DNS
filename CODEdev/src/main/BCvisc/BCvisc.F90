!> \file: SimulationSetup.F90
!> \author: Sayop Kim

MODULE BCvisc_m
   USE Parameters_m, ONLY: wp

   IMPLICIT NONE
   REAL(KIND=wp) :: wallTemp

CONTAINS

!-----------------------------------------------------------------------------!
  SUBROUTINE BCvisc()
!-----------------------------------------------------------------------------!
    USE SimulationVars_m, ONLY: imin, jmin, imax, jmax, V, cgamma, mach_ref
    USE GridJacobian_m, ONLY: PXPI, PYPI, PXPJ, PYPJ
    INTEGER :: i, j
    REAL(KIND=wp) :: coef

    j = jmin
    DO i = imin, imax
      V(2,i,j) = 0.0_wp
      V(3,i,j) = 0.0_wp

      ! Set pressure boundary condition
      coef = PXPI(i,j) ** 2 + PYPI(i,j) ** 2
      coef = (PXPI(i,j) * PXPJ(i,j) + PYPI(i,j) * PYPJ(i,j)) / coef
      IF(i .eq. imin) THEN
        V(4,i,j) = V(4,i,j+1) - coef * ( V(4,i+1,j+1) - V(4,i,j+1) )
      ELSE IF(i .eq. imax) THEN
        V(4,i,j) = V(4,i,j+1) - coef * ( V(4,i,j+1) - V(4,i-1,j+1) )
      ELSE
        V(4,i,j) = V(4,i,j+1) - coef * 0.5_wp * ( &
                   V(4,i+1,j+1) - V(4,i-1,j+1) )
      END IF

      ! Set temperature boundary condition
      IF(wallTemp .GE. 0.0_wp) THEN
        ! Set isothermal BC
        V(5,i,j) = wallTemp
      ELSE
        V(5,i,j) = (4.0_wp * V(5,i,j+1) - V(5,i,j+2)) / 3.0_wp
      END IF

      ! Set density boundary condition
      V(1,i,j) = cgamma * V(4,i,j) * mach_ref ** 2 / V(5,i,j)
    END DO


  END SUBROUTINE BCvisc
END MODULE BCvisc_m
