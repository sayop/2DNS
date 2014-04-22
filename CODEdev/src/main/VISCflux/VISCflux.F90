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

  FUNCTION SutherlandLaw(temp,isu) RESULT(SU)
    !Return dimensional viscosity and thermal conductivity
    INTEGER :: i, j, isu
    REAL(KIND=wp) :: SU, C1, C2, C3, C4, temp

    C1 = 1.458E-06_wp
    C2 = 110.4_wp
    C3 = 2.495E-03_wp
    C4 = 194.0_wp 

    IF(isu .EQ. 1) THEN
      ! Set nondimensional viscosity
      SU = C1 * temp ** 1.5 / (temp + C2)
    ELSE IF(isu .EQ. 2) THEN
      ! Set thermal conductivity
      SU = C3 * temp ** 1.5 / (temp + C4)
    END IF
  END FUNCTION SutherlandLaw

  FUNCTION ShearStress(i,j,i1,i2) RESULT(TAU)
    USE SimulationVars_m, ONLY: V, imin, imax, jmin, jmax, &
                                Temp, RE_REF, MU_REF
    USE GridJacobian_m, ONLY: PIPX, PIPY, PJPX, PJPY

    INTEGER :: i, j, i1, i2
    REAL :: DUDX, DUDY, DVDX, DVDY, &
            DUDI, DUDJ, DVDI, DVDJ, MU, &
            TAU

    IF(i1 .EQ. 1 .AND. i2 .EQ. 1) THEN
      !Set DUDX
      DUDI = V(2,i+1,j) - V(2,i,j)
      DUDJ = 0.25_wp * ( V(2,i+1,j+1) + V(2,i,j+1) - &
                         V(2,i+1,j-1) - V(2,i,j-1) )
      DUDX = 0.5_wp * ( (PIPX(i+1,j) + PIPX(i,j)) * DUDI + &
                        (PJPX(i+1,j) + PJPX(i,j)) * DUDJ )
      !Set DVDY
      DVDI = V(3,i+1,j) - V(3,i,j)
      DVDJ = 0.25_wp * ( V(3,i+1,j+1) + V(3,i,j+1) - &
                         V(3,i+1,j-1) - V(3,i,j-1) )
      DVDY = 0.5_wp * ( (PIPY(i+1,j) + PIPY(i,j)) * DVDI + &
                        (PJPY(i+1,j) + PJPY(i,j)) * DVDJ )
      !Set TAU_xx
      TAU = 2.0_wp * MU_REF / (3.0_wp * RE_REF) * ( &
            2.0_wp * DUDX - DVDY)
    ELSE IF(i1 .EQ. 1 .AND. i2 .EQ. 2) THEN
      !Set DUDY

    ELSE IF(i1 .EQ. 2 .AND. i2 .EQ. 1) THEN
      !Set DVDX

    ELSE IF(i1 .EQ. 2 .AND. i2 .EQ. 2) THEN
      !Set DVDY

    END IF
    MU = SutherlandLaw(Temp(i,j),1) / MU_REF
    

  END FUNCTION ShearStress
END MODULE VISCflux_m
