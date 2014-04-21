!> \file: GridJacobian.F90
!> \author: Sayop Kim

MODULE GridJacobian_m
   USE Parameters_m, ONLY: wp
   IMPLICIT NONE

   PUBLIC :: InverseGridMetricsArrays, GridMetricsArrays, &
             PIPX, PJPX, PIPY, PJPY, JACOBIAN

   REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:) :: PXPI, PXPJ, PYPI, PYPJ, &
                                                 PIPX, PJPX, PIPY, PJPY, &
                                                 JACOBIAN
CONTAINS

!-----------------------------------------------------------------------------!
  SUBROUTINE InverseGridMetricsArrays()
!-----------------------------------------------------------------------------!
    USE SimulationVars_m, ONLY: IMIN, JMIN, IMAX, JMAX, INCELL, JNCELL, &
                                XP
    IMPLICIT NONE
    INTEGER :: I, J

    ! PXPI: x_xi
    ! PXPJ: x_eta
    ! PYPI: y_xi
    ! PYPJ: y_eta
    ! JACOBIAN: J = 1 / (x_xi * y_eta - x_eta * y_xi)
    ALLOCATE(PXPI(INCELL, JNCELL))
    ALLOCATE(PXPJ(INCELL, JNCELL))
    ALLOCATE(PYPI(INCELL, JNCELL))
    ALLOCATE(PYPJ(INCELL, JNCELL))
    ALLOCATE(JACOBIAN(INCELL, JNCELL))
    PXPI = 0.0_wp
    PXPJ = 0.0_wp
    PYPI = 0.0_wp
    PYPJ = 0.0_wp

    DO I = IMIN, IMAX
      DO J = JMIN, JMAX
        IF(J .EQ. JMIN) THEN
          !Forward difference (2nd order)
          PXPJ(I,J) = 0.5_wp * (-3.0_wp*XP(1,I,J) + 4.0_wp*XP(1,I,J+1) &
                                - XP(1,I,J+2))
          PYPJ(I,J) = 0.5_wp * (-3.0_wp*XP(2,I,J) + 4.0_wp*XP(2,I,J+1) &
                                - XP(2,I,J+2))
        ELSE IF(J .EQ. JMAX) THEN
          !Backward difference (2nd order)
          PXPJ(I,J) = 0.5_wp * (3.0_wp*XP(1,I,J) - 4.0_wp*XP(1,I,J-1) &
                                + XP(1,I,J-2))
          PYPJ(I,J) = 0.5_wp * (3.0_wp*XP(2,I,J) - 4.0_wp*XP(2,I,J-1) &
                                + XP(2,I,J-2))
        ELSE
          !Central difference (2nd order)
          PXPJ(I,J) = 0.5_wp * (XP(1,I,J+1) - XP(1,I,J-1))
          PYPJ(I,J) = 0.5_wp * (XP(2,I,J+1) - XP(2,I,J-1))
        END IF
      END DO
    END DO

    DO J = JMIN, JMAX
      DO I = IMIN, IMAX
        IF(I .EQ. IMIN) THEN
          !Forward difference (2nd order)
          PXPI(I,J) = 0.5_wp * (-3.0_wp*XP(1,I,J) + 4.0_wp*XP(1,I+1,J) &
                                - XP(1,I+2,J))
          PYPI(I,J) = 0.5_wp * (-3.0_wp*XP(2,I,J) + 4.0_wp*XP(2,I+1,J) &
                                - XP(2,I+2,J))
        ELSE IF(I .EQ. IMAX) THEN
          !Backward difference (2nd order)
          PXPI(I,J) = 0.5_wp * (3.0_wp*XP(1,I,J) - 4.0_wp*XP(1,I-1,J) &
                                + XP(1,I-2,J))
          PYPI(I,J) = 0.5_wp * (3.0_wp*XP(2,I,J) - 4.0_wp*XP(2,I-1,J) &
                                + XP(2,I-2,J))
        ELSE
          !Central difference (2nd order)
          PXPI(I,J) = 0.5_wp * (XP(1,I+1,J) - XP(1,I-1,J))
          PYPI(I,J) = 0.5_wp * (XP(2,I+1,J) - XP(2,I-1,J))
        END IF
      END DO
    END DO

    DO I = IMIN, IMAX
      DO J = JMIN, JMAX
        ! Calculate Grid Jacobian: J
        JACOBIAN(I,J) = 1.0_wp / (PXPI(I,J) * PYPJ(I,J) - &
                                  PXPJ(I,J) * PYPI(I,J))
      END DO
    END DO
!    DO I = IMIN, IMAX
!      DO J = JMIN, JMAX
!        IF(I .EQ. IMIN) THEN
!          ! Forward difference in i-direction
!          PXPI(I,J) = 0.5_wp * (-3.0_wp*XP(1,I,J) + 4.0_wp*XP(1,I+1,J) &
!                                - XP(1,I+2,J))
!          PYPI(I,J) = 0.5_wp * (-3.0_wp*XP(2,I,J) + 4.0_wp*XP(2,I+1,J) &
!                                - XP(2,I+2,J))
!        ELSEIF(I .EQ. IMAX) THEN
!          ! Backward difference in i-direction
!          PXPI(I,J) = 0.5_wp * (3.0_wp*XP(1,I,J) - 4.0_wp*XP(1,I-1,J) &
!                                + XP(1,I-2,J))
!          PYPI(I,J) = 0.5_wp * (3.0_wp*XP(2,I,J) - 4.0_wp*XP(2,I-1,J) &
!                                + XP(2,I-2,J))
!        ELSE
!          ! Central difference in i-direction
!          PXPI(I,J) = 0.5_wp * (XP(1,I+1,J) - XP(1,I-1,J))
!          PYPI(I,J) = 0.5_wp * (XP(2,I+1,J) - XP(2,I-1,J))
!        ENDIF
!
!        IF(J .EQ. JMIN) THEN
!          ! Forward difference in j-direction
!          PXPJ(I,J) = 0.5_wp * (-3.0_wp*XP(1,I,J) + 4.0_wp*XP(1,I,J+1) &
!                                - XP(1,I,J+2))
!          PYPJ(I,J) = 0.5_wp * (-3.0_wp*XP(2,I,J) + 4.0_wp*XP(2,I,J+1) &
!                                - XP(2,I,J+2))
!        ELSEIF(J .EQ. JMAX) THEN
!          ! Backward difference in j-direction
!          PXPJ(I,J) = 0.5_wp * (3.0_wp*XP(1,I,J) - 4.0_wp*XP(1,I,J-1) &
!                                + XP(1,I,J-2))
!          PYPJ(I,J) = 0.5_wp * (3.0_wp*XP(2,I,J) - 4.0_wp*XP(2,I,J-1) &
!                                + XP(2,I,J-2))
!        ELSE
!          ! Central difference in j-direction
!          PXPJ(I,J) = 0.5_wp * (XP(1,I,J+1) - XP(1,I,J-1))
!          PYPJ(I,J) = 0.5_wp * (XP(2,I,J+1) - XP(2,I,J-1))
!        ENDIF
!        ! Calculate Grid Jacobian: J
!        JACOBIAN(I,J) = 1.0_wp / (PXPI(I,J) * PYPJ(I,J) - &
!                                  PXPJ(I,J) * PYPI(I,J))
!      END DO
!    END DO
  END SUBROUTINE InverseGridMetricsArrays

!-----------------------------------------------------------------------------!
  SUBROUTINE GridMetricsArrays()
!-----------------------------------------------------------------------------!
    USE SimulationVars_m, ONLY: IMIN, JMIN, IMAX, JMAX, INCELL, JNCELL, &
                                XP
    IMPLICIT NONE
    INTEGER :: I, J

    ! PIPX: xi_x
    ! PJPX: eta_x
    ! PIPY: xi_y
    ! PJPY: eta_y
    ALLOCATE(PIPX(INCELL, JNCELL))
    ALLOCATE(PJPX(INCELL, JNCELL))
    ALLOCATE(PIPY(INCELL, JNCELL))
    ALLOCATE(PJPY(INCELL, JNCELL))
    PIPX = 0.0_wp
    PJPX = 0.0_wp
    PIPY = 0.0_wp
    PJPY = 0.0_wp

    DO I = IMIN, IMAX
      DO J = JMIN, JMAX
        PIPX(I,J) =  PYPJ(I,J) * JACOBIAN(I,J)
        PJPX(I,J) = -PYPI(I,J) * JACOBIAN(I,J)
        PIPY(I,J) = -PXPJ(I,J) * JACOBIAN(I,J)
        PJPY(I,J) =  PXPI(I,J) * JACOBIAN(I,J)
      END DO
    END DO
  END SUBROUTINE GridMetricsArrays
END MODULE GridJacobian_m
