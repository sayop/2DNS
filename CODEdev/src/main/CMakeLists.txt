add_subdirectory(AUSMPWplus)
add_subdirectory(VISCflux)
add_subdirectory(BCvisc)
set(MAIN_SRC_FILES
   ${CMAKE_CURRENT_SOURCE_DIR}/main.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/MainLoop.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/GridJacobian.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/SimulationSetup.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/SimulationVars.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/Parameters.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/TimeIntegration.F90
   ${AUSMPWplus_SRC_FILES}
   ${VISCflux_SRC_FILES}
   ${BCVISC_SRC_FILES}
   CACHE INTERNAL "" FORCE)

