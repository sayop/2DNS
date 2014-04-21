cmake_minimum_required(VERSION 2.6)

project(CFD)

enable_language(Fortran)

#
# add sub-directories defined for each certain purpose
#
add_subdirectory(main)
add_subdirectory(io)

#
# set executable file name
#
set(CFD_EXE_NAME cfd.x CACHE STRING "CFD executable name")

#
# set source files
#
set(CFD_SRC_FILES ${MAIN_SRC_FILES}
                  ${IO_SRC_FILES})
                  
#
# define executable
#
add_executable(${CFD_EXE_NAME} ${CFD_SRC_FILES})
