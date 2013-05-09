cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

if (DOUBLE_PRECISION)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDOUBLE")
else ()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSINGLE")
endif ()

if (DO_GZIP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DGZIPSUPPORT")
  find_library(Z_LIB z)
endif ()

if (USE_CHEBY)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_CHEBY")
endif()

find_library(M_LIB m)

include_directories(${CMAKE_SOURCE_DIR}/src-flat/include)
include_directories(/usr/include)

if (NOT DO_IN_SITU_VIZ)
  add_executable(
    CoMD
      pmd.c
      cheby.c
      docompute.c
      cmdLineParser.c
      domains.c
      eam.c
      ic_fcc.c
      ljForce.c
      mycommand.c
      read.c
      write.c
      utility.c
      timestep.c
      yamlOutput.c)
  set_target_properties(CoMD PROPERTIES COMPILE_FLAGS "${OPT}") 
  target_link_libraries(CoMD ${Z_LIB} ${M_LIB})
endif ()

