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

find_library(M_LIB m)

if (DO_MPI)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDOMPI")

  # Find MPI library
  find_package(MPI REQUIRED)
  set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
  set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
  include_directories(${MPI_INCLUDE_PATH})

  add_executable(
    CoMDParallel
    ${EXECUTABLE_FLAG}
    pmd.c
    docompute.c
    domains.c
    geometry.c
    ic_fcc.c
    mycommand.c
    timestep.c
    parallel.c
    ../src-flat/cmdLineParser.c
    ../src-flat/eam.c
    ../src-flat/ljForce.c
    ../src-flat/read.c
    ../src-flat/write.c
    ../src-flat/utility.c
    )

  set_target_properties(CoMDParallel PROPERTIES COMPILE_FLAGS "${OPT} -I ${CMAKE_SOURCE_DIR}/src-parallel/include -I ${CMAKE_SOURCE_DIR}/src-flat/include -I /usr/include -DDO_MPI") 
  target_link_libraries(CoMDParallel ${MPI_LIBRARIES} ${Z_LIB} ${M_LIB})
endif ()
