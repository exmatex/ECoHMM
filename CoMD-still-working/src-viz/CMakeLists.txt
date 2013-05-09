cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

IF(APPLE)
#  SET(EXECUTABLE_FLAG MACOSX_BUNDLE)
ENDIF(APPLE)

if (DOUBLE_PRECISION)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDOUBLE")
else ()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSINGLE")
endif ()

if (APPROX_CENTROSYM)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DAPPROX_CENTROSYM")
endif ()

if (DO_GZIP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DGZIPSUPPORT")
  find_library(Z_LIB z)
endif ()

if (USE_VTK)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_VTK")
endif ()

if (VIZ_ANIMATION)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DVIZ_ANIMATION")
endif ()

find_library(M_LIB m)

if (DO_IN_SITU_VIZ)
  if (USE_VTK)
    find_package(VTK)
#    include(${VTK_USE_FILE})
  endif ()
  set(CMAKE_INSTALL_RPATH ${VTK_DIR}/bin)
  add_executable(
      CoMDViz
      ${EXECUTABLE_FLAG}
      AtomVisualize.c
      quaternion.c
      pmdViz.c
      ../src-flat/docompute.c
      ../src-flat/ic_fcc.c
      ../src-flat/cmdLineParser.c
      ../src-flat/domains.c
      ../src-flat/eam.c
      ../src-flat/ljForce.c
      ../src-flat/mycommand.c
      ../src-flat/read.c
      ../src-flat/write.c
      ../src-flat/utility.c
      ../src-flat/timestep.c
      ../src-flat/cheby.c
      )
  target_link_libraries(CoMDViz ${EXTRA_LIBS})
  if (USE_VTK)
    target_link_libraries(CoMDViz vtkRendering vtkFiltering vtkCommon)
  endif ()

  if (APPLE)
    target_link_libraries(CoMDViz ${Z_LIB} ${M_LIB} stdc++ pthread)
  else ()
    target_link_libraries(CoMDViz ${Z_LIB} ${M_LIB} stdc++ pthread glut GLU)
  endif ()

  set_target_properties(CoMDViz PROPERTIES COMPILE_FLAGS "${OPT} ${EXTRA_FLAGS} -I ${GLUT_INCLUDE_DIR} -I ${CMAKE_SOURCE_DIR}/src-flat/include -I ${CMAKE_SOURCE_DIR}/src-viz/include -I /usr/include -I include -DUSE_IN_SITU_VIZ -x c++")


  #  set_target_properties(CoMDViz PROPERTIES COMPILE_FLAGS "${OPT} ${CMAKE_C_FLAGS} ${EXTRA_FLAGS} -I ${GLUT_INCLUDE_DIR} -I ${CMAKE_SOURCE_DIR}/src-flat/include -I ${CMAKE_SOURCE_DIR}/src-viz/include -I /usr/include -I include -DUSE_IN_SITU_VIZ -x c++")
  
endif ()

