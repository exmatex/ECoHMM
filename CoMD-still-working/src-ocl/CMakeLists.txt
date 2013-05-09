cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

if (DOUBLE_PRECISION)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDOUBLE")
   set(PREC_FLAG "-DDOUBLE")
else ()
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSINGLE")
   set(PREC_FLAG "-DSINGLE")
endif ()

if (DO_GZIP)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DGZIPSUPPORT")
   find_library(Z_LIB z)
endif ()

if (USE_CHEBY)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_CHEBY")
endif()

if (DO_INTEROP_VIZ)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DINTEROP_VIZ")
endif ()

find_library(M_LIB m)

if (USE_VTK)
   find_package(VTK)
   include(${VTK_USE_FILE})
   set(CMAKE_INSTALL_RPATH ${VTK_DIR}/bin)
endif ()

include_directories(/usr/include)
include_directories(${OPENCL_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src-flat/include)
include_directories(${CMAKE_SOURCE_DIR}/src-ocl/include)

if (DO_OCL)
   find_package(OPENCL)
   add_executable(
            CoMDOCL
            ${EXECUTABLE_FLAG}
            cl_utils.c
            helpers.c
            pmdOCL.c
            docompute.c
            quaternion.c
	    yamlOutput.c
            ../src-flat/cmdLineParser.c
            ../src-flat/domains.c
            ../src-flat/eam.c
            ../src-flat/ljForce.c
            ../src-flat/ic_fcc.c
            ../src-flat/mycommand.c
            ../src-flat/read.c
            ../src-flat/cheby.c
            ../src-flat/write.c
            ../src-flat/utility.c
            ../src-flat/timestep.c
            )
   target_link_libraries(CoMDOCL ${EXTRA_LIBS})
   if (APPLE)
      if (DO_INTEROP_VIZ)
         target_link_libraries(CoMDOCL stdc++ ${Z_LIB} ${CMAKE_SOURCE_DIR}/src-ocl/libGLEW.a)
      else (DO_INTEROP_VIZ)
         target_link_libraries(CoMDOCL stdc++ ${Z_LIB} ${M_LIB})
      endif (DO_INTEROP_VIZ)
   else ()
      if (DO_INTEROP_VIZ)
         target_link_libraries(CoMDOCL stdc++ ${Z_LIB} ${OPENCL_LIBRARIES} glut GLU GLEW)
      else (DO_INTEROP_VIZ)
         target_link_libraries(CoMDOCL ${OPENCL_LIBRARIES} ${Z_LIB} ${M_LIB})
      endif (DO_INTEROP_VIZ)
   endif ()

   set_target_properties(CoMDOCL PROPERTIES COMPILE_FLAGS " ${OPT} -DUSE_OCL")

   if (APPLE)
      set_target_properties(CoMDOCL PROPERTIES LINK_FLAGS "-framework OpenCL")
   endif ()

endif ()

