cmake_minimum_required(VERSION 2.8)
project(Incompact3D Fortran)

site_name(hostname)
exec_program( ${CMAKE_Fortran_COMPILER} ARGS --version OUTPUT_VARIABLE _COMPILER_output)
string(REGEX REPLACE ".* ([0-9]*\\.[0-9]*\\.[0-9]*).*" "\\1" COMPILER_VERSION ${_COMPILER_output})
set(CMAKE_ECLIPSE_VERSION 4.2.2)


### Definition des options de compilation ###

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release)
endif()
option(openmp  "Enable OpenMP" OFF)
option(profile "Enable instrumentation of the code for profiling purpose" OFF)
option(report  "Report optimizations instead of compilation" OFF)
option(static  "Enable static linking" ON)
option(lto     "Enable link-time optimizations" ON)



### Definition des flags de compilation suivant le compilateur ###

# ifort compiler (LINUX)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS         "-fPIC -g  -traceback -init=zero")
  set(CMAKE_Fortran_FLAGS_RELEASE "-xHOST -O3 -no-prec-div -fp-model fast=2") # -prof-use -prof-dir=prof_dir/  -prof-gen
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -C -debug all -warn all,nounused -fpe0 -ftrapuv -fstack-protector-all -WB -fp-stack-check -u -gen-interfaces -stand f08 -warn nodeclarations") #  
  set(DOUBLE_FLAGS "-r8")
  set(OPTRPT_FLAGS "-qopt-report=5 -guide -diag-enable=openmp,vec,par")
  set(LTO_FLAGS    "-ipo -ipo-jobs=4")
  set(OPENMP_FLAGS "-qopenmp")
endif ()

# gfortran compiler (LINUX)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS         "-fPIC -ffree-line-length-none -g -fbacktrace  -fimplicit-none -fall-intrinsics -finit-real=zero")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -ffast-math -funsafe-math-optimizations") #
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3 -Wall -ffpe-trap=zero,overflow,invalid -Wextra  -Wimplicit-interface  -fcheck=all -Warray-temporaries -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure -pedantic -Wno-unused-parameter -Wno-unused -Wno-unused-dummy-argument -std=f2008ts")
#  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3 -Wall -ffpe-trap=zero,overflow,invalid -Wextra -finit-real=snan -fcheck=all -pedantic -std=f2008ts -Wno-unused-parameter -Wno-unused -Wno-unused-dummy-argument ")
  set(DOUBLE_FLAGS "-fdefault-real-8 -fdefault-double-8")
  set(OPTRPT_FLAGS "-ftree-vectorizer-verbose")
  set(LTO_FLAGS    "-flto ") # -combine -flto=jobserver
  set(OPENMP_FLAGS "-fopenmp")
endif ()

set(PROFILE_FLAGS "-p -pg")
set(STATIC_FLAGS "-static")


### Definitions de quelques contraintes ###

if(openmp AND CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND COMPILER_VERSION VERSION_LESS "13.0.0")
    message(STATUS "Insufficient IFort version for OpenMP\n   Deactivating it for you")
    set(openmp OFF)
    set(OPENMP_FLAGS "")
endif()
if(openmp AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND COMPILER_VERSION VERSION_LESS "4.9.0")
    message(STATUS "Insufficient GFortran version for OpenMP\n   Deactivating it for you")
    set(openmp OFF)
    set(OPENMP_FLAGS "")
endif()
if (static AND hostname STREQUAL "hulk")
    message(STATUS "Static build impossible on Hulk with IFort\n   Deactivating it for you")
    set(static OFF)
endif()

### On rajoutes les flags des options  ###

message( STATUS "")
message( STATUS "=======================================================================================")
message( STATUS "")
message( STATUS "Compilation de Incompact3D sur : " ${hostname})
message( STATUS "Fortran compiler: " ${CMAKE_Fortran_COMPILER} " - " ${CMAKE_Fortran_COMPILER_ID} " - " ${COMPILER_VERSION})
find_library(FFTW_LIBRARY NAMES fftw3)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${DOUBLE_FLAGS}")

if( openmp )
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OPENMP_FLAGS}")
  message( STATUS "With      OpenMP")
  find_library(FFTW_thread_LIBRARY NAMES fftw3_omp)
else()
  message( STATUS "Without   OpenMP")
  find_library(FFTW_thread_LIBRARY NAMES fftw3_threads)
  find_library(pthread_LIBRARY NAMES pthread)
endif()
if( profile )
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${PROFILE_FLAGS}")
  message( STATUS "With      Profiling informations")
endif()
if( report )
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OPTRPT_FLAGS}")
  message( STATUS "With      Optimization reports")
endif()
if( lto )
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${LTO_FLAGS}")
  message( STATUS "With      LTO")
endif()
if( static )
#  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${STATIC_FLAGS}")
  message( STATUS "Static    build")
else()
  message( STATUS "Dynamic   build")
endif()
if (${CMAKE_BUILD_TYPE} MATCHES "Release")
  message( STATUS "Optimized build")
  set(CMAKE_EXE_LINK_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}")
elseif (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  message( STATUS "Debugging build")
  set(CMAKE_EXE_LINK_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
endif()
message( STATUS "")
message( STATUS "=======================================================================================")
message( STATUS "")


find_path(FFTW_INCLUDES   NAMES "fftw3.f03"  PATHS ${PKG_FFTW_INCLUDE_DIRS} ${INCLUDE_INSTALL_DIR})
include_directories(${FFTW_INCLUDES})

### Définitions des dépendances  ###

if( static )
  set(BUILD_SHARED_LIBS OFF)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(LINK_SEARCH_START_STATIC ON)
  set(LINK_SEARCH_END_STATIC ON)
endif()

FILE(GLOB SRC_FILES "*.[fF]90" )
list(REMOVE_ITEM SRC_FILES "${CMAKE_CURRENT_LIST_DIR}/grille.f90")

add_executable(incompact3D ${SRC_FILES})
if( openmp )
  target_link_libraries(incompact3D ${FFTW_LIBRARY} ${FFTW_thread_LIBRARY})
else()
  target_link_libraries(incompact3D ${FFTW_LIBRARY} ${FFTW_thread_LIBRARY} ${pthread_LIBRARY})
endif()

### Autre (LTO et STATIC)  ###

set(POSITION_INDEPENDENT_CODE ON)
if( static )
  set_target_properties(incompact3D PROPERTIES LINK_SEARCH_START_STATIC ON)
  set_target_properties(incompact3D PROPERTIES LINK_SEARCH_END_STATIC   ON)
  set_property(TARGET   incompact3D PROPERTY   LINK_SEARCH_START_STATIC ON)
  set_property(TARGET   incompact3D PROPERTY   LINK_SEARCH_END_STATIC   ON)
  set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS   "")
  set(CMAKE_SHARED_LINK_DYNAMIC_Fortran_FLAGS   "")
  set(CMAKE_EXE_EXPORTS_Fortran_FLAG            "")
  set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "")
endif()
if( lto AND ${CMAKE_BUILD_TYPE} MATCHES "Release")
    set_target_properties(incompact3D PROPERTIES  INTERPROCEDURAL_OPTIMIZATION ON)
    set_property(TARGET   incompact3D PROPERTY    INTERPROCEDURAL_OPTIMIZATION ON)
endif()
