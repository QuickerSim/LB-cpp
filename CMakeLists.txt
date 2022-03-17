######################################################################
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(LatticeBoltzmann LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

list(APPEND CMAKE_MODULE_PATH "${PRJ_SOURCE_DIR}")

message(STATUS "-------------------------------------------------------------------------------" )
message(STATUS "Project name :   ${PROJECT_NAME}")
message(STATUS "-------------------------------------------------------------------------------" )

message(STATUS)
message(STATUS "-------------------------------------------------------------------------------" )
message(STATUS "Prepering build for System : ${CMAKE_SYSTEM_NAME}")

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

######################################################################
############################ DOXYGEN #################################
######################################################################
if(CMAKE_BUILD_TYPE MATCHES Release)
	# first we can indicate the documentation build as an option and set it to ON by default
	option(BUILD_DOC "Build documentation" ON)

	# check if Doxygen is installed
	find_package(Doxygen)
	if (DOXYGEN_FOUND)
		# set input and output files
		set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
		set(DOXYGEN_OUT ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.out)

		# request to configure the file
		configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
		message("Doxygen build started")

		# note the option ALL which allows to build the docs together with the application
		add_custom_target( doc_doxygen ALL
			COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc/
			COMMENT "Generating API documentation with Doxygen"
			VERBATIM )
	else (DOXYGEN_FOUND)
		message("Doxygen need to be installed to generate the doxygen documentation")
	endif (DOXYGEN_FOUND)

endif(CMAKE_BUILD_TYPE MATCHES Release)
######################################################################
############################## GIT ###################################
######################################################################
find_package(Git)
if(GIT_FOUND)
	message("Git found: ${GIT_EXECUTABLE}")
endif()

######################################################################
########################### BUILD TYPE ###############################
######################################################################

message(STATUS "Build Type is: ${CMAKE_BUILD_TYPE}")

########################### PROFILER ####################################
if(CMAKE_BUILD_TYPE MATCHES Profiler)
	find_package(easy_profiler REQUIRED)

	message(STATUS "Adding Profiler flags")
	set(CMAKE_CXX_FLAGS_OPTIM "-O3 -pipe -funroll-loops")
	set(CMAKE_C_FLAGS_OPTIM "-O3 -pipe")
	set(CMAKE_EXE_LINKER_FLAGS_OPTIM "-O3")
	add_definitions(-DPROFILER)
endif(CMAKE_BUILD_TYPE MATCHES Profiler)

########################### DEBUG ####################################
if(CMAKE_BUILD_TYPE MATCHES Debug)
	message(STATUS "Adding Debug flags")
	add_definitions(-DDEBUG)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

########################### RELEASE ####################################
if(CMAKE_BUILD_TYPE MATCHES Release)
	# OpenMP
	find_package(OpenMP)

	if(OpenMP_FOUND)
		message(STATUS "Adding OpenMP flags")

		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	endif()

	message(STATUS "Adding Release flags")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pipe -funroll-loops")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -pipe")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O3")
endif(CMAKE_BUILD_TYPE MATCHES Release)

######################################################################
#####################    PROJECT   OPTIONS   #########################
######################################################################

set(PRJ_SOURCE_DIR ${PROJECT_SOURCE_DIR})
set(PRJ_BINARY_DIR ${PROJECT_BINARY_DIR})

find_path(EXTLIB_DIR
	NAMES "Metis" "Eigen" "CGNS"
	PATHS "/home/konrad/Projekty/"
)

######################################################################
############   project instalation options   #########################
######################################################################
#SET ( INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
#SET ( INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/${Flow2_INSTALL_INCLUDE_DIR})
#SET ( INSTALL_LIB_DIR ${CMAKE_INSTALL_PREFIX}/${Flow2_INSTALL_LIB_DIR})
#SET ( INSTALL_BIN_DIR ${CMAKE_INSTALL_PREFIX}/${Flow2_INSTALL_RUNTIME_DIR})

add_subdirectory("src")

######################################################################
# print info
######################################################################

string(TOUPPER "${CMAKE_BUILD_TYPE}" MY_CMAKE_BUILD_TYPE)

message(STATUS)
message(STATUS "-------------------------------------------------------------------------------")
message(STATUS "CMAKE version       = ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}")
message(STATUS "CMAKE_SYSTEM        = ${CMAKE_SYSTEM}")
message(STATUS "CMAKE_GENERATOR     = ${CMAKE_GENERATOR}")
message(STATUS "CMAKE_BUILD_TYPE    = ${CMAKE_BUILD_TYPE}")
message(STATUS "-------------------------------------------------------------------------------")
message(STATUS "CMAKE_C_COMPILER    = ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER  = ${CMAKE_CXX_COMPILER}")
message(STATUS "-------------------------------------------------------------------------------")
message(STATUS "FLAGS :")
message(STATUS)
message(STATUS "CMAKE_C_FLAGS          = ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS        = ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_EXE_LINKER_FLAGS = ${CMAKE_EXE_LINKER_FLAGS}")
message(STATUS "-------------------------------------------------------------------------------")

#message( STATUS "Installation")
#message( STATUS "Installation_DIR = ${INSTALL_DIR}")
#message( STATUS "Installation_Include = ${INSTALL_INCLUDE_DIR}")
#message( STATUS "Installation_lib = ${INSTALL_LIB_DIR}")
#message( STATUS "Installation_bin = ${INSTALL_BIN_DIR}")

message(STATUS "-------------------------------------------------------------------------------")
