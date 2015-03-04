CMAKE_MINIMUM_REQUIRED( VERSION 2.8.4 )

project (popot) 

set(POPOT_VERSION "2.12")
set(PACKAGE_VERSION ${POPOT_VERSION})

#set( CMAKE_BUILD_TYPE Debug)
set( CMAKE_BUILD_TYPE Release)

###################################
# Get required and optional library
###################################

FIND_PACKAGE(PkgConfig REQUIRED)
pkg_check_modules(PYTHON python)
# In order to get the path to the site-packages dir to install the wrapper
execute_process ( COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
MESSAGE(${PYTHON_SITE_PACKAGES})

###################################
# Some prefix definitions
###################################

# Define the prefix where libs are going to be installed
SET (LIB_SUFFIX "" CACHE STRING "Define suffix of directory name (32/64)" )
SET (LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}" CACHE STRING "Library directory name")

# Define the path to install the headers
SET(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/popot/" CACHE PATH "The subdirectory to the header prefix")

# Define the path to install the documentation
SET(DOC_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/share/doc/popot" CACHE PATH "The subdirectory to the documentation prefix")

# Define the path to install the pc pkg-config file
SET(PKGCONFIG_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}/pkgconfig/" CACHE STRING "Base directory for pkgconfig files")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++0x" )

###################################
# Pkg-config stuff
###################################
# Let's replace the macros of the pc file
# The @ONLY means replacing only the @VAR@ variables, not the for example ${VAR}
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/popot.pc.cmake 
               ${CMAKE_CURRENT_BINARY_DIR}/popot.pc @ONLY)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/popot.pc
    	DESTINATION ${PKGCONFIG_INSTALL_DIR})


##################################
#  Subdirectories
###################################

add_subdirectory (src)
add_subdirectory (examples)
add_subdirectory (doc)

###################################
#  Custom command for source tgz
###################################

SET(DIST_DIR "${CMAKE_PROJECT_NAME}-${PACKAGE_VERSION}")
ADD_CUSTOM_TARGET(dist 
		  COMMAND rm -rf ${DIST_DIR}
		  COMMAND mkdir  ${DIST_DIR}
		  COMMAND cp -r ${CMAKE_SOURCE_DIR}/* ${DIST_DIR} || true 
		  COMMAND rm -rf ${DIST_DIR}/build
		  COMMAND mkdir ${DIST_DIR}/build
  		  COMMAND tar --exclude="*~" --exclude="._*" --exclude=".svn*" -zcvf ${DIST_DIR}.tar.gz ${DIST_DIR}
    	          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

###################################
# For building packages
###################################

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "C++ POPulation based Optimization Toolbox")
SET(CPACK_PACKAGE_VENDOR "Jeremy Fix")
SET(CPACK_PACKAGE_CONTACT "Jeremy.Fix@gmail.com")
SET(CPACK_PACKAGE_VERSION ${PACKAGE_VERSION})
SET(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_SOURCE_DIR}/LICENSE)
SET(CPACK_RESOURCE_FILE_README ${CMAKE_SOURCE_DIR}/README)
SET(CPACK_RESOURCE_FILE_WELCOME ${CMAKE_SOURCE_DIR}/README)
SET(CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION}-${CMAKE_SYSTEM}")

SET(CPACK_GENERATOR "RPM")
 
INCLUDE(CPack)
