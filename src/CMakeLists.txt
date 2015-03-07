# Generate the list of source files
# using GLOB does not RECURSE
file(
	GLOB
	lib_files
	*.cpp
)
list(REMOVE_ITEM lib_files "${CMAKE_CURRENT_SOURCE_DIR}/wrapper.cpp")

file(
	GLOB
	headers_files
	*.h
)

# Make sure the compiler can find include files from our library.
include_directories (${CMAKE_CURRENT_SOURCE_DIR})


# Define our library 
add_library(
	popot
	SHARED
	${lib_files}
)
set_target_properties(popot PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(popot PROPERTIES COMPILE_FLAGS "${PROJECT_ALL_CFLAGS}" LINK_FLAGS "${PROJECT_ALL_DEP_LDFLAGS}")

# Installation targets
install (TARGETS popot
	 DESTINATION ${LIB_INSTALL_DIR}
         COMPONENT binary)

install (FILES ${headers_files}
	 DESTINATION ${INCLUDE_INSTALL_DIR}
	 COMPONENT devel)

# Define the library for the python wrapper
include_directories(${PYTHON_INCLUDE_DIRS})
add_library(
  PyPopot
  SHARED
  wrapper.cpp ${headers_files}
)
set_target_properties(PyPopot PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(PyPopot PROPERTIES COMPILE_FLAGS "${PROJECT_ALL_CFLAGS}" LINK_FLAGS "${PROJECT_ALL_DEP_LDFLAGS}")

target_link_libraries (PyPopot ${PYTHON_LIBS} "-lboost_python")
install (TARGETS PyPopot
         DESTINATION ${PYTHON_SITE_PACKAGES}
	 COMPONENT binary)