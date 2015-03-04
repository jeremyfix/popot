# Generate the list of source files
# using GLOB does not RECURSE
file(
	GLOB
	headers_files
	*.h
)

# Define our library 
add_library(
	popot
	SHARED
	main.cpp ${headers_files}
)

# Installation targets
install (TARGETS popot
	 DESTINATION ${LIB_INSTALL_DIR})
install (FILES ${headers_files} 
	 DESTINATION ${INCLUDE_INSTALL_DIR})


# Define the library for the python wrapper
include_directories(${PYTHON_INCLUDE_DIRS})
add_library(
	PyPopot
	SHARED
	wrapper.cpp ${headers_files}
)
target_link_libraries (PyPopot ${PYTHON_LIBS} "-lboost_python")
install (TARGETS PyPopot
	 DESTINATION ${PYTHON_SITE_PACKAGES})
