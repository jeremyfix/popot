# Make sure the compiler can find include files from our library.
include_directories (${CMAKE_SOURCE_DIR}/src)

# Define our examples to compile
file(
	GLOB 
	usage_examples
	*.cc
)

file(
	GLOB 
	pyusage_examples
	*.py
)

# loop over the list to compile them
foreach(f ${usage_examples})
    get_filename_component(exampleName ${f} NAME_WE) 
    add_executable (${exampleName} ${f}) 
    set_target_properties(${exampleName} PROPERTIES LINKER_LANGUAGE CXX)
    set_target_properties(${exampleName} PROPERTIES COMPILE_FLAGS "${PROJECT_ALL_CFLAGS}" LINK_FLAGS "${PROJECT_ALL_LDFLAGS}")

    # If the project generates libs on which the examples depend
    # we must ensure that the libs are compiled before the examples (especially on multiple cores)
    add_dependencies(${exampleName} popot)	 

    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${exampleName} 
            DESTINATION bin
            RENAME ${CMAKE_PROJECT_NAME}-${exampleName}
	    COMPONENT binary)
endforeach(f)

# Define the examples for the wrapper

foreach(f ${pyusage_examples})
    get_filename_component(fname ${f} NAME) 
    configure_file(${fname} ${CMAKE_CURRENT_BINARY_DIR}/${fname} COPYONLY)
    install(FILES ${fname}
            DESTINATION bin
            RENAME ${CMAKE_PROJECT_NAME}-${fname}
	    COMPONENT binary)
endforeach(f) 
