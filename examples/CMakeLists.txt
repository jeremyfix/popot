# Make sure the compiler can find include files from our library.
include_directories (${CMAKE_SOURCE_DIR}/src)

# Define our examples to compile
set(USAGE_EXAMPLES example-000 example-001 example-002 example-003)


# loop over the list to compile them
foreach(f ${USAGE_EXAMPLES})
    add_executable (${f} ${f}.cc)
    target_link_libraries (${f} ${LIBS_FOR_EXAMPLE})	
endforeach(f) 

foreach(f ${USAGE_EXAMPLES})
    install(PROGRAMS ${CMAKE_BINARY_DIR}/examples/${f} 
            DESTINATION bin
            RENAME popot-${f})
endforeach(f) 

# Define the examples for the wrapper
set(PYTHON_USAGE_EXAMPLES example-000.py)

foreach(f ${PYTHON_USAGE_EXAMPLES})
    configure_file(${f} ${CMAKE_CURRENT_BINARY_DIR}/${f} COPYONLY)
    install(FILES ${CMAKE_BINARY_DIR}/examples/${f}
            DESTINATION bin
            RENAME popot-${f})
endforeach(f) 
