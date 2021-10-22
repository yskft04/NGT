if( ${UNIX} )
	include_directories("${PROJECT_BINARY_DIR}/lib")
        include_directories("${PROJECT_SOURCE_DIR}/lib")
        link_directories("${PROJECT_SOURCE_DIR}/lib/NGT")

	add_executable(neurips21 neurips21.cpp)
	add_dependencies(neurips21 ngt)
	target_link_libraries(neurips21 ngt pthread)

	install(TARGETS neurips21 RUNTIME DESTINATION bin)
endif()
