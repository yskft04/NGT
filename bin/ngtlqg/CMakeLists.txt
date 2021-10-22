if( ${UNIX} )
	include_directories("${PROJECT_BINARY_DIR}/lib")
        include_directories("${PROJECT_SOURCE_DIR}/lib")
        link_directories("${PROJECT_SOURCE_DIR}/lib/NGT")

	add_executable(ngtlqg_exe ngtlqg.cpp)
	add_dependencies(ngtlqg_exe ngt)
        set_target_properties(ngtlqg_exe PROPERTIES OUTPUT_NAME ngtlqg) 
	target_link_libraries(ngtlqg_exe ngt pthread)

	install(TARGETS ngtlqg_exe RUNTIME DESTINATION bin)

endif()
