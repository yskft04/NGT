if( ${UNIX} )
	include_directories("${PROJECT_BINARY_DIR}/lib")
        include_directories("${PROJECT_SOURCE_DIR}/lib")
        link_directories("${PROJECT_SOURCE_DIR}/lib/NGT")

	add_executable(opq_exe opq.cpp)
	add_dependencies(opq_exe ngt)
        set_target_properties(opq_exe PROPERTIES OUTPUT_NAME opq) 
	target_link_libraries(opq_exe ngt pthread lapack)

	install(TARGETS opq_exe RUNTIME DESTINATION bin)

endif()
