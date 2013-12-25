include_directories(${THUSAMLIB_SOURCE_DIR}/samlib)
link_directories(${THUSAMLIB_BINARY_DIR}/samlib)

#core libraries
set(CORELIBS bam m pthread)

#boost libraries
set(BOOSTLIBS boost_filesystem)

#BAM header parser executable example
add_executable(test_bam_header_parser bam_file.cpp bam_header_parser.cpp)
target_link_libraries(test_bam_header_parser ${CORELIBS} ${BOOSTLIBS}) 
