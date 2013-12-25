include_directories(${THUSAMLIB_SOURCE_DIR}/samlib)
link_directories(${THUSAMLIB_BINARY_DIR}/samlib)

#core libraries
set(CORELIBS bam m pthread)

#building just a library
add_library(bam_file bam_file.cpp)

#BAM header parser executable example
add_executable(test_bam_header_parser bam_header_parser.cpp)
target_link_libraries(test_bam_header_parser ${CORELIBS} bam_file) 
