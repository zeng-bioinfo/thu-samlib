include_directories(${THUSAMLIB_SOURCE_DIR}/samlib)
link_directories(${THUSAMLIB_BINARY_DIR}/samlib)

#samtools libraries 
include_directories(${THUSAMLIB_SOURCE_DIR}/samtools)
link_directories(${THUSAMLIB_SOURCE_DIR}/samtools)
#boost libraries
find_package(Boost REQUIRED COMPONENTS system filesystem random)
include_directories(${Boost_INCLUDE_DIR})
add_definitions("-DHAS_BOOST")

set(EXAMPLES bam_header_parser bam_random_retrieve)
foreach(example ${EXAMPLES})
    add_executable(${example} bam_file.cpp ${example}.cpp)
    target_link_libraries(${example} bam z pthread ${Boost_LIBRARIES})
endforeach(example)
#add_executable(test_bam_header_parser bam_file.cpp bam_header_parser.cpp)
#target_link_libraries(test_bam_header_parser bam z pthread ${Boost_LIBRARIES})

#calDepth executable example
add_executable(calDepth calDepth.cpp)
target_link_libraries(calDepth bam z pthread) 
