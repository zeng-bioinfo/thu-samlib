include_directories(${THUSAMLIB_SOURCE_DIR}/samlib)
link_directories(${THUSAMLIB_BINARY_DIR}/samlib)

# Samtools libraries
include_directories(${THUSAMLIB_SOURCE_DIR}/samtools)
link_directories(${THUSAMLIB_SOURCE_DIR}/samtools)

# Compressed suffix tree libraries
include_directories(${THUSAMLIB_SOURCE_DIR}/cst)
aux_source_directory(${THUSAMLIB_SOURCE_DIR}/cst CST_SRC)

# Boost libraries
find_package(Boost REQUIRED COMPONENTS system filesystem random program_options regex)
include_directories(${Boost_INCLUDE_DIR})
add_definitions("-DHAS_BOOST")

# Program Option Layout
set(PROGRAM_OPTION OptionPrinter.cpp CustomOptionDescription.cpp)

# Executable examples
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

# The executable program to extract the alignments from user's BAM file to be the training data.
add_executable(extract_train_data extract_train_data.cpp bam_file.cpp ${PROGRAM_OPTION})
target_link_libraries(extract_train_data bam z pthread ${Boost_LIBRARIES})

# The executable program to compute the statistics from alignments
set(STATSRC nucl_align.cpp nucl_align_method.cpp semi_homo_align.cpp)
add_executable(align_stat align_stat.cpp ${STATSRC} ${PROGRAM_OPTION} ${CST_SRC})
target_link_libraries(align_stat ${Boost_LIBRARIES})
