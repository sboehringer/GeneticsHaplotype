cmake_minimum_required(VERSION 2.6)
project(haplotypecpp)


set(CMAKE_CXX_FLAGS "-std=c++0x")
include_directories(/usr/local/lib64/R-3.0.1/library/Rcpp/include /home/pingu/tmp/foreign/R-3.0.1/include)
link_directories(/usr/local/lib64/R-3.0.1/lib /usr/local/lib64/R-3.0.1/library/Rcpp/lib)

add_executable(haplotypecpp helpers.cpp buffer.cpp diplotypereconstruction.cpp main.cpp pedigree.cpp)
#add_library(haplotype_lib SHARED pedigree.cpp buffer.cpp diplotypereconstruction.cpp)
target_link_libraries(haplotypecpp R Rcpp)
#target_link_libraries(haplotype_lib R Rcpp)

install(TARGETS haplotypecpp RUNTIME DESTINATION bin)


