cmake_minimum_required(VERSION 2.6)
project(haplotypecpp)


set(CMAKE_CXX_FLAGS "-std=c++11")
#include_directories(/usr/local/lib64/R-3.0.1/library/Rcpp/include /home/pingu/tmp/foreign/R-3.0.1/include)
#include_directories($ENV{HOME}/R/x86_64-unknown-linux-gnu-library/3.1/Rcpp/include $ENV{HOME}/lib64/R-3.1.3/include)
#link_directories(/usr/local/lib64/R-3.0.1/lib /usr/local/lib64/R-3.0.1/library/Rcpp/lib)
include_directories($ENV{HOME}/lib64/R-3.1.3/include $ENV{HOME}/lib64/R-3.1.3/library/Rcpp/include)
link_directories($ENV{HOME}/lib64/R-3.1.3/lib  $ENV{HOME}/lib64/R-3.1.3/library/Rcpp/libs)

add_executable(haplotypecpp valarray_ext.cpp helpers.cpp buffer.cpp diplotypereconstruction.cpp main.cpp pedigree.cpp)
#add_library(haplotype_lib SHARED pedigree.cpp buffer.cpp diplotypereconstruction.cpp)
target_link_libraries(haplotypecpp R Rcpp)
#target_link_libraries(haplotype_lib R Rcpp)

install(TARGETS haplotypecpp RUNTIME DESTINATION bin)


