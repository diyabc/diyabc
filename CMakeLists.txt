cmake_minimum_required(VERSION 2.8.6)

project (diyabc)

add_subdirectory (dcmt0.6.1)

link_directories (${dcmt_BINARY_DIR})

add_subdirectory (src-JMC-C++)





