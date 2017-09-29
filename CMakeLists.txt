
cmake_minimum_required(VERSION 3.5.1)
project(example)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DBOOST_TEST_MAIN")
#-DBOOST_TEST_ALTERNATIVE_INIT_API

ADD_LIBRARY(boost_unit_test_framework STATIC IMPORTED)
SET_TARGET_PROPERTIES(boost_unit_test_framework PROPERTIES
    IMPORTED_LOCATION /usr/lib/x86_64-linux-gnu/libboost_unit_test_framework.a)


add_executable(simple-test unittest.cpp)
target_link_libraries(simple-test boost_unit_test_framework)