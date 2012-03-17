cmake_minimum_required(VERSION 2.8)

project(toy_system_gen)

#for rapid xml
INCLUDE_DIRECTORIES("..")

if(CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

add_executable(toy_system_gen main.cpp)

