CMAKE_MINIMUM_REQUIRED (VERSION 2.8)

PROJECT(tau_leaping)

FIND_PACKAGE(CUDA REQUIRED)

SET(CUDA_NVCC_FLAGS "-arch=sm_20")

FIND_CUDA_HELPER_LIBS(cusparse)

SET(HEADERS
	cl_parser.h
	CSCHostMatrixData.h
	CSCMatrix.h
	CSRMatrix.h
	CSRMatrixPacked.h
	defines.h
	functions.h
	model.h
	ModelConst.h
	ModelState.h
	my_assert.h
	MyTime.h
	realtype.h
	statistics.h
	StochKitReader.h
	CSCHostMatrixData.cuh
	d_rng.cuh
)

SET(HOST_FILES 
	functions.cpp
	model.cpp
	statistics.cpp
	StochKitReader.cpp
	cl_parser.cpp
	CSCHostMatrixData.cpp
	tau_leaping.cpp
)

SET(CUDA_FILES
	model.cu
	ModelConst.cu
	ModelState.cu
	#test.cu
	CSRMatrixPacked.cu
	CSRMatrix.cu
	CSCMatrix.cu
)

#for making life in MSVC a bit better...
SOURCE_GROUP(CUDA FILES
	${CUDA_FILES}
)

SOURCE_GROUP("Header Files" FILES
	${HEADERS}
)

#path to rapid xml root dir
include_directories ("..")

SET(CUDA_HOST_COMPILATION_CPP OFF)

IF(CMAKE_COMPILER_IS_GNUCC)
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
ENDIF(CMAKE_COMPILER_IS_GNUCC)



#CUDA_COMPILE(CUDA_OBJ ${CUDA_FILES})

CUDA_ADD_EXECUTABLE(tau_leaping
	${HEADERS}
	${CUDA_FILES}
	${HOST_FILES}
)

TARGET_LINK_LIBRARIES (tau_leaping 
	${CUDA_cusparse_LIBRARY}
)
