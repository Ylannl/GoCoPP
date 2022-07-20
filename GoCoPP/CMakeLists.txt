set(EXECUTABLE_OUTPUT_PATH ../bin/${CMAKE_BUILD_TYPE})

find_package(Boost REQUIRED COMPONENTS filesystem)
if(${Boost_FOUND})
	set(Boost_USE_STATIC_LIBS OFF)
	set(Boost_USE_MULTITHREAD OFF)
	include_directories(${Boost_INCLUDE_DIRS})
	link_directories(${Boost_LIBRARY_DIRS})
else()
	message(STATUS "NOTICE: This demo requires Boost and will not be compiled.")
endif()

find_package(CGAL)

if(${CGAL_FOUND})
	include(${CGAL_USE_FILE})
else()
	message(STATUS "NOTICE: This demo requires CGAL and will not be compiled.")
endif()



# Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})


add_definitions(-D_CONSOLE -D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS)

file(GLOB_RECURSE 
	psd_HDRS
	include/*)
	
file(GLOB_RECURSE 
	psd_SRCS
	src/*)

include_directories(include)


if (CGAL_FOUND)

	include(AddFileDependencies)
    add_executable(GoCoPP ${psd_HDRS} ${psd_SRCS})
	target_link_libraries(GoCoPP ${CGAL_LIBRARIES})
	
	

else()
	set(OT_MISSING_DEPS "")

   
    

    if(NOT CGAL_FOUND)
        set(OT_MISSING_DEPS "the CGAL library, ${OT_MISSING_DEPS}")
    endif()

    message(STATUS "NOTICE: This demo requires ${OT_MISSING_DEPS} and will not be compiled.")
endif()
