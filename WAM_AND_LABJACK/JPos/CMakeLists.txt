cmake_minimum_required(VERSION 2.6)
project(libbarrett_examples)

########################################################################
## JACK
# set the path to the library folder - here the installer apparently installed the .so files
link_directories(/usr/local/lib)
add_executable(POS POS.cpp)
target_link_libraries(POS ljacklm)
########################################################################



# By default, build programs with compiler optimizations and debugging symbols
if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()


## libbarrett
find_package(Barrett REQUIRED)
include_directories(${BARRETT_INCLUDE_DIRS})
link_directories(${BARRETT_LIBRARY_DIRS})
add_definitions(${BARRETT_DEFINITIONS})

## curses
find_package(Curses REQUIRED)
include_directories(${CURSES_INCLUDE_DIR})


function(add_programs)
	foreach (prog ${ARGV})
		add_executable(${prog} ${prog}.cpp)
		target_link_libraries(${prog} ${BARRETT_LIBRARIES})
	endforeach()
endfunction()

target_link_libraries(POS ${BARRETT_LIBRARIES})






# Most new programs will use the ncurses library
#target_link_libraries(Basic1 ${CURSES_LIBRARIES})

