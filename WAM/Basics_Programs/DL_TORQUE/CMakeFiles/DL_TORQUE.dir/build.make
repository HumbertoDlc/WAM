# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/robot/crmlWAM/Humberto_Dlc/DL_TORQUE

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/robot/crmlWAM/Humberto_Dlc/DL_TORQUE

# Include any dependencies generated for this target.
include CMakeFiles/DL_TORQUE.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/DL_TORQUE.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/DL_TORQUE.dir/flags.make

CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.o: CMakeFiles/DL_TORQUE.dir/flags.make
CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.o: DL_TORQUE.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/robot/crmlWAM/Humberto_Dlc/DL_TORQUE/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.o -c /home/robot/crmlWAM/Humberto_Dlc/DL_TORQUE/DL_TORQUE.cpp

CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/robot/crmlWAM/Humberto_Dlc/DL_TORQUE/DL_TORQUE.cpp > CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.i

CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/robot/crmlWAM/Humberto_Dlc/DL_TORQUE/DL_TORQUE.cpp -o CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.s

CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.o.requires:
.PHONY : CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.o.requires

CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.o.provides: CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.o.requires
	$(MAKE) -f CMakeFiles/DL_TORQUE.dir/build.make CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.o.provides.build
.PHONY : CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.o.provides

CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.o.provides.build: CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.o

# Object files for target DL_TORQUE
DL_TORQUE_OBJECTS = \
"CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.o"

# External object files for target DL_TORQUE
DL_TORQUE_EXTERNAL_OBJECTS =

DL_TORQUE: CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.o
DL_TORQUE: CMakeFiles/DL_TORQUE.dir/build.make
DL_TORQUE: /usr/lib/x86_64-linux-gnu/libboost_system.so
DL_TORQUE: /usr/lib/x86_64-linux-gnu/libboost_thread.so
DL_TORQUE: /usr/lib/x86_64-linux-gnu/libboost_python.so
DL_TORQUE: /usr/lib/x86_64-linux-gnu/libpthread.so
DL_TORQUE: /usr/lib/libnative.so
DL_TORQUE: /usr/lib/libxenomai.so
DL_TORQUE: /usr/lib/librtdm.so
DL_TORQUE: /usr/lib/x86_64-linux-gnu/libpython2.7.so
DL_TORQUE: /usr/lib/x86_64-linux-gnu/libcurses.so
DL_TORQUE: /usr/lib/x86_64-linux-gnu/libform.so
DL_TORQUE: CMakeFiles/DL_TORQUE.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable DL_TORQUE"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/DL_TORQUE.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/DL_TORQUE.dir/build: DL_TORQUE
.PHONY : CMakeFiles/DL_TORQUE.dir/build

CMakeFiles/DL_TORQUE.dir/requires: CMakeFiles/DL_TORQUE.dir/DL_TORQUE.cpp.o.requires
.PHONY : CMakeFiles/DL_TORQUE.dir/requires

CMakeFiles/DL_TORQUE.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/DL_TORQUE.dir/cmake_clean.cmake
.PHONY : CMakeFiles/DL_TORQUE.dir/clean

CMakeFiles/DL_TORQUE.dir/depend:
	cd /home/robot/crmlWAM/Humberto_Dlc/DL_TORQUE && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/robot/crmlWAM/Humberto_Dlc/DL_TORQUE /home/robot/crmlWAM/Humberto_Dlc/DL_TORQUE /home/robot/crmlWAM/Humberto_Dlc/DL_TORQUE /home/robot/crmlWAM/Humberto_Dlc/DL_TORQUE /home/robot/crmlWAM/Humberto_Dlc/DL_TORQUE/CMakeFiles/DL_TORQUE.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/DL_TORQUE.dir/depend

