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
CMAKE_SOURCE_DIR = /home/robot/crmlWAM/GROUP/JPos

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/robot/crmlWAM/GROUP/JPos

# Include any dependencies generated for this target.
include CMakeFiles/NN.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/NN.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/NN.dir/flags.make

CMakeFiles/NN.dir/NN.cpp.o: CMakeFiles/NN.dir/flags.make
CMakeFiles/NN.dir/NN.cpp.o: NN.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/robot/crmlWAM/GROUP/JPos/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/NN.dir/NN.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/NN.dir/NN.cpp.o -c /home/robot/crmlWAM/GROUP/JPos/NN.cpp

CMakeFiles/NN.dir/NN.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NN.dir/NN.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/robot/crmlWAM/GROUP/JPos/NN.cpp > CMakeFiles/NN.dir/NN.cpp.i

CMakeFiles/NN.dir/NN.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NN.dir/NN.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/robot/crmlWAM/GROUP/JPos/NN.cpp -o CMakeFiles/NN.dir/NN.cpp.s

CMakeFiles/NN.dir/NN.cpp.o.requires:
.PHONY : CMakeFiles/NN.dir/NN.cpp.o.requires

CMakeFiles/NN.dir/NN.cpp.o.provides: CMakeFiles/NN.dir/NN.cpp.o.requires
	$(MAKE) -f CMakeFiles/NN.dir/build.make CMakeFiles/NN.dir/NN.cpp.o.provides.build
.PHONY : CMakeFiles/NN.dir/NN.cpp.o.provides

CMakeFiles/NN.dir/NN.cpp.o.provides.build: CMakeFiles/NN.dir/NN.cpp.o

# Object files for target NN
NN_OBJECTS = \
"CMakeFiles/NN.dir/NN.cpp.o"

# External object files for target NN
NN_EXTERNAL_OBJECTS =

NN: CMakeFiles/NN.dir/NN.cpp.o
NN: CMakeFiles/NN.dir/build.make
NN: /usr/lib/x86_64-linux-gnu/libboost_system.so
NN: /usr/lib/x86_64-linux-gnu/libboost_thread.so
NN: /usr/lib/x86_64-linux-gnu/libboost_python.so
NN: /usr/lib/x86_64-linux-gnu/libpthread.so
NN: /usr/lib/libnative.so
NN: /usr/lib/libxenomai.so
NN: /usr/lib/librtdm.so
NN: /usr/lib/x86_64-linux-gnu/libpython2.7.so
NN: CMakeFiles/NN.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable NN"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/NN.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/NN.dir/build: NN
.PHONY : CMakeFiles/NN.dir/build

CMakeFiles/NN.dir/requires: CMakeFiles/NN.dir/NN.cpp.o.requires
.PHONY : CMakeFiles/NN.dir/requires

CMakeFiles/NN.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/NN.dir/cmake_clean.cmake
.PHONY : CMakeFiles/NN.dir/clean

CMakeFiles/NN.dir/depend:
	cd /home/robot/crmlWAM/GROUP/JPos && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/robot/crmlWAM/GROUP/JPos /home/robot/crmlWAM/GROUP/JPos /home/robot/crmlWAM/GROUP/JPos /home/robot/crmlWAM/GROUP/JPos /home/robot/crmlWAM/GROUP/JPos/CMakeFiles/NN.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/NN.dir/depend
