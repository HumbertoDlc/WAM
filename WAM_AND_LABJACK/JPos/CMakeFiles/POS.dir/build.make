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
CMAKE_SOURCE_DIR = /home/robot/crmlWAM/Humberto_Dlc/WAM_AND_LABJACK/JPos

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/robot/crmlWAM/Humberto_Dlc/WAM_AND_LABJACK/JPos

# Include any dependencies generated for this target.
include CMakeFiles/POS.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/POS.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/POS.dir/flags.make

CMakeFiles/POS.dir/POS.cpp.o: CMakeFiles/POS.dir/flags.make
CMakeFiles/POS.dir/POS.cpp.o: POS.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/robot/crmlWAM/Humberto_Dlc/WAM_AND_LABJACK/JPos/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/POS.dir/POS.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/POS.dir/POS.cpp.o -c /home/robot/crmlWAM/Humberto_Dlc/WAM_AND_LABJACK/JPos/POS.cpp

CMakeFiles/POS.dir/POS.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/POS.dir/POS.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/robot/crmlWAM/Humberto_Dlc/WAM_AND_LABJACK/JPos/POS.cpp > CMakeFiles/POS.dir/POS.cpp.i

CMakeFiles/POS.dir/POS.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/POS.dir/POS.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/robot/crmlWAM/Humberto_Dlc/WAM_AND_LABJACK/JPos/POS.cpp -o CMakeFiles/POS.dir/POS.cpp.s

CMakeFiles/POS.dir/POS.cpp.o.requires:
.PHONY : CMakeFiles/POS.dir/POS.cpp.o.requires

CMakeFiles/POS.dir/POS.cpp.o.provides: CMakeFiles/POS.dir/POS.cpp.o.requires
	$(MAKE) -f CMakeFiles/POS.dir/build.make CMakeFiles/POS.dir/POS.cpp.o.provides.build
.PHONY : CMakeFiles/POS.dir/POS.cpp.o.provides

CMakeFiles/POS.dir/POS.cpp.o.provides.build: CMakeFiles/POS.dir/POS.cpp.o

# Object files for target POS
POS_OBJECTS = \
"CMakeFiles/POS.dir/POS.cpp.o"

# External object files for target POS
POS_EXTERNAL_OBJECTS =

POS: CMakeFiles/POS.dir/POS.cpp.o
POS: CMakeFiles/POS.dir/build.make
POS: /usr/lib/x86_64-linux-gnu/libboost_system.so
POS: /usr/lib/x86_64-linux-gnu/libboost_thread.so
POS: /usr/lib/x86_64-linux-gnu/libboost_python.so
POS: /usr/lib/x86_64-linux-gnu/libpthread.so
POS: /usr/lib/libnative.so
POS: /usr/lib/libxenomai.so
POS: /usr/lib/librtdm.so
POS: /usr/lib/x86_64-linux-gnu/libpython2.7.so
POS: CMakeFiles/POS.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable POS"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/POS.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/POS.dir/build: POS
.PHONY : CMakeFiles/POS.dir/build

CMakeFiles/POS.dir/requires: CMakeFiles/POS.dir/POS.cpp.o.requires
.PHONY : CMakeFiles/POS.dir/requires

CMakeFiles/POS.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/POS.dir/cmake_clean.cmake
.PHONY : CMakeFiles/POS.dir/clean

CMakeFiles/POS.dir/depend:
	cd /home/robot/crmlWAM/Humberto_Dlc/WAM_AND_LABJACK/JPos && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/robot/crmlWAM/Humberto_Dlc/WAM_AND_LABJACK/JPos /home/robot/crmlWAM/Humberto_Dlc/WAM_AND_LABJACK/JPos /home/robot/crmlWAM/Humberto_Dlc/WAM_AND_LABJACK/JPos /home/robot/crmlWAM/Humberto_Dlc/WAM_AND_LABJACK/JPos /home/robot/crmlWAM/Humberto_Dlc/WAM_AND_LABJACK/JPos/CMakeFiles/POS.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/POS.dir/depend

