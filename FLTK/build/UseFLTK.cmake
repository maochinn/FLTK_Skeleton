#
# "$Id$"
#
# UseFLTK.CMake - FLTK CMake environment configuration file for external projects.
#
# This file is deprecated and will be removed in FLTK 1.4
#
# automatically generated - do not edit
#

include_directories("D:/maochinn/Lib/FLTK/build;D:/maochinn/Lib/FLTK/fltk-1.3.5-source.tar/fltk-1.3.5-source/fltk-1.3.5")

message(STATUS "* Notice: UseFLTK.CMake is deprecated and will be removed in FLTK 1.4")
message(STATUS "*   Please use 'include_directories(\${FLTK_INCLUDE_DIRS})' or")
message(STATUS "*   'target_include_directories(<target> PUBLIC|PRIVATE \${FLTK_INCLUDE_DIRS})'")
message(STATUS "*   instead of 'include(\${FLTK_USE_FILE})'.")

#
# End of "$Id$".
#
