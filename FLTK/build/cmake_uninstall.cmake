if(NOT EXISTS "D:/maochinn/Lib/FLTK/build/install_manifest.txt")
   message(FATAL_ERROR
      "Cannot find install manifest: \"D:/maochinn/Lib/FLTK/build/install_manifest.txt\"")
endif(NOT EXISTS "D:/maochinn/Lib/FLTK/build/install_manifest.txt")

file(READ "D:/maochinn/Lib/FLTK/build/install_manifest.txt" files)
string(REGEX REPLACE "\n" ";" files "${files}")

foreach(file ${files})
message(STATUS "Uninstalling \"$ENV{DESTDIR}${file}\"")
   exec_program("D:/install/CMake/bin/cmake.exe"
      ARGS "-E remove -f \"$ENV{DESTDIR}${file}\""
      OUTPUT_VARIABLE rm_out
      RETURN_VALUE rm_retval
   )
   if(NOT "${rm_retval}" STREQUAL 0)
      message(FATAL_ERROR "Problem when removing \"$ENV{DESTDIR}${file}\"")
   endif(NOT "${rm_retval}" STREQUAL 0)
endforeach(file)
