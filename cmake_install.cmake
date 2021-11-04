# Install script for directory: /home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/iqtree" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/iqtree")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/iqtree"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/iqtree")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/iqtree" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/iqtree")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/iqtree"
         OLD_RPATH "/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/home/jm17923/Software/arm/compiler/20.2/gcc-9.3.0_Generic-AArch64_SLES-12_aarch64-linux/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/iqtree")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/." TYPE FILE FILES "/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/example/models.nex")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/." TYPE FILE FILES "/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/example/example.phy")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/." TYPE FILE FILES "/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/example/example.nex")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/." TYPE FILE FILES "/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/example/example.cf")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/pll/cmake_install.cmake")
  include("/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/ncl/cmake_install.cmake")
  include("/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/nclextra/cmake_install.cmake")
  include("/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/utils/cmake_install.cmake")
  include("/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/pda/cmake_install.cmake")
  include("/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/lbfgsb/cmake_install.cmake")
  include("/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/whtest/cmake_install.cmake")
  include("/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/sprng/cmake_install.cmake")
  include("/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/vectorclass/cmake_install.cmake")
  include("/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/model/cmake_install.cmake")
  include("/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/gsl/cmake_install.cmake")
  include("/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/alignment/cmake_install.cmake")
  include("/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/tree/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/jm17923/acrc_internship/summer20/iqtree/IQ-TREE/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
