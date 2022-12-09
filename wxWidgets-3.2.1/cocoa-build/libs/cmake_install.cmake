# Install script for directory: /Users/f006pfk/source/SecondVoicePc/wxWidgets-3.2.1/build/cmake/lib

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
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
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

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Library/Developer/CommandLineTools/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/f006pfk/source/SecondVoicePc/wxWidgets-3.2.1/cocoa-build/libs/base/cmake_install.cmake")
  include("/Users/f006pfk/source/SecondVoicePc/wxWidgets-3.2.1/cocoa-build/libs/net/cmake_install.cmake")
  include("/Users/f006pfk/source/SecondVoicePc/wxWidgets-3.2.1/cocoa-build/libs/core/cmake_install.cmake")
  include("/Users/f006pfk/source/SecondVoicePc/wxWidgets-3.2.1/cocoa-build/libs/adv/cmake_install.cmake")
  include("/Users/f006pfk/source/SecondVoicePc/wxWidgets-3.2.1/cocoa-build/libs/aui/cmake_install.cmake")
  include("/Users/f006pfk/source/SecondVoicePc/wxWidgets-3.2.1/cocoa-build/libs/html/cmake_install.cmake")
  include("/Users/f006pfk/source/SecondVoicePc/wxWidgets-3.2.1/cocoa-build/libs/propgrid/cmake_install.cmake")
  include("/Users/f006pfk/source/SecondVoicePc/wxWidgets-3.2.1/cocoa-build/libs/ribbon/cmake_install.cmake")
  include("/Users/f006pfk/source/SecondVoicePc/wxWidgets-3.2.1/cocoa-build/libs/richtext/cmake_install.cmake")
  include("/Users/f006pfk/source/SecondVoicePc/wxWidgets-3.2.1/cocoa-build/libs/webview/cmake_install.cmake")
  include("/Users/f006pfk/source/SecondVoicePc/wxWidgets-3.2.1/cocoa-build/libs/stc/cmake_install.cmake")
  include("/Users/f006pfk/source/SecondVoicePc/wxWidgets-3.2.1/cocoa-build/libs/xrc/cmake_install.cmake")
  include("/Users/f006pfk/source/SecondVoicePc/wxWidgets-3.2.1/cocoa-build/libs/media/cmake_install.cmake")
  include("/Users/f006pfk/source/SecondVoicePc/wxWidgets-3.2.1/cocoa-build/libs/gl/cmake_install.cmake")
  include("/Users/f006pfk/source/SecondVoicePc/wxWidgets-3.2.1/cocoa-build/libs/qa/cmake_install.cmake")
  include("/Users/f006pfk/source/SecondVoicePc/wxWidgets-3.2.1/cocoa-build/libs/xml/cmake_install.cmake")

endif()

