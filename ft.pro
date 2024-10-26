TEMPLATE  = app
LANGUAGE  = C++
CONFIG	 += thread c++17
QT       += core gui widgets

QMAKE_LFLAGS += -lboost_system

#QMAKE_CC              = i686-pc-mingw32-gcc
#QMAKE_CXX             = i686-pc-mingw32-g++
#QMAKE_LINK            = i686-pc-mingw32-g++
#QMAKE_LIB             = i686-pc-mingw32-ar -ru
#QMAKE_RC              = i686-pc-mingw32-windres

#QMAKE_RPATHDIR           += /opt/V-PlaySDK/5.6/gcc_64/lib/
#QMAKE_LIBDIR             += /usr/lib
#QMAKE_LFLAGS             += -static-libgcc -static-libg++
QMAKE_CXXFLAGS           += -march=native -O3

#QMAKE_CFLAGS_RELEASE += /MT
#QMAKE_CXXFLAGS_RELEASE += /MT
#QMAKE_CFLAGS_DEBUG += /MTd
#QMAKE_CXXFLAGS_DEBUG += /MTd

#QMAKE_LFLAGS_RELEASE += /FORCE:MULTIPLE
#QMAKE_LFLAGS_DEBUG += /FORCE:MULTIPLE

#QMAKE_CXXFLAGS_RELEASE += /Gy
#QMAKE_LFLAGS_RELEASE += /OPT:REF

TARGET    = ft

HEADERS	+= main.h
SOURCES	+= main.cpp
