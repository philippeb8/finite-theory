# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Phil Bouchard <phil@fornux.com>

cmake_minimum_required(VERSION 3.16)
project(ftsim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)

find_package(Qt6 REQUIRED COMPONENTS Widgets)

add_executable(ftsim main.cpp)
target_link_libraries(ftsim PRIVATE Qt6::Widgets)

# Qt5 fallback:
#   find_package(Qt5 REQUIRED COMPONENTS Widgets)
#   target_link_libraries(ftsim PRIVATE Qt5::Widgets)
