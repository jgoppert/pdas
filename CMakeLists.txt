cmake_minimum_required(VERSION 2.8)

project(pdas)

enable_language(Fortran)

add_executable(a5022wgs 2wgs/a5022wgs.f90)
add_executable(wb2wgs 2wgs/wb2wgs.f90)
add_executable(wd2wgs 2wgs/wd2wgs.f90)

add_executable(3view 3view/3view.f90)

add_library(atmos76 atmos/atmos76.f90)

add_executable(hlp hlp/hlp.f90)

add_executable(panair panair/panair.f90)

add_executable(panin panin/panin.f90)

add_executable(wgs2wrl vrml/wgs2wrl.f90)
