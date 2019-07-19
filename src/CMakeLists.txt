add_executable(md_water
  md_water_ewald.cpp
  function.cpp
  )

target_compile_features(md_water
  PUBLIC cxx_std_14
  )

target_compile_options(md_water
  PUBLIC
    -Wall
    -Wextra
    -Wpedantic
  )

target_include_directories(md_water
  PUBLIC ${PROJECT_SOURCE_DIR}/include
  )
