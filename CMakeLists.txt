cmake_minimum_required(VERSION 3.27)
project(coal_smash)

set(CMAKE_CXX_STANDARD 20)

add_executable(coal_smash
        include/config_parser.h
        include/particle.h
        include/read_config.h
        include/read_particle_data.h
        src/read_config.cpp
        src/read_particle_data.cpp
        src/coal_deutron.cpp
        include/coal_deutron.h
        src/main.cpp
        include/output_file.h
        src/output_file.cpp
        src/particle.cpp
        include/coal_alpha.h
        src/coal_alpha.cpp
        include/general_function.h
        src/general_function.cpp
        include/coal_Be.h
        src/coal_Be.cpp
        include/run.h
        src/run.cpp
)
