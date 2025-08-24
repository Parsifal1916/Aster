#pragma once

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <sstream>

using namespace Aster;
using namespace Catch;

TEST_CASE("Simulation step benchmarking", "[building]") {
    SECTION("2d Barnes"){
        auto* simulation = bake(BARNES_HUT);

        cosmic_web(simulation, 1e3, 1e4);
        simulation -> load();

        BENCHMARK("Steps 2d"){
            return simulation -> step();
        };

        free(simulation);
    }

    SECTION("3d Barnes"){
        auto* simulation = bake3d(BARNES_HUT);
        
        cosmic_web(simulation, 1e3, 1e4);
        simulation -> load();

        BENCHMARK("Steps 3d"){
            return simulation -> step();
        };

        free(simulation);
    }
}
/*
TEST_CASE("GPU Simulation step benchmarking", "[building]") {
    SECTION("GPU 2d Barnes"){
        auto* simulation = new Barnes::BHG <vec2>;

        cosmic_web(simulation, 1e3, 1e4);
        simulation -> load();

        BENCHMARK("Steps 2d"){
            return simulation -> step();
        };

        free(simulation);
    }

    SECTION("GPU 3d Barnes"){
        auto* simulation = new Barnes::BHG <vec2>;
        
        cosmic_web(simulation, 1e3, 1e4);
        simulation -> load();

        BENCHMARK("Steps 3d"){
            return simulation -> step();
        };

        free(simulation);
    }
}*/