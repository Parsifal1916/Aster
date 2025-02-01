#pragma once

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <sstream>

using namespace Aster;
using namespace Catch;

TEST_CASE("Simulation step benchmarking", "[building]") {
    SECTION("2d simulations"){
        auto* simulation = bake(BARNES_HUT);

        cosmic_web(simulation, 1e3, 1e4);

        BENCHMARK("Steps 2d"){
            return simulation -> step();
        };

        free(simulation);
    }

    SECTION("3d simulations"){
        auto* simulation = bake3d(BARNES_HUT);

        cosmic_web(simulation, 1e3, 1e4);

        BENCHMARK("Steps 3d"){
            return simulation -> step();
        };

        free(simulation);
    }
}