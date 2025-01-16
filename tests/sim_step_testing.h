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
        auto type = GENERATE(LIGHT, HEAVY, BARNES_HUT, BH_termal);
        auto* simulation = bake(type);

        presets::cosmic_web(simulation, 1e3, 1e4);

        BENCHMARK("Steps 2d"){
            return simulation -> step();
        };

        free(simulation);
    }

    SECTION("3d simulations"){
        auto type = GENERATE(LIGHT, HEAVY, BARNES_HUT);
        auto* simulation = bake3d(type);

        presets::cosmic_web3d(simulation, 1e3, 1e4);

        BENCHMARK("Steps 3d"){
            return simulation -> step();
        };

        free(simulation);
    }
}