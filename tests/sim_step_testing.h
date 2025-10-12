#pragma once

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <sstream>

using namespace Aster;
using namespace Catch;

TEST_CASE("Simulation step benchmarking", "[building]") {

    SECTION("3d Barnes"){
        auto* simulation = new Simulation();
        
        cosmic_web(simulation, 1e3, 1e4);
        simulation -> load();

        BENCHMARK("Steps 3d"){
            return simulation -> step();
        };

        free(simulation);
    }
}

TEST_CASE("update methods testing", "[parametric]") {

    for (int i = 0; i< 10; ++i) {
        SECTION("update methods testing") {
            auto* simulation = new Simulation();

            cosmic_web(simulation, 100, 10e6);
            simulation -> integrator = SABA;
            simulation-> integrator_order = i; 
            simulation -> load();

            simulation -> integrate(10);

            REQUIRE(simulation -> is_fine());

            free(simulation);
        }
    }
}