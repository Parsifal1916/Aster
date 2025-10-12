#pragma once

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <sstream>

using namespace Aster;
using namespace Catch;

TEST_CASE("Simulation inits and destructors test", "[building]") {
    SECTION("simulations"){
        auto type = GENERATE(SINGLE_THREAD, PARALLEL, BARNES_HUT);
        auto* simulation = new Simulation();

        simulation -> gravity_solver = type; 

        REQUIRE(simulation);

        simulation
            -> set_dt(.1)
            -> set_max_frames(2)
        ;

        REQUIRE(simulation -> get_dt() == Approx(.1));

        free(simulation);
    }

}
