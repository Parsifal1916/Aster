#pragma once

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <sstream>

using namespace Aster;
using namespace Catch;

TEST_CASE("Simulation inits and destructors test", "[building]") {
    SECTION("2d simulations"){
        auto type = GENERATE(LIGHT, HEAVY, BARNES_HUT, BH_termal);
        auto* simulation = bake(type);

        REQUIRE(simulation);

        simulation
            -> set_dt(.1)
            -> set_max_frames(2)
        ;

        REQUIRE(simulation -> get_dt() == Approx(.1));

        free(simulation);
    }

    SECTION("3d simulations"){
        auto type = GENERATE(LIGHT, HEAVY, BARNES_HUT);
        auto* simulation = bake3d(type);

        REQUIRE(simulation);

        simulation
            -> set_dt(.1)
            -> set_max_frames(2)
        ;

        REQUIRE(simulation -> get_dt() == Approx(.1));

        free(simulation);
    }
}
