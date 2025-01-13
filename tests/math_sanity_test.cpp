#define CATCH_CONFIG_MAIN

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <Aster.hpp>

using namespace Aster;
using namespace Catch;

TEST_CASE("Vector magnitude", "[vec2]") {
    CHECK(vec2(4, 3).magnitude() == Approx(5));
    CHECK(vec2(5, 12).magnitude() == Approx(13));

    CHECK(vec2(-4, 3).magnitude() == Approx(5));
    CHECK(vec2(4, -3).magnitude() == Approx(5));
    CHECK(vec2(-4, -3).magnitude() == Approx(5));
}

TEST_CASE("Vector magnitude", "[vec3]") {
    CHECK(vec3(1, 2, 2).magnitude() == Approx(3));
    CHECK(vec3(1, 2, 2).magnitude() == Approx(3));

    CHECK(vec3(-1, 2, 2).magnitude() == Approx(3));
    CHECK(vec3(1, -2, 2).magnitude() == Approx(3));
    CHECK(vec3(-1, -2, 2).magnitude() == Approx(3));
}