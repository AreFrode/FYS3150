#include <iostream>

#include "../cpp_scripts/planet.hpp"
#include "../cpp_scripts/solarsystem.hpp"

void test_energy_angular_momentum();

int main(int argc, char *argv[]) {
    test_energy_angular_momentum();
    std::cout << "Total energy test passed" << std::endl;
    std::cout << "Angular momentum test passed" << std::endl;

    return 0;
}

/**
* This unit test tests both the concsertvation of
* total energy, as well as conservation of angular
* momentum. Uses a modified velocityVerlet loop,
* that voids writing to std.out.
* Still, as the solver has to be initialized with
* filename String, they are replaced by "dummy"-strings
*/
void test_energy_angular_momentum() {
    Planet earth(0.000003, 1., 0., 0., 0.0, 2*M_PI, 0.);
    Planet sun(1., 0., 0., 0., 0., 0., 0.);

    std::vector<Planet> planets = {earth, sun};
    SolarSystem test(planets, 3, "foo", "bar", "foobar");

    test.test_velocityVerlet(10000, 50, 1);

}


