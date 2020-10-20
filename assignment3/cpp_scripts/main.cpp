#include <string>

#include "planet.hpp"
#include "solarsystem.hpp"

using namespace std;

int main(int argc, char* argv[]) {

    if (argc < 4) {
        cout << "missing command line input!" << endl;
        cout << "Supply mesh_points, time and program when prompted" << endl;
        exit(1);
    }

    int mesh_points = atoi(argv[1]);
    int final_time = atof(argv[2]);
    string prog = string(argv[3]);

    if (prog == "binary") {

        Planet earth(0.000003, 1., 0., 0., 0., 6.3, 0.);
        Planet sun(1., 0., 0., 0., 0., 0., 0.);

        vector<Planet> planets = {earth, sun};

        SolarSystem milkyroad(planets, 3);
        milkyroad.velocityVerlet(mesh_points, final_time, 1);
        //milkyroad.forwardEuler(mesh_points, final_time);

        cout << earth.get_pos(0) << endl;
        cout << earth.get_pos(1) << endl;
        cout << earth.get_pos(2) << endl;
    }

    return 0;
}
