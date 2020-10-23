#include <iostream>
#include <ctime>
#include <iomanip>

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
    string fname = prog + "_N_" + to_string(mesh_points) + "_t_" + to_string(final_time) + ".txt";
    string fname_E = prog + "_N_" + to_string(mesh_points) + "_t_" + to_string(final_time) + "_energy.txt";
    string fname_M = prog + "_N_" + to_string(mesh_points) + "_t_" + to_string(final_time) + "_moment.txt";
    clock_t start, end;

    if (prog == "earth-vv") {

        Planet earth(0.000003, 1., 0., 0., 0.0, 2*M_PI, 0.);
        Planet sun(1., 0., 0., 0., 0., 0., 0.);

        vector<Planet> planets = {earth, sun};

        SolarSystem binary(planets, 3, fname, fname_E, fname_M);

        start = clock();
        binary.velocityVerlet(mesh_points, final_time, planets.size()-1);
        end = clock();
        cout << setw(15) << setprecision(8) << "Time taken velocity Verlet Earth-Sun: " << 1000.*((end - start)/(double)CLOCKS_PER_SEC) << "ms" << endl;

    } else if (prog == "earth-fe") {

        Planet earth(0.000003, 1., 0., 0., 0., 2*M_PI, 0.);
        Planet sun(1., 0., 0., 0., 0., 0., 0.);

        vector<Planet> planets = {earth, sun};

        SolarSystem binary(planets, 3, fname, fname_E, fname_M);
        start = clock();
        binary.forwardEuler(mesh_points, final_time, planets.size()-1);
        end = clock();
        cout << setw(15) << setprecision(8) << "Time taken forward Euler Earth-Sun: " << 1000.*((end - start)/(double)CLOCKS_PER_SEC) << "ms" << endl;

    } else if (prog == "ellipse") {

        Planet earth(0.000003, 1., 0., 0., 0., 5.0, 0.);
        Planet sun(1., 0., 0., 0., 0., 0., 0.);

        vector<Planet> planets = {earth, sun};

        SolarSystem binary(planets, 3, fname, fname_E, fname_M);
        start = clock();
        binary.velocityVerlet(mesh_points, final_time, planets.size()-1);
        end = clock();
    }

    return 0;
}
