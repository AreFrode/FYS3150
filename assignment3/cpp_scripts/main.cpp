#include <iostream>
#include <iomanip>
#include <cstdio>

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
    double runtime;

    if (prog == "earth-vv") {

        Planet earth(0.000003, 1., 0., 0., 0.0, 2*M_PI, 0.);
        Planet sun(1., 0., 0., 0., 0., 0., 0.);

        vector<Planet> planets = {earth, sun};

        SolarSystem binary(planets, 3, fname, fname_E, fname_M);

        runtime = binary.velocityVerlet(mesh_points, final_time, planets.size()-1);

        cout << setw(15) << setprecision(8) << "Total runtime for " << mesh_points << " iterations with velocity Verlet Earth-Sun: " << 1000.*runtime << "ms" << endl;

    } else if (prog == "earth-fe") {

        Planet earth(0.000003, 1., 0., 0., 0., 2*M_PI, 0.);
        Planet sun(1., 0., 0., 0., 0., 0., 0.);

        vector<Planet> planets = {earth, sun};

        SolarSystem binary(planets, 3, fname, fname_E, fname_M);

        runtime = binary.forwardEuler(mesh_points, final_time, planets.size()-1);

        cout << setw(15) << setprecision(8) << "Total runtime for " << mesh_points << " iterations with forward Euler Earth-Sun: " << 1000.*runtime << "ms" << endl;

    } else if (prog == "ellipse") {

        Planet earth(0.000003, 1., 0., 0., 0., 5.0, 0.);
        Planet sun(1., 0., 0., 0., 0., 0., 0.);

        vector<Planet> planets = {earth, sun};

        SolarSystem binary(planets, 3, fname, fname_E, fname_M);

        binary.velocityVerlet(mesh_points, final_time, planets.size()-1);

    } else if (prog == "beta-3.5-circ") {

        Planet earth(0.000003, 1., 0., 0., 0., 2.*M_PI, 0.);
        Planet sun(1., 0., 0., 0., 0., 0., 0.);

        vector<Planet> planets = {earth, sun};

        SolarSystem binary(planets, 3, fname, fname_E, fname_M);

        binary.velocityVerlet(mesh_points, final_time, planets.size()-1, 3.5);

    } else if (prog == "beta-3.9-circ") {

        Planet earth(0.000003, 1., 0., 0., 0., 2.*M_PI, 0.);
        Planet sun(1., 0., 0., 0., 0., 0., 0.);

        vector<Planet> planets = {earth, sun};

        SolarSystem binary(planets, 3, fname, fname_E, fname_M);

        binary.velocityVerlet(mesh_points, final_time, planets.size()-1, 3.9);

    } else if (prog == "beta-4.0-circ") {

        Planet earth(0.000003, 1., 0., 0., 0., 2.*M_PI, 0.);
        Planet sun(1., 0., 0., 0., 0., 0., 0.);

        vector<Planet> planets = {earth, sun};

        SolarSystem binary(planets, 3, fname, fname_E, fname_M);

        binary.velocityVerlet(mesh_points, final_time, planets.size()-1, 4.0);

    } else if (prog == "beta-3.5-ellipse") {

        Planet earth(0.000003, 1., 0., 0., 0., 5., 0.);
        Planet sun(1., 0., 0., 0., 0., 0., 0.);

        vector<Planet> planets = {earth, sun};

        SolarSystem binary(planets, 3, fname, fname_E, fname_M);

        binary.velocityVerlet(mesh_points, final_time, planets.size()-1, 3.5);

    } else if (prog == "beta-3.9-ellipse") {

        Planet earth(0.000003, 1., 0., 0., 0., 5., 0.);
        Planet sun(1., 0., 0., 0., 0., 0., 0.);

        vector<Planet> planets = {earth, sun};

        SolarSystem binary(planets, 3, fname, fname_E, fname_M);

        binary.velocityVerlet(mesh_points, final_time, planets.size()-1, 3.9);

    } else if (prog == "beta-4.0-ellipse") {

        Planet earth(0.000003, 1., 0., 0., 0., 5., 0.);
        Planet sun(1., 0., 0., 0., 0., 0., 0.);

        vector<Planet> planets = {earth, sun};

        SolarSystem binary(planets, 3, fname, fname_E, fname_M);

        binary.velocityVerlet(mesh_points, final_time, planets.size()-1, 4.0);

    } else if (prog == "escape") {

        Planet earth(0.000003, 1., 0., 0., 0., sqrt(8.0*M_PI*M_PI), 0.);
        Planet sun(1., 0., 0., 0., 0., 0., 0.);

        vector<Planet> planets = {earth, sun};

        SolarSystem binary(planets, 3, fname, fname_E, fname_M);

        binary.velocityVerlet(mesh_points, final_time, planets.size()-1);

    } else if (prog == "earth-jupiter") {

        Planet earth(0.000003, 1., 0., 0., 0., 2*M_PI, 0.);
        Planet jupiter(0.00095, 5.20, 0., 0., 0., M_PI, 0.);
        Planet sun(1.0, 0., 0., 0., 0., 0., 0.);

        vector<Planet> planets = {earth, jupiter, sun};

        SolarSystem three(planets, 3, fname, fname_E, fname_M);

        three.velocityVerlet(mesh_points, final_time, planets.size()-1);

    } else if (prog == "earth-10jupiter") {

        Planet earth(0.000003, 1., 0., 0., 0., 2*M_PI, 0.);
        Planet jupiter(0.0095, 5.20, 0., 0., 0., M_PI, 0.);
        Planet sun(1.0, 0., 0., 0., 0., 0., 0.);

        vector<Planet> planets = {earth, jupiter, sun};

        SolarSystem three(planets, 3, fname, fname_E, fname_M);

        three.velocityVerlet(mesh_points, final_time, planets.size()-1);

    } else if (prog == "earth-1000jupiter") {


        Planet earth(0.000003, 1., 0., 0., 0., 2*M_PI, 0.);
        Planet jupiter(0.95, 5.20, 0., 0., 0., M_PI, 0.);
        Planet sun(1.0, 0., 0., 0., 0., 0., 0.);

        vector<Planet> planets = {earth, jupiter, sun};

        SolarSystem three(planets, 3, fname, fname_E, fname_M);

        three.velocityVerlet(mesh_points, final_time, planets.size()-1);

    } else if (prog == "three-body") {

        Planet earth(0.000003, 1. - 0.004943, 0., 0., 0., 2*M_PI, 0.);
        Planet jupiter(0.95, 5.20 - 0.004943, 0., 0., 0., M_PI, 0.);
        Planet sun(1.0, -0.004943, 0., 0., 0., -0.25*M_PI, 0.);

        vector<Planet> planets = {earth, jupiter, sun};

        SolarSystem three(planets, 3, fname, fname_E, fname_M);

        three.velocityVerlet(mesh_points, final_time, planets.size());

    } else if (prog == "full-system") {

        double *x, *y, *z, *vx, *vy, *vz;
        double *mass;
        int N = 10;

        x = new double[10];
        y = new double[10];
        z = new double[10];
        vx = new double[10];
        vy = new double[10];
        vz = new double[10];
        mass = new double[10];

        const char* filename_pos_and_vel = "positions_and_vel.txt";
        const char* filename_mass = "masses.txt";

        FILE *fp_init = fopen(filename_pos_and_vel, "r");
        FILE *fp_mass = fopen(filename_mass, "r");

        int foo; // g++ complains about ignoring return of fscanf

        for (int i = 0; i < N; i++) {
            foo = fscanf(fp_init, "%lf %lf %lf %lf %lf %lf", &x[i], &y[i], &z[i], &vx[i], &vy[i], &vz[i]);
            foo = fscanf(fp_mass, "%lf", &mass[i]);
        }

        fclose(fp_init);
        fclose(fp_mass);

        double com_posx, com_posy, com_posz;
        double com_velx, com_vely, com_velz;
        double sum_mass;
        com_posx = com_posy = com_posz = com_velx = com_vely = com_velz = sum_mass = 0.0;

        for (int i = 0; i < N; i++) {
            double cmass = mass[i];
            com_posx += (cmass*x[i]);
            com_posy += (cmass*y[i]);
            com_posz += (cmass*z[i]);
            com_velx += (cmass*vx[i]);
            com_vely += (cmass*vy[i]);
            com_velz += (cmass*vz[i]);
            sum_mass += cmass;
        }

        com_posx /= sum_mass;
        com_posy /= sum_mass;
        com_posz /= sum_mass;
        com_velx /= sum_mass;
        com_vely /= sum_mass;
        com_velz /= sum_mass;

        Planet sun(mass[0]/mass[0], x[0] - com_posx, y[0] - com_posy, z[0] - com_posz, vx[0] - com_velx, vy[0] - com_vely, vz[0] - com_velz);
        Planet mercury(mass[1]/mass[0], x[1] - com_posx, y[1] - com_posy, z[1] - com_posz, vx[1] - com_velx, vy[1] - com_vely, vz[1] - com_velz);
        Planet venus(mass[2]/mass[0], x[2] - com_posx, y[2] - com_posy, z[2] - com_posz, vx[2] - com_velx, vy[2] - com_vely, vz[2] - com_velz);
        Planet earth(mass[3]/mass[0], x[3] - com_posx, y[3] - com_posy, z[3] - com_posz, vx[3] - com_velx, vy[3] - com_vely, vz[3] - com_velz);
        Planet mars(mass[4]/mass[0], x[4] - com_posx, y[4] - com_posy, z[4] - com_posz, vx[4] - com_velx, vy[4] - com_vely, vz[4] - com_velz);
        Planet jupiter(mass[5]/mass[0], x[5] - com_posx, y[5] - com_posy, z[5] - com_posz, vx[5] - com_velx, vy[5] - com_vely, vz[5] - com_velz);
        Planet saturn(mass[6]/mass[0], x[6] - com_posx, y[6] - com_posy, z[6] - com_posz, vx[6] - com_velx, vy[6] - com_vely, vz[6] - com_velz);
        Planet uranus(mass[7]/mass[0], x[7] - com_posx, y[7] - com_posy, z[7] - com_posz, vx[7] - com_velx, vy[7] - com_vely, vz[7] - com_velz);
        Planet neptun(mass[8]/mass[0], x[8] - com_posx, y[8] - com_posy, z[8] - com_posz, vx[8] - com_velx, vy[8] - com_vely, vz[8] - com_velz);
        Planet pluto(mass[9]/mass[0], x[9] - com_posx, y[9] - com_posy, z[9] - com_posz, vx[9] - com_velx, vy[9] - com_vely, vz[9] - com_velz);

        vector<Planet> planets = {sun, mercury, venus, earth, mars, jupiter, saturn, uranus, neptun, pluto};

        SolarSystem mw(planets, 3, fname, fname_E, fname_M);

        mw.velocityVerlet(mesh_points, final_time, planets.size());

    }

    return 0;
}
