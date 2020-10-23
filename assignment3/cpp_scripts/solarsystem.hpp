/**
* solarsystem.hpp: header file for SolarSystem
*
* Author: Are Frode kvanum
*
* Completion Date: 22.10.2020
*/

#ifndef SolarSystem_hpp
#define SolarSystem_hpp
#include <vector>
#include <fstream>
#include <string>

#include "planet.hpp"

class SolarSystem {

private:
    std::string m_fname, m_fname_E, m_fname_M;
    std::vector<Planet> m_planets;
    int m_dimension;
    double m_G, m_step;
    void advance(Planet &current, int idx, double **a, double **a_new, double &Fx, double &Fy, double &Fz, double &Fx_new, double &Fy_new, double &Fz_new);
    void gravitational_force(Planet &current, int idx, double &Fx, double &Fy, double &Fz);
    void kinetic_energy();
    void potential_energy();
    void print_pos(std::ofstream &output, double time, int number);
    void print_energy(std::ofstream &output, double time);
    void print_angular(std::ofstream &output, double time, int number);
    double **new_matrix(int i, int j);

public:
    SolarSystem(std::vector<Planet> planets, int dimension, std::string fname, std::string fname_E, std::string fname_M);
    void forwardEuler(int mesh_points, double time, int size);
    void velocityVerlet(int mesh_points, double time, int size);
};

#endif
