/**
* solarsystem.hpp: header file for SolarSystem
*
* Author: Are Frode kvanum
*
* Completion Date: 19.10.2020
*/

#ifndef SolarSystem_hpp
#define SolarSystem_hpp
#include <iostream>
#include <vector>
#include <fstream>

#include "planet.hpp"

class SolarSystem {
private:
    std::vector<Planet> m_planets;
    int m_dimension;
    double m_G;

public:
    SolarSystem(std::vector<Planet> planets, int dimension);
    void forwardEuler(int mesh_points, double time);
    void velocityVerlet(int mesh_points, double time, int size);
    double **new_matrix(int i, int j);
    void gravitational_force(Planet &current, int idx, double &Fx, double &Fy, double &Fz);
    void print_pos(std::ofstream &output, double time, int number);
};

#endif
