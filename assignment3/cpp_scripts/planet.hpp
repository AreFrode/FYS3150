/**
* planet.hpp: heaper file for Planet
*
* Author: Are Frode Kvanum
*
* Completion Date: 22.10.2020
*/

#ifndef Planet_hpp
#define Planet_hpp
#include <cmath>

class Planet {
    friend class SolarSystem;

private:
    double m_M, m_kinetic, m_potential;
    double *m_angular;
    double *m_position;
    double *m_velocity;
    double distance(Planet othPlanet);
    void kinetic();
    void potential(Planet othPlanet, double G);
    void angular_momentum();
    void reset_energy();
    double get_M();
    double get_v(int i);
    void set_v(int i, double value);
    double get_pos(int i);
    void set_pos(int i, double value);
    double get_kinetic();
    double get_potential();
    double get_angular_momentum(int i);

public:
    Planet(double M, double x, double y, double z, double vx, double vy, double vz);
};

#endif
