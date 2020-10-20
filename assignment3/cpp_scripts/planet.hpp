/**
* planet.hpp: heaper file for Planet
*
* Author: Are Frode Kvanum
*
* Completion Date: 19.10.2020
*/

#ifndef Planet_hpp
#define Planet_hpp
#include <cmath>

class Planet {
private:
    double m_M;
    double *m_position;
    double *m_velocity;

public:
    Planet(double M, double x, double y, double z, double vx, double vy, double vz);
    double distance(Planet othPlanet);
    double kinetic();
    double get_M();
    double get_v(int i);
    void set_v(int i, double value);
    double get_pos(int i);
    void set_pos(int i, double value);
};

#endif
