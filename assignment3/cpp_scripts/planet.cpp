/**
* planet.cpp: implementation file for Planet
*
* Author: Are Frode Kvanum
*
* Completion date: 19.10.2020
*/

#include "planet.hpp"

// CONSTRUCTOR

/**
* Default constructor
* Sets the mass, as well as initial position
* and velocity for the Planet
*
* @param M mass of the planet as scaled by the Sun's mass
* @param x inital x-position for the planet
* @param y initial y-position for the planet
* @param z inital z-position for the planet
* @param vx inital x-velocity for the planet
* @param vy initial y-velocity for the planet
* @param vz initial z-velocity for the planet
*/
Planet::Planet(double M, double x, double y, double z, double vx, double vy, double vz) {
    m_M = M;
    m_position = new double[3];
    m_velocity = new double[3];
    m_position[0] = x;
    m_position[1] = y;
    m_position[2] = z;
    m_velocity[0] = vx;
    m_velocity[1] = vy;
    m_velocity[2] = vz;
}

// PUBLIC MEMBER FUNCTIONS

/**
* distance: calculates the distance between two Planets
*
* @param othPlanet the other Planet
* @return The distance between the two planets
*/
double Planet::distance(Planet othPlanet) {
    double x, y, z;
    x = m_position[0] - othPlanet.get_pos(0);
    y = m_position[1] - othPlanet.get_pos(1);
    z = m_position[2] - othPlanet.get_pos(2);

    return sqrt(x*x + y*y + z*z);
}

double Planet::kinetic() {
    double v2 = (m_velocity[0]*m_velocity[0] + m_velocity[1]*m_velocity[1] + m_velocity[2]*m_velocity[2]);
    return 0.5*m_M*v2;
}

double Planet::get_v(int i) {
    return m_velocity[i];
}

void Planet::set_v(int i, double value) {
    m_velocity[i] = value;
}

double Planet::get_pos(int i) {
    return m_position[i];
}

void Planet::set_pos(int i, double value) {
    m_position[i] = value;
}

double Planet::get_M() {
    return m_M;
}

