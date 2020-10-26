/**
* planet.cpp: implementation file for Planet
*
* Author: Are Frode Kvanum
*
* Completion date: 25.10.2020
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
    m_angular = new double[3];
    m_position[0] = x;
    m_position[1] = y;
    m_position[2] = z;
    m_velocity[0] = vx;
    m_velocity[1] = vy;
    m_velocity[2] = vz;
}

// PRIVATE MEMBER FUNCTIONS

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

/**
* kinetic: calculates the kinetic energy of the planet
*
* @return value is directly stored in private variable kinetic
*/
void Planet::kinetic() {
    double v2 = (m_velocity[0]*m_velocity[0] + m_velocity[1]*m_velocity[1] + m_velocity[2]*m_velocity[2]);
    m_kinetic = 0.5*m_M*v2;
}

/**
* potential: calculates the potential energy between two planets
*
* @param othPlanet the other planet in the calculation
* @param G gravitational constant
* @return value is directly stored in private variable potential
*/
void Planet::potential(Planet othPlanet, double G) {
    m_potential = (-G*m_M*othPlanet.get_M())/distance(othPlanet);
}

/**
* angular_momentum: calculates the angular momentum for the Planet
*                   Note: angular momentum is not scaled by Planet's mass
*
* @return values are directly stored in angular vector
*/
void Planet::angular_momentum() {
    m_angular[0] = (m_position[1]*m_velocity[2] - m_position[2]*m_velocity[1]);
    m_angular[1] = -(m_position[0]*m_velocity[2] - m_position[2]*m_velocity[0]);
    m_angular[2] = (m_position[0]*m_velocity[1] - m_position[1]*m_velocity[0]);
}

/**
* reset_energy: sets the value of both energy variables to 0
*
* @return value is stored directly in variable kinetic and potential
*/
void Planet::reset_energy() {
    m_kinetic = m_potential = 0.0;
}

/**
* get_v: getter function for private variable velocity
*
* @param i index determining dimensionsnality (0=x,1=y,2=z)
* @return velocity[i] where i is the index
*/
double Planet::get_v(int i) {
    return m_velocity[i];
}

/**
* set_v: setter function for private variable velocity
*
* @param i index determining dimensionsnality (0=x,1=y,2=z)
* @param value new value to be stored in velocity[i]
* @return value is stored dirctly in velocity[i]
*/
void Planet::set_v(int i, double value) {
    m_velocity[i] = value;
}

/**
* get_pos: getter function for private variable position
*
* @param i index determining dimensionsnality (0=x,1=y,2=z)
* @return position[i] where i is the index
*/
double Planet::get_pos(int i) {
    return m_position[i];
}

/**
* set_pos: setter function for private variable position
*
* @param i index determining dimensionsnality (0=x,1=y,2=z)
* @param value new value to be stored in positon[i]
* @return value is stored dirctly in position[i]
*/
void Planet::set_pos(int i, double value) {
    m_position[i] = value;
}

/**
* get_M: getter function for private variable M
*
* @return the mass of the planet
*/
double Planet::get_M() {
    return m_M;
}

/**
* get_kinetic: getter function for private variable
*
* @return the kinetic energy of the planet
*/
double Planet::get_kinetic() {
    return m_kinetic;
}

/**
* get_potential: getter function for private variable potential
*
* @return the potential energy of the planet
*/
double Planet::get_potential() {
    return m_potential;
}

double Planet::get_angular_momentum(int i) {
    return m_angular[i];
}
