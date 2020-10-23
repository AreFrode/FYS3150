/**
* solarsystem.cpp: Implementation file for SolarSystem
*
* Author: Are Frode Kvanum
*
* Completion date: 22.10.2020
*/
#include <iostream>
#include "solarsystem.hpp"

// CONSTRUCTOR

/**
* Default constructor
* Initializes a solar system with a list of planets
*
* @param planets list of planets to initialize the solar system
* @param dimension the dimensionality of the system
* @param fname name of the file that store positional output
* @param fname_E name of the file that store energy output
* @param fname_M name of the file that store angular momentum output
*/
SolarSystem::SolarSystem(std::vector<Planet> planets, int dimension, std::string fname, std::string fname_E, std::string fname_M) {
    m_planets = planets;
    m_dimension = dimension;
    m_fname = fname;
    m_fname_E = fname_E;
    m_fname_M = fname_M;
    m_G = 4.0*M_PI*M_PI;
}

// PUBLIC MEMBER FUNCTIONS

/**
* forwardEuler: advances the solar system by a given time using Euler's forward algorithm
*
* @param mesh_points: number of points to perform calculations
* @param time: the final time of the advancement
* @param size: the number of planets to be considered
*              (clarification: if less than size of private
*               variable planets, Planets at index > size
*               will not be advanced in the simulation)
* @return nothing is returned
*/
void SolarSystem::forwardEuler(int mesh_points, double time, int size) {
    double step = time /((double) mesh_points);
    double current_time = 0.0;

    std::ofstream output_file(m_fname);
    std::ofstream output_file_E(m_fname_E);

    double **a = new_matrix(m_planets.size(), m_dimension);
    double Fx, Fy, Fz;

    print_pos(output_file, current_time, size);
    print_energy(output_file_E, current_time);

    current_time += step;
    while (current_time < time) {
        for (int i = 0; i < size; i++) {
            Planet &current = m_planets[i];
            Fx = Fy = Fz = 0.0;

            gravitational_force(current, i, Fx, Fy, Fz);

            a[i][0] = Fx / current.get_M();
            a[i][1] = Fy / current.get_M();
            a[i][2] = Fz / current.get_M();

            for (int j = 0; j < m_dimension; j++) {
                double new_pos = current.get_pos(j) + step*current.get_v(j);
                current.set_pos(j, new_pos);
                double new_vel = current.get_v(j) + step*a[i][j];
                current.set_v(j, new_vel);
            }
        }

        print_pos(output_file, current_time, size);
        print_energy(output_file_E, current_time);
        current_time += step;
    }

    output_file.close();
    output_file_E.close();

    for (int i = 0; i < m_planets.size(); i++) {
        delete[] a[i];
    }
    delete [] a;
}

/**
* velocityVerlet: advances the solar system by a given time using the velocity Verlet algorithm
*
* @param mesh_points: number of points to perform calculations
* @param time: the final time of the advancement
* @param size: the number of planets to be considered
*              (clarification: if less than size of private
*               variable planets, Planets at index > size
*               will not be advanced in the simulation)
* @return nothing is returned
*/
void SolarSystem::velocityVerlet(int mesh_points, double time, int size) {
    double step = time / ((double) mesh_points);
    double current_time = 0.0;

    std::ofstream output_file(m_fname);
    std::ofstream output_file_E(m_fname_E);
    std::ofstream output_file_M(m_fname_M);

    double **a = new_matrix(m_planets.size(), m_dimension);
    double **a_new = new_matrix(m_planets.size(), m_dimension);

    double Fx, Fy, Fz, Fx_new, Fy_new, Fz_new;

    print_pos(output_file, current_time, size);
    print_energy(output_file_E, current_time);
    print_angular(output_file_M, current_time, size);

    current_time += step;
    while (current_time < time) {
        for (int i = 0; i < size; i++) {
            Planet &current = m_planets[i];

            Fx = Fy = Fz = Fx_new = Fy_new = Fz_new = 0.0;

            gravitational_force(current, i, Fx, Fy, Fz);

            a[i][0] = Fx / current.get_M();
            a[i][1] = Fy / current.get_M();
            a[i][2] = Fz / current.get_M();

            for (int j = 0; j < m_dimension; j++) {
                double new_pos = current.get_pos(j) + current.get_v(j)*step + 0.5*step*step*a[i][j];
                current.set_pos(j, new_pos);
            }

            gravitational_force(current, i, Fx_new, Fy_new, Fz_new);
            a_new[i][0] = Fx_new / current.get_M();
            a_new[i][1] = Fy_new / current.get_M();
            a_new[i][2] = Fz_new / current.get_M();

            for (int j = 0; j < m_dimension; j++) {
                double new_vel = current.get_v(j) + 0.5*step*(a_new[i][j] + a[i][j]);
                current.set_v(j, new_vel);
            }
        }

        print_pos(output_file, current_time, size);
        print_energy(output_file_E, current_time);
        print_angular(output_file_M, current_time, size);

        current_time += step;
    }

    output_file.close();
    output_file_E.close();

    for (int i = 0; i < m_planets.size(); i++) {
        delete[] a[i];
        delete[] a_new[i];
    }
    delete [] a;
    delete [] a_new;
}


// PRIVATE MEMBER FUNCTIONS

/**
* gravitational_force: calculates the gravitational force between one planet and the others
*
* @param current planet to calculate force for
* @param idx index for current Planet in vector planets
* @param Fx force in x-direction
* @param Fy force in y-direction
* @param Fz force in z-direction
* @return values are stores appropriately in Fx, Fy and Fz
*/
void SolarSystem::gravitational_force(Planet &current, int idx, double &Fx, double &Fy, double &Fz) {
    for (int j = 0; j < m_planets.size(); j++) {
        if (j != idx) {
            Planet &other = m_planets[j];
            double rel_dist[3];
            for (int k = 0; k < m_dimension; k++) {
                rel_dist[k] = current.get_pos(k) - other.get_pos(k);
            }
            double r = current.distance(other);
            Fx += other.get_M()*rel_dist[0]/(r*r*r);
            Fy += other.get_M()*rel_dist[1]/(r*r*r);
            Fz += other.get_M()*rel_dist[2]/(r*r*r);
        }
    }
    Fx *= -m_G*current.get_M();
    Fy *= -m_G*current.get_M();
    Fz *= -m_G*current.get_M();
}

/**
* kinetic_energy: calculates the kinetic energy for all Planets
*
* @return nothing is returned
*/
void SolarSystem::kinetic_energy() {
    for (int i = 0; i < m_planets.size(); i++) {
        Planet &current = m_planets[i];
        current.kinetic();
    }
}

/**
* potential_energy: calculates the potential energy for all Planets
*
* @returns nothing is returned
*/
void SolarSystem::potential_energy() {
    for (int i = 0; i < m_planets.size(); i++) {
        Planet &current = m_planets[i];
        for (int j = 0; j < m_planets.size(); j++) {
            if (j != i) {
                Planet &other = m_planets[j];
                current.potential(other, m_G);
            }
        }
    }
}

/**
* print_pos: writes positional data for current time to file
*
* @param output: outputstream for the outfile
* @param time: current timestep
* @param number: amount of planets being considered
* @return nothing is returned
*/
void SolarSystem::print_pos(std::ofstream &output, double time, int number) {
    for (int i = 0; i < number; i++) {
        Planet &current = m_planets[i];
        output << time << "\t" << i+1 << "\t" << current.get_M();
        for (int j = 0; j < m_dimension; j++)
            output << "\t" << current.get_pos(j);

        for (int j = 0; j < m_dimension; j++) {
            output << "\t" << current.get_v(j);
        }
        output << std::endl;
    }
}

/**
* print_energy: writes energy data for current time to file
*
* @param output: outputstream for the outfile
* @param time: current timestep
* @return nothing is returned
*/
void SolarSystem::print_energy(std::ofstream &output, double time) {
    kinetic_energy();
    potential_energy();
    int counter = 0;
    for (Planet &planet : m_planets) {
        output << time << "\t" << counter << "\t";
        output << planet.get_kinetic() << "\t" << planet.get_potential() << "\t";
        output << planet.get_kinetic() + planet.get_potential() << std::endl;
        counter++;
    }
}

/**
*
*/
void SolarSystem::print_angular(std::ofstream &output, double time, int number) {
    int counter = 0;
    for (Planet &planet : m_planets) {
        planet.angular_momentum();
        output << time << "\t" << counter << "\t";
        for (int i = 0; i < m_dimension; i++) {
            output << planet.get_angular_momentum(i) << "\t";
        }
        output << std::endl;
        counter++;
    }
}

/**
* new_matrix: constructs a matrix to store acceleration values
*
* @param n height of the matrix
* @param m width of the matrix
* @return calloced matrix
*/
double **SolarSystem::new_matrix(int n, int m) {
    double **mat;
    mat = new double*[n];

    for (int i = 0; i < n; i++)
        mat[i] = new double[m];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            mat[i][j] = 0.0;
        }
    }

    return mat;
}
