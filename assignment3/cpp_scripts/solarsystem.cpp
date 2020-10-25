/**
* solarsystem.cpp: Implementation file for SolarSystem
*
* Author: Are Frode Kvanum
*
* Completion date: 25.10.2020
*/
#include <ctime>

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
* @return the total runtime of all advance_fe instances is returned
*/
double SolarSystem::forwardEuler(int mesh_points, double time, int size) {
    double step = time /((double) mesh_points);
    double current_time = 0.0;
    clock_t start, end;
    double runtime;

    std::ofstream output_file(m_fname);
    std::ofstream output_file_E(m_fname_E);

    double **a = new_matrix(m_planets.size(), m_dimension);
    double Fx, Fy, Fz;
    Fx = Fy = Fz = 0.0;

    print_pos(output_file, current_time, size);
    print_energy(output_file_E, current_time);

    current_time += step;
    while (current_time < time) {
        for (int i = 0; i < size; i++) {
            Planet &current = m_planets[i];
            start = clock();
            advance_fe(current, i, step, a, Fx, Fy, Fz);
            end = clock();
            runtime += (double)(end - start)/((double)CLOCKS_PER_SEC);
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

    return runtime;
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
* @param beta: exponential factor for radius in force (default=3.0)
* @return the total runtime of all advance_vv instances is returned
*/
double SolarSystem::velocityVerlet(int mesh_points, double time, int size, double beta) {
    double step = time / ((double) mesh_points);
    double current_time = 0.0;
    clock_t start, end;
    double runtime;

    std::ofstream output_file(m_fname);
    std::ofstream output_file_E(m_fname_E);
    std::ofstream output_file_M(m_fname_M);

    double **a = new_matrix(m_planets.size(), m_dimension);
    double **a_new = new_matrix(m_planets.size(), m_dimension);

    double Fx, Fy, Fz, Fx_new, Fy_new, Fz_new;
    Fx = Fy = Fz = Fx_new = Fy_new = Fz_new = 0.0;

    print_pos(output_file, current_time, size);
    print_energy(output_file_E, current_time);
    print_angular(output_file_M, current_time, size);

    current_time += step;
    while (current_time < time) {
        for (int i = 0; i < size; i++) {
            Planet &current = m_planets[i];
            start = clock();
            advance_vv(current, i, step, a, a_new, Fx, Fy, Fz, Fx_new, Fy_new, Fz_new, beta);
            end = clock();
            runtime += (double)(end - start)/((double)CLOCKS_PER_SEC);
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

    return runtime;
}

/**
* test_velocityVerlet: method for unit testing total energy
*                      and conservation of angular momentum
*
*/
void SolarSystem::test_velocityVerlet(int mesh_points, double time, int size) {
    double step = time / ((double) mesh_points);
    double current_time = 0.0;
    double energy[2];
    double momentum[2];
    double tol = 1e-8;

    double **a = new_matrix(m_planets.size(), m_dimension);
    double **a_new = new_matrix(m_planets.size(), m_dimension);

    double Fx, Fy, Fz, Fx_new, Fy_new, Fz_new;
    Fx = Fy = Fz = Fx_new = Fy_new = Fz_new = 0.0;

    kinetic_energy();
    potential_energy();

    energy[0] = m_planets[0].get_kinetic() + m_planets[0].get_potential();

    m_planets[0].angular_momentum();
    momentum[0] = m_planets[0].get_angular_momentum(2);

    current_time += step;
    while (current_time < time) {
        for (int i = 0; i < size; i++) {
            Planet &current = m_planets[i];
            advance_vv(current, i, step, a, a_new, Fx, Fy, Fz, Fx_new, Fy_new, Fz_new);
        }

        kinetic_energy();
        potential_energy();
        energy[1] = m_planets[0].get_kinetic() + m_planets[0].get_potential();
        if (fabs(energy[0] - energy[1]) > tol) {
            std::cout << "Energy at time " << current_time << ": " << energy[1] << " not equal to energy at time " << current_time - step <<": " << energy[0] << std::endl;
            exit(1);
        }

        m_planets[0].angular_momentum();
        momentum[1] = m_planets[0].get_angular_momentum(2);
        if (fabs(momentum[0] - momentum[1]) > tol) {
            std::cout << "Angular momentum at time " << current_time << ": " << momentum[1] << " not equal to angular momentum at time " << current_time - step <<": " << momentum[0] << std::endl;
            exit(1);
        }
        energy[0] = energy[1];
        momentum[0] = momentum[1];

        current_time += step;
    }

    for (int i = 0; i < m_planets.size(); i++) {
        delete[] a[i];
        delete[] a_new[i];
    }
    delete [] a;
    delete [] a_new;
}

// PRIVATE MEMBER FUNCTIONS

/**
* advance_fe: advances the solar system by one timestep using forward Euler
*
* @param current: current planet to advance
* @param idx: index of current in planets vector
* @param step: the step length
* @param a: acceleration matrix
* @param Fx: force in x-direction
* @param Fy: force in y-direction
* @param Fz: force in z-direction
* @return nothing is returned
*/
void SolarSystem::advance_fe(Planet &current, int idx, double step, double **a, double &Fx, double &Fy, double &Fz) {
    Fx = Fy = Fz = 0.0;

    gravitational_force(current, idx, Fx, Fy, Fz);

    a[idx][0] = Fx / current.get_M();
    a[idx][1] = Fy / current.get_M();
    a[idx][2] = Fz / current.get_M();

    for (int j = 0; j < m_dimension; j++) {
        double new_pos = current.get_pos(j) + step*current.get_v(j);
        current.set_pos(j, new_pos);
        double new_vel = current.get_v(j) + step*a[idx][j];
        current.set_v(j, new_vel);
    }
}

/**
* advance_vv: advances the solar system by one timestep using velocity Verlet
*
* @param current: current planet to advance
* @param idx: index of current in planets vector
* @param step: the step length
* @param a: acceleration matrix
* @param a_new: acceleration matrix next timestep
* @param Fx: force in x-direction
* @param Fy: force in y-direction
* @param Fz: force in z-direction
* @param Fx_new: force in x-direction next-timestep
* @param Fy_new: force in y-direction next-timestep
* @param Fz_new: force in z-direction next-timestep
* @param beta: exponential factor for radius in force (default=3.0)
* @return nothing is returned
*/
void SolarSystem::advance_vv(Planet &current, int idx, double step, double **a, double **a_new, double &Fx, double &Fy, double &Fz, double &Fx_new, double &Fy_new, double &Fz_new, double beta) {
    Fx = Fy = Fz = Fx_new = Fy_new = Fz_new = 0.0;

    gravitational_force(current, idx, Fx, Fy, Fz, beta);

    a[idx][0] = Fx / current.get_M();
    a[idx][1] = Fy / current.get_M();
    a[idx][2] = Fz / current.get_M();

    for (int j = 0; j < m_dimension; j++) {
        double new_pos = current.get_pos(j) + current.get_v(j)*step + 0.5*step*step*a[idx][j];
        current.set_pos(j, new_pos);
    }

    gravitational_force(current, idx, Fx_new, Fy_new, Fz_new, beta);
    a_new[idx][0] = Fx_new / current.get_M();
    a_new[idx][1] = Fy_new / current.get_M();
    a_new[idx][2] = Fz_new / current.get_M();

    for (int j = 0; j < m_dimension; j++) {
        double new_vel = current.get_v(j) + 0.5*step*(a_new[idx][j] + a[idx][j]);
        current.set_v(j, new_vel);
    }
}

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
void SolarSystem::gravitational_force(Planet &current, int idx, double &Fx, double &Fy, double &Fz, double beta) {
    for (int j = 0; j < m_planets.size(); j++) {
        if (j != idx) {
            Planet &other = m_planets[j];
            double rel_dist[3];
            for (int k = 0; k < m_dimension; k++) {
                rel_dist[k] = current.get_pos(k) - other.get_pos(k);
            }
            double r = current.distance(other);
            Fx += other.get_M()*rel_dist[0]/(pow(r,beta));
            Fy += other.get_M()*rel_dist[1]/(pow(r,beta));
            Fz += other.get_M()*rel_dist[2]/(pow(r,beta));
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
