/**
* solarsystem.cpp: Implementation file for SolarSystem
*
* Author: Are Frode Kvanum
*
* Completion date: 19.10.2020
*/

#include "solarsystem.hpp"

using namespace std;

// CONSTRUCTOR

/**
* Default constructor
* Initializes a solar system with a list of planets
*
* @param planets list of planets to initialize the solar system
* @param dimension the dimensionality of the system
*/
SolarSystem::SolarSystem(vector<Planet> planets, int dimension) {
    m_planets = planets;
    m_dimension = dimension;
    m_G = 4.0*M_PI*M_PI;
}

void SolarSystem::velocityVerlet(int mesh_points, double time, int size) {
    double step = time / ((double) mesh_points);
    double current_time = 0.0;

    char *filename = new char[1000];
    sprintf(filename, "PlanetsVV_%lu_%.3f.txt", m_planets.size(), step);
    ofstream output_file(filename);

    double **a = new_matrix(m_planets.size(), m_dimension);
    double **a_new = new_matrix(m_planets.size(), m_dimension);

    double Fx, Fy, Fz, Fx_new, Fy_new, Fz_new;

    print_pos(output_file, current_time, 1);

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

        print_pos(output_file, current_time, 1);
        current_time += step;
    }

    output_file.close();

    for (int i = 0; i < m_planets.size(); i++) {
        delete[] a[i];
        delete[] a_new[i];
    }
    delete [] a;
    delete [] a_new;
}

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

void SolarSystem::forwardEuler(int mesh_points, double time) {
    double step = time /((double) mesh_points);
    double current_time = 0.0;

    double **a = new_matrix(m_planets.size(), m_dimension);
    double Fx, Fy, Fz;

    current_time += step;
    while (current_time < time) {
        for (int i = 0; i < m_planets.size(); i++) {
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
        current_time += step;
    }

    for (int i = 0; i < m_planets.size(); i++) {
        delete[] a[i];
    }
    delete [] a;
}

void SolarSystem::print_pos(ofstream &output, double time, int number) {
    for (int i = 0; i < number; i++) {
        Planet &current = m_planets[i];
        output << time << "\t" << i+1 << "\t" << current.get_M();
        for (int j = 0; j < m_dimension; j++)
            output << "\t" << current.get_pos(j);

        for (int j = 0; j < m_dimension; j++) {
            output << "\t" << current.get_v(j);
        }
        output << endl;
    }
}


