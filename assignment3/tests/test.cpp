#include "../cpp_scripts/planet.hpp"
#include "../cpp_scripts/solarsystem.hpp"

void test_energy();
void test_momentum();

int main(int argc, char *argv[]) {


    return 0;
}

void test_energy() {
        Planet earth(0.000003, 1., 0., 0., 0.0, 2*M_PI, 0.);
        Planet sun(1., 0., 0., 0., 0., 0., 0.);

        std::vector<Planet> planets = {earth, sun};

        double current_time = 0.0;
        double time = 50.;
        double mesh = 10000;
        double step = time /(double)mesh;
        double Fx, Fy, Fz, Fx_new, Fy_new, Fz_new;

        current_time += step;
        while (current_time < time) {
            for (int i = 0; i < 1; i++) {
                Planet &current = planets[i];

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

}
