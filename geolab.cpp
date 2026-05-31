#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

using Grid = std::vector<std::vector<double>>;

struct Params {
    // Corrected geolab dimensions: 225 cm x 635 cm
    double Lx = 2.25;
    double Ly = 6.35;

    double dx = 0.10;
    int Nx;
    int Ny;

    double dt = 10.0;
    double tmax = 1800.0;  // 30 minutes, enough for 1, 5, 10, 20, 30 min maps
    int nmax;

    // Effective parameters recalculated after geometry correction.
    // D is an effective diffusivity: it includes unresolved air mixing/convection.
    double D = 0.0013013117895950114;
    double h = 0.08103161457851407;
    double h_w = 0.08103161457851407;

    // Temperatures in Kelvin
    double T_inf = 276.15;   // outside/environment temperature assumption
    double T_init = 290.15;  // initial indoor temperature

    // Heater source strength. This stays approximately unchanged under the same source-area model.
    double Smax = 5.766483445358526;

    int kmax = 100;
    double TOL = 1e-8;

    // Sensor location used for optional point comparison.
    // Keep/check this against the real geolab sensor position if available.
    int ic;
    int jc;

    Params() {
        Nx = static_cast<int>(std::round(Lx / dx)) + 1;
        Ny = static_cast<int>(std::round(Ly / dx)) + 1;
        nmax = static_cast<int>(std::round(tmax / dt));

        double x_sensor = 1.0;
        double y_sensor = Ly;

        ic = static_cast<int>(std::round(x_sensor / dx));
        jc = static_cast<int>(std::round(y_sensor / dx));

        ic = std::clamp(ic, 0, Nx - 1);
        jc = std::clamp(jc, 0, Ny - 1);
    }
};

double average_temperature(const Grid& T) {
    double sum = 0.0;
    int count = 0;

    for (const auto& row : T) {
        for (double value : row) {
            sum += value;
            ++count;
        }
    }

    return sum / static_cast<double>(count);
}

void apply_all_boundaries(const Params& p, Grid& T) {
    int Nx = p.Nx;
    int Ny = p.Ny;

    double D = p.D;
    double dx = p.dx;
    double Tinf = p.T_inf;
    double Troom = p.T_init;

    // Convection wall coefficients from: -D dT/dn = h(T - T_inf)
    double A = (p.h * dx) / (D + p.h * dx);
    double B = D / (D + p.h * dx);

    // ==================== INTERNAL WALLS ====================
// WEST (x=0) stays at room temperature
    for (int j = 1; j < Ny - 1; ++j) {
        T[0][j] = Troom;
    }

// EAST (x=Lx) stays at room temperature
    for (int j = 1; j < Ny - 1; ++j) {
        T[Nx - 1][j] = Troom;
    }

// ==================== EXTERNAL WALLS ====================
// SOUTH wall: convective heat loss
    for (int i = 1; i < Nx - 1; ++i) {
        T[i][0] = A * Tinf + B * T[i][1];
    }

// NORTH wall: convective heat loss
    for (int i = 1; i < Nx - 1; ++i) {
        T[i][Ny - 1] = A * Tinf + B * T[i][Ny - 2];
    }

    // ==================== CORNERS ============================
    T[0][0] = Troom;
    T[0][Ny - 1] = Troom;
    T[Nx - 1][Ny - 1] = Troom;
    T[Nx - 1][0] = Troom;
}

void compute_R(const Params& p,
               const Grid& Tn,
               const Grid& S,
               Grid& R,
               double w_old) {
    int Nx = p.Nx;
    int Ny = p.Ny;

    double D = p.D;
    double dt = p.dt;
    double dx = p.dx;

    double alpha = D * dt / (2.0 * dx * dx);
    double beta = dt / 2.0;

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            R[i][j] = 0.0;
        }
    }

    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            double laplace_T =
                    Tn[i + 1][j] + Tn[i - 1][j]
                    + Tn[i][j + 1] + Tn[i][j - 1]
                    - 4.0 * Tn[i][j];

            R[i][j] = Tn[i][j]
                      + alpha * laplace_T
                      + beta * w_old * S[i][j];
        }
    }
}

void gauss_seidel_CN_step(const Params& p,
                          const Grid& R,
                          const Grid& S,
                          Grid& T,
                          double w) {
    int Nx = p.Nx;
    int Ny = p.Ny;

    double D = p.D;
    double dt = p.dt;
    double dx = p.dx;

    double alpha = D * dt / (2.0 * dx * dx);
    double beta = dt / 2.0;
    double denom = 1.0 + 4.0 * alpha;

    for (int k = 0; k < p.kmax; ++k) {
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                double Tnew =
                        (R[i][j]
                         + alpha * (T[i - 1][j] + T[i + 1][j]
                                    + T[i][j - 1] + T[i][j + 1])
                         + beta * w * S[i][j])
                        / denom;

                T[i][j] = Tnew;
            }
        }

        apply_all_boundaries(p, T);

        double residual = 0.0;

        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                double Lap =
                        T[i + 1][j] + T[i - 1][j]
                        + T[i][j + 1] + T[i][j - 1]
                        - 4.0 * T[i][j];

                double Lij = T[i][j]
                             - alpha * Lap
                             - beta * w * S[i][j];

                double diff = Lij - R[i][j];
                residual += diff * diff * dx * dx;
            }
        }

        if (std::sqrt(residual) < p.TOL) {
            break;
        }
    }
}

void init_temperature(const Params& p, Grid& T) {
    for (int i = 0; i < p.Nx; ++i) {
        for (int j = 0; j < p.Ny; ++j) {
            T[i][j] = p.T_init;
        }
    }
}

void init_source(const Params& p, Grid& S) {
    double dx = p.dx;

    // Heater physical location in meters.
    // Update these coordinates if you know the exact heater position in the corrected geolab plan.
    double hx = 1.0;
    double hy = 4.8;

    double heater_radius = 0.20;
    double r2 = heater_radius * heater_radius;

    for (int i = 0; i < p.Nx; ++i) {
        double x = i * dx;

        for (int j = 0; j < p.Ny; ++j) {
            double y = j * dx;

            double dx2 = (x - hx) * (x - hx);
            double dy2 = (y - hy) * (y - hy);

            if (dx2 + dy2 <= r2) {
                S[i][j] = p.Smax;
            } else {
                S[i][j] = 0.0;
            }
        }
    }
}

void save_temperature_field(const std::string& dir,
                            int n,
                            const Grid& T,
                            const Params& p) {
    std::ostringstream fname;
    fname << dir << "/T_" << std::setw(4) << std::setfill('0') << n << ".dat";

    std::ofstream file(fname.str());

    for (int i = 0; i < p.Nx; ++i) {
        for (int j = 0; j < p.Ny; ++j) {
            file << T[i][j];
            if (j < p.Ny - 1) {
                file << " ";
            }
        }
        file << "\n";
    }
}

void CN_algorithm(const std::string& dir,
                  const Params& p,
                  const std::vector<double>& save_times) {
    Grid T(p.Nx, std::vector<double>(p.Ny));
    Grid Tn(p.Nx, std::vector<double>(p.Ny));
    Grid S(p.Nx, std::vector<double>(p.Ny));
    Grid R(p.Nx, std::vector<double>(p.Ny));

    init_temperature(p, T);
    init_source(p, S);
    apply_all_boundaries(p, T);

    std::vector<int> save_steps;
    save_steps.reserve(save_times.size());
    for (double save_time : save_times) {
        save_steps.push_back(static_cast<int>(std::round(save_time / p.dt)));
    }

    std::ofstream avg_file(dir + "/average_temperature.dat");
    avg_file << "# time_s average_T_K sensor_T_K\n";
    avg_file << 0.0 << " " << average_temperature(T) << " " << T[p.ic][p.jc] << "\n";

    for (int n = 1; n <= p.nmax; ++n) {
        Tn = T;

        // Heater is always ON in this simplified model.
        int w = 1;
        int w_old = 1;

        compute_R(p, Tn, S, R, w_old);
        gauss_seidel_CN_step(p, R, S, T, w);

        double time = n * p.dt;
        avg_file << time << " " << average_temperature(T) << " " << T[p.ic][p.jc] << "\n";

        if (std::find(save_steps.begin(), save_steps.end(), n) != save_steps.end()) {
            save_temperature_field(dir, n, T, p);
        }
    }
}

int main() {
    std::filesystem::create_directory("geolab");

    Params p;

    std::vector<double> save_times = {
            60.0,    // 1 min  -> T_0006.dat
            300.0,   // 5 min  -> T_0030.dat
            600.0,   // 10 min -> T_0060.dat
            1200.0,  // 20 min -> T_0120.dat
            1800.0   // 30 min -> T_0180.dat
    };

    CN_algorithm("geolab", p, save_times);

    std::cout << "Geolab simulation completed.\n";
    std::cout << "Grid: " << p.Nx << " x " << p.Ny << "\n";
    std::cout << "Output: geolab/T_*.dat and geolab/average_temperature.dat\n";

    return 0;
}
