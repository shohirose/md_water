#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>

#include "atom.hpp"
#include "function.hpp"
#include "system.hpp"

using namespace std;

int main(int argc, char** argv) {
  using std::cout;

  if (argc < 4) {
    cout << "usage: ./a.out trjout psf pdb config\n";
    return 1;
  }

  srand((unsigned int)time(NULL));

  ofstream fo(argv[1]);
  ifstream fs(argv[2]);
  ifstream fp(argv[3]);
  ifstream fc(argv[4]);
  System system;
  if (!load_config(fc, system)) return 1;
  vector<Atom> atomVector;
  if (!load_psf(fs, system, atomVector)) return 1;

  const int natom = system.natom;
  const int nwat = natom / 3;
  const int nfree = natom * 3 - nwat * 3;

  if (!load_pdb(fp, atomVector)) return 1;

  int print_energy_step = 1;
  int print_trj_step = 1;
  int nstep = system.nstep;

  const double ps2asu = 20.455;
  const double dt_ps = system.dt * 0.001;  // in ps
  const double gamma_ps = system.gamma;    // in ps-1
  const double dt = dt_ps * ps2asu;
  const double gamma = gamma_ps / ps2asu;
  const double dt_div2 = dt * 0.5;
  const double dtdt_div2 = dt * dt * 0.5;
  const double kb = 0.001987191;  // kcal/mol/K
  const double T = system.T;      // K
  const double A = 1. - gamma * dt * 0.5;
  const double B = 1. + gamma * dt * 0.5;
  const double inB = 1. / B;

  cout << "REMARK Number of atoms " << natom << '\n';
  for (int i = 0; i < natom; i++) {
    Atom& at = atomVector[i];
    at.R = sqrt(2. * kb * T * gamma * at.mass / dt);
  }

  cout << "REMARK Number of water molecules " << nwat << '\n';
  cout << "REMARK Degrees of freedom " << nfree << '\n';
  cout << "REMARK dt[fs] " << dt_ps * 1e3 << '\n';
  cout << "REMARK gamma[ps-1] " << gamma_ps << '\n';
  cout << "REMARK T[K] " << T << '\n';
  cout << "REMARK Cutoff[ang.] " << system.cutoff << '\n';

  const double rOH = 0.9572;
  const double aHOH = 104.52 / 180. * acos(-1.0);

  // make reciprocal vectors
  vector<Vector3> g;
  int kmax = system.kmax;
  int sqkmax = kmax - 1;
  sqkmax *= sqkmax;
  double dum = 0;
  for (int k = 0; k < kmax; k++)
    for (int i = -kmax + 1; i < kmax; i++)
      for (int j = -kmax + 1; j < kmax; j++) {
        if (k * k + j * j + i * i > sqkmax or (i == 0 && j == 0 && k == 0))
          continue;
        Vector3 vtmp(i * system.g1.x, j * system.g2.y, k * system.g3.z);
        if (vabs(vtmp) > dum) dum = vabs(vtmp);
        g.push_back(vtmp);
      }
  cout << "REMARK gmax = " << dum << '\n';

  double ew_self = 0;
  for (int i = 0; i < atomVector.size(); i++) {
    ew_self += atomVector[i].charge * atomVector[i].charge;
  }
  ew_self *= -system.ewcoeff / sqrt(M_PI) * 332.0636;

  int icnt = 0;
  for (int i = 0; i < natom; i++) {
    Atom& at = atomVector[i];
    at.velocity.x = gauss() * sqrt(kb * T / at.mass);
    at.velocity.y = gauss() * sqrt(kb * T / at.mass);
    at.velocity.z = gauss() * sqrt(kb * T / at.mass);
  }

  output(fo, atomVector, system);

  vector<int> lj_pair_list, el_pair_list, shake_list;
  make_lj_pair(atomVector, lj_pair_list);
  make_el_pair(atomVector, el_pair_list);
  make_shake_pair(atomVector, shake_list);

  calc_pot(atomVector, lj_pair_list, el_pair_list, system, g);
  double Ktmp = calc_kin(atomVector);
  system.es += ew_self;
  print_ene(0, system, Ktmp, nfree);

  calc_frc(atomVector, lj_pair_list, el_pair_list, system, g);
  for (int i = 0; i < atomVector.size(); i++) {
    Atom& at = atomVector[i];
    at.fold = at.fnew;
    at.rold = at.position;
    at.position = at.rold + dt * at.velocity + dtdt_div2 / at.mass * at.fold;
  }

  double accum = 0;

  for (int istep = 1; istep <= nstep; istep++) {
    calc_frc(atomVector, lj_pair_list, el_pair_list, system, g);
    for (int i = 0; i < atomVector.size(); i++) {
      Vector3 noise(gauss(), gauss(), gauss());
      Atom& at = atomVector[i];
      at.rnew = 2. * at.position - at.rold + gamma * dt_div2 * at.rold +
                dt * dt / at.mass * (at.fnew + at.R * noise);
      at.rnew = at.rnew * inB;
    }

    for (int ishake = 0; ishake < 100; ishake++) {
      if (shake(atomVector, shake_list)) break;
      if (ishake == 99) {
        cerr << "error: shake does not converged.\n";
        return 1;
      }
    }

    for (int i = 0; i < atomVector.size(); i++) {
      Atom& at = atomVector[i];
      at.vnew = 0.5 / dt * (at.rnew - at.rold);
    }

    double K = calc_kin(atomVector);
    calc_pot(atomVector, lj_pair_list, el_pair_list, system, g);
    system.es += ew_self;
    if (istep % print_energy_step == 0) print_ene(istep, system, K, nfree);

    if (istep % print_trj_step == 0) {
      output(fo, atomVector, system);
    }

    for (int i = 0; i < atomVector.size(); i++) {
      Atom& at = atomVector[i];
      at.rold = at.position;
      at.position = at.rnew;
      at.velocity = at.vnew;
    }
  }
}
