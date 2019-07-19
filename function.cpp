#include "function.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iostream>
#include <sstream>

void print_ene(int istep, System& system, double K, int nfree) {
  static double kb = 0.001987191;  // kcal/mol/K
  double totpot = system.es + system.lj;
  using std::setw;
  std::cout << std::setprecision(4) << std::fixed << setw(12) << istep
            << setw(16) << system.lj << setw(16) << system.es << setw(16)
            << totpot << setw(16) << K << setw(16) << totpot + K << setw(16)
            << K * 2. / nfree / kb << '\n';
}

bool shake(std::vector<Atom>& atoms, std::vector<int>& shake_list) {
  static const double eps = 1e-6;
  static const double eps2 = eps * eps;
  static const double rOH = 0.9572;
  static const double aHOH = 104.52 / 180. * acos(-1.0);
  static const double dOH = rOH * rOH;
  static const double rHH = rOH * sin(aHOH / 2.) * 2.;
  static const double dHH = rHH * rHH;
  for (int i = 0; i < shake_list.size() / 2; i++) {
    Atom& at1 = atoms[shake_list[2 * i]];
    Atom& at2 = atoms[shake_list[2 * i + 1]];
    double gamma;
    if (at1.atomname == "O" && at2.atomname == "H") {
      gamma = (dOH - norm(at1.rnew - at2.rnew)) /
              (2. * (1. / at1.mass + 1. / at2.mass) *
               ((at1.position - at2.position) * (at1.rnew - at2.rnew)));
    } else {
      gamma = (dHH - norm(at1.rnew - at2.rnew)) /
              (2. * (1. / at1.mass + 1. / at2.mass) *
               ((at1.position - at2.position) * (at1.rnew - at2.rnew)));
    }
    at1.rnew = at1.rnew + gamma * (at1.position - at2.position) / at1.mass;
    at2.rnew = at2.rnew + gamma * (at2.position - at1.position) / at2.mass;
  }

  for (int i = 0; i < shake_list.size() / 2; i++) {
    Atom& at1 = atoms[shake_list[2 * i]];
    Atom& at2 = atoms[shake_list[2 * i + 1]];
    double r = vabs(at1.rnew - at2.rnew);
    double error;
    if (at1.atomname == "O" && at2.atomname == "H") {
      error = abs(r - rOH);
    } else {
      error = abs(r - rHH);
    }
    if (error > eps) return false;
  }

  return true;
}

double calc_kin(std::vector<Atom>& atoms) {
  double k = 0;
  for (int i = 0; i < atoms.size(); i++) {
    Atom& at = atoms[i];
    k += at.mass * norm(at.vnew);
  }
  return k * 0.5;
}

void calc_frc(std::vector<Atom>& atoms, std::vector<int>& lj_pair_list,
              std::vector<int>& el_pair_list, System& system,
              std::vector<Vector3>& g) {
  static double cutoff = system.cutoff;
  static double cutoff2 = cutoff * cutoff;
  static double ewcoeff = system.ewcoeff;
  static double ewcoeff2 = ewcoeff * ewcoeff;
  static double boxsize = system.boxsize;
  static double factor = 1. / (4. * ewcoeff * ewcoeff);
  static double const_intra = -2. * ewcoeff / sqrt(acos(-1.));
  static double const_recipro = 4. * acos(-1.) / system.volume;
  static double A = 582. * 1e3;
  static double A_2 = A * 2.;
  static double B = 595.0;
  static double C = 332.0636;
  for (int i = 0; i < atoms.size(); i++) {
    atoms[i].fnew.x = 0.;
    atoms[i].fnew.y = 0.;
    atoms[i].fnew.z = 0.;
  }

  for (int i = 0; i < lj_pair_list.size(); i++) {
    Atom& at1 = atoms[lj_pair_list[i]];

    for (int j = i + 1; j < lj_pair_list.size(); j++) {
      Atom& at2 = atoms[lj_pair_list[j]];
      Vector3 del = at1.position - at2.position;
      del.x -= boxsize * floor(del.x / boxsize + 0.5);
      del.y -= boxsize * floor(del.y / boxsize + 0.5);
      del.z -= boxsize * floor(del.z / boxsize + 0.5);
      double r2 = norm(del);
      if (r2 > cutoff2) continue;
      double r6 = r2 * r2 * r2;
      double r8 = r6 * r2;
      Vector3 f12 = 6. / r8 * (B - A_2 / r6) * del;
      at1.fnew -= f12;
      at2.fnew += f12;
    }
    //	if (!i) print(f12);
  }
  //	print(atoms[0].fnew);

  // ewald intra force
  for (int i = 0; i < atoms.size() / 3; i++) {
    Vector3 frc(0., 0., 0.);
    int j = 3 * i;
    Atom& at1 = atoms[j];
    Atom& at2 = atoms[j + 1];
    Atom& at3 = atoms[j + 2];
    Vector3 del1 = at1.position - at2.position;
    del1.x -= boxsize * floor(del1.x / boxsize + 0.5);
    del1.y -= boxsize * floor(del1.y / boxsize + 0.5);
    del1.z -= boxsize * floor(del1.z / boxsize + 0.5);
    double adel1 = vabs(del1);
    double sqadel1 = adel1 * adel1;
    frc = at1.charge * at2.charge * C *
          (const_intra * exp(-ewcoeff2 * sqadel1) +
           erfl(ewcoeff * adel1) / adel1) /
          sqadel1 * del1;
    at1.fnew -= frc;
    at2.fnew += frc;
    Vector3 del2 = at1.position - at3.position;
    del2.x -= boxsize * floor(del2.x / boxsize + 0.5);
    del2.y -= boxsize * floor(del2.y / boxsize + 0.5);
    del2.z -= boxsize * floor(del2.z / boxsize + 0.5);
    double adel2 = vabs(del2);
    double sqadel2 = adel2 * adel2;
    frc = at1.charge * at3.charge * C *
          (const_intra * exp(-ewcoeff2 * sqadel2) +
           erfl(ewcoeff * adel2) / adel2) /
          sqadel2 * del2;
    at1.fnew -= frc;
    at3.fnew += frc;
    Vector3 del3 = at3.position - at2.position;
    del3.x -= boxsize * floor(del3.x / boxsize + 0.5);
    del3.y -= boxsize * floor(del3.y / boxsize + 0.5);
    del3.z -= boxsize * floor(del3.z / boxsize + 0.5);
    double adel3 = vabs(del3);
    double sqadel3 = adel3 * adel3;
    frc = at3.charge * at2.charge * C *
          (const_intra * exp(-ewcoeff2 * sqadel3) +
           erfl(ewcoeff * adel3) / adel3) /
          sqadel3 * del3;
    at3.fnew -= frc;
    at2.fnew += frc;
  }

  // ewald direct space force summation
  for (int i = 0; i < el_pair_list.size() / 2; i++) {
    Atom& at1 = atoms[el_pair_list[2 * i]];
    Atom& at2 = atoms[el_pair_list[2 * i + 1]];

    Vector3 del = at1.position - at2.position;
    del.x -= boxsize * floor(del.x / boxsize + 0.5);
    del.y -= boxsize * floor(del.y / boxsize + 0.5);
    del.z -= boxsize * floor(del.z / boxsize + 0.5);
    double r2 = norm(del);
    if (r2 > cutoff2) continue;
    double r = sqrt(r2);
    Vector3 frc(0., 0., 0.);
    frc = at1.charge * at2.charge * C *
          (-const_intra * exp(-ewcoeff2 * r2) + erfc(ewcoeff * r) / r) / r2 *
          del;
    at1.fnew += frc;
    at2.fnew -= frc;
  }

  // ewald reciprocal space summation
  for (int ii = 0; ii < atoms.size(); ii++) {
    Atom& iat = atoms[ii];
    Vector3 frc(0., 0., 0.);
    for (int i = 0; i < g.size(); i++) {
      double ag = vabs(g[i]);
      double agag = ag * ag;
      double pre, dtmp;
      pre = 0;
      for (int j = 0; j < atoms.size(); j++) {
        Atom& jat = atoms[j];
        Vector3 del = iat.position - jat.position;
        del.x -= boxsize * floor(del.x / boxsize + 0.5);
        del.y -= boxsize * floor(del.y / boxsize + 0.5);
        del.z -= boxsize * floor(del.z / boxsize + 0.5);
        double dot = g[i] * del;

        pre += jat.charge * sin(dot);
      }
      g[i].z ? dtmp = 1. : dtmp = 0.5;
      frc = frc + (dtmp * exp(-agag * factor) / agag * pre) * g[i];
    }
    iat.fnew += frc * const_recipro * C * iat.charge;
    //		cout << vabs(frc) << '\n';
  }
  //	print(atoms[0].fnew);
}

void calc_pot(std::vector<Atom>& atoms, std::vector<int>& lj_pair_list,
              std::vector<int>& el_pair_list, System& system,
              std::vector<Vector3>& g) {
  static double cutoff = system.cutoff;
  static double cutoff2 = cutoff * cutoff;
  static double boxsize = system.boxsize;
  static double ewcoeff = system.ewcoeff;
  static double factor = 1. / (4. * ewcoeff * ewcoeff);
  static double const_recipro = 4. * acos(-1.) / system.volume;
  static double A = 582. * 1e3;
  static double B = 595.0;
  static double qO = -0.834;
  static double qH = 0.417;
  static double C = 332.0636;

  system.lj = 0;
  system.es = 0;
  double ew_direct = 0;
  double ew_recipro = 0;
  double ew_intra = 0;
  double U = 0;

  for (int i = 0; i < lj_pair_list.size(); i++) {
    Atom& iat = atoms[lj_pair_list[i]];
    for (int j = i + 1; j < lj_pair_list.size(); j++) {
      Atom& jat = atoms[lj_pair_list[j]];
      Vector3 del = iat.position - jat.position;
      del.x -= boxsize * floor(del.x / boxsize + 0.5);
      del.y -= boxsize * floor(del.y / boxsize + 0.5);
      del.z -= boxsize * floor(del.z / boxsize + 0.5);
      double roo2 = norm(del);
      if (roo2 > cutoff2) continue;
      double roo6 = roo2 * roo2 * roo2;
      system.lj += A / (roo6 * roo6) - B / roo6;
    }
  }

  // ewald intra energy
  for (int i = 0; i < atoms.size() / 3; i++) {
    int j = 3 * i;
    Atom& at1 = atoms[j];
    Atom& at2 = atoms[j + 1];
    Atom& at3 = atoms[j + 2];
    Vector3 del1 = at1.position - at2.position;
    del1.x -= boxsize * floor(del1.x / boxsize + 0.5);
    del1.y -= boxsize * floor(del1.y / boxsize + 0.5);
    del1.z -= boxsize * floor(del1.z / boxsize + 0.5);
    double adel1 = vabs(del1);
    Vector3 del2 = at1.position - at3.position;
    del2.x -= boxsize * floor(del2.x / boxsize + 0.5);
    del2.y -= boxsize * floor(del2.y / boxsize + 0.5);
    del2.z -= boxsize * floor(del2.z / boxsize + 0.5);
    double adel2 = vabs(del2);
    Vector3 del3 = at3.position - at2.position;
    del3.x -= boxsize * floor(del3.x / boxsize + 0.5);
    del3.y -= boxsize * floor(del3.y / boxsize + 0.5);
    del3.z -= boxsize * floor(del3.z / boxsize + 0.5);
    double adel3 = vabs(del3);
    ew_intra += at1.charge * at2.charge * erfl(ewcoeff * adel1) / adel1;
    ew_intra += at1.charge * at3.charge * erfl(ewcoeff * adel2) / adel2;
    ew_intra += at3.charge * at2.charge * erfl(ewcoeff * adel3) / adel3;
  }

  // ewald direct space summation
  for (int i = 0; i < el_pair_list.size() / 2; i++) {
    Atom& at1 = atoms[el_pair_list[2 * i]];
    Atom& at2 = atoms[el_pair_list[2 * i + 1]];

    Vector3 del = at1.position - at2.position;
    del.x -= boxsize * floor(del.x / boxsize + 0.5);
    del.y -= boxsize * floor(del.y / boxsize + 0.5);
    del.z -= boxsize * floor(del.z / boxsize + 0.5);
    double r = vabs(del);
    if (r > cutoff) continue;
    ew_direct += at1.charge * at2.charge * erfc(ewcoeff * r) / r;
  }

  // ewald reciprocal space summation
  for (int i = 0; i < g.size(); i++) {
    double ag = vabs(g[i]);
    double agag = ag * ag;
    double re, im, dtmp;
    re = 0;
    im = 0;
    for (int j = 0; j < atoms.size(); j++) {
      Vector3 del = atoms[j].position;
      del.x -= boxsize * floor(del.x / boxsize + 0.5);
      del.y -= boxsize * floor(del.y / boxsize + 0.5);
      del.z -= boxsize * floor(del.z / boxsize + 0.5);
      double dot = g[i] * del;

      re += atoms[j].charge * cos(dot);
      im += atoms[j].charge * sin(dot);
    }
    g[i].z ? dtmp = 1. : dtmp = 0.5;
    ew_recipro += dtmp * exp(-agag * factor) / agag * (re * re + im * im);
  }

  ew_recipro *= const_recipro;

  //	cout << "el_direct  " << el_direct * 332.0636 << '\n';
  //	cout << "el_intra    " << el_intra * 332.0636 << '\n';
  //	cout << "el_recipro " << el_recipro * 332.0636<< '\n';

  system.es += ew_direct + ew_recipro - ew_intra;
  system.es *= C;
}

void output(std::ofstream& fo, std::vector<Atom>& atoms, System& system) {
  static const double boxsize = system.boxsize;
  fo << atoms.size() << '\n' << '\n';
  fo << std::setprecision(6) << std::fixed;
  for (int i = 0; i < atoms.size() / 3; i++) {
    int j = 3 * i;
    Atom& at1 = atoms[j];
    Atom& at2 = atoms[j + 1];
    Atom& at3 = atoms[j + 2];
    Vector3 p1 = at1.position;
    Vector3 p2 = at2.position;
    Vector3 p3 = at3.position;
    Vector3 vcom = (p1 + p2 + p3) / 3.;

    double d = floor(0.5 + vcom.x / boxsize);
    p1.x -= boxsize * d;
    p2.x -= boxsize * d;
    p3.x -= boxsize * d;

    d = floor(0.5 + vcom.y / boxsize);
    p1.y -= boxsize * d;
    p2.y -= boxsize * d;
    p3.y -= boxsize * d;

    d = floor(0.5 + vcom.z / boxsize);
    p1.z -= boxsize * d;
    p2.z -= boxsize * d;
    p3.z -= boxsize * d;

    using std::setw;
    fo << " " << at1.atomname << setw(20) << p1.x << setw(20) << p1.y
       << setw(20) << p1.z << '\n';
    fo << " " << at2.atomname << setw(20) << p2.x << setw(20) << p2.y
       << setw(20) << p2.z << '\n';
    fo << " " << at3.atomname << setw(20) << p3.x << setw(20) << p3.y
       << setw(20) << p3.z << '\n';
  }
}

double gauss() {
  const double A1 = 3.949846138, A3 = 0.252408784;
  const double A5 = 0.076542912, A7 = 0.008355968;
  const double A9 = 0.029899776;
  double sum = 0.;
  for (int i = 0; i < 12; i++) {
    sum += rand() / static_cast<double>(RAND_MAX);
  }

  double R = (sum - 6.) / 4.;
  double R2 = R * R;
  return ((((A9 * R2 + A7) * R2 + A5) * R2 + A3) * R2 + A1) * R;
}

void make_lj_pair(std::vector<Atom>& atoms, std::vector<int>& lj_pair_list) {
  for (int i = 0; i < atoms.size(); i++) {
    if (atoms[i].atomname == "O") lj_pair_list.push_back(i);
  }
}

void make_el_pair(std::vector<Atom>& atoms, std::vector<int>& el_pair_list) {
  for (int i = 0; i < atoms.size(); i++) {
    Atom& at1 = atoms[i];

    for (int j = i + 1; j < atoms.size(); j++) {
      Atom& at2 = atoms[j];

      if (i / 3 == j / 3) continue;

      el_pair_list.push_back(i);
      el_pair_list.push_back(j);
    }
  }
}

void make_shake_pair(std::vector<Atom>& atoms, std::vector<int>& shake_list) {
  for (int i = 0; i < atoms.size() / 3; i++) {
    int j = 3 * i;
    shake_list.push_back(j);
    shake_list.push_back(j + 1);
    shake_list.push_back(j);
    shake_list.push_back(j + 2);
    shake_list.push_back(j + 1);
    shake_list.push_back(j + 2);
  }
}

bool load_config(std::ifstream& fc, System& system) {
  std::string s;
  while (getline(fc, s)) {
    std::istringstream is(s);
    std::string stmp;
    if (s.find("#", 0) != std::string::npos) {
      // found comment line.
      continue;
    } else if (s.find("dt", 0) != std::string::npos) {
      is >> stmp >> system.dt;
    } else if (s.find("temperature", 0) != std::string::npos) {
      is >> stmp >> system.T;
    } else if (s.find("gamma", 0) != std::string::npos) {
      is >> stmp >> system.gamma;
    } else if (s.find("numstep", 0) != std::string::npos) {
      is >> stmp >> system.nstep;
    } else if (s.find("boxsize", 0) != std::string::npos) {
      is >> stmp >> system.boxsize;
      system.setup_box();
    } else if (s.find("origin", 0) != std::string::npos) {
      is >> stmp >> system.origin.x >> system.origin.y >> system.origin.z;
    } else if (s.find("cutoff", 0) != std::string::npos) {
      is >> stmp >> system.cutoff;
    } else if (s.find("ewcoeff", 0) != std::string::npos) {
      is >> stmp >> system.ewcoeff;
    } else if (s.find("kmax", 0) != std::string::npos) {
      is >> stmp >> system.kmax;
    } else {
      is >> stmp;
      std::cerr << "unknown config param \"" << stmp << "\" found.\n";
      return false;
    }
  }

  return true;
}

bool load_psf(std::ifstream& fs, System& system, std::vector<Atom>& atoms) {
  std::string s;
  while (getline(fs, s)) {
    if (s.find("NATOM", 10) != std::string::npos) {
      std::istringstream is(s);
      is >> system.natom;
      for (int i = 0; i < system.natom; i++) {
        Atom at(0., 0., 0.);
        getline(fs, s);
        std::istringstream iss(s);
        int itmp;
        std::string stmp, atname;
        iss >> itmp >> stmp >> itmp >> stmp >> atname >> stmp >> at.charge >>
            at.mass;
        // cout << atname << '\n';
        // cout << s << '\n';
        at.atomname = atname.substr(0, 1);
        atoms.push_back(at);
      }
      return true;
    }
  }
  return false;
}

bool load_pdb(std::ifstream& fp, std::vector<Atom>& atoms) {
  std::string s;
  int icnt = 0;
  while (getline(fp, s)) {
    if (s.find("ATOM", 0) != std::string::npos) {
      Atom& at = atoms[icnt++];
      std::istringstream is(s.substr(30, 24));
      is >> at.position.x >> at.position.y >> at.position.z;
    }
  }
  return atoms.size() == icnt ? true : false;
}
