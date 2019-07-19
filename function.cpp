#include "function.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iostream>
#include <sstream>

void print_ene(int istep, System& system, double K, int nfree) {
  constexpr auto kb = 0.001987191;  // kcal/mol/K
  double totpot = system.es + system.lj;
  using std::setw;
  std::cout << std::setprecision(4) << std::fixed  //
            << setw(12) << istep                   //
            << setw(16) << system.lj               //
            << setw(16) << system.es               //
            << setw(16) << totpot                  //
            << setw(16) << K                       //
            << setw(16) << totpot + K              //
            << setw(16) << K * 2. / nfree / kb << std::endl;
}

bool shake(std::vector<Atom>& atoms, const std::vector<int>& shake_list) {
  constexpr auto eps = 1e-6;
  // constexpr auto eps2 = eps * eps;
  constexpr auto rOH = 0.9572;
  constexpr auto aHOH = 104.52 / 180. * acos(-1.0);
  constexpr auto dOH = rOH * rOH;
  constexpr auto rHH = rOH * sin(aHOH / 2.) * 2.;
  constexpr auto dHH = rHH * rHH;

  for (size_t i = 0; i < shake_list.size() / 2; i++) {
    auto& at1 = atoms[shake_list[2 * i]];
    auto& at2 = atoms[shake_list[2 * i + 1]];

    const auto gamma =
        (at1.atomname == "O" && at2.atomname == "H")
            ? (dOH - norm(at1.rnew - at2.rnew)) /
                  (2. * (1. / at1.mass + 1. / at2.mass) *
                   ((at1.position - at2.position) * (at1.rnew - at2.rnew)))
            : (dHH - norm(at1.rnew - at2.rnew)) /
                  (2. * (1. / at1.mass + 1. / at2.mass) *
                   ((at1.position - at2.position) * (at1.rnew - at2.rnew)));

    at1.rnew = at1.rnew + gamma * (at1.position - at2.position) / at1.mass;
    at2.rnew = at2.rnew + gamma * (at2.position - at1.position) / at2.mass;
  }

  for (size_t i = 0; i < shake_list.size() / 2; i++) {
    auto& at1 = atoms[shake_list[2 * i]];
    auto& at2 = atoms[shake_list[2 * i + 1]];

    const auto r = vabs(at1.rnew - at2.rnew);
    const auto error = (at1.atomname == "O" && at2.atomname == "H")
                           ? abs(r - rOH)
                           : abs(r - rHH);
    if (error > eps) return false;
  }

  return true;
}

double calc_kin(const std::vector<Atom>& atoms) {
  double k = 0;
  for (size_t i = 0; i < atoms.size(); i++) {
    const auto& at = atoms[i];
    k += at.mass * norm(at.vnew);
  }
  return k * 0.5;
}

void calc_frc(std::vector<Atom>& atoms, const std::vector<int>& lj_pair_list,
              const std::vector<int>& el_pair_list, System& system,
              std::vector<Vector3>& g) {
  const auto cutoff = system.cutoff;
  const auto cutoff2 = cutoff * cutoff;
  const auto ewcoeff = system.ewcoeff;
  const auto ewcoeff2 = ewcoeff * ewcoeff;
  const auto boxsize = system.boxsize;
  const auto factor = 1. / (4. * ewcoeff * ewcoeff);
  const auto const_intra = -2. * ewcoeff / sqrt(acos(-1.));
  const auto const_recipro = 4. * acos(-1.) / system.volume;

  constexpr auto A = 582. * 1e3;
  constexpr auto A_2 = A * 2.;
  constexpr auto B = 595.0;
  constexpr auto C = 332.0636;

  for (size_t i = 0; i < atoms.size(); i++) {
    atoms[i].fnew.x = 0.;
    atoms[i].fnew.y = 0.;
    atoms[i].fnew.z = 0.;
  }

  for (size_t i = 0; i < lj_pair_list.size(); i++) {
    auto& at1 = atoms[lj_pair_list[i]];

    for (size_t j = i + 1; j < lj_pair_list.size(); j++) {
      auto& at2 = atoms[lj_pair_list[j]];
      Vector3 del = at1.position - at2.position;
      del.x -= boxsize * floor(del.x / boxsize + 0.5);
      del.y -= boxsize * floor(del.y / boxsize + 0.5);
      del.z -= boxsize * floor(del.z / boxsize + 0.5);
      const double r2 = norm(del);

      if (r2 > cutoff2) continue;

      const auto r6 = r2 * r2 * r2;
      const auto r8 = r6 * r2;
      const Vector3 f12 = 6. / r8 * (B - A_2 / r6) * del;
      at1.fnew -= f12;
      at2.fnew += f12;
    }
    //	if (!i) print(f12);
  }
  //	print(atoms[0].fnew);

  // ewald intra force
  for (size_t i = 0; i < atoms.size() / 3; i++) {
    Vector3 frc(0., 0., 0.);

    const auto j = 3 * i;
    auto& at1 = atoms[j];
    auto& at2 = atoms[j + 1];
    auto& at3 = atoms[j + 2];

    Vector3 del1 = at1.position - at2.position;
    del1.x -= boxsize * floor(del1.x / boxsize + 0.5);
    del1.y -= boxsize * floor(del1.y / boxsize + 0.5);
    del1.z -= boxsize * floor(del1.z / boxsize + 0.5);

    const double adel1 = vabs(del1);
    const double sqadel1 = adel1 * adel1;
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

    const double adel2 = vabs(del2);
    const double sqadel2 = adel2 * adel2;
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

    const double adel3 = vabs(del3);
    const double sqadel3 = adel3 * adel3;
    frc = at3.charge * at2.charge * C *
          (const_intra * exp(-ewcoeff2 * sqadel3) +
           erfl(ewcoeff * adel3) / adel3) /
          sqadel3 * del3;
    at3.fnew -= frc;
    at2.fnew += frc;
  }

  // ewald direct space force summation
  for (size_t i = 0; i < el_pair_list.size() / 2; i++) {
    auto& at1 = atoms[el_pair_list[2 * i]];
    auto& at2 = atoms[el_pair_list[2 * i + 1]];

    Vector3 del = at1.position - at2.position;
    del.x -= boxsize * floor(del.x / boxsize + 0.5);
    del.y -= boxsize * floor(del.y / boxsize + 0.5);
    del.z -= boxsize * floor(del.z / boxsize + 0.5);

    const double r2 = norm(del);
    if (r2 > cutoff2) continue;

    const double r = sqrt(r2);

    Vector3 frc(0., 0., 0.);
    frc = at1.charge * at2.charge * C *
          (-const_intra * exp(-ewcoeff2 * r2) + erfc(ewcoeff * r) / r) / r2 *
          del;
    at1.fnew += frc;
    at2.fnew -= frc;
  }

  // ewald reciprocal space summation
  for (size_t ii = 0; ii < atoms.size(); ii++) {
    auto& iat = atoms[ii];
    Vector3 frc(0., 0., 0.);

    for (size_t i = 0; i < g.size(); i++) {
      double ag = vabs(g[i]);
      double agag = ag * ag;
      double pre, dtmp;
      pre = 0;

      for (size_t j = 0; j < atoms.size(); j++) {
        auto& jat = atoms[j];
        Vector3 del = iat.position - jat.position;
        del.x -= boxsize * floor(del.x / boxsize + 0.5);
        del.y -= boxsize * floor(del.y / boxsize + 0.5);
        del.z -= boxsize * floor(del.z / boxsize + 0.5);
        double dot = g[i] * del;

        pre += jat.charge * sin(dot);
      }
      dtmp = g[i].z ? 1. : 0.5;
      frc = frc + (dtmp * exp(-agag * factor) / agag * pre) * g[i];
    }
    iat.fnew += frc * const_recipro * C * iat.charge;
    //		cout << vabs(frc) << '\n';
  }
  //	print(atoms[0].fnew);
}

void calc_pot(std::vector<Atom>& atoms, const std::vector<int>& lj_pair_list,
              const std::vector<int>& el_pair_list, System& system,
              std::vector<Vector3>& g) {
  const double cutoff = system.cutoff;
  const double cutoff2 = cutoff * cutoff;
  const double boxsize = system.boxsize;
  const double ewcoeff = system.ewcoeff;
  const double factor = 1. / (4. * ewcoeff * ewcoeff);
  const double const_recipro = 4. * acos(-1.) / system.volume;

  constexpr auto A = 582. * 1e3;
  constexpr auto B = 595.0;
  // constexpr auto qO = -0.834;
  // constexpr auto qH = 0.417;
  constexpr auto C = 332.0636;

  system.lj = 0;
  system.es = 0;
  double ew_direct = 0;
  double ew_recipro = 0;
  double ew_intra = 0;
  // double U = 0;

  for (size_t i = 0; i < lj_pair_list.size(); i++) {
    auto& iat = atoms[lj_pair_list[i]];

    for (size_t j = i + 1; j < lj_pair_list.size(); j++) {
      auto& jat = atoms[lj_pair_list[j]];

      Vector3 del = iat.position - jat.position;
      del.x -= boxsize * floor(del.x / boxsize + 0.5);
      del.y -= boxsize * floor(del.y / boxsize + 0.5);
      del.z -= boxsize * floor(del.z / boxsize + 0.5);

      const auto roo2 = norm(del);
      if (roo2 > cutoff2) continue;

      const auto roo6 = roo2 * roo2 * roo2;
      system.lj += A / (roo6 * roo6) - B / roo6;
    }
  }

  // ewald intra energy
  for (size_t i = 0; i < atoms.size() / 3; i++) {
    const auto j = 3 * i;
    auto& at1 = atoms[j];
    auto& at2 = atoms[j + 1];
    auto& at3 = atoms[j + 2];

    Vector3 del1 = at1.position - at2.position;
    del1.x -= boxsize * floor(del1.x / boxsize + 0.5);
    del1.y -= boxsize * floor(del1.y / boxsize + 0.5);
    del1.z -= boxsize * floor(del1.z / boxsize + 0.5);
    const auto adel1 = vabs(del1);

    Vector3 del2 = at1.position - at3.position;
    del2.x -= boxsize * floor(del2.x / boxsize + 0.5);
    del2.y -= boxsize * floor(del2.y / boxsize + 0.5);
    del2.z -= boxsize * floor(del2.z / boxsize + 0.5);
    const auto adel2 = vabs(del2);

    Vector3 del3 = at3.position - at2.position;
    del3.x -= boxsize * floor(del3.x / boxsize + 0.5);
    del3.y -= boxsize * floor(del3.y / boxsize + 0.5);
    del3.z -= boxsize * floor(del3.z / boxsize + 0.5);
    const auto adel3 = vabs(del3);

    ew_intra += at1.charge * at2.charge * erfl(ewcoeff * adel1) / adel1;
    ew_intra += at1.charge * at3.charge * erfl(ewcoeff * adel2) / adel2;
    ew_intra += at3.charge * at2.charge * erfl(ewcoeff * adel3) / adel3;
  }

  // ewald direct space summation
  for (size_t i = 0; i < el_pair_list.size() / 2; i++) {
    auto& at1 = atoms[el_pair_list[2 * i]];
    auto& at2 = atoms[el_pair_list[2 * i + 1]];

    Vector3 del = at1.position - at2.position;
    del.x -= boxsize * floor(del.x / boxsize + 0.5);
    del.y -= boxsize * floor(del.y / boxsize + 0.5);
    del.z -= boxsize * floor(del.z / boxsize + 0.5);
    const auto r = vabs(del);
    if (r > cutoff) continue;
    ew_direct += at1.charge * at2.charge * erfc(ewcoeff * r) / r;
  }

  // ewald reciprocal space summation
  for (size_t i = 0; i < g.size(); i++) {
    const auto ag = vabs(g[i]);
    const auto agag = ag * ag;
    double re, im, dtmp;
    re = 0;
    im = 0;

    for (size_t j = 0; j < atoms.size(); j++) {
      Vector3 del = atoms[j].position;
      del.x -= boxsize * floor(del.x / boxsize + 0.5);
      del.y -= boxsize * floor(del.y / boxsize + 0.5);
      del.z -= boxsize * floor(del.z / boxsize + 0.5);
      const auto dot = g[i] * del;

      re += atoms[j].charge * cos(dot);
      im += atoms[j].charge * sin(dot);
    }
    dtmp = g[i].z ? 1. : 0.5;
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
  const double boxsize = system.boxsize;

  fo << atoms.size() << '\n' << '\n';
  fo << std::setprecision(6) << std::fixed;

  for (size_t i = 0; i < atoms.size() / 3; i++) {
    const auto j = 3 * i;
    auto& at1 = atoms[j];
    auto& at2 = atoms[j + 1];
    auto& at3 = atoms[j + 2];

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
    fo << " " << at1.atomname  //
       << setw(20) << p1.x     //
       << setw(20) << p1.y     //
       << setw(20) << p1.z << '\n'
       << " " << at2.atomname  //
       << setw(20) << p2.x     //
       << setw(20) << p2.y     //
       << setw(20) << p2.z << '\n'
       << " " << at3.atomname  //
       << setw(20) << p3.x     //
       << setw(20) << p3.y     //
       << setw(20) << p3.z << '\n';
  }
}

double gauss() {
  constexpr auto A1 = 3.949846138;
  constexpr auto A3 = 0.252408784;
  constexpr auto A5 = 0.076542912;
  constexpr auto A7 = 0.008355968;
  constexpr auto A9 = 0.029899776;
  double sum = 0.;
  for (int i = 0; i < 12; i++) {
    sum += rand() / static_cast<double>(RAND_MAX);
  }

  const auto R = (sum - 6.) / 4.;
  const auto R2 = R * R;
  return ((((A9 * R2 + A7) * R2 + A5) * R2 + A3) * R2 + A1) * R;
}

void make_lj_pair(const std::vector<Atom>& atoms,
                  std::vector<int>& lj_pair_list) {
  for (size_t i = 0; i < atoms.size(); i++) {
    if (atoms[i].atomname == "O") lj_pair_list.push_back(i);
  }
}

void make_el_pair(const std::vector<Atom>& atoms,
                  std::vector<int>& el_pair_list) {
  for (size_t i = 0; i < atoms.size(); i++) {
    // const auto& at1 = atoms[i];
    for (size_t j = i + 1; j < atoms.size(); j++) {
      // const auto& at2 = atoms[j];
      if (i / 3 == j / 3) continue;

      el_pair_list.push_back(i);
      el_pair_list.push_back(j);
    }
  }
}

void make_shake_pair(const std::vector<Atom>& atoms,
                     std::vector<int>& shake_list) {
  for (size_t i = 0; i < atoms.size() / 3; i++) {
    const auto j = 3 * i;
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
  size_t icnt = 0;

  while (getline(fp, s)) {
    if (s.find("ATOM", 0) != std::string::npos) {
      Atom& at = atoms[icnt++];
      std::istringstream is(s.substr(30, 24));
      is >> at.position.x >> at.position.y >> at.position.z;
    }
  }
  return atoms.size() == icnt;
}
