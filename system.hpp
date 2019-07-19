#ifndef MD_WATER_SYSTEM_HPP
#define MD_WATER_SYSTEM_HPP

#include <cmath>
#include "vector3.h"

struct System {
  int nstep, natom, kmax;
  double dt, T, gamma, boxsize, cutoff, volume, ewcoeff, lj, es;
  Vector3 origin, a1, a2, a3, g1, g2, g3, b1, b2, b3;

  System()
      : nstep{},
        natom{},
        kmax{},
        dt{9999},
        T{9999},
        gamma{9999},
        boxsize{1},
        cutoff{1},
        volume{},
        ewcoeff{},
        lj{},
        es{},
        origin{},
        a1{},
        a2{},
        a3{},
        g1{},
        g2{},
        g3{},
        b1{},
        b2{},
        b3{} {}

  void setup_box() {
    volume = boxsize * boxsize * boxsize;
    a1.x = boxsize;
    a1.y = 0;
    a1.z = 0;
    a2.x = 0;
    a2.y = boxsize;
    a2.z = 0;
    a3.x = 0;
    a3.y = 0;
    a3.z = boxsize;
    const auto rvol2pi = 2. * M_PI / volume;
    g1 = (a2 % a3) * rvol2pi;
    g2 = (a3 % a1) * rvol2pi;
    g3 = (a1 % a2) * rvol2pi;
  }
};

#endif
