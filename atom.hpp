#ifndef MD_WATER_ATOM_HPP
#define MD_WATER_ATOM_HPP

#include <string>
#include "vector3.h"

struct Atom {
  Vector3 position;
  Vector3 velocity;
  Vector3 fold, fnew, vnew, rold, rnew;
  std::string atomname;
  double charge, mass, R;

  Atom(double x, double y, double z)
      : position(x, y, z),
        velocity{},
        fold{},
        fnew{},
        vnew{},
        rold{},
        rnew{},
        atomname{},
        charge{},
        mass{},
        R{} {}
};

#endif
