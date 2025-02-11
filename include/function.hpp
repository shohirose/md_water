#ifndef MD_WATER_FUNCTION_HPP
#define MD_WATER_FUNCTION_HPP

#include <fstream>
#include <vector>

#include "atom.hpp"
#include "system.hpp"

double gauss();

void print_ene(int istep, System& system, double K, int nfree);

bool shake(std::vector<Atom>& atoms, const std::vector<int>& shake_list);

void calc_frc(std::vector<Atom>& atoms, const std::vector<int>& lj_pair_list,
              const std::vector<int>& el_pair_list, System& system,
              std::vector<Vector3>& g);

double calc_kin(const std::vector<Atom>& atoms);

void calc_pot(std::vector<Atom>& atoms, const std::vector<int>& lj_pair_list,
              const std::vector<int>& el_pair_list, System& system,
              std::vector<Vector3>& g);

void output(std::ofstream& fo, std::vector<Atom>& atoms, System& system);

void make_lj_pair(const std::vector<Atom>& atoms,
                  std::vector<int>& lj_pair_list);

void make_el_pair(const std::vector<Atom>& atoms,
                  std::vector<int>& el_pair_list);

void make_shake_pair(const std::vector<Atom>& atoms,
                     std::vector<int>& shake_list);

bool load_config(std::ifstream& fi, System& system);

bool load_psf(std::ifstream& fs, System& system, std::vector<Atom>& atoms);

bool load_pdb(std::ifstream& fp, std::vector<Atom>& atoms);

#endif
