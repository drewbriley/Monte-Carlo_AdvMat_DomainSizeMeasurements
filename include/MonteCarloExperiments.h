#pragma once
#include "Exciton.h"
#include "Lattice.h"
#include <algorithm>
#include <bits/stdc++.h>
#include <iostream>
#include <omp.h>
#include <random>
#include <unistd.h>
namespace monteCarloExperiments {

double single_exciton_diffusion_3D(int, double, std::string,
                                   monteCarlo::Lattice &);
void setup_multiple_Exciton_PL_3D(double, double, double, double, double,double, double, int, int, int, double, double,double,int, std::string);
std::tuple<double, double, double> run_multiple_Exciton_PL_3D(double, double, double, std::string,monteCarlo::Lattice &, int, int, bool,std::default_random_engine &);
void test(monteCarlo::Lattice &, std::string);
void saveDiffusionCoeffData(std::vector<double>, monteCarlo::Lattice &,std::string, int, int);
void measureDiffusionLength(int, monteCarlo::Lattice &, std::string);
void saveExcitonDiffusionData_3D(std::vector<double>, double, std::vector<double>, std::string);
void domainSizeExp(int, double, double, double, double, double, double,double, int, double, double, double, std::string);
std::string createDirectory();
} // namespace monteCarloExperiments
