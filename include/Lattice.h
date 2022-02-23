#pragma once
#include <bits/stdc++.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#define Kb 8.617333262145e-5 // eV/k	//1.38064852e-23//J/K
#define elementaryCharge 1.60217662e-19
namespace monteCarlo {
class Lattice {
private:
  int xsize, ysize, zsize;
  int yBCs, xBCs, zBCs; // boundary conditions for y,x,z,
                        // axis:0-periodic,1-reflecting,2-extracting
  int numOfStreams;
  double averageE, disorder, temperature, escapeFreq, dx;
  std::vector<double> energies;
  std::vector<std::vector<int>> siteOccupation;
  std::vector<std::vector<double>> hoppingRates;
  double calculateEnergy();
  double gaussian(double, double);
  void initalizeEnergiesandOccupation(double *);
  void inializeLattice(int, int, int, double, double, double, double, double,
                       double, std::string);
  // void initalizeRates();
  double calculateRate(int, int);
  void updateHoppingRates(int, std::vector<int> &);

  enum BC { PERIODIC, REFLECTING, EXTRACTING };

public:
  Lattice(int, int, int, double, double, double, double, double, int, int, int,
          int, std::string);
  int onBorder_legacy(int);
  double getDisorder();
  double getAverageE();
  double getTemperature();
  int getXsize();
  int getYsize();
  int getZsize();
  int getXBCs();
  int getYBCs();
  int getZBCs();
  int getNumOfStreams();
  double getEscapeFreq();
  double getdx();
  double getHoppingRate(int, int);
  int getLatticeSize();
  double getSiteEnergy(int);
  int getSiteOccupation(int, int);
  void setSiteOccupation(int, int, int);
  void resetOccupation(int);
  int flattenPosition(int, int, int);
  double calculateDistance(
      int, int); // Returns distance bewteen the two inducies, in units of dx
  std::tuple<int, int, int> openPosition(int);
  std::vector<int> getNeighbours_legacy(int);
  std::vector<int> getNeighbours(int);
};

} // namespace monteCarlo
