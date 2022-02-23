#include "MonteCarloExperiments.h"
#include <iostream>

int main() 
/**Doc String:
*This program measures the diffusion length as a function of disorder
*/
{
  int dimension, xBC, yBC, zBC;
  double dt;
  std::string folderName = monteCarloExperiments::createDirectory();
  double disorderStart = 0.001, disorderStop = 500, disorderStep = 0.5;

  dimension = 100, xBC = 1, yBC = 1, zBC = 1, dt = 1e-12;
  std::vector<double> DisorderVector;
  std::vector<double> Ld;
  for (double temp = 400; temp >= 50; temp -= 50) {
    DisorderVector.clear();
    Ld.clear();
    for (double disorder = disorderStart; disorder <= disorderStop;
         disorder *= pow(10, disorderStep)) {
      std::cout << "Disorder: " << disorder << std::endl;
      double D = 1e-7;
      double dx = sqrt(6 * D * dt);

      monteCarlo::Lattice *L =
          new monteCarlo::Lattice(dimension, dimension, dimension, dx, disorder,
                                  0, 300, 1 / dt, 0, xBC, yBC, zBC, folderName);
      Ld.push_back(monteCarloExperiments::single_exciton_diffusion_3D(
          1000, D, folderName, *L));
      std::cout << "Ld: " << Ld[Ld.size() - 1] << std::endl;
      delete L;
      DisorderVector.push_back(disorder);
    }

    monteCarloExperiments::saveExcitonDiffusionData_3D(DisorderVector, temp, Ld,
                                                       folderName);
  }
}
