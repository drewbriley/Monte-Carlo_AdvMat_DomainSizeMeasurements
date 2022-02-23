#include "Lattice.h"
namespace monteCarlo {

double Lattice::gaussian(double sigma, double mean) 
/**Doc String:
*This function uses Marsaglia Polar Method to calculate a value within a gaussian distribution
*Input:
*	standard deviation of gaussian distribution
*	average of gaussian distribution
*Return:
*	a random value with the gaussian distribution
*/
{
  double w1, w2, R2;
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0, 1.0);
  do {
    w1 = distribution(generator) * 2 - 1;
    w2 = distribution(generator) * 2 - 1;
    R2 = w1 * w1 + w2 * w2;
  } while ((R2 > 1.0) || (R2 == 0));
  return sqrt(-2 * log(R2) / R2) * w1 * sigma + mean;
}

double Lattice::calculateEnergy() 
/**Doc String:
*This function calculates the energy for a lattice point
*/
{ return gaussian(disorder, averageE); }

void Lattice::inializeLattice(int x_passed,               //
                              int y_passed,               //
                              int z_passed,               //
                              double dx_passed,           //
                              double d_passed,            //
                              double averageE_passed,     //
                              double temp_passed,         //
                              double escapeFreq_passed,   //
                              double numOfStreams_passed, //
                              std::string folderName) 
/**Doc String:
*This function initalizes the lattice by acting as the constructor
*Input:
*	X/Y/Z sizes
*	Spacial step size
*	Disorder
*	Average energy of lattice points
*	Temperature of the lattice
*	The escape frequency of the lattice (1/dt)
*	The total number of streams used in the simulation
*/							  
{
  xsize = x_passed;
  ysize = y_passed;
  zsize = z_passed;
  disorder = d_passed;
  averageE = averageE_passed;
  temperature = temp_passed;
  escapeFreq = escapeFreq_passed;
  dx = dx_passed;
  numOfStreams = numOfStreams_passed;
  double n = (double)x_passed * (double)y_passed * (double)z_passed;
  std::cout << n << " " << std::endl;

  std::vector<double> temp = {-1, -1, -1, -1, -1, -1, -1};
  std::vector<int> tempSiteOccupations;
  for (int i = 0; i < n; i++) {
    tempSiteOccupations.push_back(0);
    energies.push_back(calculateEnergy());
    hoppingRates.push_back(temp);
  }
  for (int i = 0; i <= numOfStreams; i++) {
    siteOccupation.push_back(tempSiteOccupations);
  }
}

// Calculate the hopping rate between site i and j
double Lattice::calculateRate(int i, int j) 
/**Doc String:
*This function takes in two lattice points and calcualtes the Miller Abrahams hopping rate between them
*Input:
*	Two lattice points
*Return:
*	The Miller Abrahams hopping rate between the two points
*/
{
  double E = energies[i], Eprime = energies[j], deltaE = Eprime - E, f;
  if (deltaE < 0) {
    f = 1.0;
  } else {
    f = exp(-deltaE / (Kb * (temperature)));
  }
  return (escapeFreq)*exp(-2.0) * f;
}

void Lattice::updateHoppingRates(int index, std::vector<int> &neighbours) 
/**Doc String:
*This function updates and saves the hopping rates between a lattice site and it's nearest neighbours
*Input:
*	Index of the lattice site in question
*	Vector containing nearest neighbours of lattice site in question
*/
{
  double sum = 0, rate[neighbours.size() + 1];
  std::vector<double> temp;
  for (int i = 1; i < neighbours.size() + 1; i++) {
    rate[i] = calculateRate(index, neighbours[i - 1]);
    sum += rate[i];
  }

  rate[0] = sum;
  for (int i = 0; i < 7; i++) {
    temp.push_back(rate[i]);
  }
  hoppingRates[index] = temp;
}

std::tuple<int, int> nn_1d_periodic(int c, int size) 
/**Doc String:
*This function finds the inducies of the nearest neighbours in one dimension of a site for periodic boundary conditions
*Input:
*	Index of site
*	Size of lattice in the dimension in question
*Return:
*	(tuple) index indicating the 'minus' and 'plus' indicies for nearest neighbours
*/
{
  int max = size - 1;
  int cm = (c > 0) ? (c - 1) : max;
  int cp = (c < max) ? (c + 1) : 0;

  return std::make_tuple(cm, cp);
}

std::tuple<int, int> nn_1d_reflecting(int c, int size) 
/**Doc String:
*This function finds the inducies of the nearest neighbours in one dimension of a site for reflecting boundary conditions
*Input:
*	Index of site
*	Size of lattice in the dimension in question
*Return:
*	(tuple) index indicating the 'minus' and 'plus' indicies for nearest neighbours
*/
{
  int max = size - 1;
  int cm = (c > 0) ? (c - 1) : 0;
  int cp = (c < max) ? (c + 1) : max;
  return std::make_tuple(cm, cp);
}

/*Calculates the indicies of neighbouring cells*/
/*returns array with inducies of movement {X--,X++,Y--,Y++,Z--,Z++}*/
std::vector<int> Lattice::getNeighbours(int index) 
/**Doc String:
*This function takes in a index of a lattice site and returns the inducies of the nearest neighbours
*Input:
*	Index of interest
*Return:
*	vector of the nearest neighbours inducies of type {X--,X++,Y--,Y++,Z--,Z++}
*/
{
  int x, y, z;
  std::vector<int> neighbours;
  std::tie(x, y, z) = openPosition(index);
  int xm, xp;
  if (xBCs == BC::PERIODIC)
    std::tie(xm, xp) = nn_1d_periodic(x, xsize);
  else
    std::tie(xm, xp) = nn_1d_reflecting(x, xsize);

  int ym, yp;
  if (yBCs == BC::PERIODIC)
    std::tie(ym, yp) = nn_1d_periodic(y, ysize);
  else
    std::tie(ym, yp) = nn_1d_reflecting(y, ysize);

  int zm, zp;
  if (zBCs == BC::PERIODIC)
    std::tie(zm, zp) = nn_1d_periodic(z, zsize);
  else
    std::tie(zm, zp) = nn_1d_reflecting(z, zsize);

  // clang-format off
  neighbours.push_back(flattenPosition(xm, y , z ));
  neighbours.push_back(flattenPosition(xp, y , z ));
  neighbours.push_back(flattenPosition(x , ym, z ));
  neighbours.push_back(flattenPosition(x , yp, z ));
  neighbours.push_back(flattenPosition(x , y , zm));
  neighbours.push_back(flattenPosition(x , y , zp));
  // clang-format on
  return neighbours;
}

/*Calculates the indicies of neighbouring cells*/
/*returns array with inducies of movement {X--,X++,Y--,Y++,Z--,Z++}*/
std::vector<int> Lattice::getNeighbours_legacy(int index) 
/**Doc String:
*This function takes in a index of a lattice site and returns the inducies of the nearest neighbours, this is an old function and not needed.
*Input:
*	Index of interest
*Return:
*	vector of the nearest neighbours inducies of type {X--,X++,Y--,Y++,Z--,Z++}
*/
{
  int xPos, yPos, zPos;
  std::vector<int> neighbours;
  std::tie(xPos, yPos, zPos) = openPosition(index);
  int borderCheck = onBorder_legacy(index);

  if (borderCheck == 0) {
    // calculate inducies around location
    neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
    neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
    neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
    neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));

  }
  /*Create list of neighbours using BCs*/
  else if (borderCheck == 1) { // lower left corner on bottom face
    if (xBCs == 0) {
      neighbours.push_back(flattenPosition(xsize - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    }
    if (yBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, ysize - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    }
    if (zBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos, zsize - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
    }
  } else if (borderCheck == 2) { // upper left corner on bottom face
    if (xBCs == 0) {
      neighbours.push_back(flattenPosition(xsize - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    }
    if (yBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, 0, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
    if (zBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos, zsize - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
    }
  } else if (borderCheck == 3) { // lower right corner on bottom face
    if (xBCs == 0) {
      neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(0, yPos, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
    if (yBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, ysize - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    }
    if (zBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos, zsize - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
    }
  } else if (borderCheck == 4) { // upper right corner on bottom face
    if (xBCs == 0) {
      neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(0, yPos, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
    if (yBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, 0, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
    if (zBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos, zsize - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
    }
  } else if (borderCheck == 5) { // lower left corner on top face
    if (xBCs == 0) {
      neighbours.push_back(flattenPosition(xsize - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    }
    if (yBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, ysize - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    }
    if (zBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, 0));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
  } else if (borderCheck == 6) { // upper left corner on top face
    if (xBCs == 0) {
      neighbours.push_back(flattenPosition(xsize - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    }
    if (yBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, 0, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
    if (zBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, 0));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
  } else if (borderCheck == 7) { // lower right corner on top face
    if (xBCs == 0) {
      neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(0, yPos, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
    if (yBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, ysize - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    }
    if (zBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, 0));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
  } else if (borderCheck == 8) { // upper right corner on top face
    if (xBCs == 0) {
      neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(0, yPos, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
    if (yBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, 0, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
    if (zBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, 0));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
  } else if (borderCheck == 9) { // lower edge on bottom face
    neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
    neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    if (yBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, ysize - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    }
    if (zBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos, zsize - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
    }
  } else if (borderCheck == 10) { // upper edge on bottom face
    neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
    neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    if (yBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, 0, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
    if (zBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos, zsize - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
    }
  } else if (borderCheck == 11) { // left edge on bottom face
    if (xBCs == 0) {
      neighbours.push_back(flattenPosition(xsize - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    }
    neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
    neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    if (zBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos, zsize - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
    }
  } else if (borderCheck == 12) { // right edge on bottom face
    if (xBCs == 0) {
      neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(0, yPos, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
    neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
    neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    if (zBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos, zsize - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
    }
  } else if (borderCheck == 13) { // lower edge on top face
    neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
    neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    if (yBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, ysize - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    }
    if (zBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, 0));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
  } else if (borderCheck == 14) { // upper edge on top face
    neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
    neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    if (yBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, 0, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
    if (zBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, 0));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
  } else if (borderCheck == 15) { // left edge on top face
    if (xBCs == 0) {
      neighbours.push_back(flattenPosition(xsize - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    }
    neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
    neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    if (zBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, 0));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
  } else if (borderCheck == 16) { // right edge on top face
    if (xBCs == 0) {
      neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(0, yPos, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
    neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
    neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    if (zBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, 0));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
  } else if (borderCheck == 17) { // lower left edge
    if (xBCs == 0) {
      neighbours.push_back(flattenPosition(xsize - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    }
    if (yBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, ysize - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    }
    neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
    neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
  } else if (borderCheck == 18) { // lower right edge
    if (xBCs == 0) {
      neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(0, yPos, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
    if (yBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, ysize - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    }
    neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
    neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
  } else if (borderCheck == 19) { // upper left edge
    if (xBCs == 0) {
      neighbours.push_back(flattenPosition(xsize - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    }
    if (yBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, 0, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
    neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
    neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
  } else if (borderCheck == 20) { // upper right edge
    if (xBCs == 0) {
      neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(0, yPos, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
    if (yBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, 0, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
    neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
    neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
  } else if (borderCheck == 21) { // left face
    if (xBCs == 0) {
      neighbours.push_back(flattenPosition(xsize - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    }
    neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
    neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
    neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
  } else if (borderCheck == 22) { // lower face
    neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
    neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    if (yBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, ysize - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    }
    neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
    neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
  } else if (borderCheck == 23) { // right face
    if (xBCs == 0) {
      neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(0, yPos, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
    neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
    neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
    neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
  } else if (borderCheck == 24) { // upper face
    neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
    neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    if (yBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, 0, zPos));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
    neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
    neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
  } else if (borderCheck == 25) { // bottom face
    neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
    neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
    neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    if (zBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos, zsize - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos + 1));
    }
  } else if (borderCheck == 26) { // top face
    neighbours.push_back(flattenPosition(xPos - 1, yPos, zPos));
    neighbours.push_back(flattenPosition(xPos + 1, yPos, zPos));
    neighbours.push_back(flattenPosition(xPos, yPos - 1, zPos));
    neighbours.push_back(flattenPosition(xPos, yPos + 1, zPos));
    if (zBCs == 0) {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, 0));
    } else {
      neighbours.push_back(flattenPosition(xPos, yPos, zPos - 1));
      neighbours.push_back(flattenPosition(xPos, yPos, zPos));
    }
  }

  return neighbours;
}

/*checks if the index of interest is on a border(s) and returns the appropriate
 * int stating it's location*/
int Lattice::onBorder_legacy(int index) 
/**Doc String
*This function takes in an index and returns if this is on any border. This function no longer neeeded.
*Input:
*	Index of interest
*Returns:
*	(int) an integer relating to a specific border (corner, side, etc) 
*/{
  int xPos, yPos, zPos;
  std::tie(xPos, yPos, zPos) = openPosition(index);
  if (zPos == 0) {
    if (xPos == 0) {
      if (yPos == 0) {
        return 1;
      } // lower left corner on bottom face
      else if (yPos == ysize - 1) {
        return 2;
      } // upper left corner on bottom face
      else {
        return 11;
      } // left edge on bottom face
    } else if (xPos == xsize - 1) {
      if (yPos == 0) {
        return 3;
      } // lower right corner on bottom face
      else if (yPos == ysize - 1) {
        return 4;
      } // upper right corner on bottom face
      else {
        return 12;
      } // right edge on bottom face
    } else if (yPos == 0) {
      return 9;
    } // lower edge on bottom face
    else if (yPos == ysize - 1) {
      return 10;
    } // upper edge on bottom face
    else {
      return 25;
    } // bottom face
  } else if (zPos == zsize - 1) {
    if (xPos == 0) {
      if (yPos == 0) {
        return 5;
      } // lower left corner on top face
      else if (yPos == ysize - 1) {
        return 6;
      } // upper left corner on top face
      else {
        return 15;
      } // left edge on bottom face
    } else if (xPos == xsize - 1) {
      if (yPos == 0) {
        return 7;
      } // lower right corner on top face
      else if (yPos == ysize - 1) {
        return 8;
      } // upper right corner on top face
      else {
        return 16;
      } // right edge on top face
    } else if (yPos == 0) {
      return 13;
    } // lower edge on top face
    else if (yPos == ysize - 1) {
      return 14;
    } // upper edge on top face
    else {
      return 26;
    } // top face
  } else {
    if (xPos == 0) {
      if (yPos == 0) {
        return 17;
      } // lower left edge
      else if (yPos == ysize - 1) {
        return 19;
      } // upper left edge
      else {
        return 21;
      } // left face
    } else if (xPos == xsize - 1) {
      if (yPos == 0) {
        return 18;
      } // lower right edge
      else if (yPos == ysize - 1) {
        return 20;
      } // upper right edge
      else {
        return 23;
      } // right face
    } else if (yPos == 0) {
      return 22;
    } // lower face
    else if (yPos == ysize - 1) {
      return 24;
    } // upper face
    else {
      return 0;
    } // not on a border
  }
}

/*Constructor*/
Lattice::Lattice(int x_passed,             //
                 int y_passed,             //
                 int z_passed,             //
                 double dx_passed,         //
                 double d_passed,          //
                 double averageE_passed,   //
                 double temp_passed,       //
                 double escapeFreq_passed, //
                 int numOfStreams_passed,  //
                 int xBC_passed,           //
                 int yBC_passed,           //
                 int zBC_passed,           //
                 std::string folderName) 
/**Doc String:
*This constructs the lattice
*/
{
  std::cout << "3D constructor" << std::endl;
  xBCs = xBC_passed;
  yBCs = yBC_passed;
  zBCs = zBC_passed;

  inializeLattice(x_passed, y_passed, z_passed, dx_passed, d_passed,
                  averageE_passed, temp_passed, escapeFreq_passed,
                  numOfStreams_passed, folderName);
}

/*Public Methods*/
double Lattice::getDisorder() { return disorder; }																		//Returns the lattice disorder
double Lattice::getAverageE() { return averageE; }																		//Returns the average energy of the lattice points
double Lattice::getTemperature() { return temperature; }																//Returns the temperature of the lattice
int Lattice::getXsize() { return xsize; }																				//Returns the x size of the lattice
int Lattice::getYsize() { return ysize; }																				//Returns the y size of the lattice
int Lattice::getZsize() { return zsize; }																				//Returns the z size of the lattice
int Lattice::getXBCs() { return xBCs; }																					//Returns the boundary conditions on the x-axis
int Lattice::getYBCs() { return yBCs; }																					//Returns the boundary conditions on the y-axis
int Lattice::getZBCs() { return zBCs; }																					//Returns the boundary conditions on the z-axis
int Lattice::getNumOfStreams() { return numOfStreams; }																	//Returns the number of streams created on the lattice
double Lattice::getEscapeFreq() { return escapeFreq; }																	//Returns the escape frequency of the lattice
double Lattice::getdx() { return dx; }																					//Returns the spacial step size of the lattice
int Lattice::getLatticeSize() { return energies.size(); }																//Returns the size of the lattice
int Lattice::getSiteOccupation(int index, int stream){return siteOccupation[stream][index];}							//Returns the number of excitons on the lattice site of interest
double Lattice::getSiteEnergy(int index) { return energies[index]; }													//Returns the site energy of site of interest
void Lattice::setSiteOccupation(int index, int stream, int occupation) {siteOccupation[stream][index] += occupation;}	//Changes the occupation of the site on the stream of interest by the number occupation 
double Lattice::getHoppingRate(int index, int count) 
/**Doc String:
*This function takes in two inducies and returns the hopping rates between them.
*Input:
*	Index of the occupied site
*	Index of the site to hop to
*Return:
*	(double) hopping rate from occupied site to site of interest
*/
{					
  if (hoppingRates[index][0] == -1) {
    std::vector<int> n = getNeighbours_legacy(index);
    updateHoppingRates(index, n);
  }

  return hoppingRates[index][count];
}


void Lattice::resetOccupation(int stream) 
/**Doc String:
*This function sets the occupation of all sites on the lattice in the stream of interest to zero.
*Input:
*	Stream of interest
*/{
  for (int i = 0; i < siteOccupation[stream].size(); i++) {
    siteOccupation[stream][i] = 0;
  }
}

int Lattice::flattenPosition(int x, int y, int z) 
/**Doc String:
*This function takes in an x,y, and z coordinate of a lattice point and returns the index of this point. 
*Input:
*	x,y,z positions
*Return:
*	(int) index of position of interest
*/
{
  return x + (xsize) * (y + ysize * z);
}

std::tuple<int, int, int> Lattice::openPosition(int index) 
/**Doc String:
*This function takes in an index from the lattice and returns the x,y, and z coordinates of the index.
*Input:
*	the index of interest
*Returns:
*	(tuple) the x,y,z coordinates of the lattice index in question
*/
{

  // clang-format off
  int x =  index % xsize;
  int y = (index / xsize) % ysize;
  int z = (index / (xsize * ysize));
  // clang-format on

  return std::make_tuple(x, y, z);
}

double Lattice::calculateDistance(int i, int j) 
/**Doc String:
*This function takes in two inducies from the lattice and calculates the absolute value of the distance between the two points.
*Input:
*	Two inducies
*Returns:
*	The absolute distance between the two inducies
*/
{
  int x1, x2, y1, y2, z1, z2;
  std::tie(x1, y1, z1) = openPosition(i);
  std::tie(x2, y2, z2) = openPosition(j);
  return sqrt((x2 - x1) * (x2 - x1) + //
              (y2 - y1) * (y2 - y1) + //
              (z2 - z1) * (z2 - z1));
}

} // namespace monteCarlo
