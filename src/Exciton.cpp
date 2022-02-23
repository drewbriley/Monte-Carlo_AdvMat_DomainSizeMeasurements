#include "Exciton.h"
#include <bits/stdc++.h>
namespace monteCarlo {
/*Private methods*/
/*moves Exciton to neirest neighbours*/
void Exciton::moveExciton(monteCarlo::Lattice &L,
                          std::default_random_engine &generator) 
/**Doc String:
*This function implements the Gillipse algorithm to move the exciton to an adjacent site.
*Input:
*	Lattice
*	Random Number Generator
*/						  
{
  std::uniform_real_distribution<double> distribution(0, 1.0);
  int tempX = xPos, tempY = yPos, tempZ = zPos;
  L.setSiteOccupation(index, excitonStream, -1);
  // generate random number from 0-1 and check gillipse algorithum
  double r = distribution(generator);
  bool done = false; //(double)rand()/RAND_MAX;bool done = false;
  double sum = 0, i = 0, totalRate = L.getHoppingRate(index, 0);
  do {
    i++;
    sum += L.getHoppingRate(index, i) / totalRate;
  } while (sum <= r);
  i--;
  std::vector<int> n = L.getNeighbours(index);
  index = n[i];
  std::tie(xPos, yPos, zPos) = L.openPosition(index);
  L.setSiteOccupation(index, excitonStream, 1);
  // Check if hop occurred and track movement
  if (index != L.flattenPosition(tempX, tempY, tempZ)) {
    // clang-format off
    if      (i == 0) xhops--; //
    else if (i == 1) xhops++; //
    else if (i == 2) yhops--; //
    else if (i == 3) yhops++; //
    else if (i == 4) zhops--; //
    else if (i == 5) zhops++; //
    // clang-format on
  }
}

double Exciton::exponentialRand(double k,
                                std::default_random_engine &generator) 
/**Doc String:
*This returns a number from an exponential distribution with variance of k
*Input:
*	Variance of the distribution
*Returns:
*	(double) random number on distribution
*/								
{
  std::uniform_real_distribution<double> distribution(0, 1.0);
  double R = distribution(generator); //(double)rand()/RAND_MAX;
  return -1 / k * log(1 - R);
}

double Exciton::calculateDwellTime(double gamma,std::default_random_engine &generator) 
/**Doc String:
*This function calculates the dwell time for an exciton given the total rate of transfer away from the lattice site.
*Input:
*	The total rate of transfer away from the lattice site
*	Random number generator
*Returns:
*	(double) dwell time of exciton
*/
{
  double T = exponentialRand(1.0, generator);
  double DT = T / gamma;
  if (std::isinf(DT)) {
    std::cout << "ERROR, DwellTime T: " << T << "/" << gamma << "\t"
              << T / gamma << std::endl;
  }
  return T / gamma;
}

/*Saving methods*/
void Exciton::saveExcitonRoute(std::string folderName, monteCarlo::Lattice &L,
                               int excitonNum) 
/**Doc String:
*This function saves the route an exction took throughout the simultion
*Input:
*	folder name for data saving
*	Lattice
*	Total number of excitons
*/							   
{
  if (folderName != "NAN") {
    double dis = L.getDisorder();
    std::stringstream strs;
    strs << excitonNum;
    std::string str = strs.str();
    std::stringstream strs1;
    strs1 << dis;
    std::string str1 = strs1.str();
    std::string fileName = folderName + "\\ExcitonRoutes\\ExcitonRoute_" +
                           str1 + "_" + str + ".txt";
    std::ofstream output;
    output.open(fileName);
    output << "Lattice Parameters:\n"
           << "\t xsize: " << L.getXsize() << "\n"
           << "\t ysize: " << L.getYsize() << "\n"
           << "\t Average E: " << L.getAverageE() << "eV\n"
           << "\t Disorder: " << L.getDisorder() << "eV\n"
           << "\t Temperature: " << L.getTemperature() << "K\n"
           << "Exciton Parameters:\n"
           << "\tStarting Position: (" << startIndex << ")\n"
           << "\tLifeTime: " << getTimeOfRecombination() << " s\n"

           << "Data:\n"
           << "Index\tx\ty\tz\tEnergy (eV)\t DwellTime (s)\n";
    int ind, X, Y, Z;
    for (int i = 0; i < ExcitonRoute.size(); i++) {
      ind = ExcitonRoute[i];
      std::tie(X, Y, Z) = L.openPosition(ind);
      output << ind << "\t" << X << "\t" << Y << "\t" << Z << "\t"
             << L.getSiteEnergy(ind) << "\t" << dwellTimes[i] << "\t"
             << "\t" << xhops << "\t" << yhops << "\t" << zhops << "\t"
             << getLengthMoved() << std::endl;
    }
    output.close();
  }
}

void Exciton::construct(monteCarlo::Lattice &L, double ensembleLifetime,
                        int stream, bool Running_Extraction_exp,
                        std::default_random_engine &generator) 
/**Doc String:
*This function acts as the constructor for each exciton
*Input:
*	Lattice
*	Ensemble singlet lifetime
*	Stream number for 
*	Bool to check if running extraction check, makes sure the exciton is generated in the domain
*	Random number generator 
*/
{
  excitonStream = stream;
  decayed = false;
  stationaryTime = 0;

  std::uniform_int_distribution<int> distribution(
      0, L.getXsize() * L.getYsize() * L.getZsize());
  if (Running_Extraction_exp == true) {
    double radius = double(L.getXsize() / 2.0), distance;
    int centerIndex = L.flattenPosition(floor(L.getXsize() / 2),
                                        floor(L.getYsize() / 2), 0),
        position;
    do {
      index = distribution(generator); // rand()%(L.getXsize()*L.getYsize()*L.getZsize());
      std::tie(xPos, yPos, zPos) = L.openPosition(index);
      position = L.flattenPosition(xPos, yPos, 0);
      distance = L.calculateDistance(centerIndex, position);
    } while (distance > radius);
  } else {
    index = distribution(generator);
  } // rand()%(L.getXsize()*L.getYsize()*L.getZsize());}
  L.setSiteOccupation(index, excitonStream, 1);
  std::tie(xPos, yPos, zPos) = L.openPosition(index);
  startIndex = index;
  ExcitonRoute.push_back(index);
  LatticeEnergies.push_back(L.getSiteEnergy(index));
  averageLifeTime = ensembleLifetime;
  double totalRate = L.getHoppingRate(index, 0);
  dwellTimes.push_back(calculateDwellTime(totalRate, generator));
  xhops = 0;
  yhops = 0;
  zhops = 0;
  typeOfRecombination = -1;
}

/*Constructor*/
Exciton::Exciton(monteCarlo::Lattice &L, double ensembleLifetime, int stream,
                 bool Running_Extraction_exp,
                 std::default_random_engine &generator) 
/**Doc String:
*Calls the constructor function
*/				 
{
  construct(L, ensembleLifetime, stream, Running_Extraction_exp, generator);
}

/*Public Methods returning exciton parameters*/
int Exciton::getXpos() { return xPos; } 													//Returns the excitons x grid position
int Exciton::getYpos() { return yPos; } 													//Return the exctions y grid position
int Exciton::getZpos() { return zPos; }														//Return the exctions z grid position
int Exciton::getIndex() {return index;}														//Returns the excitons index on lattice L
int Exciton::getStream() { return excitonStream; }											//Returns the stream that the exciton is on
int Exciton::getExcitonRouteLength() {return ExcitonRoute.size();}							//Returns the number of hops the exciton has taken.
int Exciton::getStartIndex() {return ExcitonRoute[0];}										//Returns the excitons starting index on lattice L
double Exciton::getLifeTime() { return averageLifeTime; }									//Returns the average lifetime of the exciton ensemble
double Exciton::getTimeOfRecombination() { return timeOfRecombination; }					//Returns the time of recombination undergone by the exciton
int Exciton::getTypeOfRecombination() { return typeOfRecombination; }						//Returns the type of recombination undergone by the exciton

bool Exciton::getDecayed() { return decayed; }
double Exciton::getLengthMoved() {return xhops * xhops + yhops * yhops + zhops * zhops;}	//Returns the square of total number of hops made by the exciton

void Exciton::resetExciton(monteCarlo::Lattice &L,
                           std::default_random_engine &generator) 
/**Doc String:
*This function resets the exciton to its original position on the lattice
*Input:
*	Lattice
*	Random number generator
*/						   
{
  dwellTimes.clear();
  LatticeEnergies.clear();
  ExcitonRoute.clear();
  L.setSiteOccupation(index, excitonStream, -1);
  index = startIndex;
  L.setSiteOccupation(index, excitonStream, 1);
  std::tie(xPos, yPos, zPos) = L.openPosition(index);
  ExcitonRoute.push_back(index);
  LatticeEnergies.push_back(L.getSiteEnergy(index));
  decayed = false;
  stationaryTime = 0;
  double totalRate = L.getHoppingRate(L.flattenPosition(xPos, yPos, zPos),0);
  dwellTimes.push_back(calculateDwellTime(totalRate, generator));
  xhops = 0;
  yhops = 0;
  zhops = 0;
}

void Exciton::EEA(monteCarlo::Lattice &L, double timeStamp) 
/**Doc String:
*This function sets an exciton as decayed via EEA (type 1, non-radiative) at the passed in timestamp
*Input:
*	Lattice
*	Time of simulation
*/
{
  decayed = true;
  timeOfRecombination = timeStamp;
  typeOfRecombination = 1;
  L.setSiteOccupation(L.flattenPosition(xPos, yPos, zPos), excitonStream, -1);
}

void Exciton::Extraction(monteCarlo::Lattice &L, double timeStamp) 
/**Doc String:
*This function sets an exciton as decayed via extraction (type 2, non-radiative) at the passed in timestamp
*Input:
*	Lattice
*	Time of simulation
*/
{
  decayed = true;
  timeOfRecombination = timeStamp;
  typeOfRecombination = 2;
  L.setSiteOccupation(L.flattenPosition(xPos, yPos, zPos), excitonStream, -1);
}

void Exciton::setDecayed(bool set) 
/**Doc String:
*This function sets the exciton state to decayed and if it did not decay in during the simulation it notes this.
*Input:
*	The state to set the exciton to true-decayed,false-not decayed.
*/
{
  decayed = set;
  if (not decayed) {
    typeOfRecombination = -1;
    timeOfRecombination = -1;
  }
}

/*Advances Exciton*/
void Exciton::advanceExciton(monteCarlo::Lattice &L, double timeStamp,
                             std::default_random_engine &generator) 
/**Doc String:
*This function advances an exciton based on the possibilties:
*The exciton decays radiatively
*The exciton moves to an adjacent lattice site
*OR
*The exciton dwells at it's current location
*Input:
*	Lattice
*	Time of the simulation
*	Random number generator
*/							 
{
  if (decayed == false) {
    std::uniform_real_distribution<double> distribution(0, 1.0);
    double event = distribution(generator); //(double)rand()/RAND_MAX;
    double decayProb = 1.0 / (L.getEscapeFreq() * (averageLifeTime));
    if (event < decayProb) {
      decayed = true;
      timeOfRecombination = timeStamp;
      typeOfRecombination = 0;
      L.setSiteOccupation(L.flattenPosition(xPos, yPos, zPos), excitonStream,
                          -1);
    } else if (dwellTimes[dwellTimes.size() - 1] > stationaryTime) {
      stationaryTime = 0;
      double totalRate =
      L.getHoppingRate(L.flattenPosition(xPos, yPos, zPos), 0);
      moveExciton(L, generator);
      LatticeEnergies.push_back(L.getSiteEnergy(index));
      ExcitonRoute.push_back(index);
      dwellTimes.push_back(calculateDwellTime(totalRate, generator));
    } else {
      stationaryTime += 1.0 / L.getEscapeFreq();
    }
  }
}
} // namespace monteCarlo
