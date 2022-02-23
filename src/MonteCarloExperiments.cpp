#include "MonteCarloExperiments.h"
#include <sys/stat.h>
#include <sys/types.h>

namespace monteCarloExperiments {

std::string createDirectory()
/**Doc String:
*Creates directory named as 'day-month-hour-minute-second' in the running directory
*Returns:
*	the name of the directory created.
*/
{
  /*Create folder to save data*/
  time_t *rawtime = new time_t;
  struct tm *timeinfo;
  time(rawtime);
  timeinfo = localtime(rawtime);
  char buffer[80];
  strftime(buffer, sizeof(buffer), "%d-%m-%Y-%H-%M-%S", timeinfo);
  std::string folderName(buffer);
  const char *folderPath = folderName.c_str();
  // CreateDirectory(folderPath,NULL);
  mkdir(folderPath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  std::string subFolder = folderName + "/ExcitonRoutes";
  const char *fp = subFolder.c_str();
  mkdir(fp, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  delete rawtime;

  return folderName;
}

void savePLData(std::vector<monteCarlo::Exciton> excitonArray,
                monteCarlo::Lattice &L, std::string folderName, double rho,
                double numberOfExcitons, int OMPthreadID) 
/**Doc String:				
*Saves PL and EEA counts as a function of time. FILE: folderName/PLData_DISORDER_DENSITY.txt
*Parameters: 
*	vector of excitons
*	lattice
*	folder name
*	density of excitation
*	number of excitons used in simulation
*	thread number that simualtion was run on
*
*/			
{
  std::cout << "\n\t\tStarting PL Data Save" << std::endl;
  double Vol = L.getXsize() * L.getYsize() * L.getZsize() * L.getdx() *
               L.getdx() * L.getdx();
  std::stringstream strs;
  strs << rho;
  std::string str = strs.str();
  std::stringstream strs1;
  strs1 << L.getDisorder();
  std::string str1 = strs1.str();
  std::string fileName = folderName + "/PLData_" + str1 + "_" + str + ".txt";
  std::cout << "\t\t\tFile Name: " << fileName << std::endl;
  std::ofstream output;
  output.open(fileName);
  output << "Lattice Parameters:\n"
         << "\t xsize: " << L.getXsize() << "\n"
         << "\t ysize: " << L.getYsize() << "\n"
         << "\t zsize: " << L.getZsize() << "\n"
         << "\t dx: " << L.getdx() << "\n"
         << "\t Volume: " << Vol << "\n"
         << "\t Average E: " << L.getAverageE() << "eV\n"
         << "\t Disorder: " << L.getDisorder() << "eV\n"
         << "\t Temperature: " << L.getTemperature() << "K\n"
         << "\t dt: " << 1.0 / L.getEscapeFreq() << " (s)\n"
         << "Simulation Parameters:\n"
         << "\tDiffusion Coefficient: "
         << L.getdx() * L.getdx() * L.getEscapeFreq() / 6 << " (m^2/s)\n"
         << "\tSiglet Lifetime: " << excitonArray[0].getLifeTime() << " (s)\n"
         << "\tDensity: " << rho << " (m^-3)\n"
         << "\tNumber of Exitons: " << excitonArray.size() << "\n"
         << "\tRunning simulation on stream: " << excitonArray[0].getStream()
         << "\n"
         << "\tRunning simulation on OMP Thread: " << OMPthreadID << "\n";
  int PLcount, EEAcount, ExtractedCount, totalCount = 0, numOfEEA = 0,
                                         numOfPL = 0, numOfExtraced = 0;
  bool done = false;
  double dt = 1.0 / L.getEscapeFreq(), time = 0;
  for (int i = 0; i < excitonArray.size(); i++) {
    // std::cout<<i<<"\t"<<excitonArray[i].getTimeOfRecombination()<<"\t"<<excitonArray[i].getTypeOfRecombination()<<std::endl;
    if (excitonArray[i].getTypeOfRecombination() == 1) {
      numOfEEA++;
    } else if (excitonArray[i].getTypeOfRecombination() == 0) {
      numOfPL++;
    } else if (excitonArray[i].getTypeOfRecombination() == 2) {
      numOfExtraced++;
    }
  }
  output << "Number of EEA events: " << numOfEEA << "\n"
         << "Data:\nTime (s)\t PL events (counts)\tEEA events\tExtracted "
            "events\t Density(m^-3)\n";

  while (time < 3 * excitonArray[0].getLifeTime()) { // done == false){
    output << time << "\t";
    PLcount = 0;
    EEAcount = 0;
    ExtractedCount = 0;
    for (int i = 0; i < excitonArray.size(); i++) {
      if ((excitonArray[i].getTimeOfRecombination() >= time) &&
          (excitonArray[i].getTimeOfRecombination() < time + dt)) {
        if (excitonArray[i].getTypeOfRecombination() == 0) {
          PLcount++;
          totalCount++;
        } else if (excitonArray[i].getTypeOfRecombination() == 1) {
          EEAcount++;
          totalCount++;
        } else if (excitonArray[i].getTypeOfRecombination() == 2) {
          ExtractedCount++;
          totalCount++;
        }
      }
    }
    output << PLcount << "\t" << EEAcount << "\t" << ExtractedCount << "\t"
           << (numberOfExcitons - totalCount) / Vol << "\n";
    if (totalCount == excitonArray.size()) {
      done = true;
    }
    time = time + dt;
    // std::cout<<totalCount<<"\t"<<excitonArray.size()<<std::endl;
  }
  output.close();
  std::cout << "\t\tCompleted PL Data Save" << std::endl;
}

void saveExcitonRouteData(std::vector<monteCarlo::Exciton> excitonArray,
                          monteCarlo::Lattice &L, std::string folderName,
                          double rho, double numberOfExcitons,
                          int OMPthreadID) 
/**Doc String:
*Saves length of exciton route, total route length, start and end position, and recombination type of each excitons. FILE: folderName/ExcitonRouteData_DISORDER_DENSITY.txt
*Inputs:
*	vector of excitons
*	lattice
*	folder name to save to
*	density of simulation
*	total number of excitons in simulation
*	thread number the simulation was run on
*/						  
{
  std::cout << "\n\t\tStarting Exciton Route Data Save" << std::endl;
  double Vol = L.getXsize() * L.getYsize() * L.getZsize() * L.getdx() *
               L.getdx() * L.getdx();
  std::stringstream strs;
  strs << rho;
  std::string str = strs.str();
  std::stringstream strs1;
  strs1 << L.getDisorder();
  std::string str1 = strs1.str();
  std::string fileName =
      folderName + "\\ExcitonRouteData_" + str1 + "_" + str + ".txt";
  std::cout << "\t\t	File Name: " << fileName << std::endl;
  std::ofstream output;
  output.open(fileName);
  output << "Lattice Parameters:\n"
         << "\t xsize: " << L.getXsize() << "\n"
         << "\t ysize: " << L.getYsize() << "\n"
         << "\t zsize: " << L.getZsize() << "\n"
         << "\t dx: " << L.getdx() << "\n"
         << "\t Volume: " << Vol << "\n"
         << "\t Average E: " << L.getAverageE() << "eV\n"
         << "\t Disorder: " << L.getDisorder() << "eV\n"
         << "\t Temperature: " << L.getTemperature() << "K\n"
         << "\t dt: " << 1.0 / L.getEscapeFreq() << " (s)\n"
         << "Simulation Parameters:\n"
         << "\tDiffusion Coefficient: "
         << L.getdx() * L.getdx() * L.getEscapeFreq() / 6 << " (m^2/s)\n"
         << "\tDensity: " << rho << " (m^-3)\n"
         << "\tNumber of Exitons: " << excitonArray.size() << "\n"
         << "\tRunning simulation on stream: " << excitonArray[0].getStream()
         << "\n"
         << "\tRunning simulation on OMP Thread: " << OMPthreadID << "\n";

  output << "Exciton\tTotal Hops\tSquared Distance Travelled (au)\t Start X\t "
            "Start Y\t Start Z\tEnd X\tEnd Y\tEnd Z\tExciton Stream\t Type of "
            "Recombination (0-Rad,1-nonRad)\n";
  for (int i = 0; i < excitonArray.size(); i++) {
    auto start = L.openPosition(excitonArray[i].getStartIndex());
    auto end = L.openPosition(excitonArray[i].getIndex());
    output << i << "\t"                                       //
           << excitonArray[i].getExcitonRouteLength() << "\t" //
           << excitonArray[i].getLengthMoved() << "\t"        //
           << std::get<0>(start) << "\t"                      //
           << std::get<1>(start) << "\t"                      //
           << std::get<2>(start) << "\t"                      //
           << std::get<0>(end) << "\t"                        //
           << std::get<1>(end) << "\t"                        //
           << std::get<2>(end) << "\t"                        //
           << excitonArray[i].getStream() << "\t"             //
           << excitonArray[i].getTypeOfRecombination() << "\n";
  }
  output.close();
  std::cout << "\t\tCompleted Exciton Route Data Save" << std::endl;
}

void saveLatticeInfo(monteCarlo::Lattice &L, double disorder,
                     std::string folderName) 
/**Doc String:
*Saves the distribution of energies of lattice points. File: folderName/LatticeSiteEnergy_DISORDER.txt
*Inputs:
*	Lattice
*	disorder (double)
*	folder to save to
*/					 
{
  std::vector<double> data;
  for (int i = 0; i < L.getXsize() * L.getYsize() * L.getZsize(); i++) {
    data.push_back(L.getSiteEnergy(i));
  }

  std::vector<double> range;
  std::vector<double>::iterator result;
  result = std::max_element(data.begin(), data.end());
  int index = std::distance(data.begin(), result);
  double max = data[index];
  result = std::min_element(data.begin(), data.end());
  index = std::distance(data.begin(), result);
  double min = data[index];
  double delta = disorder / 10;
  double BN = (max - min) / (delta) + 1;
  int binNum = static_cast<int>(BN);
  for (int i = 0; i < binNum + 1; i++) {
    range.push_back(min + i * delta);
  }

  std::vector<int> binValues;
  if (range.size() == 0) {
    binValues.push_back(data.size());
    range.push_back(0);
    binNum = 1;
  } else {
    for (int j = 0; j < range.size() + 1; j++) {
      binValues.push_back(0);
    }
    int j;
    for (int i = 0; i < data.size(); i++) {
      j = 0;
      while (data[i] > range[j]) {
        j++;
      }

      binValues[j]++;
    }
  }
  std::stringstream strs;
  double dis = L.getDisorder();
  strs << dis;
  std::string str = strs.str();
  std::cout << "\n\t\tStarting Lattice Data Save" << std::endl;
  std::string fileName = folderName + "/LatticeSiteEnergy_" + str + ".txt";
  std::cout << "\t\t	File Name: " << fileName << std::endl;
  std::ofstream output;
  output.open(fileName);
  output << "Lattice Parameters:\n"
         << "\t xsize: " << L.getXsize() << "\n"
         << "\t ysize: " << L.getYsize() << "\n"
         << "\t zsize: " << L.getZsize() << "\n"
         << "\t dx: " << L.getdx() << " (m)\n"
         << "\t Volume: "
         << L.getXsize() * L.getXsize() * L.getZsize() * L.getdx() * L.getdx() *
                L.getdx()
         << " (m^3)\n"
         << "\t Average E: " << L.getAverageE() << " (eV)\n"
         << "\t Disorder: " << L.getDisorder() << " (eV)\n"
         << "\t Temperature: " << L.getTemperature() << " (K)\n"
         << "\t dt: " << 1.0 / L.getEscapeFreq() << " (s)\n"
         << "Analysis Parameters:\n"
         << "\tMin Energy: " << min << " (meV)\n"
         << "\tMax Energy: " << max << " (meV)\n"
         << "\tNumber of Bins: " << binNum << "\n"
         << "\tStep Size: " << delta << " \n"

         << "Energy Value (eV)\tNumber of Sites\n";
  for (int i = 0; i < binNum; i++) {
    output << range[i] << "\t" << binValues[i] << std::endl;
  }
  output.close();
  std::cout << "\t\tCompleted Lattice Data Save\n" << std::endl;
}

void saveDomainSizeExp_PLQYvsRho(double rho[], double PLQY[], double EEA[],
                                 double Ext[], double numOfExcitonsVector[],
                                 int OMPThreads[], double domainSize,
                                 int zSize,double injectionBarrier,
                                 std::string folderName, double dx,
                                 double disorder, double temp, double averages,
                                 double rhoStart, double rhoEnd, double rhoStep,
                                 double lifetime, double dt) 
/**Doc String:
*Saves the PLQY as a function of density. FILE:folderName/RhoVsPLQY_DOMAINSIZE_INJECTIONBARRRIER.txt
*Inputs:
*	array containing densities
*	array containing PLQY 
*	array containing number of EEA events as %
*	array conatining number of extraction events as %
*	array containing the total number of excitons for each density
*	array containing the thread number for each density
*	domain size used in simulation
*	injection barrier used for re-injection of trapped excitons
*	folder name to save to
*	spacial stepsize
*	disorder
*	temperature of lattice
*	number of averages used in simulation
*	starting density of the simulation
*	ending density of the simulation
*	lifetime of singlet excitons
*	temporal stepsize
*/								 
{

  bool done = false;
  int array_length = 0;
  do {
    if (rho[array_length] < 1) {
      done = true;
    }
    array_length++;
  } while (done == false);
  array_length--;

std::cout<<"SAVE: "<<domainSize*dx<<std::endl;

  std::stringstream strs;
  strs << domainSize*dx;
  std::string str = strs.str();
  std::stringstream strs2;
  strs2 << injectionBarrier;
  std::string str2 = strs2.str();
  std::string fileName = folderName + "/RhoVsPLQY_" + str + "_" + str2 + ".txt";
  std::ofstream output;
  output.open(fileName);
  output << "Lattice Parameters:\n"
         << "\t xsize: " << domainSize*dx*1e9 << " nm\n"
         << "\t ysize: " << domainSize*dx*1e9 << " nm\n"
         << "\t zsize: " << zSize*dx*1e9 << " nm\n"
         << "\t dx: " << dx << "\n"
         << "\t Volume: " << domainSize * domainSize * zSize * dx * dx * dx
         << "\n"
         << "\t Average E: " << 0 << " eV\n"
         << "\t Injection Barrier: " << injectionBarrier << " eV\n"
         << "\t Disorder: " << disorder << "eV\n"
         << "\t Temperature: " << temp << "K\n"
         << "\t dt: " << dt << " (s)\n"
         << "Simulation Parameters:\n"
         << "\tDiffusion Coefficient: " << dx * dx / (6 * dt) << " (m^2/s)\n"
         << "\tSiglet Lifetime: " << lifetime << " (s)\n"
         << "\tExciton Diffusion Length (3D): "
         << sqrt(6 * dx * dx / (6 * dt) * lifetime) << " (s)\n"
         << "\tStart Density: " << rhoStart << " (m^-3)\n"
         << "\tEnd Density: " << rhoEnd << " (m^-3)\n"
         << "\tDensity Step: " << rhoStep << " (m^-3)\n"
         << "\tNumber of Averages: " << averages << "\n"
         << "Data\nDensity (m^-3)\tPLQY (au)\tEEA events (au)\tExtraction "
            "Events (au)\tNumber of Excitons\tNumber of Averages\tOMP thread "
            "number of simulation\n";
  for (int i = 0; i < array_length; i++) {
    output << rho[i] << "\t" << PLQY[i] << "\t" << EEA[i] << "\t" << Ext[i]
           << "\t" << floor(numOfExcitonsVector[i]) << "\t" << OMPThreads[i]
           << "\n";
  }
  output.close();
}

void saveExcitonDiffusionData_3D(std::vector<double> disorder, double temp,
                                 std::vector<double> Ld,
                                 std::string folderName) 
/**Doc String:
*Saves the diffusion length as a function of disorder. FILE: folderName/LdvsDisorder_TEMPERATURE.txt.
*Input:
*	vector containing disorders in simulation
*	temperature of the lattice
*	vector containing diffusion lengths
*	folder name to save to
*/								 
{
  std::stringstream strs;
  strs << temp;
  std::string str = strs.str();
  std::string fileName = folderName + "/LdvsDisorder_" + str + ".txt";
  std::ofstream output;
  output.open(fileName);

  for (int i = 0; i < disorder.size(); i++) {
    output << disorder[i] << "\t" << Ld[i] << "\n";
  }
  output.close();
}

void EEA(monteCarlo::Exciton &Exciton, monteCarlo::Lattice &L,
         double timeStamp) 
/**Doc String:
*This function takes in a exciton and the lattice it occupies and checks if the posision or it's nearest neighbours contains additional excitons. 
*If either condition is true, the exciton is decayed non-radiatively through EEA
*Inputs:
*	Exciton to check and possibly decay
*	Lattice of the simulation
*	time step of the simulation.
*/		 
{
  std::vector<int> neighbours = L.getNeighbours_legacy(L.flattenPosition(
      Exciton.getXpos(), Exciton.getYpos(), Exciton.getZpos()));
  int count = 0;
  bool done = false;
  if (L.getSiteOccupation(Exciton.getIndex(), Exciton.getStream()) > 1) {
    done = true;
    Exciton.EEA(L, timeStamp);
  }
  while ((count < neighbours.size()) && (done == false)) {
    if ((L.getSiteOccupation(neighbours[count], Exciton.getStream()) > 0) &&
        (neighbours[count] != Exciton.getIndex())) {
      done = true;
      Exciton.EEA(L, timeStamp);
    }
    count++;
  }
}

void Extraction(monteCarlo::Exciton &Exciton, monteCarlo::Lattice &L,
                double timeStamp) 
/**Doc String:
*This function takes in an exciton and checks if it is outside the defined domain size.
*If this is true, the exciton is counted as extracted.
*Inputs:
*	Exciton to check
*	Lattice of simulation
*	time step of the simulation
*/				
{
  int centerIndex =
      L.flattenPosition(floor(L.getXsize() / 2), floor(L.getYsize() / 2), 0);
  double ExcitonIndex =
      L.flattenPosition(Exciton.getXpos(), Exciton.getYpos(), 0);
  double distance = L.calculateDistance(centerIndex, ExcitonIndex);
  double radius = double(L.getXsize()) / 2.0;
  if (distance >= radius) { //(L.onBorder(Exciton.getIndex())!=0){
    Exciton.Extraction(L, timeStamp);
  }
}

void Re_injection(monteCarlo::Exciton &Exciton, monteCarlo::Lattice &L,
                 double injectionBarrier,
                 std::default_random_engine &generator) 
/**Doc String:
*This function takes in an exciton and uses a Miller-Abrahams rate to see the possiblity of re-injection into the lattice. 
*Inputs:
*	Exciton to test for re-injection
*	Lattice used in simulations
*	Injection barrier from simulation
*	random number generator
*/				 
{
  std::uniform_real_distribution<double> distribution(0, 1.0);
  double extractionRate = 1.0;
  double injectionRate =
      1.0 * exp(-injectionBarrier / (Kb * L.getTemperature()));
  double r = distribution(generator); //(double)rand()/RAND_MAX;
  double prob = injectionRate / (extractionRate + injectionRate);
  double radius = double(L.getXsize()) / 2.0, distance;
  int centerIndex =
      L.flattenPosition(floor(L.getXsize() / 2), floor(L.getYsize() / 2), 0);
  double ExcitonIndex;
  if (r < prob) {
    do {
      Exciton.moveExciton(L, generator);
      ExcitonIndex = L.flattenPosition(Exciton.getXpos(), Exciton.getYpos(), 0);
      distance = L.calculateDistance(centerIndex, ExcitonIndex);
    } while (distance > radius);
    Exciton.setDecayed(false);
  }
}

double single_exciton_diffusion_3D(int numberOfTrials, double D,
                                   std::string folderName,
                                   monteCarlo::Lattice &L) 
/**Doc String:
*This function calculates the diffusion length of an ensemble of size 'numberOfTrials'
*Each exciton starts from the same location on the lattice. 
*Inputs:
*	number of trials to average over
*	Not needed
*	folder to save data to
*	Lattice
*Returns:
*	The diffusion length in nm
*/
{
  saveLatticeInfo(L, L.getDisorder(), folderName);
  double averageLifeTime = 300e-12, //
      squaredDisTravelled = 0,      //
      dt = 1.0 / L.getEscapeFreq(), //
      LT = 0,                       //
      sum = 0;
  std::default_random_engine generator;
  monteCarlo::Exciton E =
      monteCarlo::Exciton(L, averageLifeTime, 0, false, generator);
  int timeIndex;
  for (int i = 0; i < numberOfTrials; i++) {
    timeIndex = 0;
    E.resetExciton(L, generator);
    if (i == 0) {
      std::cout << "Trial: " << i << " xpos: " << E.getXpos()
                << " ypos: " << E.getYpos() << " zPos: " << E.getZpos()
                << std::endl;
    }
    while (E.getDecayed() == false) {
      E.advanceExciton(L, (double)timeIndex / L.getEscapeFreq(), generator);
      timeIndex++;
    }
    // E.saveExcitonRoute(folderName,L,i);
    squaredDisTravelled += E.getLengthMoved();
    LT += E.getLifeTime();
  }
  squaredDisTravelled /= numberOfTrials;
  double LD = sqrt(squaredDisTravelled) * L.getdx();
  double averageTime = E.getLifeTime();
  return LD;
}

int findXYSize(double xSize, double zSize, double dx, double density,
               double upperlim, double lowerlim) 
/**Doc String:
*This function calculates the appropriate X and Y sizes of a square lattice for a given z-dimension to accomodate a particular denisty.
*Input:
*	starting x-size
*	fixed z-size
*	spacial step size
*	input density
*	upper limit to number of excitons used
*	lower limit to the number of excitons
*Returns:
*	(int) the size to use for x and y-dimension (arbitrary units)
*/
{
  double vol = xSize * xSize * zSize * dx * dx * dx;
  double num = (vol * density);
  while ((num > upperlim) || (num < lowerlim)) {
    if (num == 0) { /*std::cout<<"ERROR: num==0: "<<(num==0)<<std::endl;*/
    }
    if (num > upperlim) {
      xSize--;
      vol = xSize * xSize * zSize * dx * dx * dx;
      num = ((vol * density));
      // std::cout<<"Too big"<<num<<"\t"<<xSize<<"\t"<<density<<std::endl;
    } else if (num < lowerlim) {
      xSize++;
      vol = xSize * xSize * zSize * dx * dx * dx;
      num = ((vol * density));
      // std::cout<<"Too small"<<num<<"\t"<<xSize<<"\t"<<density<<std::endl;
    }
  }
  return xSize;
}

int findZSize_Cylindical(double xSize, double zSize, double dx, double density,
                         double upperlim, double lowerlim) 
/**Doc String:
*This function finds the appropriate z-dimension for a given denisty and x and y-dimension in a cylindrical geometry.
*Inputs:
*	Fixed x-dimension
*	starting z-dimension
*	spacial step size
*	density in simulation
*	upper limit to total number of excitons
*	lower limit to total number of excitons.
*Returns:
*	(int) z-dimension
*/
{
  double vol = zSize * M_PI * xSize * xSize * dx * dx * dx;
  double num = (vol * density);
  while ((num > upperlim) || (num < lowerlim)) {
    if (num == 0) { /*std::cout<<"ERROR: num==0: "<<(num==0)<<std::endl;*/
    }
    if (num > upperlim) {
      zSize--;
      vol = zSize * M_PI * xSize * xSize * dx * dx * dx;
      num = ((vol * density));
      // std::cout<<"Too big"<<num<<"\t"<<xSize<<"\t"<<density<<std::endl;
    } else if (num < lowerlim) {
      zSize++;
      vol = zSize * M_PI * xSize * xSize * dx * dx * dx;
      num = ((vol * density));
      // std::cout<<"Too small"<<num<<"\t"<<xSize<<"\t"<<density<<std::endl;
    }
  }
  return zSize;
}

void setup_multiple_Exciton_PL_3D(double averageLifeTime, double diffusionCoeff,
                                  double rhoStart, double rhoStop,
                                  double rhoStep, double disorder, double dt,
                                  int xBC, int yBC, int zBC, double upperlim,
                                  double lowerlim, double zDimension, int nProcessors,
                                  std::string folderName) 
/**Doc String:
*This function sets up the required lattice and parallelization to complete the EEA simulations on a semi-infinite square lattice.
*Input:
*	singlet exciton lifetime
*	diffusion coefficient of excitons
*	starting density of simulation
*	ending density of simulation
*	step size of densities
*	disorder
*	temporal step size
*	X/Y/Z boundary conditions
*	upper limit of excitons on lattice
*	lower limit of excitons on lattice
*	folder name for saveing data in
*/								  
{
	
  int xsize = 10;
  double numberOfExcitons_previous, numberOfExcitons;
  double dx = sqrt(6 * diffusionCoeff * dt);
  int zsize = int(floor(zDimension / dx));

  int totalNumofStreams = (int)(std::log10(rhoStop) - std::log10(rhoStart)) / rhoStep;
  // initalize lattice to be right size for density
  bool sameLattice = false;
  xsize = findXYSize(xsize, zsize, dx, rhoStart, upperlim, lowerlim);
  std::cout << "Creating first Lattice, xsize: " << xsize << " zSize: " << zsize
            << " Num: " << rhoStart * xsize * xsize * zsize * dx * dx * dx
            << ", disorder: " << disorder
            << " totalNumofStreams: " << totalNumofStreams << std::endl;

  monteCarlo::Lattice *L =
      new monteCarlo::Lattice(xsize, xsize, zsize, dx, disorder, 0, 300, 1 / dt,
                              totalNumofStreams, xBC, yBC, zBC, folderName);

  std::cout << "Lattice Created" << std::endl;
  // for (double rho = rhoStart;rho<=rhoStop;rho*=std::pow(10,rhoStep)){

  if (nProcessors ==0){nProcessors = omp_get_max_threads();}//set processors number to max if input is 0
  std::cout << "Running with: " << nProcessors << " processors." << std::endl;
  omp_set_dynamic(0);
  omp_set_num_threads(nProcessors);
  saveLatticeInfo(*L, disorder, folderName);
#pragma omp parallel for
  for (int stream = 0; stream <= totalNumofStreams; stream++) {
    std::default_random_engine generator; //random number generator to pass to each stream of the simulation
    int OMPthreadID = omp_get_thread_num();
    double rho = rhoStart * std::pow(10, rhoStep * stream);
    usleep(50000 * OMPthreadID);
    std::cout << "\tDensity: " << rho << " Disorder: " << disorder
              << " Stream: " << stream
              << " Num: " << rho * xsize * xsize * zsize * dx * dx * dx
              << " Thread: " << OMPthreadID << std::endl;

    // run diffusion with EEA for given density and lattice.
    run_multiple_Exciton_PL_3D(rho, averageLifeTime, 0, folderName, *L, stream,
                               OMPthreadID, false, generator);
    numberOfExcitons_previous = rho * xsize * xsize * zsize * dx * dx * dx;
    numberOfExcitons =
        (rho * std::pow(10, rhoStep)) * xsize * xsize * zsize * dx * dx * dx;
  }
  delete L; // delete lattice from heap
}

std::tuple<double, double, double> run_multiple_Exciton_PL_3D(
    double rho, double averageLifeTime, double injectionBarrierForDomainSizeExp,
    std::string folderName, monteCarlo::Lattice &L, int stream, int OMPthreadID,
    bool ExtractionCheck, std::default_random_engine &generator) 
/**Doc String:
*This function runs a simulation of multiple excitons on a lattice including (if properly informed) EEA and extraction
*Input:
*	Density of the simulation
*	injection barrier (in case of extraction)
*	folder name to save data to 
*	Lattice
*	lattice stream to run the simulation on 
*	OMP parallelization thread number
*	boolean to indicate if extraction happens (true-yes extraction)
*	random number generator
*Returns:
*	(tuple<double,double,double>) with:
*				Total number of PL events
*				Total number of EEA events
*				Total number of extraction events
*/	
{

  // create array of excitons on the grid with deisty given
  std::vector<monteCarlo::Exciton> excitonArray;
  double vol, numberOfExcitons, numberOfTimeSteps, timeIndex,
      numberOfExcitons_previous = 0;
  if (ExtractionCheck == true) {
    double radius = (L.getXsize() * L.getdx() / 2.0);
    vol = M_PI * radius * radius * L.getZsize() * L.getdx();
  } else {
    vol = (double)L.getXsize() * (double)L.getYsize() * (double)L.getZsize() *
          L.getdx() * L.getdx() * L.getdx();
  }
  int decayedCheck, check, excitonCount;
  bool liveExcitons = true, EEAcheck = true;
  numberOfExcitons = rho * vol;
  for (int i = numberOfExcitons_previous; i < numberOfExcitons - 1; i++) {
    excitonArray.push_back(monteCarlo::Exciton(L, averageLifeTime, stream,
                                               ExtractionCheck, generator));
  }
  std::cout << "\t" << rho << " Excitons Initialized: " << excitonArray.size()
            << " on thread " << OMPthreadID << std::endl;
  decayedCheck = 0;
  timeIndex = 0;
  while (liveExcitons == true) {
    for (int i = 0; i < excitonArray.size(); i++) {
      if (excitonArray[i].getDecayed() == false) {
        if (EEAcheck == true) {
          EEA(excitonArray[i], L, (double)timeIndex / L.getEscapeFreq());
        }
        if (ExtractionCheck == true) {
          Extraction(excitonArray[i], L, (double)timeIndex / L.getEscapeFreq());
        }
        if (timeIndex != 0) {
          excitonArray[i].advanceExciton(
              L, (double)timeIndex / L.getEscapeFreq(), generator);
        }
      } else if ((ExtractionCheck == true) and
                 (excitonArray[i].getTypeOfRecombination() == 2)) {
        Re_injection(excitonArray[i], L, injectionBarrierForDomainSizeExp,
                     generator);
      }
    }
    while ((excitonArray[decayedCheck].getDecayed() == true) &&
           (decayedCheck < excitonArray.size())) {
      decayedCheck++;
    }
    if (decayedCheck == excitonArray.size()) {
      liveExcitons = false;
    }

    timeIndex++;
  }
  numberOfExcitons_previous = numberOfExcitons;

  double x, xsquared = 0, sum = 0;
  for (int i = 0; i < excitonArray.size(); i++) {
    x = L.calculateDistance(excitonArray[i].getStartIndex(),
                            excitonArray[i].getIndex());
    xsquared += x * x;
    sum += excitonArray[i].getLengthMoved();
  }
  xsquared /= excitonArray.size();
  sum /= excitonArray.size();

  if (ExtractionCheck == true) {
    double PLcounts = 0, EEAcounts = 0, Extcounts = 0, noDecay = 0;
    double dt = 1 / L.getEscapeFreq();
    for (int i = 0; i < excitonArray.size(); i++) {

      if ((excitonArray[i].getTimeOfRecombination() !=
           -1)) { // and(excitonArray[i].getTypeOfRecombination()==2)){}
        if (excitonArray[i].getTypeOfRecombination() == 0) {
          PLcounts++;
        } else if (excitonArray[i].getTypeOfRecombination() == 1) {
          EEAcounts++;
        } else if (excitonArray[i].getTypeOfRecombination() == 2) {
          Extcounts++;
        } else if (excitonArray[i].getTypeOfRecombination() == -1) {
          noDecay++;
        }
      }
    }
    L.resetOccupation(stream);
    return std::make_tuple(PLcounts, EEAcounts, Extcounts);
  } else {
    savePLData(excitonArray, L, folderName, rho, numberOfExcitons, OMPthreadID);
    saveExcitonRouteData(excitonArray, L, folderName, rho, numberOfExcitons,
                         OMPthreadID);
    return {0, 0, 0};
  }
}

void test(monteCarlo::Lattice &L, std::string folderName) 
/**Doc String:
*This is a test function to test various functions of the lattice and exciton interactions
*/
{
  std::default_random_engine generator;
  monteCarlo::Exciton E = monteCarlo::Exciton(L, 300e-12, 0, false, generator);

  std::cout << "Start Position:\t" << E.getIndex() << "\t" << E.getXpos()
            << "\t" << E.getYpos() << "\t" << E.getZpos() << std::endl;
  for (int i = 0; i < 6; i++) {
    E.advanceExciton(L, 0, generator);

    std::cout << "End Position:\t" << E.getIndex() << "\t" << E.getXpos()
              << "\t" << E.getYpos() << "\t" << E.getZpos() << std::endl;
    E.construct(L, 300e-12, 0, true, generator);
  }

  std::cout << L.getXBCs() << "\t" << L.getYBCs() << "\t" << L.getZBCs() << "\t"
            << std::endl;

  int xPos = 2, yPos = 2, zPos = 4;
  std::cout << "next index: " << L.flattenPosition(xPos, yPos, zPos)
            << std::endl;
}

void measureDiffusionLength(int numberOfTrials, monteCarlo::Lattice &L,
                            std::string folderName) 
/**Doc String:
*This function measures the diffusion length of an ensemble of excitons on the passed lattice.
*Inputs:
*	number of trials to test over
*	Lattice
*	folder name to save data to
*/						
{
  std::vector<double> DiffusionCoeff;
  int numberOfStarts = 100;
  double LT;
  std::default_random_engine generator;
  for (int j = 0; j < numberOfStarts; j++) {
    double averageLifeTime = 300e-12, x, xSquared = 0, xave = 0,
           dt = 1.0 / L.getEscapeFreq();
    LT = 0;
    monteCarlo::Exciton E =
        monteCarlo::Exciton(L, 300e-12, 0, false, generator);
    int timeIndex;
    for (int i = 0; i < numberOfTrials; i++) {
      timeIndex = 0;
      E.resetExciton(L, generator);
      while ((E.getXpos() < 20) || (E.getXpos() > L.getXsize() - 20) ||
             (E.getYpos() < 20) || (E.getYpos() > L.getYsize() - 20) ||
             (E.getZpos() < 20) || (E.getZpos() > L.getZsize() - 20)) {
        E.construct(L, averageLifeTime, 0, true, generator);
      }
      if (i == 0) {
        std::cout << "Start: " << j << " xpos: " << E.getXpos()
                  << " ypos: " << E.getYpos() << " zPos: " << E.getZpos()
                  << std::endl;
      }

      while (E.getDecayed() == false) {
        E.advanceExciton(L, (double)timeIndex / L.getEscapeFreq(), generator);
        timeIndex++;
      }
      x = L.calculateDistance(E.getStartIndex(), E.getIndex());
      xSquared += x * x;
      xave += x;
      LT += E.getTimeOfRecombination();
      // std::cout<<"\ti:
      // "<<i<<"\t"<<timeIndex<<"\t"<<timeIndex*dt<<"\t"<<E.getTimeOfRecombination()<<"\t"<<LT<<std::endl;
    }
    // std::cout<<"done "<<timeIndex<<"\t";
    std::cout << numberOfTrials << "\t" << E.getTimeOfRecombination() << "\t"
              << LT << "\t" << LT / numberOfTrials << std::endl;
    // std::cout<<"dx: "<<x*L.getdx()<<std::endl;
    E.saveExcitonRoute(folderName, L, 0);
    xSquared /= numberOfTrials;
    xave /= numberOfTrials;
    double Ld = sqrt(xSquared);
    double averageTime = LT / ((double)numberOfTrials);
    // std::cout<<"\t"<<xave<<"\t"<<L.getdx()<<"\t"<<numberOfTrials<<std::endl;
    std::cout << "Ld: " << Ld << std::endl;
    // DiffusionCoeff.push_back(Diff/(2*averageTime));
  }
  saveDiffusionCoeffData(DiffusionCoeff, L, folderName, numberOfTrials,
                         numberOfStarts);
}

void saveDiffusionCoeffData(std::vector<double> DiffusionCoeff,
                            monteCarlo::Lattice &L, std::string folderName,
                            int numberOfTrials, int numberOfStarts) 
/**Doc String:
*This function takes in data from measureDiffusionLength function and saves it to FILE: folderName/DiffusionCoeffData_DISORDER.txt
*Inputs:
*	vector containing diffusion coefficients from simulation
*	Lattice
*	folder name to save data to
*	number of trials to average over
*	numer of places to start the exciton at
*/							
{

  std::stringstream strs;
  strs << L.getDisorder();
  std::string str = strs.str();
  std::string fileName = folderName + "\\DiffusionCoeffData_" + str + ".txt";
  std::ofstream output;
  output.open(fileName);
  output << "Lattice Parameters:\n"
         << "\t xsize: " << L.getXsize() << "\n"
         << "\t ysize: " << L.getYsize() << "\n"
         << "\t Average E: " << L.getAverageE() << "eV\n"
         << "\t Disorder: " << L.getDisorder() << "eV\n"
         << "\t Temperature: " << L.getTemperature() << "K\n"
         << "\t dt: " << 1.0 / L.getEscapeFreq() << " (s)\n"
         << "Simulation Parameters:\n"
         << "\tNumber of Trials: " << numberOfTrials << "\n"
         << "\tNumber of Starts: " << numberOfStarts << "\n"
         << "Data:\nTrial\tDiffusionCoef\n";
  for (int i = 0; i < DiffusionCoeff.size(); i++) {
    output << i << "\t" << DiffusionCoeff[i] << "\n";
  }
  output.close();
}

int domainSizeStoppingCondition(double PLQY, double PLQYCuttOff) 
/**Doc String:
*This function determines the stopping condition for the parallelized while loop in WHILEdomainSizeExp
*Inputs:
*	The PLQY value in the simulation
*	The defined PLQY cuttoff value
*Returns:
*	The result of the condition
*/
{
  int condition;
  if (PLQY < PLQYCuttOff) {
    condition = 1;
  } else {
    condition = 0;
  }
  return condition;
}

void WHILEdomainSizeExp(int domainSize, double densityStart, double rhoStep,
                        double dx, double disorder, double dt, int zBC,
                        double lowerExcitonLim, double injectionBarrier,
                        double PLQYCuttOff, std::string folderName) 
/**Doc String:
*This program runs the run_multiple_Exciton_PL_3D function in a parallelized while loop. This function is uncompleted.
*/						
{
  double upperExcitonLim = 12000;
  int xBC = 2, yBC = 2, totalNumofStreams = 4;
  std::cout << "PLQY Cuttoff: " << PLQYCuttOff << std::endl;

  int zsize = findZSize_Cylindical(domainSize / 2.0, 1, dx, densityStart,
                                   upperExcitonLim, lowerExcitonLim);
  double num = densityStart * zsize * M_PI * domainSize / 2.0 * domainSize /
               2.0 * dx * dx * dx;

  monteCarlo::Lattice *L = new monteCarlo::Lattice(
      domainSize, domainSize, zsize, dx, disorder, 0, 300, 1 / dt,
      totalNumofStreams, xBC, yBC, zBC, folderName);

  double rho, tempPL, tempEEA, tempExt, tau = 300e-12;
  std::vector<double> PLQYVector, EEAVector, EXTVector, rhoVector,
      numOfExcitonsVector;
  double maxArrayLength = 6 / rhoStep;
  double rhoArray[(int)ceil(maxArrayLength)],
      PLQYArray[(int)ceil(maxArrayLength)], EEAArray[(int)ceil(maxArrayLength)],
      EXTArray[(int)ceil(maxArrayLength)],
      numOfExcitonsArray[(int)ceil(maxArrayLength)];
  int OMPthreadArray[(int)ceil(maxArrayLength)];

  /*Create Paralell openMP While loop over densities to end on condition that
   * the PLQY is blow the threshold specified*/
  /*Shared and private variables*/
  int i_shared, stop, condition_shared, thread_number, i_private,
      condition_private;
  omp_set_num_threads(totalNumofStreams);

  i_shared = -1;
  condition_shared = 0;
#pragma omp parallel private(thread_number, i_private, condition_private, rho)
  {
    thread_number = omp_get_thread_num();
    while (!condition_shared) {
#pragma omp critical
      {
        i_shared++;
        i_private = i_shared;
        rho = densityStart * std::pow(10, rhoStep * i_private);
      }
      std::default_random_engine generator;
      usleep(10000 * thread_number * thread_number * thread_number);
      num = floor(rho * zsize * M_PI * domainSize / 2.0 * domainSize / 2.0 *
                  dx * dx * dx);
      std::cout << "Start " << rho << " " << domainSize << " " << zsize << " "
                << num << std::endl;

      std::tie(tempPL, tempEEA, tempExt) = run_multiple_Exciton_PL_3D(
          rho, tau, injectionBarrier, folderName, *L, thread_number,
          thread_number, true, generator);
      PLQYArray[i_private] = tempPL / num;
      EEAArray[i_private] = tempEEA / num;
      EXTArray[i_private] = tempExt / num;
      numOfExcitonsArray[i_private] = num;
      rhoArray[i_private] = rho;

      std::cout << "\t" << rho << ": " << tempPL << " " << tempEEA << " "
                << tempExt << " "
                << " " << num << " " << tempPL / num << " " << tempExt / num
                << std::endl;
      condition_private = domainSizeStoppingCondition(rho, 1e26);
      if (condition_private == 0) {
        condition_shared = 1;
      }
    }
  }

  saveDomainSizeExp_PLQYvsRho(rhoArray, PLQYArray, EEAArray, EXTArray,
                              numOfExcitonsArray, OMPthreadArray, zsize, domainSize,
                              injectionBarrier, folderName, dx, disorder, 300,
                              1, densityStart, rho, rhoStep, tau, dt);
}

void domainSizeExp(int domainSize, double densityStart, double densityStop,
                   double rhoStep, double dx, double disorder, double dt,double tau,
                   int zBC, double lowerExcitonLim, double injectionBarrier,
                   double nProcessors, std::string folderName) 
/**Doc String:
*This function runs the run_multiple_Exciton_PL_3D function in a parallelized for loop over a set denisty range, in a cylindrical domain.
*Inputs:
*	The diameter of the cylindrical domain.
*	The lowest density in the simulation
*	The highest density in the simulation
*	The step size over densities used in the simulation
*	The spacial step size
*	The disorder
*	The temporal step size
*	The boundary condition on the z-axis
*	The lower exciton limit in the simualtion
*	Injection barrier
*	The PLQY cuttoff defined in the simulation (not used)
*	Folder name for saving data
*/ 
{
std::cout<<"HERE: "<<tau<<std::endl;
  double upperExcitonLim = 12000;
  int xBC = 2, yBC = 2;

  int totalNumofStreams =
      (int)(std::log10(densityStop) - std::log10(densityStart)) / rhoStep;

  int zsize = findZSize_Cylindical(domainSize / 2.0, 1, dx, densityStart,
                                   upperExcitonLim, lowerExcitonLim);

  monteCarlo::Lattice *L = new monteCarlo::Lattice(
      domainSize, domainSize, zsize, dx, disorder, 0, 300, 1 / dt,
      totalNumofStreams, xBC, yBC, zBC, folderName);
  std::cout << "Lattice Created, Streams: " << totalNumofStreams + 1
            << std::endl;
  if (nProcessors ==0){
	nProcessors = omp_get_max_threads();
  }
  double rho, tempPL, tempEEA, tempExt;

  double rhoArray[totalNumofStreams], PLQYArray[totalNumofStreams],
      EEAArray[totalNumofStreams], EXTArray[totalNumofStreams],
      numOfExcitonsArray[totalNumofStreams];
  int OMPthreadArray[totalNumofStreams];
  omp_set_dynamic(0);
  omp_set_num_threads(nProcessors);
  saveLatticeInfo(*L, disorder, folderName);

  std::cout << "Processors: " << nProcessors << std::endl;
#pragma omp parallel for
  for (int stream = 0; stream <= totalNumofStreams; stream++) {
    std::default_random_engine generator;
    int OMPthreadID = omp_get_thread_num();
    double rho = densityStart * std::pow(10, rhoStep * stream);
    usleep(50000 * OMPthreadID);
    double num = floor(rho * zsize * M_PI * domainSize / 2.0 * domainSize /
                       2.0 * dx * dx * dx);
    std::cout << "Start " << stream << " " << OMPthreadID << " " << rho << " "
              << domainSize << " " << zsize << " " << num << std::endl;

    std::tie(tempPL, tempEEA, tempExt) =
        run_multiple_Exciton_PL_3D(rho, tau, injectionBarrier, folderName, *L,
                                   stream, OMPthreadID, true, generator);
    PLQYArray[stream] = tempPL / num;
    EEAArray[stream] = tempEEA / num;
    EXTArray[stream] = tempExt / num;
    numOfExcitonsArray[stream] = num;
    rhoArray[stream] = rho;
    OMPthreadArray[stream] = OMPthreadID;

    std::cout << "End " << stream << " " << rho << ": " << tempPL / num << " "
              << tempExt / num << std::endl;
  }

  saveDomainSizeExp_PLQYvsRho(rhoArray, PLQYArray, EEAArray, EXTArray,
                              numOfExcitonsArray, OMPthreadArray, domainSize, zsize,
                              injectionBarrier, folderName, dx, disorder, 300,
                              1, densityStart, densityStop, rhoStep, tau, dt);
  delete L;
  std::cout << "Lattice Deleted" << std::endl;
}

} // namespace monteCarloExperiments
