#include "MonteCarloExperiments.h"
#include "tools.h"
#include <iostream>
#include <experimental/filesystem>

void DomainSizeExperiment(std::string folderName,std::string parametersFilePath, ParameterMap constants) 
/**Doc String:
*This function runs the domainSizeExp function using the parameters given in the file: parameters/DomainSizeMeasurements.txt
*Input:
*	Folder name for saving data
*	file path to parameters file to copy to data folder
*	ParametersMap of constants found from parameters file
*/
{
  std::string path = folderName + "\\"+parametersFilePath.substr(parametersFilePath.find('/')+1,parametersFilePath.length());//Parameters_DomainSizeMeasurement.txt";
  const char *paramsPath = path.c_str();

  std::cout << constants["domainSizeStart"] << " " //
            << constants["domainSizeStop"] << " "   //
            << constants["domainSizeStep"] << std::endl;

  const std::experimental::filesystem::path ParametersSource = parametersFilePath;
  const std::experimental::filesystem::path ParametersTarget = folderName + "/ParametersFile.txt";
  std::experimental::filesystem::copy(ParametersSource, ParametersTarget, std::experimental::filesystem::copy_options::overwrite_existing);

  for (double domainSize = constants["domainSizeStart"];domainSize >= constants["domainSizeStop"];domainSize -= constants["domainSizeStep"]) {
    std::cout << domainSize << std::endl;
    std::cout << "Starting Domainsize: " << domainSize << std::endl;
	
	
    monteCarloExperiments::domainSizeExp(domainSize,                     //
                                         constants["densityStart"],      //
                                         constants["densityStop"],       //
                                         constants["densityStep"],       //
                                         constants["dx"],                //
                                         constants["disorder"],          //
                                         constants["dt"],                //
                                         constants["SingletLifetime"],	//
					 constants["zBC"],               //
                                         constants["lowerExcitonLimit"], //
                                         constants["injectionBarrier"],  //
                                         constants["NumerofProcessors"],       //
                                         folderName);
  }
}
int main(int argc, char **argv) {
  std::string folderName = monteCarloExperiments::createDirectory();
  auto constants = read_parameters(argv[1]);
  DomainSizeExperiment(folderName,argv[1], constants);
}
