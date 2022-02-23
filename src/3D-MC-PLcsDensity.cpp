#include "MonteCarloExperiments.h"
#include "tools.h"
#include <iostream>
#include <experimental/filesystem>

void MultipleExcitonPL(std::string folderName,std::string parametersFilePath, ParameterMap constants) 
/**Doc String:
*This program runs setup_multiple_Exciton_PL_3D given the parameters in the file: parameters/Parameters_MultipleExcitonPL.txt
*Input:
*	folder name to save data to
*	parameters file path to copy to folder name
*	ParameterMap of constants read from parameters file
*/
{

  // Using a lambda for this slightly complicated condition
  auto loop_condition = [&](double disorder, int count) {
    if (constants["disorderStep"] == 0)
      return count == 0;
    else
      return disorder <= constants["DisorderEnd"];
  };

  int count = 0;
  for (double disorder = constants["DisorderStart"]; //
       loop_condition(disorder, count);              //
       disorder *= std::pow(10, constants["DisorderStep"]), count++) {

    std::cout << "Disorder Step: " << constants["DisorderStep"] << std::endl;
	
	const std::experimental::filesystem::path ParametersSource = parametersFilePath;
	const std::experimental::filesystem::path ParametersTarget = folderName + "/ParametersFile.txt";
	std::experimental::filesystem::copy(ParametersSource, ParametersTarget, std::experimental::filesystem::copy_options::overwrite_existing);
	
    monteCarloExperiments::setup_multiple_Exciton_PL_3D( //
        constants["SigletLifetime"],                    //
        constants["DiffusionCoeff"],                     //
        constants["DensityStart"],                       //
        constants["DensityEnd"],                         //
        constants["DensityStep"],                        //
        disorder,                                        //
        constants["dt"],                                 //
        constants["xBC"],                                //
        constants["yBC"],                                //
        constants["zBC"],                                //
        constants["UpperExcitonLimit"],                  //
        constants["LowerExcitonLimit"],                  //
        constants["zDimension"],                         //      
		constants["nProcessors"],						 //
        folderName);
  }
}

int main(int argc, char **argv) {
	std::cout<<argc<<" "<<argv[1]<<std::endl;
  std::string folderName = monteCarloExperiments::createDirectory();
  std::cout<<"Outside"<<std::endl;
  auto constants = read_parameters(argv[1]);
  MultipleExcitonPL(folderName,argv[1], constants);
}
