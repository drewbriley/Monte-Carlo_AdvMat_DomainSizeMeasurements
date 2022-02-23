#include "tools.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

ParameterMap read_parameters(std::string filePath) {
	
  ParameterMap out;
  std::cout << "Reading file " << filePath << std::endl;
  auto file = std::ifstream(filePath);
  if (not file){
	std::cout<<"No Input File found"<<std::endl;
    throw std::runtime_error("Error opening file!");
  }
  for (std::string line; getline(file, line);) {
    std::string str;
    double value;
    std::istringstream iss(line);
    if (not((iss >> str) && (iss >> value))) {
      std::cout << "Breaking at line: " << line;
      // continuing until the first read failure.
      break;
    }
    out[str] = value;
    std::cout << std::setw(20) << str << ':' << std::setw(10) << value << std::endl;
  }
  return out;
}
