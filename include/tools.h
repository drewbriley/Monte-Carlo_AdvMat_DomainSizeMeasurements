#pragma once
#include <map>
#include <string>
using ParameterMap = std::map<std::string, double>;

ParameterMap read_parameters(std::string filePath);
