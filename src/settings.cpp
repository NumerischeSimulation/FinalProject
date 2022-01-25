#include "settings.h"
#include <fstream>   // for file operations
#include <iostream>  // for cout

void Settings::printSettings()
{
  std::cout << "Settings: " << std::endl
    << "  physicalSize: " << physicalSize[0] << " x " << physicalSize[1] << ", nCells: " << nCells[0] << " x " << nCells[1] << std::endl
    << "  endTime: " << endTime << " s, re: " << re << ", g: (" << g[0] << "," << g[1] << "), tau: " << tau << ", maximum dt: " << maximumDt << std::endl
    << "  dirichletBC: bottom: (" << dirichletBcBottom[0] << "," << dirichletBcBottom[1]  << ")"
    << ", top: ("  << dirichletBcTop[0] << "," << dirichletBcTop[1]  << ")"
    << ", left: ("  << dirichletBcLeft[0] << "," << dirichletBcLeft[1] << ")"
    << ", right: ("  << dirichletBcRight[0] << "," << dirichletBcRight[1] << ")" << std::endl
    << "  useDonorCell: " << std::boolalpha << useDonorCell << ", alpha: " << alpha << std::endl
    << "  pressureSolver: " << pressureSolver << ", omega: " << omega << ", epsilon: " << epsilon << ", maximumNumberOfIterations: " << maximumNumberOfIterations << std::endl;
}

void Settings::loadFromFile(std::string filename)
{
  // open file
  std::ifstream file(filename.c_str(), std::ios::in);

  // check if file is open
  if (!file.is_open())
  {
    std::cout << "Could not open parameter file \"" << filename << "\"." << std::endl;
    return;
  }

  // loop over lines of file
  for (int lineNo = 0;; lineNo++)
  {
    // read line
    std::string line;
    getline(file, line);

    // at the end of the file break for loop
    if (file.eof())
      break;

    int first_character_idx = line.find_first_not_of(" \t\n");

    // std::cout << line << ":" << lineNo << std::endl;

    // check whether to ignore the line
    if (line[0] == '#'){ // commented out
      continue;
    } else if (!(first_character_idx != std::string::npos)) { // no content
      continue;
    }

    // find equality
    int eq_idx = line.find('=');
    if (!(eq_idx != std::string::npos)){
      std::cout << "There is a badly defined line\"" << lineNo << " of " << filename << "\"." << std::endl;
      std::cout << line << std::endl;
      return;
    }

    // left side of equality - parameter name
    int last_character_idx = line.find_first_of(" \t", first_character_idx);
    // std::cout << first_character_idx << " " << last_character_idx << " " << eq_idx << std::endl;
    int end_idx = 0;
    if (last_character_idx >= eq_idx){
      end_idx = eq_idx - first_character_idx;
    } else {
      end_idx = last_character_idx - first_character_idx;
    }
    std::string parameterName = line.substr(first_character_idx, end_idx);

    // right side of equality - parameter value
    int first_value_idx = line.find_first_not_of(" \t", eq_idx+1);
    int last_value_idx = line.find_first_of(" \t\n#", first_value_idx);
    std::string valueString = line.substr(first_value_idx, (last_value_idx-first_value_idx));
    // std::cout << first_value_idx << " " << last_value_idx << std::endl;

    // std::cout << parameterName << "=" << valueString << ":" << eq_idx << std::endl;

    // mapping of parameter to read value - all the if's
    // ordering according to struct
    if (parameterName == "nCellsX"){
      nCells[0] = int(atof(valueString.c_str())); // int
    }
    else if (parameterName == "nCellsY"){
      nCells[1] = int(atof(valueString.c_str()));  // int
    } 
    else if (parameterName == "physicalSizeX"){
      physicalSize[0] = atof(valueString.c_str()); // double
    } 
    else if (parameterName == "physicalSizeY"){
      physicalSize[1] = atof(valueString.c_str()); // double
    }
    else if (parameterName == "re"){
      re = atof(valueString.c_str()); // double
    }
    else if (parameterName == "endTime"){
      endTime = atof(valueString.c_str()); // double
    }
    else if (parameterName == "tau"){
      tau = atof(valueString.c_str()); // double
    }
    else if (parameterName == "maximumDt"){
      maximumDt = atof(valueString.c_str()); // double
    }
    else if (parameterName == "gX"){
      g[0] = atof(valueString.c_str()); // double
    }
    else if (parameterName == "gY"){
      g[1] = atof(valueString.c_str()); // double
    }
    else if (parameterName == "useDonorCell"){
      if (valueString == "true"){
        useDonorCell = true;
      } else if (valueString == "false"){
        useDonorCell = false;
      } else {
        std::cout << "The donor cell parameter is not understood: " << parameterName << std::endl;
        return;
      }
      useDonorCell = (valueString == "true"); // bool
    }
    else if (parameterName == "alpha"){
      alpha = atof(valueString.c_str()); // double
    }
    else if (parameterName == "dirichletBottomX"){
      dirichletBcBottom[0] = atof(valueString.c_str()); // double
    }
    else if (parameterName == "dirichletBottomY"){
      dirichletBcBottom[1] = atof(valueString.c_str()); // double
    } 
    else if (parameterName == "dirichletTopX"){
      dirichletBcTop[0] = atof(valueString.c_str()); // double
    }
    else if (parameterName == "dirichletTopY"){
      dirichletBcTop[1] = atof(valueString.c_str()); // double
    }
    else if (parameterName == "dirichletLeftX"){
      dirichletBcLeft[0] = atof(valueString.c_str()); // double
    }
    else if (parameterName == "dirichletLeftY"){
      dirichletBcLeft[1] = atof(valueString.c_str()); // double
    }
    else if (parameterName == "dirichletRightX"){
      dirichletBcRight[0] = atof(valueString.c_str()); // double
    }
    else if (parameterName == "dirichletRightY"){
      dirichletBcRight[1] = atof(valueString.c_str()); // double
    }
    else if (parameterName == "pressureSolver"){
      pressureSolver = valueString; // std:string
    }
    else if (parameterName == "omega"){
      omega = atof(valueString.c_str()); // double
    }
    else if (parameterName == "epsilon"){
      epsilon = atof(valueString.c_str()); // double
    }
    else if (parameterName == "maximumNumberOfIterations"){
      maximumNumberOfIterations = int(atof(valueString.c_str())); // int, but else the scientific notation is not understood
    }
    else {
      std::cout << "The parameter " << parameterName << " in this line " << lineNo << " could not be found in " << filename << "\"." << std::endl;
      std::cout << line << std::endl;
      return;
    }

/*     if (value == 0.0 && (line_right.find("0.0") == std::string::npos)){ // check whether there is truely no int conversion
      // maybe its a boolean
      if (line_right == 'true' || line_right == 'True' || line_right == 'TRUE'){
          bool value = true;
      }
      else if (line_right == 'false' || line_right == 'False' || line_right == 'FALSE'){
          bool value = false;
      }
        
      // then it must be a string
      std::string value = line_right;

      std::cout << "No valid conversation could be performed in" << lineNo << " of " << filename << "\"." << std::endl;
      return;
    } */

/*     // set corresponding parameter
    if (line_left == "nCellsX"){
      // TODO: I would like to check/cast the type
      Settings.nCells[0] = value;
    } */
/*     // print line
    std::cout << "line " << lineNo << ": " << line << std::endl; */
  }
  std::cout << "Reading is complete." << std::endl;
}

