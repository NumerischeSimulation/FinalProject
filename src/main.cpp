#include <iostream>
#include <cstdlib>

#include "settings.h"
#include "computation.h"

int main(int argc, char *argv[])
{
  // if the number of given command line arguments is only 1 (= the program name), print out usage information and exit
  if (argc == 1)
  {
    std::cout << "usage: " << argv[0] << " <filename>" << std::endl;

    return EXIT_FAILURE;
  }
  std::cout << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "Parsing paramter file..." << std::endl;
  std::cout << std::endl;
  
  // construct computation obj: parses parameter file and prints settings
  Computation computation = Computation();
  computation.initialize(argc, argv);
  
  std::cout << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
  std::cout << "Starting simulation ..."  << std::endl;
  std::cout << std::endl;
  
  // // start simulating
  computation.runSimulation();

  // test some stuff
  //computation.runTest();

  std::cout << "-------------------------------------------------" << std::endl;
  
  return EXIT_SUCCESS;
}
