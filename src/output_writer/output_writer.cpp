#include "output_writer/output_writer.h"
#include "discretization/1_discretization.h"

#include <iostream>

OutputWriter::OutputWriter(std::shared_ptr<Discretization> discretization, std::string uID)
 : discretization_(discretization), fileNo_(0)
{
  // create "out" subdirectory if it does not yet exist
  std::string command = "mkdir -p ";
  std::string folder_name = "out_"  + uID;
  int returnValue = system((command + folder_name).c_str());
  if (returnValue != 0)
    std::cout << "Could not create subdirectory \"out_"  << uID << "\"." << std::endl;
  uID_ = uID;
}
