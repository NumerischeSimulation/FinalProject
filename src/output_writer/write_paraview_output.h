//! Write a file out/output_<fileNo>.vti to be visualized in ParaView.
//! It contains 10x10 nodes with an artifical pressure field.
//! This method is only for demonstration purpose and does nothing useful.
//! However, we will provide similar files, e.g. "output_writer_paraview.h", to be used in the submission code.

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include <cstdlib>
#include <iostream>
#include <string>

void writeParaviewOutput(int fileNo, std::string uID);
