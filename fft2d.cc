// Distributed two-dimensional Discrete FFT transform
// YOUR NAME HERE
// ECE8893 Project 1


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <signal.h>
#include <math.h>
#include <mpi.h>

#include "Complex.h"
#include "InputImage.h"

using namespace std;


void Transform2D(const char* inputFN) 
{ 

  // Do the 2D transform here.
  
  // 1) Use the InputImage object to read in the Tower.txt file and
  //    find the width/height of the input image.
  InputImage image(inputFN);
  int imgHeight, imgWidth;
  Complex *imgData;

  imgHeight = image.GetHeight();
  imgWidth = image.GetWidth();
  imgData = image.GetImageData();

  printf("\nImage height:%d\tWidth:%d\n", imgHeight, imgWidth);

  // 2) Use MPI to find how many CPUs in total, and which one
  //    this process is
  int numtasks, rank;

  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  printf("\nTasks:%d\tRank:%d\n", numtasks, rank);


  // 3) Allocate an array of Complex object of sufficient size to
  //    hold the 2d DFT results (size is width * height)
  // 4) Obtain a pointer to the Complex 1d array of input data
  // 5) Do the individual 1D transforms on the rows assigned to your CPU
  // 6) Send the resultant transformed values to the appropriate
  //    other processors for the next phase.
  // 6a) To send and receive columns, you might need a separate
  //     Complex array of the correct size.
  // 7) Receive messages from other processes to collect your columns
  // 8) When all columns received, do the 1D transforms on the columns
  // 9) Send final answers to CPU 0 (unless you are CPU 0)
  //   9a) If you are CPU 0, collect all values from other processors
  //       and print out with SaveImageData().

  //InputImage image(inputFN);  // Create the helper object for reading the image
  // Step (1) in the comments is the line above.
  // Your code here, steps 2-9
}

void Transform1D(Complex* h, int w, Complex* H)
{
  // Implement a simple 1-d DFT using the double summation equation
  // given in the assignment handout.  h is the time-domain input
  // data, w is the width (N), and H is the output array.
}

int main(int argc, char** argv)
{

  string fn("Tower.txt");   // default file name
  if (argc > 1) {
    fn = string(argv[1]);   // if name specified on cmd line
  }  

  // MPI initialization here
  int rc;

  rc = MPI_Init(&argc, &argv);
  if(rc != MPI_SUCCESS){
    printf("Error starting MPI Program. Terminating\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  Transform2D(fn.c_str()); // Perform the transform.

  // Finalize MPI here
  MPI_Finalize();

}  
  

  
