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

void Transform1D(Complex* h, int w, Complex* H);



//gather:flag=0  scatter:flag=1
void ScatterGather(Complex* H, int numtasks, int rank, int size, int flag)
{
  int rc;
  
  // Reciever at 0
  if(flag == 0){
    //Send the computed values from all the other ranks to rank 0.
    if(rank != 0)
    {
      Complex* sendBuff = H;
      MPI_Request request;
      rc = MPI_Isend(sendBuff, size*2/numtasks, MPI_COMPLEX, 0, 0, MPI_COMM_WORLD, &request);
    }

    //Rank 0 will receive blocks of 1D Transform data from other ranks.
    if(rank == 0)
    {
      for(int i = 1; i<numtasks; i++)
      {
        Complex* recBuff = H+size*i/numtasks;
        MPI_Status status;
        rc = MPI_Recv(recBuff, size*2/numtasks, MPI_COMPLEX, i, 0, MPI_COMM_WORLD, &status);
      }
    }
  }
  
  //Transmitter at 0
  if(flag == 1){
    //Send values to all other ranks from rank0
    if(rank == 0)
    {
      for(int i = 1; i<numtasks; i++)
      {
        Complex* sendBuff = H;  
        MPI_Request request;
        rc = MPI_Isend(sendBuff, size*2, MPI_COMPLEX, i, 0, MPI_COMM_WORLD, &request);
      }
    }
    //All ranks recieve values from rank0
    if(rank != 0)
    {
      Complex* recBuff = H;
      MPI_Status status;
      rc = MPI_Recv(recBuff, size*2, MPI_COMPLEX, 0, 0, MPI_COMM_WORLD, &status);
    }
  }
}


void transpose(Complex* in, Complex* out, int w, int h, int rank)
{
  if(rank == 0)
  {
    //Calculating the transpose.
    int pos = 0;
    for(int i = 0; i<h; i++)
    {
      for(int j = 0; j<w; j++)
      {
        out[pos] = in[i + j * w];
        pos++;
      }
    }
  }
}



void Inverser1D(Complex* h, int w, Complex* H)
{

  //cout<< "in inverse1d" << endl;

  Complex sum = Complex(0,0);
  Complex W;

  for(int n = 0; n<w; n++){
    for(int k = 0; k<w; k++){
      W = Complex( +cos(2 * M_PI * n * k / w), +sin(2 * M_PI * n * k / w) );
      sum = sum + W * h[k];
    }
    H[n] = sum;
    H[n].real = H[n].real/w;
    H[n].imag = H[n].imag/w;
    sum = Complex(0,0);
    //cout <<"\n" << H[n].Mag();
  }
}



void comp1d(Complex* data, Complex* H, int imgWidth, int imgHeight, int rank, int numtasks, int flag)
{
  int startRow, offset;
  startRow = imgHeight * rank / numtasks;

  //COMPUTING 1D TRANSFORM
  for (int i=0; i < imgHeight / numtasks; i++){
    offset = imgWidth * (startRow + i);
    if(flag == 0){
      Transform1D(data + offset, imgWidth, H + imgWidth * i);
    }
    else{
      Inverser1D(data + offset, imgWidth, H + imgWidth * i);
    }
  }
}



void t1d(Complex* data, Complex* H, int imgWidth, int imgHeight, int rank, int numtasks, int flag)
{
  //scatter
  ScatterGather(data, numtasks, rank, imgWidth * imgHeight, 1);

  comp1d(data, H, imgWidth, imgHeight, rank, numtasks, flag);

  // GATEHER
  ScatterGather(H, numtasks, rank, imgWidth*imgHeight, 0);
}



void t2d(Complex* H, Complex* out, int imgWidth, int imgHeight, int rank, int numtasks, int flag)
{
  Complex tpose[imgHeight * imgWidth];

  transpose(H, tpose, imgWidth, imgHeight, rank);
  t1d(tpose, H, imgWidth, imgHeight, rank, numtasks, flag);
  transpose(H, out, imgWidth, imgHeight, rank);
}





void Transform2D(const char* inputFN) 
{ 

  // Do the 2D transform here.
  // 1) Use the InputImage object to read in the Tower.txt file and
  //    find the width/height of the input image.
  // 2) Use MPI to find how many CPUs in total, and which one
  //    this process is
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
  InputImage image(inputFN);  // Create the helper object for reading the image
  // Step (1) in the comments is the line above.
  // Your code here, steps 2-9

  //InputImage image(inputFN);
  int imgHeight, imgWidth;
  Complex *imgData;

  imgHeight = image.GetHeight();
  imgWidth = image.GetWidth();
  imgData = image.GetImageData();

  int numtasks, rank;

  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  Complex *H;
  H = new Complex[imgWidth * imgHeight];

  // COMPUTING 1D TRANSFORM
  t1d(imgData, H, imgWidth, imgHeight, rank, numtasks, 0);


  if(rank == 0)
  {
    cout<<"Generating Image File MyAfter1d.txt"<<endl;
    image.SaveImageData("MyAfter1d.txt", H, imgWidth, imgHeight);
  }

  // COMPUTING 2D TRANFORM

  Complex* H2D = new Complex[imgHeight * imgWidth];

  t2d(H, H2D, imgWidth, imgHeight, rank, numtasks, 0);
  
  if(rank == 0)
  {
    cout<<"Generating Image File MyAfter2d.txt"<<endl;
    image.SaveImageData("MyAfter2d.txt", H2D, imgWidth, imgHeight);
  }

  // COMPUTING INVERSE TRANFORM

  Complex* IH = new Complex[imgHeight * imgWidth];
  
  // flag set to 1 for inverse transform
  t2d(H2D, IH, imgWidth, imgHeight, rank, numtasks, 1);
  t1d(IH, H, imgWidth, imgHeight, rank, numtasks, 1);

  //t1d(H2D, IH, imgWidth, imgHeight, rank, numtasks, 1);
  //t2d(IH, H, imgWidth, imgHeight, rank, numtasks, 1);

  if(rank == 0)
  {
    cout<<"Generating Image File MyAfterInverse.txt"<<endl;
    image.SaveImageData("MyAfterInverse.txt", H, imgWidth, imgHeight);
  }

}



void Transform1D(Complex* h, int w, Complex* H)
{
  // Implement a simple 1-d DFT using the double summation equation
  // given in the assignment handout.  h is the time-domain input
  // data, w is the width (N), and H is the output array.
  //cout<<"\nimage data: "<< h <<"\n";

  Complex sum = Complex(0,0);
  Complex W;

  for(int n = 0; n<w; n++){
    for(int k = 0; k<w; k++){
      W = Complex( +cos(2 * M_PI * n * k / w), -sin(2 * M_PI * n * k / w) );
      sum = sum + W * h[k];
    }
    H[n] = sum;
    sum = Complex(0,0);
    //cout <<"\n" << H[n].Mag();
  }
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
  

  
