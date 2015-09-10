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

void TxRx(Complex* H, int numtasks, int rank, int size, int flag)
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

  //printf("\nTasks:%d\tRank:%d\n", numtasks, rank);


  // 3) Allocate an array of Complex object of sufficient size to
  //    hold the 2d DFT results (size is width * height)
  Complex *H;

  H = new Complex[imgWidth * imgHeight];


  // 4) Obtain a pointer to the Complex 1d array of input data
  int startRow;
  startRow = imgHeight * rank / numtasks;
  //printf("\nstartrow:%d", startRow);
  int offset;

  // 5) Do the individual 1D transforms on the rows assigned to your CPU

  for (int i=0; i<imgHeight/numtasks; i++){
    offset = imgWidth*(startRow+i);
    //printf("\noffset:%d", offset);
    Transform1D(imgData + offset, imgWidth, H + imgWidth * i);
  }

  //Now, Master CPU will gather data from the rest of the CPUs.
  TxRx(H, numtasks, rank, imgWidth*imgHeight, 0);

  if(rank == 0)
  {
    cout<<"Generating Image File MyAfter1d.txt"<<endl;
    image.SaveImageData("MyAfter1d.txt", H, imgWidth, imgHeight);
  }

  // 6) Send the resultant transformed values to the appropriate
  //    other processors for the next phase.
  // 6a) To send and receive columns, you might need a separate
  //     Complex array of the correct size.

  Complex tpose[imgHeight * imgWidth];
  transpose(H, tpose, imgWidth, imgHeight, rank);


  TxRx(tpose, numtasks, rank, imgWidth*imgHeight, 1);

  Complex* after2D = new Complex[imgHeight*imgWidth];
  startRow = imgHeight*rank/numtasks;
  
  for(int i=0; i<imgHeight/numtasks; i++)
  {
    offset = imgWidth*(startRow + i);
    Transform1D(tpose + offset, imgWidth, after2D + imgWidth*i);
  }

  //Now, Master CPU will gather data from the rest of the CPUs.
  TxRx(after2D, numtasks, rank, imgHeight*imgWidth, 0);

  Complex tpose2D[imgWidth*imgHeight];

  transpose(after2D, tpose2D, imgWidth, imgHeight, rank);
  
  if(rank == 0)
  {
    cout<<"Generating Image File MyAfter2d.txt"<<endl;
    image.SaveImageData("MyAfter2d.txt", tpose2D, imgWidth, imgHeight);
  }
  



  // 7) Receive messages from other processes to collect your columns
  // 8) When all columns received, do the 1D transforms on the columns
  // 9) Send final answers to CPU 0 (unless you are CPU 0)
  //   9a) If you are CPU 0, collect all values from other processors
  //       and print out with SaveImageData().

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
    cout <<"\n" << H[n].Mag();
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
  

  
