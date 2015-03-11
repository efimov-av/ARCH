#include"stdafx.h"
#include <mpi.h>
#include <iostream>
#include<stdio.h>
#include<stdlib.h>
#include<conio.h>
#include<math.h>
#include<time.h>
<<<<<<< HEAD
#include <iostream>
int ProcNum;
int ProcRank;
int flag = 0;
int Size;
double *A;  double *B; double *C;

//***********************************************************
void PrintMatrix(double* pMatrix, int Size) {
	for (int i = 0; i<Size; i++) {
		printf("\n");
		for (int j = 0; j<Size; j++)
			printf("%7.4f ", pMatrix[i*Size + j]);
	}
}
//------------------------------------------------------------
void RandInit(double* pMatrix, int Size) {
	srand(100);
	for (int i = 0; i<Size; i++) {
		for (int j = 0; j<Size; j++)  pMatrix[i*Size + j] = rand() / double(1000);
	}
}
//-------------------------------------------------
void InitProcess(double* &A, double* &B, double* &C, int &Size) {
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	if (ProcRank == 0) {
		do {
			printf("\n--Square matrix multiplication--");
			printf("\nPlease, enter matrix size: ");
			std::cout << "Enter size n: "; std::cin >> Size;
			if (Size< ProcNum) printf("Matrix size is less than the number of processes! \n");
			if (Size%ProcNum != 0) printf("Matrix size should be dividable by the number of processes! \n");
		} while ((Size< ProcNum) || (Size%ProcNum != 0));
	}
	if (Size<10) flag = 1;
	MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (ProcRank == 0) {
		A = new double[Size*Size];
		B = new double[Size*Size];
		C = new double[Size*Size];
		RandInit(A, Size); RandInit(B, Size);
	}
}
//-------------------------------------------------
void Flip(double *&B, int dim) {
	double temp = 0.0;
	for (int i = 0; i<dim; i++){
		for (int j = i + 1; j<dim; j++){
			temp = B[i*dim + j]; B[i*dim + j] = B[j*dim + i]; B[j*dim + i] = temp;
=======

using namespace std;
int ProcNum;  
int ProcRank; 
int flag=0;
int Size; 
double *A;  double *B; double *C;

//***********************************************************
void PrintMatrix(double* pMatrix,int Size) {
	for (int i=0;i<Size;i++) {printf("\n");
	for (int j=0;j<Size;j++)
		printf("%7.4f ",pMatrix[i*Size+j]);
	}
}
//------------------------------------------------------------
void RandInit (double* pMatrix, int Size) {
	srand(100);
	for (int i=0; i<Size; i++) {
		for (int j=0;j<Size;j++)  pMatrix[i*Size+j]=rand()/double(1000);
	}
}
//-------------------------------------------------
void InitProcess (double* &A,double* &B,double* &C ,int &Size) {
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	printf("ProcRank %7.4f",ProcRank);
	if (ProcRank == 0) {
		do {printf("\n--Square matrix multiplication--");
		
		printf("\nPlease, enter matrix size: "); scanf_s("%d", &Size); 
		Size=4096;
		if (Size< ProcNum) printf("Matrix size is less than the number of processes! \n");
		if (Size%ProcNum!=0) printf("Matrix size should be dividable by the number of processes! \n");
		}
		while ((Size< ProcNum)||(Size%ProcNum!=0));
	}
	if (Size<10) flag=1;
	MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (ProcRank == 0) {
		A = new double [Size*Size];
		B = new double [Size*Size];
		C = new double [Size*Size];
		RandInit (A, Size); RandInit (B, Size);
	}
}
//-------------------------------------------------
void Flip (double *&B, int dim) { 
	double temp=0.0; 
	for (int i=0; i<dim; i++){
		for (int j=i+1; j<dim; j++){
			temp=B[i*dim+j];B[i*dim+j]=B[j*dim+i];B[j*dim+i]=temp;
>>>>>>> origin/master
		}
	}
}
//-------------------------------------------------
void MatrixMultiplicationMPI(double *&A, double *&B, double *&C, int &Size) {
<<<<<<< HEAD
	int dim = Size;
	int i, j, k, p, ind;
	double temp;
	MPI_Status Status;
	int ProcPartSize = dim / ProcNum;
	int ProcPartElem = ProcPartSize*dim;
	double* bufA = new double[dim*ProcPartSize];
	double* bufB = new double[dim*ProcPartSize];
	double* bufC = new double[dim*ProcPartSize];
	int ProcPart = dim / ProcNum, part = ProcPart*dim;
	if (ProcRank == 0) {
		Flip(B, Size);
	}

	MPI_Scatter(A, part, MPI_DOUBLE, bufA, part, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(B, part, MPI_DOUBLE, bufB, part, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	temp = 0.0;
	for (i = 0; i<ProcPartSize; i++){
		for (j = 0; j<ProcPartSize; j++){
			for (k = 0; k<dim; k++) temp += bufA[i*dim + k] * bufB[j*dim + k];
			bufC[i*dim + j + ProcPartSize*ProcRank] = temp; temp = 0.0;
		}
	}

	int NextProc; int PrevProc;
	for (p = 1; p<ProcNum; p++) {
		NextProc = ProcRank + 1;
		if (ProcRank == ProcNum - 1) NextProc = 0;
		PrevProc = ProcRank - 1;
		if (ProcRank == 0) PrevProc = ProcNum - 1;
		MPI_Sendrecv_replace(bufB, part, MPI_DOUBLE, NextProc, 0, PrevProc, 0, MPI_COMM_WORLD, &Status);
		temp = 0.0;
		for (i = 0; i<ProcPartSize; i++) {
			for (j = 0; j<ProcPartSize; j++) {
				for (k = 0; k<dim; k++){
					temp += bufA[i*dim + k] * bufB[j*dim + k];
				}
				if (ProcRank - p >= 0)
					ind = ProcRank - p;
				else ind = (ProcNum - p + ProcRank);
				bufC[i*dim + j + ind*ProcPartSize] = temp;
				temp = 0.0;
			}
		}
	}
	MPI_Gather(bufC, ProcPartElem, MPI_DOUBLE, C, ProcPartElem, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	delete[]bufA;
	delete[]bufB;
	delete[]bufC;
=======
	int dim =Size;
	int i, j, k, p, ind;
	double temp;
	MPI_Status Status;
	int ProcPartSize=dim/ProcNum; 
	int ProcPartElem=ProcPartSize*dim; 
	double* bufA=new double[dim*ProcPartSize];
	double* bufB=new double[dim*ProcPartSize];
	double* bufC=new double[dim*ProcPartSize];
	int ProcPart = dim/ProcNum, part = ProcPart*dim;
	if (ProcRank == 0) {
		Flip(B,Size);
	}

	MPI_Scatter(A, part, MPI_DOUBLE, bufA, part, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(B, part, MPI_DOUBLE, bufB, part, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	temp=0.0;
	for(i=0;i<ProcPartSize;i++){
		for(j=0;j<ProcPartSize; j++){
			for(k=0;k<dim;k++) temp+=bufA[i*dim+k]*bufB[j*dim+k];
			bufC[i*dim+j+ProcPartSize*ProcRank]=temp;temp=0.0;
		}
	}

	int NextProc; int PrevProc;
	for(p=1; p<ProcNum; p++) {
		NextProc = ProcRank+1;
		if ( ProcRank ==ProcNum-1 ) NextProc = 0;
		PrevProc = ProcRank - 1;
		if ( ProcRank == 0 ) PrevProc = ProcNum-1;
		MPI_Sendrecv_replace(bufB,part, MPI_DOUBLE, NextProc, 0, PrevProc, 0, MPI_COMM_WORLD, &Status);
		temp=0.0;
		for (i=0; i<ProcPartSize; i++) {
			for(j=0; j<ProcPartSize; j++) {
				for(k=0; k<dim; k++){
					temp +=bufA[i*dim+k]*bufB[j*dim+k];}
				if (ProcRank - p >= 0 ) 
					ind=ProcRank - p;
				else ind=(ProcNum - p+ProcRank);
				bufC[i*dim+j+ind*ProcPartSize]=temp;
				temp=0.0;
			}
		}
	}
	MPI_Gather(bufC, ProcPartElem, MPI_DOUBLE, C, ProcPartElem, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	delete []bufA;
	delete []bufB;
	delete []bufC;
>>>>>>> origin/master
}

//--------------------------------------------------------

void main(int argc, char* argv[]) {
	double beg, end, serial, parallel=0;
	MPI_Init ( &argc, &argv );   
	InitProcess (A,B,C,Size);
	beg=MPI_Wtime();
	MatrixMultiplicationMPI(A,B,C,Size);
	end=MPI_Wtime(); parallel= end-beg;  
	if (ProcRank == 0) {
		printf("\n",&ProcNum);
		printf ("\n");
		printf("\nTime of execution -  Parallel calculation:\n");
		printf("%7.4f",parallel);
		if (flag) {
			printf("\nMatrix C - Parallel calculation\n");
			PrintMatrix(C,Size); 
		}    
	}
	MPI_Finalize();
<<<<<<< HEAD
	delete[] A; delete[] B; delete[] C;
}
=======
	delete [] A; delete [] B; delete [] C;   
	

}
>>>>>>> origin/master
