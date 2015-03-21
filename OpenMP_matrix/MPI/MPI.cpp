#include "stdafx.h"
#include<mpi.h>
#include<iostream>
#include<stdlib.h>
#include<conio.h>
#include<math.h>
#include<time.h>
#include <omp.h>
int ProcNum;
int ProcRank;
int flag = 0;
int Size;
double *A;  double *B; double *C; double *CS; double *CM;



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
			printf("\nPlease, enter matrix size: "); //scanf_s("%d", &Size);	
			Size = 16;
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
			temp = B[i*dim + j]; 
			B[i*dim + j] = B[j*dim + i]; 
			B[j*dim + i] = temp;
		}
	}
}
//-------------------------------------------------
void MatrixMultiplicationMPI(double *&A, double *&B, double *&C, int &Size) {
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
}

//--------------------------------------------------------

void main(int argc, char* argv[]) {
	double beg, end, parallel = 0;
	double result;
	int row, col, inner;
	bool passed = true;	
	MPI_Init(&argc, &argv);
	InitProcess(A, B, C, Size);
	
	beg = MPI_Wtime();
	MatrixMultiplicationMPI(A, B, C, Size);
	end = MPI_Wtime(); parallel = end - beg;
	if (ProcRank == 0) {		
		
		
		printf("\nTime of execution -  MPI calculation:\n");
		printf("%7.4f", parallel);

		CS = new double[Size*Size];
		
		beg = MPI_Wtime();
		for ( row = 0; row < Size; row++) {			
			for ( col = 0; col < Size; col++) {
				result = 0;
				// Multiply the row of A by the column of B to get the row, column of product.
				for ( inner = 0; inner < Size; inner++) {
					result += A[row*Size + inner] * B[col*Size + inner];
				}
				CS[row*Size + col] = result;
			}			
		}
		end = MPI_Wtime(); parallel = end - beg;
		
		for ( row = 0; row < Size; row++) {
			for ( col = 0; col < Size; col++) {
				// Multiply the row of A by the column of B to get the row, column of product.
				
				if (CS[row*Size + col] != C[row*Size + col])passed = false;
				
			}
		}
		if (passed) {
			printf("\nMatrix CS = C");
			printf("\nMatrix CS - Single calculation\n");
			printf("%7.4f", parallel);			
		}else printf("\nMatrix CS != C");

		CM = new double[Size*Size];
		
		beg = MPI_Wtime();
		#pragma omp parallel for private(row,col,result)
		for (row = 0; row < Size; row++) {
			for (col = 0; col < Size; col++) {
				result = 0;
				
				for ( inner = 0; inner < Size; inner++) {
					result += A[row*Size + inner] * B[col*Size + inner];
				}
				CM[row*Size + col] = result;
			}
		}
		end = MPI_Wtime(); parallel = end - beg;

		for (row = 0; row < Size; row++) {
			for ( col = 0; col < Size; col++) {
				// Multiply the row of A by the column of B to get the row, column of product.

				if (CM[row*Size + col] != C[row*Size + col])passed = false;

			}
		}
		if (passed) {
			printf("\nMatrix CM = C");
			printf("\nMatrix CM - Multi calculation\n");
			printf("%7.4f", parallel);
		}
		else printf("\nMatrix CM != C");
		scanf_s("%d", &Size);
	}
	MPI_Finalize();
	delete[] A; delete[] B; delete[] C;	
	
}