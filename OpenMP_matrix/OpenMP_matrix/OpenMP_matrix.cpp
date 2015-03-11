// OpenMP_matrix.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <omp.h>



#include <conio.h>

#include <intrin.h>

#include <stdlib.h>

#include <math.h>

#include <iostream>
#include <assert.h>
#include <amp.h>
using namespace std;






void MatrixInit(std::vector<float> &M, int dim)

{

	for (int i = 0; i<dim; i++)

	{

		for (int j = 0; j<dim; j++)

		{

			M[i*dim + j] = (float)((float)rand() * 100 / (float)(RAND_MAX + 1));

		}

	}

}


int _tmain(int argc, _TCHAR* argv[])

{
	

	int n;

	double sum;

	int i, j, k;

	__int64 time2;
	srand(2012);

	while (true)
	{



		cout << "Enter size n: ";

		cin >> n;

		std::vector<float> MatrixA(n*n);

		std::vector<float> MatrixB(n*n);

		std::vector<float> MatrixC(n*n);

		MatrixInit(MatrixA, n);

		MatrixInit(MatrixB, n);

		time2 = __rdtsc();
#pragma omp parallel for private(j,k,sum)
		for (i = 0; i < n; i++)
		{
			

			for (k = 0; k < n; k++)
			{
				sum = 0;

				for (j = 0; j < n; j++)
				{

					sum += MatrixA[i*n + j] * MatrixB[j*n + k];

				}
				MatrixC[i*n + k] = sum;
			}

		}
		cout << "Parallel Matrix Multiplication\n";
		time2 = __rdtsc() - time2;
		cout << "Time2: " << time2 / 1595400000 << " seconds" << "\n";		

		system("pause");
	}

	return 0;
}