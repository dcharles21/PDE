// Assignment 2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdafx.h"
#include <iostream> 
#include <fstream> 
#include <math.h> 
#include <string.h>

using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	int i, j, n, N, MaxIterations, Method, Count;
	double h, p, w, h2rho;

	cout << "Select a Method \nFor Jacobi Method press 1 \nFor Gauss-Seidel Method press 2 \nFor Successive Over-Relaxation (SOR) Method press 3 \n\nMethod= ";
	cin >> Method;

	cout << "Select the grid size from. N= ";
	cin >> N;
	cout << "N= " << N << endl;

	cout << "Select the Max Number of Iterations= ";
	cin >> MaxIterations;
	cout << "Max Iterations= " << MaxIterations << endl;

	cout << "Select a value for Rho. p= ";
	cin >> p;
	cout << "p= " << p << endl;

	double *x = NULL;
	double *y = NULL;
	double *MaxError = NULL;

	x = new double[N];
	y = new double[N];
	MaxError = new double[MaxIterations];

	double phi[100][100];
	double phiold[100][100];
	double error[100][100];

	h = 1.0 / (N - 1.0);
	Count = 0;
	MaxError[0] = 0.001;

	h2rho = h*h*p;

	// This loop sets the values of x(i) and y(j)

	for (i = 0; i < N; i++)

	{

		x[i] = i*h;
		y[i] = i*h;

	}

	// The loop below fills in the boundaries of the grid

	for (i = 0; i < N; i++)
	{
		phiold[i][0] = pow(y[0], 2) - pow(x[i], 2);
		phiold[i][N - 1] = pow(y[N - 1], 2) - pow(x[i], 2);


		phiold[0][i] = pow(y[i], 2) - pow(x[0], 2);
		phiold[N - 1][i] = pow(y[i], 2) - pow(x[N - 1], 2);
	}

	// The series of loops below fill in the rest of the grid with an inital guess of 0

	for (i = 1; i < N - 1; i++)
		for (j = 1; j < N - 1; j++)
			phiold[i][j] = 0;
	

	// Copy the values of phiold into phi

	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			phi[i][j] = phiold[i][j];

	// initializes the value of the error[i][j]

	for (i = 0; i < N; i++)
	{
		error[i][0] = 0;
		error[i][N - 1] = 0;

		error[0][i] = 0;
		error[N - 1][i] = 0;
	}

	// Jacobi Method

	if (Method == 1)

	{

		for (n = 1; n < MaxIterations + 1; n++)

	
		{
			if (n > 1)

			{

				for (i = 1; i < N - 1; i++)
					for (j = 1; j < N - 1; j++)
						phiold[i][j] = phi[i][j];
			}


			if (MaxError[n - 1] <= 0.0001) break;

			Count++;

			for (i = 1; i < N - 1; i++)


			{

				for (j = 1; j < N - 1; j++)
				{

					phi[i][j] = 0.25*(phiold[i + 1][j] + phiold[i - 1][j] + phiold[i][j + 1] + phiold[i][j - 1] - h2rho);

					error[i][j] = abs(phi[i][j] - phiold[i][j]);

					if (error[i][j] >= MaxError[n])
						MaxError[n] = error[i][j];

					if (i == j)
						phi[i][j] = 0;

				}

			}

		}

		cout << "The Maximum Error is " << MaxError[n - 1] << " after " << Count << " iterations" << endl;

		ofstream file0;
		file0.open("cw2-1.dat");

		if (file0.is_open())

			for (i = 0; i < N; i++)
				for (j = 0; j < N; j++)

					file0 << x[i] << "  " << y[j] << " " << phi[i][j] << endl;


		file0.close();

	}

	// Gauss-Seidel Method

	if (Method == 2)

	{

		for (n = 1; n < MaxIterations + 1; n++)

		{

			if (MaxError[n - 1] <= 0.0001) break;

			Count++;

			for (i = 1; i < N - 1; i++)


			{
				for (j = 1; j < N - 1; j++)

				{

					phi[i][j] = 0.25*(phiold[i + 1][j] + phi[i - 1][j] + phiold[i][j + 1] + phi[i][j - 1] - h2rho);

					error[i][j] = abs(phi[i][j] - phiold[i][j]);

					if (error[i][j] >= MaxError[n])
						MaxError[n] = error[i][j];

					phiold[i][j] = phi[i][j];

				}

			}

		}

		cout << "The Maximum Error is " << MaxError[n - 1] << " after " << Count << " iterations" << endl;

		ofstream file1;
		file1.open("cw2-2.dat");

		if (file1.is_open())

			for (i = 0; i < N; i++)
				for (j = 0; j < N; j++)

					file1 << x[i] << "  " << y[j] << " " << phiold[i][j] << endl;

	}

	// Successive Over-Relaxation (SOR) Method 

	if (Method == 3)

	{

		cout << "Enter a value for omega. w= ";
		cin >> w;
		cout << "w= " << endl << endl;

		for (n = 1; n < MaxIterations + 1; n++)

		{

			if (MaxError[n - 1] <= 0.0001) break;

			Count++;

			for (i = 1; i < N - 1; i++)


			{
				for (j = 1; j < N - 1; j++)

				{

					phi[i][j] = w*0.25*(phiold[i + 1][j] + phi[i - 1][j] + phiold[i][j + 1] + phi[i][j - 1] - h2rho) + (1.0 - w)*phiold[i][j];

					error[i][j] = abs(phi[i][j] - phiold[i][j]);

					if (error[i][j] >= MaxError[n])
						MaxError[n] = error[i][j];

					phiold[i][j] = phi[i][j];

				}

			}

		}

		cout << "The Maximum Error is " << MaxError[n - 1] << " after " << Count << " iterations" << endl;

		ofstream file3;
		file3.open("cw2-3.dat");

		if (file3.is_open())

			for (i = 0; i < N; i++)
				for (j = 0; j < N; j++)

					file3 << x[i] << "  " << y[j] << " " << phiold[i][j] << endl;
	}
	
	system("pause");
	return 0;
}

