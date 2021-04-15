// Assignment 3.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream> 
#include <fstream> 
#include <math.h> 
#include <string.h>

using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{

	int N, a, b, i, j, n, A, Scheme;
	double h, dt, dtc, dth2, mc, te, nt;
	const double PI = 3.141592653589793;

	cout << "Select a Scheme \nIf you want to plot the initial function press 0\n\nFor Basic Explicit Scheme press 1 \nFor Fully Implicit Scheme press 2 \nFor Crank-Nicholson Scheme press 3 \n\nScheme= ";
	cin >> Scheme;

	cout << "Select Grid Size N= ";
	cin >> N;
	cout << "N= " << N << endl << endl;

	cout << "Input the end time te= ";
	cin >> te;
	cout << "te= " << te << endl << endl;

	cout << "Select the time step factor mc= ";
	cin >> mc;
	cout << "mc= " << mc << endl << endl;

	a = 0;
	b = 1;

	h = (b - a) / (N - 1.0);
	dtc = (h*h)/2.0;
	dt = mc*dtc;
	nt = ceil(te / dt);

	cout << "The number of time steps is nt= " << nt << endl << endl;
	cout << "Do you wish to use this many time steps?" << endl;
	cout << "Input 0 for Yes or 1 for No" << endl;
	cin >> A;


	double *x = NULL;
	double *u = NULL;
	double *ut = NULL;
	double *ut2 = NULL;

	x = new double[N];
	u = new double[N];
	ut = new double[N];
	ut2 = new double[N];

	// Initialize the value of u[n]

	for (i = 0; i < N; i++)
		x[i] = i*h;

	for (i = 0; i < N; i++)
	{
		u[i] = 0;
		ut[i] = 0;
		ut2[i] = 0;
	}

	for (i = 1; i < N - 1; i++)
		u[i] = 5.0 * sin(PI*x[i]) + 3.0 * sin(8.0 * PI*x[i]) + sin(64.0 * PI*x[i]);

	if (te == 0) Scheme = 0;
	
	if (A == 1) return 0;

	if (A == 0)
		dt = te / nt;

	dth2 = dt / (h*h);
	
	if (Scheme == 0)

	{

		ofstream file0;
		file0.open("cw3-0.dat");

		if (file0.is_open())

			for (i = 0; i < N; i++)

				file0 << x[i] << "  " << u[i] << endl;

		file0.close();
	}

	if (Scheme == 1)

	{

		for (n = 1; n < nt + 1; n++)

		{

			if (n > 1)
			{
				for (i = 1; i < N; i++)
					u[i] = ut[i];
			}

			for (i = 1; i < N - 1; i++)
				ut[i] = u[i] + dth2*(u[i + 1] - 2 * u[i] + u[i - 1]);
			
		}


		ofstream file1;
		file1.open("cw3-1.dat");

		if (file1.is_open())

			for (i = 0; i < N; i++)

				file1 << x[i] << "  " << ut[i] << endl;

		file1.close();
	
	}

	if (Scheme == 2)

	{

		// Initializing the Tridiogonal Matrix

		double B1, B2, C, D, E;

		double *alpha = NULL;
		double *beta = NULL;
		double *F = NULL;

		alpha = new double[N];
		beta = new double[N];
		F = new double[N];

		C = -dth2;
		D = 1.0 + (2.0*dth2);
		E = C;

		B1 = u[0];
		B2 = u[N - 1];

		for (n = 1; n < nt + 1; n++)

		{

			if (n > 1)

			{

				for (i = 1; i < N; i++)
					u[i] = ut[i];

			}

			if (n > 0)

			{

				// Solve the Matrix

				for (i = 1; i < N - 1; i++)
					F[i] = u[i];

				alpha[1] = -(D / E);
				beta[1] = (F[1] - C*B1) / E;

				for (i = 2; i < N - 1; i++)

				{
					alpha[i] = (-1.0 / E)*(D + (C / alpha[i - 1]));
					beta[i] = (1.0 / E)*(F[i] + ((C*beta[i - 1]) / alpha[i - 1]));
				}

				ut[N - 1] = B2;

				for (i = N - 2; i > 0; i--)
					ut[i] = (1 / alpha[i])*(ut[i + 1] - beta[i]);

			}

			for (i = 1; i < N - 1; i++)
				ut2[i] = u[i] + dth2*(ut[i + 1] - 2 * ut[i] + ut[i - 1]);

		}


		ofstream file2;
		file2.open("cw3-2.dat");

		if (file2.is_open())

			for (i = 0; i < N; i++)

				file2 << x[i] << "  " << ut2[i] << endl;

		file2.close();

	}

	if (Scheme == 3)

	{

		// Initializing the Crank-Nicholson Tridiogonal Matrix

		double B1, B2, C, D, E, G;

		double *alpha = NULL;
		double *beta = NULL;
		double *F = NULL;

		alpha = new double[N];
		beta = new double[N];
		F = new double[N];

		C = -dth2 / 2.0;
		E = C;
		D = 1.0 + dth2;
		G = 1.0 - dth2;

		B1 = u[0];
		B2 = u[N - 1];

		for (n = 1; n < nt + 1; n++)

		{

			if (n > 1)

			{

				for (i = 1; i < N; i++)
					u[i] = ut[i];
				
			}

			if (n > 0)

			{

				// Solve the Matrix

				for (i = 1; i < N - 1; i++)
					F[i] = -C*u[i + 1] + G*u[i] + -E*u[i - 1];

				alpha[1] = -(D / E);
				beta[1] = (F[1] - C*B1) / E;

				for (i = 2; i < N - 1; i++)

				{
					alpha[i] = (-1.0 / E)*(D + (C / alpha[i - 1]));
					beta[i] = (1.0 / E)*(F[i] + ((C*beta[i - 1]) / alpha[i - 1]));
				}

				ut[N - 1] = B2;

				for (i = N - 2; i > 0; i--)
					ut[i] = (1 / alpha[i])*(ut[i + 1] - beta[i]);

			}

			for (i = 1; i < N - 1; i++)
				ut2[i] = u[i] + (dth2 / 2.0)*((ut[i + 1] - 2 * ut[i] + ut[i - 1]) + (u[i + 1] - 2 * u[i] + u[i - 1]));

			
		}


		ofstream file3;
		file3.open("cw3-3.dat");

		if (file3.is_open())

			for (i = 0; i < N; i++)

				file3 << x[i] << "  " << ut2[i] << endl;

		file3.close();

	}
	
	system("pause");
	return 0;
}