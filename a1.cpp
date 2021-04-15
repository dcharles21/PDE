// Assignment 1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream> // needed for input and output
#include <fstream> // needed for streams
#include <math.h> // needed for mathematical functions
#include <string.h>

using namespace std;

//Function for Lax-Wendroff Method

double F(double u, double c)

{
	double f;
	f = u + 0.5*c*pow(u, 2);
	return f;
}


int _tmain(int argc, _TCHAR* argv[])
{

	int i, j, n, t, N, nt, Method;
	double n1, h, dt, c;

	double *x = NULL;
	double *u = NULL;
	double *ut = NULL;
	double *up = NULL;

	cout << "Select a Method \nFor Forward Differencing press 1 \nFor Backwards Differencing press 2 \nFor Lax-Wendroff Scheme press 3 \n\nMethod= ";
	cin >> Method;

	cout << "Select the number of grind points. N= ";
	cin >> N;
	cout << "N= " << N << endl << endl;

	cout << "Set the number of time steps nt=: ";
	cin >> nt;
	cout << "nt= " << nt << endl << endl;

	cout << "Set value of time constant dt=: ";
	cin >> dt;
	cout << "dt= " << dt << endl << endl;

	x = new double[N + 1];
	u = new double[N];
	ut = new double[N];
	up = new double[N];

	h = 1.0 / N;
	n = 0;
	
	for (i = 0; i <= N; i++)

	{
		x[i] = i*h;

		if  (x[i] <= 0.25 || x[i] >= 0.75)
			u[i] = 0;

		else
			u[i] = 1.0 - 16.0*(pow((x[i] - 0.5), 2));


		u[N] = u[0];
	}

	ofstream file0;
	file0.open("cw0.dat"); //cw1-0 x values for time t=0

	if (file0.is_open())

		for (i = 0; i <= N - 1; i++)
		{
			file0 << x[i] << "  " << u[i] << endl;
		}

	file0 << x[N] << " " << u[0] << endl;
	
	file0.close();
	
	// Forward Difference Loop

	if (Method == 1)

		{

		for (n = 1; n <= nt; n++)

			for (i = 0; i <= N - 1; i++)

			{

				if (i == N - 1)

				u[i + 1] = u[0];

				ut[i] = u[i] - (dt / h)*(u[i + 1] - u[i]);

				u[i] = ut[i];

			}

		ofstream file1;
		file1.open("cw1.dat");

		if (file1.is_open())

			for (i = 0; i <= N - 1; i++)
				file1 << x[i] << " " << u[i] << endl;
			
		file1 << x[N] << " " << u[0] << endl;
		
		file1.close();
	}
	
	//Backwards Difference Loop

if (Method == 2)

	{

		for (n = 1; n <= nt; n++)

		{

			for (i = 0; i <= N - 1; i++)
			{

			ut[0] = u[0] - (dt / h)*(u[0] - u[N-1]);
			ut[i + 1] = u[i + 1] - (dt / h)*(u[i + 1] - u[i]);

			u[i] = ut[i];

			}
		}

		ofstream file2;
		file2.open("cw2.dat");

		if (file2.is_open())
		{

			for (i = 0; i <= N - 1; i++) 
				file2 << x[i] << " " << u[i] << endl;

		file2 << x[N] << " " << u[0] << endl;

		file2.close();
		}
	}


   //Lax-Wendroff Loop 

if (Method == 3)

{

	cout << "choose a constant c= ";
	cin >> c;
	cout << "c= " << c << endl;


	for (n = 1; n <= nt; n++)


	{

		for (i = 0; i <= N - 1; i++)

		{

			if (i == N - 1)
				u[i + 1] = u[0];

			up[i] = 0.5*(u[i] + u[i + 1]) - 0.5*(dt / h)*(F(u[i + 1], c) - F(u[i], c));
		}

		for (j = 0; j <= N - 1; j++)

			{
			

			ut[0] = u[0] - (dt / h)*(F(up[0], c) - F(up[N - 1], c));
			ut[j + 1] = u[j + 1] - (dt / h)*(F(up[j + 1], c) - F(up[j], c));

			u[j] = ut[j];

			}
	}
	
	ofstream file3;
	file3.open("cw3.dat");

	if (file3.is_open())
	{  

		for (i = 0; i <= N - 1; i++)
					file3 << x[i] << " " << u[i] << endl;
		file3 << x[N] << " " << u[0] << endl;

		file3.close();
	}
	
}
		
}
	system("pause");
	return 0;	
}