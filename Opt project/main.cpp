/***************************************************
Code written for the optimization exercises purposes
by Lukasz Sztangret, PhD
Department of Applied Computer Science and Modelling
AGH University of Science and Technology
***************************************************/

#include <stdlib.h>
#include <time.h>
#include"opt_alg.h"


int main()
{

	try
	{
		std::cout << "LAB NUMBER " << LAB_NO << endl;
		std::cout << "LAB PART " << LAB_PART << endl << endl;

#if LAB_NO==0
		
#elif LAB_NO==1 && LAB_PART==1
		srand(time(NULL));
	//	double m = 10, b = 2.5, k = 100, F = 10;
		double t0 = 0, tend = 50, dt = 0.1, d, Nmax;
		double alpha = rand() % 100 + 101;
		double x0 = rand() % 200 - 100;
		alpha /= 100;
		d =	1; 
		double* p = new double[2];

	//	matrix Y0 = matrix(2, new double[2]{ 0,0 });
	//	matrix P = matrix(4, new double[4]{ m,b,k,F }) ;
	//	matrix* Y = solve_ode(t0, dt, tend, Y0, &P);
	//	matrix OUT = hcat(Y[0], Y[1]);
	//	ofstream S("out_1_1.csv");
	//	S << OUT;
	//	S.close();
		matrix* ud = nullptr;
		matrix* ad = nullptr;


		ofstream plik("konspekt_1_f_3.csv");
		cout << alpha << endl;


		//for (int i = 0; i < 100; i++)
		//{
		//	//alpha = rand() % 100 + 101;
		//	x0 = rand() % 200 - 100;
		//	//alpha /= 100;

		//	p = expansion(x0, d, alpha, 100000, ud, ad);
		//	plik << x0 << ";" << p[0] << ";" << p[1] << ";" << solution::f_calls<<";";
		//	solution::f_calls = 0;

		//	solution wynik2 = fib(p[0], p[1], 0.0001);
		//	plik << wynik2.x(0) << ";" << wynik2.y(0) << ";" << solution::f_calls << ";";
		//	string lok = "";
		//	solution::f_calls = 0;

		//	if (wynik2.y(0) < -0.7)
		//		lok = "tak";
		//	else
		//		lok = "nie";

		//	plik << lok << ";";
		//	
		//	solution wynik = lag(p[0], p[1], 0.0001, 0.000001, 100000);
		//	plik << wynik.x(0) << ";" << wynik.y(0) << ";"<< solution::f_calls << ";";
		//	solution::f_calls = 0;

		//	if (wynik.y(0) < -0.7)
		//		lok = "tak";
		//	else
		//		lok = "nie";

		//	plik << lok << ";\n";
		//}
		//plik.close();
		solution::f_calls = 0;
		solution wynik2 = fib(-100, 100, 0.00001);
		cout << "fib: x: " << wynik2.x(0) << " y: "<< wynik2.y()<<"n wywolan " << solution::f_calls<<endl;
		solution::f_calls = 0;
		solution wynik = lag(-100, 100, 0.00001, 0.000001, 100000);
		cout << "lag: x: " << wynik.x(0) << " y: " << wynik.y() << "n wywolan " << solution::f_calls<<endl;

		std::cout << "p = [" << p[0] << ", " << p[1] << "]" << endl;

	/*	std::cout << "m = " << P(0) << endl;
		std::cout << "b = " << P(1) << endl;
		std::cout << "k = " << P(2) << endl;
		std::cout << "F = " << P(3) << endl;*/


#elif LAB_NO==1 && LAB_PART==2
		double m = 10, b = 2.5, k = 100, FA = 10, Ff = 2;
		double t0 = 0, tend = 50, dt = 0.1;
		matrix Y0 = matrix(2, new double[2]{ 0,0 });
		matrix P = matrix(5, new double[5]{ m,b,k,FA, Ff });
		matrix* Y = solve_ode(t0, dt, tend, Y0, &P);
		matrix OUT = hcat(Y[0], Y[1]);
		ofstream S("out_1_2.csv");
		S << OUT;
		S.close();
		cout << "m = " << P(0) << endl;
		cout << "b = " << P(1) << endl;
		cout << "k = " << P(2) << endl;
		cout << "FA = " << P(3) << endl;
		cout << "Ff = " << P(4) << endl;
		
#elif LAB_NO==2 && LAB_PART==1

		
#elif LAB_NO==2 && LAB_PART==2
		
#elif LAB_NO==2 && LAB_PART==3
		std::srand(time(NULL));
		double alpha = 1.02;
		double x0 = 0.005; // [0,0001 ; 0,01] 
		double d = 0.0001;
		double* p = new double[2];
		matrix* ad = new matrix(0.00252941);
		matrix* ud = nullptr;
		//p = expansion(x0, d, alpha, 100000, ud, ad);
		//solution::f_calls = 0;
		//solution fibonacci = fib(p[0], p[1], 0.00000000001);
		//std::cout << "y: " << fibonacci.y(0) << " x:" << fibonacci.x(0) << "  n wywolan " << solution::f_calls << endl;
		//solution::f_calls = 0;
		
		//solution lagrange = lag(p[0], p[1], 1e-8, 1e-8, 100000);
		//std::cout << "y: " << lagrange.y(0) << " x:" << lagrange.x(0) << "  n wywolan " << solution::f_calls << endl;
		matrix Y0 = matrix(3, new double[3]{ 5,1,10 });
		matrix* Y = solve_ode(0, 1, 1000, Y0, ud, ad);


#elif LAB_NO==3 && LAB_PART==1
		std::srand(time(NULL));
		double s[3] = { 0.1, 0.5, 0.7 };
		ofstream plik1("konspekt_2_f_1.csv");
		ofstream plik2("konspekt_2_f_2.csv");

		for (int i = 0; i < 100; i++)
		{
			double x01 = rand() % 200 - 100;
			double x02 = rand() % 200 - 100;
			x01 /= 100;
			x02 /= 100;
			//std::cout << endl << "punkt startowy: x=[" << x01 << ", " << x02 << "]" << endl;
			matrix x0(2, 1);
			x0(0, 0) = x01;
			x0(1, 0) = x02;
			//std::cout << x0 << endl;
			solution hjSol = HJ(x0, s[2], 0.00001, 0.00001, 10000000);
			plik1 << x01<<";"<<x02<<";"<<  hjSol.x(0, 0) << ";" << hjSol.x(1, 0) << ";" << hjSol.y[0] <<";"<< solution::f_calls << endl;
			solution::f_calls = 0;

			
			matrix s0 = matrix(2, 1, s[2]);
			solution rosSol = Rosen(x0, s0, 1.01, 0.1, 0.00001, 10000000);
			plik2 << x01 << ";" << x02 << ";" << rosSol.x(0, 0) << ";" << rosSol.x(1, 0) << ";" << rosSol.y[0] << solution::f_calls << endl;
			solution::f_calls = 0;
		}

		plik1.close();
		plik2.close();
#elif LAB_NO==3 && LAB_PART==2
		matrix x0(2, 1);
		x0(0, 0) = 0;
		x0(1, 0) = 0;
		
		solution hjSol = HJ(x0, 0.3, 0.00001, 0.00001, 10000000);
		cout <<  hjSol.x(0, 0) << ";" << hjSol.x(1, 0) << ";" << hjSol.y[0]<< solution::f_calls << endl;
		solution::f_calls = 0;

		cout << endl;

		matrix s0 = matrix(2, 1, 0.3);
		solution rosSol = Rosen(x0, s0, 1.01, 0.1, 0.00001, 10000000);
		cout  << rosSol.x(0, 0) << ";" << rosSol.x(1, 0) << ";" << rosSol.y[0] << solution::f_calls << endl;

#elif LAB_NO==3 && LAB_PART==3
	
		
		matrix x(2, 1);
		x(0, 0) = 1.5;
		x(1, 0) = 1.5;

		matrix x2(2, 1);
		x2(0, 0) = 1.5347;
		x2(1, 0) = 1.67163;

		matrix* ud = nullptr;
		matrix Y0(2, 1);

		//matrix* Y = solve_ode(0, 0.1, 100, Y0, ud, &x);
		matrix* Y = solve_ode(0, 0.1, 100, Y0, ud, &x2);

#elif LAB_NO==4 && LAB_PART==1

		matrix x0;
		double c = 1;
		double a[3] = { 4, 4.4934, 5 };
		for (int i = 0; i < 100; i++)
		{

			do
				x0 = 5 * rand_mat(2, 1) + 1;
			while (norm(x0) > a[0]); //a to s¹ alfy z konspektu

			cout << x0(0, 0) << " " << x0(1, 0)<<endl;

			solution sn = pen(x0, c = 1, 2, 0.00001, 100000);
			cout<< sn
			cout<<

		}


#elif LAB_NO==4 && LAB_PART==2
		
#elif LAB_NO==5 && LAB_PART==1
		
#elif LAB_NO==5 && LAB_PART==2
		
#elif LAB_NO==5 && LAB_PART==3
		
#elif LAB_NO==6 && LAB_PART==1
		
#elif LAB_NO==6 && LAB_PART==2
		
#elif LAB_NO==7 && LAB_PART==1
		
#elif LAB_NO==7 && LAB_PART==2
		
#endif
	}
	catch (char * EX_INFO)
	{
		std::cout << EX_INFO << endl;
	}
	system("pause");
	return 0;
}
