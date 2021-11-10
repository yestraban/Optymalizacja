#include"opt_alg.h"
#include <math.h>

#if LAB_NO>=1
double *expansion(double x0, double d, double alpha, int Nmax, matrix *ud, matrix *ad)
{
	double *p = new double[2];
	solution X0(x0), X1(x0+d);
	X0.fit_fun(ud, ad);
	X1.fit_fun(ud, ad);
	if (X0.y==X1.y)
	{
		p[0]=X0.x(0);
		p[1]=X1.x(0);
		return p;
	}
	if (X1.y>X0.y)
	{
		d *= -1;
		X1.x = X0.x+d;
		X1.fit_fun(ud, ad);
		if (X1.y>=X0.y)
		{
			p[0] = X1.x(0);
			p[1] = -X1.x(0);
			return p;
		}
	}
	solution X2;
	int i = 1;
	while (true)
	{
		X2.x = x0 + pow(alpha, i)*d; //to jest ju� dobrze
		X2.fit_fun(ud, ad);
		if (X2.y>=X1.y || solution::f_calls>Nmax)
			break;
		X0 = X1;
		X1 = X2;
		++i;
	}
	if (d > 0)
	{
		p[0] = X0.x(0);
		p[1] = X2.x(0);
	}
	else
	{
		p[0] = X2.x(0);
		p[1] = X0.x(0);
	}
	//d > 0 ? p = { X0.x(0), X2.x(0) } : p = { X2.x(0), X0.x(0) };
	return p;
}
int k(double a, double b, double epsilon)
{
	int i = 2, temp;
	int Fib[2]{ 1, 1 };
		while ((b - a) / epsilon >= Fib[1])
		{
			i++;
			temp = Fib[1];
			Fib[1] = Fib[0] + Fib[1];
			Fib[0] = Fib[1];
		}
	return i;
}
solution fib(double a, double b, double epsilon, matrix *ud, matrix *ad)
{
	int n = k(a, b, epsilon);
	int *F = new int[n] {1, 1};
	for (int i = 2; i < n; ++i)
		F[i] = F[i-1]+F[i-2];

	solution A(a), B(b), C, D;
	C.x(0) = b - ((1.0*F[n-2])/(1.0*F[n-1])) * (b - a);
	D.x(0) = a + b - C.x(0);
	C.fit_fun(ud, ad);
	D.fit_fun(ud, ad);
	//cout << b - D.x(0) << " " << C.x(0) - a << endl;
	for (int i = 0; i <= n - 3; ++i)
	{
		//cout << B.x(0) - A.x(0) << endl;
		if (C.y<D.y)
			B=D;
		else
			A=C;
		C.x = B.x(0) - ((1.0*F[n-i-2])/(1.0*F[n-i-1]))*(B.x(0)-A.x(0));
		D.x = A.x(0) + B.x(0) - C.x(0);
		C.fit_fun(ud, ad);
		D.fit_fun(ud, ad);
#if LAB_NO==2 && LAB_PART==2
		*(ud).addrow((B.x()-A.x());
#endif
	}
	return C;
}

solution lag(double a, double b, double epsilon, double gamma, int Nmax, matrix *ud, matrix *ad)
{
	solution A(a), B(b), C, D, D_old(a);
	C.x = (a+b)/2;
	A.fit_fun(ud, ad);
	B.fit_fun(ud, ad);
	C.fit_fun(ud, ad);
	double l, m;
	while (true)
	{
		//cout << B.x(0) - A.x(0) << endl;
		l = A.y(0)*(pow(B.x(0), 2) - pow(C.x(0), 2)) + B.y(0)*(pow(C.x(0), 2) - pow(A.x(0), 2)) + C.y(0)*(pow(A.x(0), 2) - pow(B.x(0), 2));
		m = A.y(0)*(B.x(0) - C.x(0)) + B.y(0)*(C.x(0) - A.x(0)) + C.y(0)*(A.x(0) - B.x(0));
		if (m<=0)
		{
			C.x = NAN;
			C.y = NAN;
			return C;
		}
		D.x = l / (2 * m);
		D.fit_fun(ud, ad);
		if (A.x(0) < D.x(0) || D.x(0)<C.x(0))
		{
			if (D.y(0)<C.y(0))
			{
				B.x = C.x;
				C.x = D.x;
			}
			else
				A.x = D.x;
		}
		else if (C.x(0)<D.x(0)||D.x(0)<B.x(0))
		{
			if (D.y(0) < C.y(0))
			{
				A.x = C.x;
				C.x = D.x;
			}
			else
				B.x = D.x;
		}
		else
		{
			C.x = NAN;
			C.y = NAN;
			return C;
		}
		
#if LAB_NO==2 && LAB_PART==2
		* (ud).addrow((B.x() - A.x());
#endif
		
		//if (B.x - A.x <epsilon || abs(D.x(0) - D_old.x(0))<gamma || solution::f_calls>Nmax)
		if (B.x - A.x < epsilon || abs(D.x(0) - D_old.x(0)) < gamma || solution::f_calls > Nmax)
			return C;
		D_old = D;
	}
}
#endif
#if LAB_NO>2
solution HJ(matrix x0, double s, double alpha, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	ofstream plik2("konspekt_2_f_wykres2.csv");
	solution XB(x0), XB_old, X;
	XB.fit_fun(ud, ad);
	while (true)
	{
		X = HJ_trial(XB, s, ud, ad);
		if (X.y < XB.y)
		{
			while (true)
			{
				XB_old = XB;
				XB = X;
#if LAB_NO==3 && LAB_PART==2
				? ? ?
#endif
				X.x = 2 * XB.x - XB_old.x;
				X.fit_fun(ud, ad);
				X = HJ_trial(X, s, ud, ad);
				if (X.y > XB.y)
					break;
				if (X.y == XB.y)
					return XB;
			}
			cout << XB.x(0, 0) << ";" << XB.x(1, 0) << endl;
			plik2 << XB.x(0, 0) << ";" << XB.x(1, 0) << endl;
			
		}
		else
			s = alpha * s;
		if (s<epsilon)
			return XB;
	}
	cout << endl << endl;
}

solution HJ_trial(solution XB, double s, matrix *ud, matrix *ad)
{
	int n=get_dim(XB);
	matrix D = ident_mat(n);		//czy to jest macierz kierunk�w? tak
	solution X;
	for (int i = 0; i < n; ++i)
	{
		X.x = XB.x + s * D[i];	//czy to dobrze
		X.fit_fun(ud, ad);
		if (X.y<XB.y)
			XB =X;
		else
		{
			X.x = XB.x - s * D[1];
			X.fit_fun(ud, ad);
			if (X.y < XB.y) 
				XB=X;
		}
	}
	return XB;
}

solution Rosen(matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	ofstream plik("konspekt_2_f_wykres.csv");
	solution X(x0), Xt;
	int n = get_dim(X);
	matrix l(n, 1), p(n, 1), s(s0), D = ident_mat(n);
	X.fit_fun(ud, ad);
	cout << X.x(0, 0) << ";" << X.x(1, 0) << endl;
	while (true)
	{
		for (int i = 0; i < n; ++i)
		{
			Xt.x = X.x+s(i)*D[i];
			Xt.fit_fun(ud, ad);
			if (Xt.y<X.y)
			{
				X = Xt;
				l(i) = l(i) + s(i);
				s(i) = alpha * s(i);  //poprawi� na okr�g�e
			}
			else
			{
				s(i) = -beta * s(i);
				p(i) = p(i) + 1;
			}
		}
#if LAB_NO==3 && LAB_PART==2
		(*ud).add_row(trans(X.x));
#endif
		bool change = true;
		for (int i = 0; i < n; ++i)
			if (l(i)==0 || p(i)==0)
			{
				change = false;
				break;
			}
		if (change)
		{
			matrix Q(n,n), v(n,1);
			for (int i = 0; i<n; ++i)
				for (int j = 0; j<=i; ++j)
					Q(i, j) = l(i);
			Q = Q * D;
			v = Q[0] / norm(Q[0]);
			D.set_col(v, 0);  //dowiedzie� sie o co cho
			for (int i = 1; i<n; ++i)  //tu te�
			{
				matrix temp(n, 1);
				for (int j = 0; j < n; ++j)
					temp = temp + (trans(Q[i]) * D[j])*D[j];  //sprawdzi to
				v = (Q[i] - temp) / norm(Q[i] - temp);
				D.set_col(v, i);
			}
			s = s0;
			l = matrix(n, 1);
			p = matrix(n, 1);
		}
		plik << X.x(0, 0) << ";" << X.x(1, 0) << endl;
		//cout << X.x(0, 0) << ";" << X.x(1, 0) << endl;
		double max_s = abs(s(0));
		for (int i = 1; i < n; ++i)
			if (max_s < abs(s(i)))
				max_s = abs(s(i));
		if (solution::f_calls>Nmax  || max_s < epsilon) 
			return X;
	}
}

//solution Rosen(matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix* ud, matrix* ad)
//{
//	solution X(x0), Xt;
//	int n = get_dim(X);
//	matrix l(n, 1), p(n, 1), s(s0), D = ident_mat(n);
//	X.fit_fun(ud, ad);
//	while (true)
//	{
//		for (int i = 0; i < n; ++i)
//		{
//			Xt.x = X.x + s(i) * D[i];
//			Xt.fit_fun(ud, ad);
//			if (Xt.y < X.y)
//			{
//				X = Xt;
//				l(i) += s(i);
//				s(i) *= alpha;
//			}
//			else
//			{
//				++p(i);
//				s(i) *= -beta;
//			}
//		}
//#if LAB_NO==3 && LAB_PART==2
//		(*ud).add_row(trans(X.x));
//#endif
//		bool change = true;
//		for (int i = 0; i < n; ++i)
//			if (l(i) == 0 || p(i) == 0)
//			{
//				change = false;
//				break;
//			}
//		if (change)
//		{
//			matrix Q(n, n), v(n, 1);
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j <= i; ++j)
//					Q(i, j) = l(i);
//			Q = D * Q;
//			v = Q[0] / norm(Q[0]);
//			D.set_col(v, 0);
//			for (int i = 1; i < n; ++i)
//			{
//				matrix temp(n, 1);
//				for (int j = 0; j < i; ++j)
//					temp = temp + trans(Q[i]) * D[j] * D[j];
//				v = (Q[i] - temp) / norm(Q[i] - temp);
//				D.set_col(v, i);
//			}
//			s = s0;
//			l = matrix(n, 1);
//			p = matrix(n, 1);
//		}
//		double max_s = abs(s(0));
//		for (int i = 1; i < n; ++i)
//			if (max_s < abs(s(i)))
//				max_s = abs(s(i));
//		if (max_s<epsilon || solution::f_calls>Nmax)
//			return X;
//	}
//}
#endif
#if LAB_NO>3
solution pen(matrix x0, double c0, double dc, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	double alpha = 1, beta = 0.5, gamma = 2, delta = 0.5, s = 0.5;
	solution X(???), X1;
	matrix c(2, new double[2]{ c0,dc });
	while (true)
	{
		X1 = sym_NM(X.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud, &c);
		if (???)
			return X1;
		???
		???
	}
}

solution sym_NM(matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	int n = get_len(x0);
	matrix D = ident_mat(n);
	int N = n + 1;
	solution *S = new solution[N];
	S[0].x = ???;
	S[0].fit_fun(ud, ad);
	for (int i = 1; i < N; ++i)
	{
		S[i].x = ???
		S[i].fit_fun(ud, ad);
	}
	solution PR, PE, PN;
	matrix pc;
	int i_min, i_max;
	while (true)
	{
		i_min = i_max = 0;
		for (int i = 1; i < N; ++i)
		{
			if (???)
				i_min = i;
			if (???)
				i_max = i;
		}
		pc = matrix(???);
		for (int i = 0; i < N; ++i)
			if (???)
				???
		pc = pc / (???);
		PR.x = ???
		PR.fit_fun(ud, ad);
		if (???)
			S[i_max] = ???
		else if (???)
		{
			PE.x = ???
			PE.fit_fun(ud, ad);
			if (???)
				S[i_max] = ???
			else
				S[i_max] = ???
		}
		else
		{
			PN.x = ???
			PN.fit_fun(ud, ad);
			if (???)
				S[i_max] = ???
			else
			{
				for (int i = 0; i < N; ++i)
					if (???)
					{
						S[i].x = ???
						S[i].fit_fun(ud, ad);
					}
			}
		}
		double max_s = ???
		for (int i = 1; i < N; ++i)
			if (max_s < ???)
				max_s = ???
		if (???)
			return S[i_min];
	}
}
#endif
#if LAB_NO>4
solution SD(matrix x0, double h0, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	int n = get_len(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n, 1), *P = new matrix[2];
	solution h;
	double *ab;
	while (true)
	{
		X.grad();
		d = ???
		if (h0<0)
		{
			P[0] = ???
			P[1] = ???
			ab = ???
			h = ???
			X1.x = ???
		}
		else
			X1.x = ???
#if LAB_NO==5 && LAB_PART==2
		???
#endif
		if (???)
		{
			X1.fit_fun(ud, ad);
			return X1;
		}
		???
	}
}

solution CG(matrix x0, double h0, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	int n = get_len(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n, 1), *P = new matrix[2];
	solution h;
	double *ab, beta;
	X.grad();
	d = ???
	while (true)
	{
		if (h0<0)
		{
			P[0] = ???
			P[1] = ???
			ab = ???
			h = ???
			X1.x = ???
		}
		else
			X1.x = ???
#if LAB_NO==5 && LAB_PART==2
		???
#endif
		if (???)
		{
			X1.fit_fun(ud);
			return X1;
		}
		X1.grad();
		beta = ???
		d = ???
		???
	}
}

solution Newton(matrix x0, double h0, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	int n = get_len(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n, 1), *P = new matrix[2];
	solution h;
	double *ab;
	while (true)
	{
		X.grad();
		X.hess();
		d = ???
		if (h0<0)
		{
			P[0] = ???
			P[1] = ???
			ab = ???
			h = ???
			X1.x = ???
		}
		else
			X1.x = ???
#if LAB_NO==5 && LAB_PART==2
		???
#endif
		if (???)
		{
			X1.fit_fun(ud);
			return X1;
		}
		???
	}
}

solution golden(double a, double b, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	double alfa = ???
	solution A, B, C, D;
	A.x = a;
	B.x = b;
	C.x = ???
	C.fit_fun(ud, ad);
	D.x = ???
	D.fit_fun(ud, ad);
	while (true)
	{
		if (???)
		{
			???
			???
			???
			C.fit_fun(ud, ad);
		}
		else
		{
			???
			???
			???
			D.fit_fun(ud, ad);
		}
		if (???)
		{
			A.x = ???
			A.fit_fun(ud, ad);
			return A;
		}
	}
}

#endif
#if LAB_NO>5
solution Powell(matrix x0, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	int n = get_len(x0);
	matrix D = ident_mat(n), *A = new matrix[2];
	solution X, P, h;
	X.x = x0;
	double *ab;
	while (true)
	{
		P = ???
		for (int i = 0; ???; ++i)
		{
			A[0] = ???
			A[1] = ???
			ab = ???
			h = ???
			P.x = ???
		}
		if (???)
		{
			P.fit_fun(ud);
			return P;
		}
		for (int i = 0; i < n - 1; ++i)
			D.set_col(???);
		D.set_col(???);
		A[0] = ???
		A[1] = ???
		ab = ???
		h = ???
		X.x = ???
	}
}
#endif
#if LAB_NO>6
solution EA(int N, matrix limits, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix *ud, matrix *ad)
{
	solution *P = new solution[mi + lambda];
	solution *Pm = new solution[mi];
	random_device rd;
	default_random_engine gen;
	gen.seed(static_cast<unsigned int>(chrono::system_clock::now().time_since_epoch().count()));
	normal_distribution<double> distr(0.0, 1.0);
	matrix IFF(mi, 1), temp(N, 2);
	double r, s, s_IFF;
	double tau = pow(2 * N, -0.5), tau1 = pow(2 * pow(N, 0.5), -0.5);
	int j_min;
	for (int i = 0; i < mi; ++i)
	{
		P[i].x = matrix(N, 2);
		for (int j = 0; j < N; ++j)
		{
			P[i].x(j, 0) = (limits(j, 1) - limits(j, 0))*rand_mat(1, 1)() + limits(j, 0);
			P[i].x(j, 1) = sigma0(j);
		}
		P[i].fit_fun(ud, ad);
		if (P[i].y < epsilon)
			return P[i];
	}
	while (true)
	{
		s_IFF = 0;
		for (int i = 0; ???; ++i)
		{
			IFF(i) = ???
			s_IFF += ???
		}
		for (int i = 0; ???; ++i)
		{
			r = ???
			s = ???
			for (int j = 0; ???; ++j)
			{
				s += ???
				if (???)
				{
					P[mi + i] = ???
					break;
				}
			}
		}
		for (int i = 0; ???; ++i)
		{
			r = ???
			for (int j = 0; ???; ++j)
			{
				???
				???
			}
		}
		for (int i = 0; ???; i += 2)
		{
			r = ???
			temp = P[mi + i].x;
			P[mi + i].x = ???
			P[mi + i + 1].x = ???
		}
		for (int i = 0; ???; ++i)
		{
			P[mi + i].fit_fun(ud, ad);
			if (???)
				return P[mi + i];
		}
		for (int i = 0; ???; ++i)
		{
			j_min = 0;
			for (int j = 1; ???; ++j)
				if (???)
					j_min = j;
			Pm[i] = ???
			P[j_min].y = ???
		}
		for (int i = 0; i < mi; ++i)
			P[i] = Pm[i];
		if (???)
			return P[0];
	}
}
#endif
