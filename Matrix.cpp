#include "stdafx.h"
#include "Matrix.h"
#include <time.h>


Matrix::Matrix()
{
}


Matrix::~Matrix()
{
}

void Matrix::show() {//打印矩阵
	int m = A.size();
	int n = A[0].size();
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			cout << A[i][j] << " ";
		}
		cout << "\n";
	}
}

void Matrix::show_diag_dual() {//打印矩阵对角以及次对角块
	int m = A.size();
	int n = A[0].size();
	if (m < n) { n = m; }
	for (int i = 1; i <n; i++) {
		if (A[i][i - 1] != 0) {
			cut(i, i + 1, i, i + 1).show();
			cout << "\n";
			i = i + 1;
			continue;
		}
		if (A[i][i - 1] == 0) {
			cut(i, i, i, i).show();
			cout << "\n";

		}
	}

	if (A[n - 1][n - 2] == 0) {
		cut(n, n, n, n).show();
	}
}

void Matrix::show_diag() {//打印矩阵对角
	int m = A.size();
	int n = A[0].size();
	if (m < n) { n = m; }
	for (int i = 1; i < n; i++) {
		cout << A[i][i] << "\t";
	}
}

double Matrix::max() {//返回第一列模最大的元素
	int n = A.size();
	int m = 0;
	for (int j = 0; j < n; j++) {
		if (abs(A[j][0])>abs(A[m][0])) {
			m = j;
		}
	}
	return A[m][0];
}

double Matrix::tr() {
	int n = A.size();
	int m = A[0].size();
	if (m < n) { n = m; }
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += A[i][i];
	}
	return sum;
}

void Matrix::construct(int n, double a1, double a2, double a3) {//--------------矩阵的构造函数
	A.clear();
	vector<double> a;
	for (int j = 0; j < n; j++) {
		a.push_back(0);
	}
	a[0] = a2;
	a[1] = a3;
	A.push_back(a);
	for (int i = 0; i < n - 2; i++) {
		for (int j = 0; j < n; j++) {
			a[j] = 0;
		}
		a[i] = a1;
		a[i + 1] = a2;
		a[i + 2] = a3;
		A.push_back(a);
	}
	for (int j = 0; j < n; j++) {
		a[j] = 0;
	}
	a[n - 2] = a1;
	a[n - 1] = a2;
	A.push_back(a);
}

void Matrix::construct_companion_matrix_of_polynomial(double *a, int n) {//友方阵构造函数
	A.clear();
	vector<double> u;
	for (int i = 0; i < n; i++) { u.push_back(-a[i]); }
	vector<double> v;
	for (int i = 0; i < n; i++) { v.push_back(0); }
	A.push_back(u);
	for (int i = 0; i < n - 1; i++) {
		v[i] = 1;
		A.push_back(v);
		v[i] = 0;
	}
}

void Matrix::construct2_3(double x) {
	zeros(4, 4);
	double a[4][4] = { { 9.1,3.0,2.6,4.0 },{ 4.2,5.3,4.7,1.6 },{ 3.2,1.7,9.4,x },{ 6.1,4.9,3.5,6.2 } };
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			A[i][j] = a[i][j];
		}
	}
}

void Matrix::zeros(int m, int n) {//全0矩阵构造函数
	A.clear();
	vector<double> v;
	for (int i = 0; i < m; i++) {
		v.clear();
		for (int j = 0; j < n; j++) { v.push_back(0); }
		A.push_back(v);
	}
}
void Matrix::ones(int m, int n) {//全1矩阵构造函数
	A.clear();
	vector<double> v;
	for (int i = 0; i < m; i++) {
		v.clear();
		for (int j = 0; j < n; j++) { v.push_back(1); }
		A.push_back(v);
	}
}
void Matrix::identity(int n) {//n阶单位矩阵构造函数
	zeros(n, n);
	for (int i = 0; i < n; i++) {
		A[i][i] = 1;
	}
}

double Matrix::power_iteration() {//幂法迭代
	double lambda = 1, miu = 0, err = 0.00001;
	int n = A[0].size();
	Matrix u;
	u.ones(n, 1);
	Matrix y = u;
	while (abs(lambda - miu)>err)
	{
		lambda = miu;
		y = Multiply(*this, u);
		miu = y.max();
		u = Scalar(1 / miu, y);

	}

	return lambda;
}


void Matrix::Givens_down(int i, int j, double c, double s) {//左乘效果为A[i]=-sA[i]+cA[j],A[j]=cA[i]+sA[j]
	int n = A[0].size();
	double tempi, tempj;
	for (int k = 0; k < n; k++) {
		tempi = A[i][k];
		tempj = A[j][k];
		A[i][k] = -s*tempi + c*tempj;
		A[j][k] = s*tempj + c*tempi;
	}
}

void Matrix::Givens_up(int i, int j, double c, double s) {//左乘效果为A[j]=-sA[i]+cA[j],A[i]=cA[i]+sA[j]
	int n = A[0].size();
	double tempi, tempj;
	for (int k = 0; k < n; k++) {
		tempi = A[i][k];
		tempj = A[j][k];
		A[j][k] = -s*tempi + c*tempj;
		A[i][k] = s*tempj + c*tempi;
	}
}

void Matrix::Givens_right(int i, int j, double c, double s) {//右乘效果为A[][i]=-sA[][i]+cA[][j],A[][j]=cA[][i]+sA[][j]
	int n = A[0].size();
	double tempi, tempj;
	for (int k = 0; k < n; k++) {
		tempi = A[k][i];
		tempj = A[k][j];
		A[k][i] = -s*tempi + c*tempj;
		A[k][j] = s*tempj + c*tempi;
	}
}
void Matrix::Givens_left(int i, int j, double c, double s) {//右乘效果为A[][j]=sA[][i]+cA[][j],A[][i]=cA[][i]-sA[][j]
	int n = A[0].size();
	double tempi, tempj;
	for (int k = 0; k < n; k++) {
		tempi = A[k][i];
		tempj = A[k][j];
		A[k][j] = +s*tempi + c*tempj;
		A[k][i] = -s*tempj + c*tempi;
	}
}

Matrix Matrix::householder() {//householder变换(已检验完全准确)
	int n = A.size();
	if (A[0].size() > 1) {
		cout << "Illegal calculation";
		exit(0);
	}

	vector<vector<double>> x;
	double yita = abs(max()), sigma = 0, beta, alpha;
	for (int i = 0; i < n; i++) {
		A[i][0] = A[i][0] / yita;
	}
	for (int i = 1; i < n; i++) {
		sigma += A[i][0] * A[i][0];
	}
	x = A;
	if (sigma == 0) {
		x[0][0] = 0;
	}
	else {
		alpha = sqrt(A[0][0] * A[0][0] + sigma);
		if (A[0][0] <= 0) {
			x[0][0] = A[0][0] - alpha;
		}
		else {
			x[0][0] = -sigma / (A[0][0] + alpha);
		}
		beta = 2 * x[0][0] * x[0][0] / (sigma + x[0][0] * x[0][0]);
		for (int i = 1; i < n; i++) {
			x[i][0] = x[i][0] / x[0][0];
		}
		x[0][0] = beta;
	}
	Matrix C;
	C.A = x;
	return C;
}



Matrix Matrix::jacobi_iteration_with_threshold() {
	int n = A.size();
	vector<vector<double>> B;
	double tau, threshold = Square_sum(*this) / (n*n), t, c, s;
	int flag, count = 0;
	while (threshold > 0.000000001) {
		flag = 0;
		for (int p = 0; p < n; p++) {
			for (int q = 0; q < p; q++) {
				if (abs(A[p][q]) > threshold) {
					B = A;
					flag = flag + 1;
					count = count + 1;
					tau = (A[q][q] - A[p][p]) / (2 * A[p][q]);
					t = sign(tau) / (abs(tau) + sqrt(1 + tau*tau));
					c = 1 / sqrt(1 + t*t);
					s = t*c;
					for (int i = 0; i < n; i++) {
						if (i != p&&i != q) {
							B[i][p] = c*A[i][p] - s*A[i][q];
							B[p][i] = c*A[i][p] - s*A[i][q];
							B[i][q] = s*A[i][p] + c*A[i][q];
							B[q][i] = s*A[i][p] + c*A[i][q];
						}
					}
					B[p][p] = c*c*A[p][p] - 2 * s*c*A[p][q] + s*s*A[q][q];
					B[q][q] = s*s*A[p][p] + 2 * s*c*A[p][q] + c*c*A[q][q];
					B[p][q] = (c*c - s*s)*A[p][q] + s*c*(A[p][p] - A[q][q]);
					B[q][p] = (c*c - s*s)*A[p][q] + s*c*(A[p][p] - A[q][q]);
					A = B;
				}
			}
		}
		if (flag == 0) {
			threshold = threshold / n;
		}
	}

	cout << "iteration times: " << count << endl;
	return *this;
}







Matrix Matrix::up_Hessenburg_decomposition() {//上hessenburg分解
	int n = A.size();
	if (n < 3) {
		cout << "nonsense calculation";
		exit(0);
	}
	Matrix V, temp;
	double beta;
	for (int k = 1; k < n - 1; k++) {
		//house (A[k+1:n,k])
		V = cut(k + 1, n, k, k).householder();
		beta = V.A[0][0];
		V.A[0][0] = 1;
		temp = cut(k + 1, n, k, n);
		//A[k+1:n,k:n] = A[k+1:n,k:n]-beta V V^T A[k+1:n,k:n]
		paste(Minus(temp, Scalar(beta, Multiply(V, Multiply(V.T(), temp)))), k + 1, n, k, n);

		temp = cut(1, n, k + 1, n);
		//A[1:n,k+1:n]=A[1:n,k+1:n](I-beta V V^T)
		paste(Minus(temp, Multiply(Scalar(beta, Multiply(temp, V)), V.T())), 1, n, k + 1, n);
	}
	for (int i = 2; i < n; i++) {
		for (int k = 0; k < i - 1; k++) {
			A[i][k] = 0;
		}
	}
}

void Matrix::check_real_or_complex(int i) {//判断(i,i+1)的2*2下次对角元非0矩阵特征值是否是实的
	Matrix temp = cut(i + 1, i + 2, i + 1, i + 2);
	double s, t, delta;
	s = temp.A[0][0] + temp.A[1][1];
	t = temp.A[0][0] * temp.A[1][1] - temp.A[0][1] * temp.A[1][0];
	delta = s*s - 4 * t;
	if (abs(delta) < 0.000001) {
		delta = 0;
	}
	if (delta < 0) { return; }
	temp.direct_QR_iteration();
	paste(temp, i + 1, i + 2, i + 1, i + 2);
}


Matrix Matrix::implicilt_QR_iteration() {//隐式QR算法
	int n = A.size();

	if (n < 3) {
		cout << "nonsense calculation";
		exit(0);
	}
	Matrix V, temp;
	double beta;
	for (int k = 1; k < n - 1; k++) {
		//house (A[k+1:n,k])
		V = cut(k + 1, n, k, k).householder();

		beta = V.A[0][0];
		V.A[0][0] = 1;
		temp = cut(k + 1, n, k, n);
		//A[k+1:n,k:n] = A[k+1:n,k:n]-beta V V^T A[k+1:n,k:n]
		paste(Minus(temp, Scalar(beta, Multiply(V, Multiply(V.T(), temp)))), k + 1, n, k, n);

		temp = cut(1, n, k + 1, n);
		//A[1:n,k+1:n]=A[1:n,k+1:n](I-beta V V^T)
		paste(Minus(temp, Multiply(Scalar(beta, Multiply(temp, V)), V.T())), 1, n, k + 1, n);
	}
	for (int i = 2; i < n; i++) {
		for (int k = 0; k < i - 1; k++) {
			A[i][k] = 0;
		}
	}
	int m = 0, l = 0, flag;
	Matrix H22, H12, H23, P;
	while (m < n - 1)//m为真实的行
	{
		for (int i = 1; i < n; i++) { if (abs(A[i][i - 1]) <= 0.000001*(abs(A[i][i]) + abs(A[i - 1][i - 1]))) { A[i][i - 1] = 0; } }
		flag = 0;
		//寻找最大的满足条件的m
		for (int i = n - 1; i > 0; i--) {
			if (A[i][i - 1] != 0 && flag == 0) {
				flag = 1;
				m = m + 1;
				continue;
			}
			if (A[i][i - 1] != 0 && flag == 1) {
				m = m - 1;
				//寻找最小的l
				for (int j = i - 1; j > 0; j--) {
					if (A[j][j - 1] == 0) {
						l = j;
					}
				}
				break;
			}
			if (A[i][i - 1] == 0 && flag == 1) {
				flag = 0;
				check_real_or_complex(i);//判断(i+1,i+2)的2*2下次对角元非0矩阵特征值是否是实的
			}
			m = m + 1;
		}
		if (m >= n - 2) {
			check_real_or_complex(0);
			for (int i = 1; i < n; i++) { if (abs(A[i][i - 1]) <= 0.0001*(abs(A[i][i]) + abs(A[i - 1][i - 1]))) { A[i][i - 1] = 0; } }
			break;
		}
		H22 = cut(l + 1, n - m, l + 1, n - m);//(n-l-m)*(n-l-m)

		P = H22.double_shift_QR_iteration();//(n-l-m)*(n-l-m)
		paste(H22, l + 1, n - m, l + 1, n - m);
		if (l != 0) {
			H12 = cut(1, l, l + 1, n - m);//l*(n-l-m)
			H12 = Multiply(H12, P);
			paste(H12, 1, l, l + 1, n - m);
		}
		if (m != 0) {
			H23 = cut(l + 1, n - m, n - m + 1, n);//(n-l-m)*m
			H23 = Multiply(P.T(), H23);
			paste(H23, l + 1, n - m, n - m + 1, n);
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (abs(A[i][j]) < 0.0000001) { A[i][j] = 0; }
		}
	}
	return *this;
}


Matrix Matrix::QR_decomposition() {//----------------------------------------------QR分解(利用算法3.3.1)
	int n = A.size();
	int m = A[0].size();//m>=n
	Matrix temp, tempQ, v, Q;
	Q.identity(n);
	double beta;
	for (int j = 1; j < n + 1; j++) {
		if (j < m) {
			temp = cut(j, m, j, n);
			v = cut(j, m, j, j).householder();
			beta = v.A[0][0];
			v.A[0][0] = 1;
			temp = Minus(temp, Multiply(Scalar(beta, v), Multiply(v.T(), temp)));
			paste(temp, j, m, j, n);
			temp.identity(n);
			tempQ.identity(m - j + 1);
			tempQ = Minus(tempQ, Multiply(Scalar(beta, v), v.T()));
			temp.paste(tempQ, j, m, j, n);
			Q = Multiply(Q, temp);
		}
	}
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < j; i++) {
			A[j][i] = 0;
		}
	}
	return Q;
}


Matrix Matrix::direct_QR_iteration() {//----------------------------------------------直接QR算法
	int n = A.size(), count = 0;
	Matrix Q, recordQ;
	recordQ.identity(n);
	while (Lower_diag_square_sum(*this)>1e-12)
	{
		Q = QR_decomposition();
		recordQ = Multiply(recordQ, Q);
		A = Multiply(*this, Q).A;
		//count = count + 1;
	}
	//cout << "iteration times:" << count << endl;
	return recordQ;
}

Matrix Matrix::double_shift_QR_iteration() {//双重步位移QR迭代
	int n = A.size();
	Matrix P, P_temp, temp, V, H;
	P.identity(n);
	temp.zeros(3, 1);
	int m = n - 1, q, r;
	double s, t, x, y, z, beta;
	s = A[m - 1][m - 1] + A[n - 1][n - 1];
	t = A[m - 1][m - 1] * A[n - 1][n - 1] - A[m - 1][n - 1] * A[n - 1][m - 1];
	x = A[0][0] * A[0][0] + A[0][1] * A[1][0] - s*A[0][0] + t;
	y = A[1][0] * (A[0][0] + A[1][1] - s);
	z = A[1][0] * A[2][1];
	for (int k = 0; k < n - 2; k++) {
		//[v,beta]=house([x,y,z]^T)
		temp.A[0][0] = x;
		temp.A[1][0] = y;
		temp.A[2][0] = z;

		V = temp.householder();

		beta = V.A[0][0];
		V.A[0][0] = 1;
		//q=max{1,k}
		q = k;
		if (1 > k) { q = 1; }
		H = cut(k + 1, k + 3, q, n);
		//H(k + 1, k + 3, q, n)=(I-beta V V^T )H(k + 1, k + 3, q, n))
		H = Minus(H, Multiply(V, Scalar(beta, Multiply(V.T(), H))));

		paste(H, k + 1, k + 3, q, n);
		r = k + 4;
		if (k + 4 > n) { r = n; }
		//H(1, r, k+1, k+3)=H(1, r, k+1, k+3)(I-beta V V^T )
		H = cut(1, r, k + 1, k + 3);
		H = Minus(H, Multiply(Multiply(H, V), Scalar(beta, V.T())));
		paste(H, 1, r, k + 1, k + 3);
		P_temp = P.cut(1, n, k + 1, k + 3);
		P_temp = Minus(P_temp, Multiply(Multiply(P_temp, V), Scalar(beta, V.T())));
		P.paste(P_temp, 1, n, k + 1, k + 3);
		x = A[k + 1][k];
		y = A[k + 2][k];
		if (k < n - 3) { z = A[k + 3][k]; }
	}
	//[v,beta]=house([x,y,z]^T)
	temp.zeros(2, 1);
	temp.A[0][0] = x;
	temp.A[1][0] = y;
	V = temp.householder();
	beta = V.A[0][0];
	V.A[0][0] = 1;

	H = cut(n - 1, n, n - 2, n);
	H = Minus(H, Multiply(V, Scalar(beta, Multiply(V.T(), H))));
	paste(H, n - 1, n, n - 2, n);
	H = cut(1, n, n - 1, n);
	H = Minus(H, Multiply(Multiply(H, V), Scalar(beta, V.T())));
	paste(H, 1, n, n - 1, n);
	//记录P
	P_temp = P.cut(1, n, n - 1, n);
	P_temp = Minus(P_temp, Multiply(Multiply(P_temp, V), Scalar(beta, V.T())));
	P.paste(P_temp, 1, n, n - 1, n);
	return P;
}


Matrix Matrix::dual_diag_with_householder(Matrix &U, Matrix &V) {//------------------------------------------二对角化(householder变换)
	int m = A.size();
	int n = A[0].size();
	double beta;
	vector<double> b, c;
	Matrix uT, u, v, temp, tempU, tempV;
	U.identity(m);
	V.identity(n);
	for (int k = 1; k < n + 1; k++) {
		v = cut(k, m, k, k).householder();
		beta = v.A[0][0];
		v.A[0][0] = 1;
		temp = cut(k, m, k, n);
		uT = Multiply(Scalar(beta, v.T()), temp);
		paste(Minus(temp, Multiply(v, uT)), k, m, k, n);

		thread r(stimulate_householderU, k, m, v, beta, ref(U));


		/*
		tempU.identity(m);
		tempU.paste(Minus(tempU.cut(k, m, k, m), Multiply(v, Scalar(beta, v.T()))),k,m,k,m);
		U = Multiply(tempU, U);
		*/


		if (k < n - 1) {
			v = cut(k, k, k + 1, n).T().householder();
			beta = v.A[0][0];
			v.A[0][0] = 1;

			thread q(stimulate_householderV, k, n, v, beta, ref(V));
			temp = cut(k, m, k + 1, n);
			u = Multiply(temp, Scalar(beta, v));
			paste(Minus(temp, Multiply(u, v.T())), k, m, k + 1, n);
			q.join();
			/*
			tempV.identity(n);
			tempV.paste(Minus(tempV.cut(k+1, n, k+1, n), Multiply(v, Scalar(beta, v.T()))), k+1, n, k+1, n);
			V = Multiply(V, tempV);
			*/

		}
		r.join();

	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (i > j || i < j - 1) {
				A[i][j] = 0;
			}

		}
	}
	return *this;
}

Matrix Matrix::all_row_zeros(int i) {
	int n = A[0].size();
	int m = A.size();
	double a, b, c, s, t;
	Matrix G;
	G.identity(m);
	if (A[i][i] != 0 || A[i][i + 1] == 0) {
		cout << "error input";
		exit(0);
	}
	for (int j = i + 1; j < n; j++) {
		a = A[i][j];
		b = A[j][j];
		if (a == 0) {
			break;
		}
		if (abs(b) > abs(a)) {
			t = a / b;
			s = 1 / sqrt(1 + t*t);
			c = s*t;
		}
		else {
			t = b / a;
			c = 1 / sqrt(1 + t*t);
			s = c*t;
		}
		Givens_down(i, j, c, s);
		G.Givens_down(i, j, c, s);
	}
	return G;
}

Matrix Matrix::small_to_zeros(double thresh) {//------------------------------------------------------把小的元素打成0
	int n = A.size();
	int m = A[0].size();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			if (abs(A[i][j]) < thresh) { A[i][j] = 0; }
		}
	}
	return *this;
}


Matrix Matrix::SVD_iteration(Matrix &P, Matrix &Q) {
	int n = A.size();
	P.identity(n);
	Q.identity(n);
	double alpha = A[n - 1][n - 1] * A[n - 1][n - 1] + A[n - 2][n - 1] * A[n - 2][n - 1], delta, beta, miu, y, z, a, b, c, i, j, r, t, s;
	double dk, gk, dk1, gk1;
	delta = (A[n - 2][n - 2] * A[n - 2][n - 2] + A[n - 3][n - 2] * A[n - 3][n - 2] - alpha) / 2;
	beta = A[n - 2][n - 2] * A[n - 2][n - 1];
	miu = alpha - beta*beta / (delta + sign(delta)*sqrt(delta*delta + beta*beta));
	y = A[0][0] * A[0][0] - miu;
	z = A[0][0] * A[0][1];
	for (int k = 0; k < n - 1; k++) {
		dk = A[k][k];
		dk1 = A[k + 1][k + 1];
		gk = A[k][k + 1];
		a = y;
		b = z;
		if (b == 0) {
			break;
		}
		if (abs(b) > abs(a)) {
			t = -a / b;
			s = 1 / sqrt(1 + t*t);
			r = -b / s;
			c = s*t;
		}
		else {
			t = -b / a;
			c = 1 / sqrt(1 + t*t);
			r = a / c;
			s = c*t;
		}


		Q.Givens_left(k, k + 1, c, s);
		y = dk*c - gk*s;
		z = -dk1*s;
		gk = s*dk + c*gk;
		dk1 = c*dk1;
		if (k > 0) {
			A[k - 1][k] = r;
		}



		a = y;
		b = z;
		a = y;
		b = z;
		if (b == 0) {
			break;
		}
		if (abs(b) > abs(a)) {
			t = -a / b;
			s = 1 / sqrt(1 + t*t);
			r = -b / s;
			c = s*t;
		}
		else {
			t = -b / a;
			c = 1 / sqrt(1 + t*t);
			r = a / c;
			s = c*t;
		}


		P.Givens_left(k, k + 1, c, s);
		dk = r;
		if (k < n - 2) {
			gk1 = A[k + 1][k + 2];
			y = c*gk - s*dk1;
			z = -s*gk1;
			dk1 = s*gk + c*dk1;
			gk1 = c*gk1;
			A[k + 1][k + 2] = gk1;
		}
		else {
			t = gk;
			gk = c*gk - s*dk1;
			dk1 = s*t + c*dk1;
		}

		A[k][k] = dk;
		A[k + 1][k + 1] = dk1;
		A[k][k + 1] = gk;
	}
	/*Matrix B;
	B.A = A;
	show();
	P.show();
	Q.show();
	Multiply(Multiply(P,B), Q.T()).small_to_zeros(1e-6).show();*/
	return *this;
}

int Matrix::check_zeros_of_diag(int i, int j) {//--------------------------------------检测i到j个对角元中有没有零
	int k = -1;
	for (int m = i - 1; m > j - 1; m--) {
		if (A[m][m] == 0) {
			return m;
		}
	}
	return k;
}
Matrix Matrix::SVD(Matrix &U, Matrix &V) {
	int m = A.size(), n = A[0].size(), q = 1, i = n - 1, j, k, count = 0;
	Matrix G, tempU, tempV, B22, P, Q;
	dual_diag_with_householder(U, V);//B=U^T A V
	V = V.T();
	double Norm = sqrt(Square_sum(*this)), epsilon = 1e-7;
	small_to_zeros((epsilon)*Norm);
	while (q < n) {
		count = count + 1;
		if (count > 100 * m*n) {
			cout << "WARNING!!:BAD Matrix, return after " << 100 * m*n << " iteration" << endl;
			break;
		}
		while (i > 0 && A[i - 1][i] == 0) {
			i = i - 1;
			q = q + 1;
		}
		if (q == n) {
			break;
		}
		j = i;
		while (j > 0 && A[j - 1][j] != 0) {
			j = j - 1;
		}
		k = check_zeros_of_diag(i, j);//找倒着数第一个0
		if (k > -1) {
			G = all_row_zeros(k);
			U = Multiply(U, G);
			continue;
		}

		B22 = cut(j + 1, i + 1, j + 1, i + 1);
		//B22.show();
		if (i - j <= 1) {
			P = Multiply(B22, B22.T()).direct_QR_iteration().T();
			Q = Multiply(B22.T(), B22).direct_QR_iteration().T();
			B22 = Multiply(Multiply(P, B22), Q.T());



		}
		else {
			B22.SVD_iteration(P, Q);
			P = P.T();
			Q = Q.T();


		}

		/*
		tempU.identity(m);
		tempU.paste(P, j + 1, i + 1, j + 1, i + 1);
		U = Multiply(tempU, U);
		tempV.identity(n);
		tempV.paste(Q, j + 1, i + 1, j + 1, i + 1);
		V = Multiply(tempV, V);
		*/
		thread t(stimulate, m, i, j, P, ref(U));
		thread s(stimulate, n, i, j, Q, ref(V));
		paste(B22, j + 1, i + 1, j + 1, i + 1);
		small_to_zeros((epsilon)*Norm);
		t.join();
		s.join();
	}
	return *this;
}
Matrix Matrix::SVD_direct_QR(Matrix &U, Matrix &V) {//---------------------------------------奇异值分解用直接QR方法简约版本
	int m = A.size(), n = A[0].size(), q = 1, i = n - 1, j, k;
	Matrix G, tempU, tempV, B22, P, Q;
	double epsilon = 1e-7;
	dual_diag_with_householder(U, V);//B=U^T A V
	V = V.T();

	while (q < n) {

		small_to_zeros(epsilon);
		while (i > 0 && A[i - 1][i] == 0) {
			i = i - 1;
			q = q + 1;
		}
		if (q == n) {
			break;
		}
		j = i;
		while (j > 0 && A[j - 1][j] != 0) {
			j = j - 1;
		}
		k = check_zeros_of_diag(i, j);//找倒着数第一个0
		if (k > -1) {
			G = all_row_zeros(k);
			U = Multiply(U, G);
			continue;
		}

		B22 = cut(j + 1, i + 1, j + 1, i + 1);
		if (j - i == 1) {
			P = Multiply(B22, B22.T()).direct_QR_iteration().T();
			Q = Multiply(B22.T(), B22).direct_QR_iteration().T();
			B22 = Multiply(Multiply(P, B22), Q.T());
			B22.small_to_zeros(epsilon);
		}
		else {
			//B22.SVD_iteration(P, Q);
			P = Multiply(B22, B22.T()).direct_QR_iteration().T();
			Q = Multiply(B22.T(), B22).direct_QR_iteration().T();
			B22 = Multiply(Multiply(P, B22), Q.T());

			B22.small_to_zeros(epsilon);
		}
		paste(B22, j + 1, i + 1, j + 1, i + 1);
		tempU.identity(m);
		tempU.paste(P, j + 1, i + 1, j + 1, i + 1);
		U = Multiply(tempU, U);
		tempV.identity(n);
		tempV.paste(Q, j + 1, i + 1, j + 1, i + 1);
		V = Multiply(tempV, V);

	}
	return *this;
}

Matrix Matrix::cut(int j, int m, int k, int n) {//真实的行列，不需要-1
	Matrix C;
	vector<double> v;
	for (int i = j - 1; i < m; i++) {
		v.clear();
		for (int l = k - 1; l < n; l++) { v.push_back(A[i][l]); }
		C.A.push_back(v);
	}
	return C;
}
void Matrix::paste(Matrix C, int j, int m, int k, int n) {//真实的行列，不需要-1
	for (int i = j - 1; i < m; i++) {
		for (int l = k - 1; l < n; l++) { A[i][l] = C.A[i - j + 1][l - k + 1]; }
	}
}

Matrix Matrix::T() {//求转置
	int m = A.size();
	int n = A[0].size();
	vector<double>  v;
	vector<vector<double>> B;
	for (int k = 0; k < n; k++) {
		for (int i = 0; i < m; i++) { v.push_back(A[i][k]); }
		B.push_back(v);
		v.clear();
	}
	Matrix C;
	C.A = B;
	return C;
}

Matrix Matrix::small_to_zeros(double thresh) {//------------------------------------------------------把小的元素打成0
	int n = A.size();
	int m = A[0].size();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			if (abs(A[i][j]) < thresh) { A[i][j] = 0; }
		}
	}
	return *this;
}


Matrix Minus(Matrix A, Matrix B) {//两矩阵相减
	int m = A.A.size();
	int n = A.A[0].size();
	if (m != B.A.size() || n != B.A[0].size()) {
		cout << "Illegal calculation!";
		return A;
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) { A.A[i][j] -= B.A[i][j]; }
	}
	return A;
}

Matrix Plus(Matrix A, Matrix B) {//两矩阵相减
	int m = A.A.size();
	int n = A.A[0].size();
	if (m != B.A.size() || n != B.A[0].size()) {
		cout << "Illegal calculation!";
		return A;
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) { A.A[i][j] += B.A[i][j]; }
	}
	return A;
}

Matrix Multiply(Matrix A, Matrix B) {//两矩阵相乘
	int m = A.A.size();
	int n = A.A[0].size();
	double sum = 0;
	Matrix C;
	vector<double> v;
	if (n != B.A.size()) {
		cout << "Illegal calculation!";
		return A;
	}
	for (int i = 0; i < m; i++) {
		v.clear();
		for (int j = 0; j < B.A[0].size(); j++) {
			sum = 0;
			for (int k = 0; k < n; k++) { sum += A.A[i][k] * B.A[k][j]; }
			v.push_back(sum);
		}
		C.A.push_back(v);
	}
	return C;
}

Matrix Scalar(double beta, Matrix A) {//数乘
	int m = A.A.size();
	int n = A.A[0].size();
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) { A.A[i][j] = beta*A.A[i][j]; }
	}
	return A;
}

double Square_sum(Matrix &A) {
	int m = A.A.size();
	int n = A.A[0].size();
	double sum = 0;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) { sum = sum + A.A[i][j] * A.A[i][j]; }
	}
	return sum;
}

double Diag_square_sum(Matrix &A) {
	int m = A.A.size();
	int n = A.A[0].size();
	if (m < n) { n = m; }
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += A.A[i][i];
	}
	return sum;
}

double Not_diag_square_sum(Matrix &A) {
	return Square_sum(A) - Diag_square_sum(A);
}

double Lower_diag_square_sum(Matrix &A) {
	int m = A.A.size();
	int n = A.A[0].size();
	double sum = 0;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < i; j++) { sum = sum + A.A[i][j] * A.A[i][j]; }
	}
	return sum;
}

double sign(double delta) {
	if (delta > 0) {
		return 1.0;
	}
	if (delta == 0) {
		return 0;
	}
	return -1.0;
}

double bisection_for_symmetrtic_triple_diagonal(vector<double> x, vector<double> y, int k) {//--------------------二分法求(三对角阵)第k个特征值
	int s = 0, n = y.size(), count = 0;
	if (n != x.size()) {
		cout << "illegal calculation" << endl;
		return 0;
	}
	double sum = 0, q, miu = 0;
	for (int i = 0; i < n; i++) {
		sum = sum + abs(y[i]) + abs(x[i]);
		y[i] = y[i] * y[i];
	}
	double upper = sum, lower = -sum;
	while (upper - lower > 0.00000001) {
		q = x[0] - miu;
		s = 0;
		for (int i = 0; i < n; i++) {
			if (q < 0) { s = s + 1; }
			if (i < n - 1) {
				if (q == 0) {
					q = abs(y[i + 1])*0.000000001;
				}
				q = x[i + 1] - miu - y[i + 1] / q;
			}
		}
		if (s >= k) {
			upper = miu;
		}
		if (s < k) {
			lower = miu;
		}
		miu = (lower + upper) / 2;
		count = count + 1;
	}
	cout << "iteration times: " << count << endl;
	return miu;
}

void stimulate(int m, int i, int j, Matrix P, Matrix &U) {
	Matrix tempU;
	tempU.identity(m);
	tempU.paste(P, j + 1, i + 1, j + 1, i + 1);
	U = Multiply(tempU, U);
}
void stimulate_householderU(int k, int m, Matrix v, double beta, Matrix &U) {
	Matrix tempU;
	tempU.identity(m);
	tempU.paste(Minus(tempU.cut(k, m, k, m), Multiply(v, Scalar(beta, v.T()))), k, m, k, m);
	U = Multiply(tempU, U);
}

void stimulate_householderV(int k, int n, Matrix v, double beta, Matrix &V) {
	Matrix tempV;
	tempV.identity(n);
	tempV.paste(Minus(tempV.cut(k + 1, n, k + 1, n), Multiply(v, Scalar(beta, v.T()))), k + 1, n, k + 1, n);
	V = Multiply(V, tempV);
}