#pragma once

#include <iostream>
#include<thread>
#include<time.h>
#include <vector>

using namespace std;

class Matrix//-------------------------------------------------------------------定义矩阵类
{
public:
	Matrix();
	~Matrix();
public:
	void construct(int n, double a1, double a2, double a3);//---------------------矩阵的构造函数
	void construct2_3(double x);//------------------------------------------------矩阵的构造函数

	void construct_companion_matrix_of_polynomial(double *a, int n);//------------构造多项式的友方阵
	void zeros(int m, int n);//---------------------------------------------------构造m*n的全0矩阵
	void ones(int m, int n);//----------------------------------------------------构造m*n的全0矩阵
	void identity(int n);//-------------------------------------------------------构造n阶单位矩阵

	Matrix cut(int j, int m, int k, int n);//------------------------------------选出A(j:m,k:n)
	void paste(Matrix C, int j, int m, int k, int n);//--------------------------把C粘贴到A(j:m,k:n)
	Matrix T();//----------------------------------------------------------------求转置
	void show();//----------------------------------------------------------------打印矩阵
	void show_diag_dual();//------------------------------------------------------打印矩阵拟对角
	void show_diag();//-----------------------------------------------------------打印矩阵对角

	double max();//---------------------------------------------------------------求最大行
	double tr();//----------------------------------------------------------------计算迹

	void Givens_down(int i, int j, double c, double s);//----------------------------Givens变换,向下化零
	void Givens_up(int i, int j, double c, double s);//---------------------------Givens变换，向上化零
	void Givens_left(int i, int j, double c, double s);//-------------------------Givens变换,向左化零
	void Givens_right(int i, int j, double c, double s);//------------------------Givens变换，向右化零

	Matrix Matrix::jacobi_iteration_with_threshold();//-------------------------过关jacobi

	Matrix householder();//------------------------------------------------------householder变换

	Matrix direct_QR_iteration();//----------------------------------------------直接QR算法
	Matrix QR_decomposition();//-------------------------------------------------QR分解

	Matrix dual_diag_with_householder(Matrix &U, Matrix &V);//-----------------二对角化(householder变换)
	Matrix small_to_zeros(double thresh);//--------------------------------------把小的元素打成0
	Matrix all_row_zeros(int i);//-----------------------------------------------把对角元为0的整行变成0
	int check_zeros_of_diag(int i, int j);//--------------------------------------检测i到j个对角元中有没有零
	Matrix SVD(Matrix &U, Matrix &V);//----------------------------------------奇异值分解
	Matrix SVD_direct_QR(Matrix &U, Matrix &V);//------------------------------奇异值分解用直接QR方法
	Matrix SVD_iteration(Matrix &P, Matrix &Q);//------------------------------算法7.6.2

	double power_iteration();//---------------------------------------------------幂法
	Matrix householder();//------------------------------------------------------householder变换
	Matrix up_Hessenburg_decomposition();//--------------------------------------上Hessenburg分解
	Matrix implicilt_QR_iteration();//-------------------------------------------隐式QR算法
	Matrix double_shift_QR_iteration();//----------------------------------------双重步位移QR迭代
	void check_real_or_complex(int i);//------------------------------------------判断(i,i+1)的2*2下次对角元非0矩阵特征值是否是实的
	Matrix small_to_zeros(double thresh);//--------------------------------------把小的元素打成0

	vector<vector<double>> A;
};

Matrix Minus(Matrix A, Matrix B);//--------------------------------------------两矩阵相减
Matrix Plus(Matrix A, Matrix B);//---------------------------------------------两矩阵相加
Matrix Multiply(Matrix A, Matrix B);//-----------------------------------------两矩阵相乘
Matrix Scalar(double beta, Matrix A);//-----------------------------------------两矩阵相乘
double sign(double delta);//---------------------------------------------------------取符号
double Square_sum(Matrix &A);//--------------------------------------------------计算所有元素的平方和
double Diag_square_sum(Matrix &A);//---------------------------------------------计算对角元的平方和
double Not_diag_square_sum(Matrix &A);//-----------------------------------------计算对角元的平方和
double Lower_diag_square_sum(Matrix &A);//---------------------------------------计算对角元以下的平方和

void stimulate(int m, int i, int j, Matrix P, Matrix &U);//-------------在新线程中积累U,V
void stimulate_householderU(int k, int m, Matrix v, double beta, Matrix &U);//---在新线程中积累U
void stimulate_householderV(int k, int n, Matrix v, double beta, Matrix &V);//---在新线程中积累V

double bisection_for_symmetrtic_triple_diagonal(vector<double> x, vector<double> y, int k);//二分法求(三对角阵)第k个特征值