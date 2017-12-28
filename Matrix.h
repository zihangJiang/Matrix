#pragma once

#include <iostream>
#include<thread>
#include<time.h>
#include <vector>

using namespace std;

class Matrix//-------------------------------------------------------------------���������
{
public:
	Matrix();
	~Matrix();
public:
	void construct(int n, double a1, double a2, double a3);//---------------------����Ĺ��캯��
	void construct2_3(double x);//------------------------------------------------����Ĺ��캯��

	void construct_companion_matrix_of_polynomial(double *a, int n);//------------�������ʽ���ѷ���
	void zeros(int m, int n);//---------------------------------------------------����m*n��ȫ0����
	void ones(int m, int n);//----------------------------------------------------����m*n��ȫ0����
	void identity(int n);//-------------------------------------------------------����n�׵�λ����

	Matrix cut(int j, int m, int k, int n);//------------------------------------ѡ��A(j:m,k:n)
	void paste(Matrix C, int j, int m, int k, int n);//--------------------------��Cճ����A(j:m,k:n)
	Matrix T();//----------------------------------------------------------------��ת��
	void show();//----------------------------------------------------------------��ӡ����
	void show_diag_dual();//------------------------------------------------------��ӡ������Խ�
	void show_diag();//-----------------------------------------------------------��ӡ����Խ�

	double max();//---------------------------------------------------------------�������
	double tr();//----------------------------------------------------------------���㼣

	void Givens_down(int i, int j, double c, double s);//----------------------------Givens�任,���»���
	void Givens_up(int i, int j, double c, double s);//---------------------------Givens�任�����ϻ���
	void Givens_left(int i, int j, double c, double s);//-------------------------Givens�任,������
	void Givens_right(int i, int j, double c, double s);//------------------------Givens�任�����һ���

	Matrix Matrix::jacobi_iteration_with_threshold();//-------------------------����jacobi

	Matrix householder();//------------------------------------------------------householder�任

	Matrix direct_QR_iteration();//----------------------------------------------ֱ��QR�㷨
	Matrix QR_decomposition();//-------------------------------------------------QR�ֽ�

	Matrix dual_diag_with_householder(Matrix &U, Matrix &V);//-----------------���Խǻ�(householder�任)
	Matrix small_to_zeros(double thresh);//--------------------------------------��С��Ԫ�ش��0
	Matrix all_row_zeros(int i);//-----------------------------------------------�ѶԽ�ԪΪ0�����б��0
	int check_zeros_of_diag(int i, int j);//--------------------------------------���i��j���Խ�Ԫ����û����
	Matrix SVD(Matrix &U, Matrix &V);//----------------------------------------����ֵ�ֽ�
	Matrix SVD_direct_QR(Matrix &U, Matrix &V);//------------------------------����ֵ�ֽ���ֱ��QR����
	Matrix SVD_iteration(Matrix &P, Matrix &Q);//------------------------------�㷨7.6.2

	double power_iteration();//---------------------------------------------------�ݷ�
	Matrix householder();//------------------------------------------------------householder�任
	Matrix up_Hessenburg_decomposition();//--------------------------------------��Hessenburg�ֽ�
	Matrix implicilt_QR_iteration();//-------------------------------------------��ʽQR�㷨
	Matrix double_shift_QR_iteration();//----------------------------------------˫�ز�λ��QR����
	void check_real_or_complex(int i);//------------------------------------------�ж�(i,i+1)��2*2�´ζԽ�Ԫ��0��������ֵ�Ƿ���ʵ��
	Matrix small_to_zeros(double thresh);//--------------------------------------��С��Ԫ�ش��0

	vector<vector<double>> A;
};

Matrix Minus(Matrix A, Matrix B);//--------------------------------------------���������
Matrix Plus(Matrix A, Matrix B);//---------------------------------------------���������
Matrix Multiply(Matrix A, Matrix B);//-----------------------------------------���������
Matrix Scalar(double beta, Matrix A);//-----------------------------------------���������
double sign(double delta);//---------------------------------------------------------ȡ����
double Square_sum(Matrix &A);//--------------------------------------------------��������Ԫ�ص�ƽ����
double Diag_square_sum(Matrix &A);//---------------------------------------------����Խ�Ԫ��ƽ����
double Not_diag_square_sum(Matrix &A);//-----------------------------------------����Խ�Ԫ��ƽ����
double Lower_diag_square_sum(Matrix &A);//---------------------------------------����Խ�Ԫ���µ�ƽ����

void stimulate(int m, int i, int j, Matrix P, Matrix &U);//-------------�����߳��л���U,V
void stimulate_householderU(int k, int m, Matrix v, double beta, Matrix &U);//---�����߳��л���U
void stimulate_householderV(int k, int n, Matrix v, double beta, Matrix &V);//---�����߳��л���V

double bisection_for_symmetrtic_triple_diagonal(vector<double> x, vector<double> y, int k);//���ַ���(���Խ���)��k������ֵ