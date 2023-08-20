//ʵ�ֹ��ܵ���Ҫ�㷨
#pragma once
#include<iostream>
#include<vector>
#include"helper.h"
#include<tuple>
#include<complex>
using namespace std;

//��һ�����Է������ֱ�ӽⷨ��Ҫ��ĺ���
void forward_subs(vector<vector<double>>& L, vector<double>& b);//ǰ����

void forward_subs1(vector<vector<double>>& L, vector<double>& b);//�Խ�ԪΪ1��ǰ����

void back_subs(vector<vector<double>>& U, vector<double>& b);//�ش���

void back_subs1(vector<vector<double>>& U, vector<double>& b);//�Խ�ԪΪ1�Ļش���

void gauss_elim(vector<vector<double>>& A);//Gauss��ȥ��

void gauss_elim_solvequa(vector<vector<double>> A, vector<double>& b);//gauss��ȥ����ⷽ��

void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v);//ȫ��ԪGauss��ȥ��

void gauss_elim_full_pivoting_solvequa(vector<vector<double>> A, vector<double>& b); //ȫ��ԪGauss��ȥ���ⷽ��

void gauss_elim_col_pivoting(vector<vector<double>>& A, vector<int>& u);//����ԪGauss��ȥ��

void gauss_col_piv_solvequa(vector<vector<double>> A, vector<int> u, vector<double>& b);//����֪A������ԪLU�ֽ������£��ⷽ��

void gauss_elim_col_pivoting_solvequa(vector<vector<double>> A, vector<double>& b);//����Ԫ��ȥ���ⷽ��


void cholesky_decomp(vector<vector<double>>& A);//�Գ��������׼Cholesky�ֽ�

void cholesky_decomp_solvequa(vector<vector<double>> A, vector<double>& b);//����Cholesky�ֽ�ⷽ��

void modified_cholesky_decomp(vector<vector<double>>& A);//�Ľ���ƽ������

void modified_cholesky_decomp_solvequa(vector<vector<double>> A, vector<double>& b);//�øĽ���ƽ�������ⷽ��

//�ڶ��¼����ľ��ȹ��ƺ͵����Ľ���Ҫ��ĺ���

double matrix_1norm(vector<vector<double>>); //���ƾ����1�������Ż�����

double matrix_inverse_infnorm(vector<vector<double>>);//���ƾ������1�������Ż�����

double inf_condition_num(vector<vector<double>>); //���ƾ����������

double equasolv_error(vector<vector<double>>, vector<double>, vector<double>);//�����ľ��ȹ���

//��������С��������Ľⷨ��Ҫ��ĺ���

tuple<vector<double>, double> house(vector<double>); //��ʹ��������һ��Ԫ��ȫ��Ϊ0��household�任

void QR(vector<vector<double>>&, vector<double>&); //QR�ֽ�

void QR_equsolve(vector<vector<double>>, vector<double>&);//����QR�ֽ�ⷽ��

vector<double> QR_least_square(vector<vector<double>>, vector<double>);//����С���˽�

//�����¹ŵ�����ⷨ��Ҫ��ĺ���

void jacobi_itr(vector<vector<double>>, vector<double>, vector<double>&, double tol = 1e-6);//jacobi ����

void GS_itr(vector<vector<double>> A, vector<double> b, vector<double>& x, double tol = 1e-6);//GS����

void SOR_itr(vector<vector<double>> A, vector<double> b, vector<double>& x, double w, double tol = 1e-6);//SOR����



//�����¹ŵ������Ҫ��ĺ���

void CG_method(vector<vector<double>> A, vector<double> b, vector<double>& x, double tol = 1e-7); //�����ݶȷ�

//�����·ǶԳ�����ֵ����ļ��㷽����Ҫ��ĺ���

double find_largest_root(const vector<double> a, int times = 1000);  //�ݷ�������

tuple<vector<vector<double>>, vector<double>> hessenberg_decomp(vector<vector<double>>& A);
//��Hessenberg�ֽ⣬Hessenberg�����A��,Hessenberg�����A�У�Householder�任��v��beta����V��b��

void two_step_displacement_qr(vector<vector<double>>& H);//˫�ز�λ�Ƶ�QR��,����H�У���ֱ�Ӷ�H12��H23���任

vector<complex<double>> implicit_qr(vector<vector<double>>& A);//��ʽQR�㷨

//�����¶Գ�����ֵ����ļ��㷽����Ҫ��ĺ���


void Jacobi_Eig(vector<vector<double>>& A, vector<vector<double>>& Q); //jocabi������������ֵ

double bisection_eig(vector<vector<double>>A, int m);//���ַ���ָ������ֵ

vector<double> inverse_power_method(vector<vector<double>> A, double mu);//���ݷ�������������

//SVD�㷨��Ҫ��ĺ���

tuple<vector<vector<double>>, vector<double>, vector<vector<double>>, vector<double>> bidiagnize(vector<vector<double>>& A);//���Խǻ�

tuple<vector<vector<double>>, vector<vector<double>>> SVD_itr(vector<double>& D, vector<double>& S);//SVD����


vector<vector<double>> modify(vector<double>& D, vector<double>& S, int ind); //���Խ�ԪΪ0���������������Ϊ0


tuple<vector<double>, vector<vector<double>>, vector<vector<double>>> SVD(vector<vector<double>> A); //SVD�㷨���հ�