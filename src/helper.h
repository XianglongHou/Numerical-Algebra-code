//����������ʵ����Ҫ�㷨
#pragma once
#include<vector>
using namespace std;
vector<double> get_col(vector<vector<double>>, int);
//�����ָ���е�������
//���룺���� �������
//�����ָ���е������� vector��ʽ

vector<double> get_col(vector<vector<double>>, int, int, int);
//��ȡ����ָ���еĲ���������

vector<double> get_row(vector<vector<double>>, int, int, int);
//��ȡ����ָ���еĲ���������

void swap_col(vector<vector<double>>&, int, int);
//�������������



double get_innerproduct(vector<double>, vector<double>);
//�������ڻ�

vector <double> vector_minus(vector<double>, vector<double>);
//�������ǰ-��



vector <double> vector_plus(vector<double>x, vector<double>y);


void scalar_mul(vector<double>&, double);
//��������


vector <double> scalar_mul(double, vector<double>);
//��������

vector <double> mat_vec_mult(vector<vector<double>> A, vector<double>& b);
//����������˻�

vector<vector<double>> martrix_mult(vector<vector<double>>, vector<vector<double>>);
//�����˻�

vector<vector<double>> martrix_minus(vector<vector<double>>, vector<vector<double>>);
//��������

vector <double> sub_vector(vector <double>, int, int);
//��������

vector<vector<double>> sub_matrix(vector<vector<double>>, int, int, int, int);
//���Ӿ���

double vector_1norm(vector<double>);
//��������1����

double vector_infnorm(vector<double>);
//�������������

double matrix_infnorm(vector<vector<double>>);
//�������������

double vector_infnorm(vector<double>, int&);
//�������������������index�����أ�

void matrix_dispaly(vector<vector<double>>);
//��ӡ����

void vector_display(vector<double>);
//��ӡ����

void square_transpose(vector<vector<double>>&);
//����ת��

double compute_error(vector<vector<double>>, vector<double>, vector<double>);
//�������

double sgn(double);
//���ź���

double bidiag_inf_norm(vector<double> D, vector<double> S); 
//���Խ�����������m>=n)

double biggest_elim(vector<vector<double>>);
//�Ҿ���ļ���Ԫ

void diagnize(vector<double> d, vector<vector<double>>&);
//��diag��d��

vector<vector<double>> eye(int n);
//n�׵�λ��
