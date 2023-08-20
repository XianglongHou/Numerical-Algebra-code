//辅助函数以实现主要算法
#pragma once
#include<vector>
using namespace std;
vector<double> get_col(vector<vector<double>>, int);
//求矩阵指定列的列向量
//输入：矩阵 所求的列
//输出：指定列的列向量 vector形式

vector<double> get_col(vector<vector<double>>, int, int, int);
//截取矩阵指定列的部分列向量

vector<double> get_row(vector<vector<double>>, int, int, int);
//截取矩阵指定列的部分列向量

void swap_col(vector<vector<double>>&, int, int);
//交换矩阵的两列



double get_innerproduct(vector<double>, vector<double>);
//求向量内积

vector <double> vector_minus(vector<double>, vector<double>);
//求向量差（前-后）



vector <double> vector_plus(vector<double>x, vector<double>y);


void scalar_mul(vector<double>&, double);
//向量数乘


vector <double> scalar_mul(double, vector<double>);
//向量数乘

vector <double> mat_vec_mult(vector<vector<double>> A, vector<double>& b);
//求矩阵向量乘积

vector<vector<double>> martrix_mult(vector<vector<double>>, vector<vector<double>>);
//求矩阵乘积

vector<vector<double>> martrix_minus(vector<vector<double>>, vector<vector<double>>);
//求矩阵减法

vector <double> sub_vector(vector <double>, int, int);
//求子向量

vector<vector<double>> sub_matrix(vector<vector<double>>, int, int, int, int);
//求子矩阵

double vector_1norm(vector<double>);
//求向量的1范数

double vector_infnorm(vector<double>);
//求向量的无穷范数

double matrix_infnorm(vector<vector<double>>);
//求向量的无穷范数

double vector_infnorm(vector<double>, int&);
//求向量的无穷范数并返回index（重载）

void matrix_dispaly(vector<vector<double>>);
//打印矩阵

void vector_display(vector<double>);
//打印向量

void square_transpose(vector<vector<double>>&);
//方阵转置

double compute_error(vector<vector<double>>, vector<double>, vector<double>);
//计算误差

double sgn(double);
//符号函数

double bidiag_inf_norm(vector<double> D, vector<double> S); 
//二对角阵的无穷范数（m>=n)

double biggest_elim(vector<vector<double>>);
//找矩阵的极大元

void diagnize(vector<double> d, vector<vector<double>>&);
//求diag（d）

vector<vector<double>> eye(int n);
//n阶单位阵
