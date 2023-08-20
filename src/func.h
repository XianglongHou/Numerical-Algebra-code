//实现功能的主要算法
#pragma once
#include<iostream>
#include<vector>
#include"helper.h"
#include<tuple>
#include<complex>
using namespace std;

//第一章线性方程组的直接解法所要求的函数
void forward_subs(vector<vector<double>>& L, vector<double>& b);//前代法

void forward_subs1(vector<vector<double>>& L, vector<double>& b);//对角元为1的前代法

void back_subs(vector<vector<double>>& U, vector<double>& b);//回代法

void back_subs1(vector<vector<double>>& U, vector<double>& b);//对角元为1的回代法

void gauss_elim(vector<vector<double>>& A);//Gauss消去法

void gauss_elim_solvequa(vector<vector<double>> A, vector<double>& b);//gauss消去法后解方程

void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v);//全主元Gauss消去法

void gauss_elim_full_pivoting_solvequa(vector<vector<double>> A, vector<double>& b); //全主元Gauss消去法解方程

void gauss_elim_col_pivoting(vector<vector<double>>& A, vector<int>& u);//列主元Gauss消去法

void gauss_col_piv_solvequa(vector<vector<double>> A, vector<int> u, vector<double>& b);//在已知A的列主元LU分解的情况下，解方程

void gauss_elim_col_pivoting_solvequa(vector<vector<double>> A, vector<double>& b);//列主元消去法解方程


void cholesky_decomp(vector<vector<double>>& A);//对称正定阵标准Cholesky分解

void cholesky_decomp_solvequa(vector<vector<double>> A, vector<double>& b);//利用Cholesky分解解方程

void modified_cholesky_decomp(vector<vector<double>>& A);//改进的平方根法

void modified_cholesky_decomp_solvequa(vector<vector<double>> A, vector<double>& b);//用改进的平方根法解方程

//第二章计算解的精度估计和迭代改进所要求的函数

double matrix_1norm(vector<vector<double>>); //估计矩阵的1范数（优化法）

double matrix_inverse_infnorm(vector<vector<double>>);//估计矩阵逆的1范数（优化法）

double inf_condition_num(vector<vector<double>>); //估计矩阵的条件数

double equasolv_error(vector<vector<double>>, vector<double>, vector<double>);//计算解的精度估计

//第三章最小二乘问题的解法所要求的函数

tuple<vector<double>, double> house(vector<double>); //求使向量除第一个元素全变为0的household变换

void QR(vector<vector<double>>&, vector<double>&); //QR分解

void QR_equsolve(vector<vector<double>>, vector<double>&);//利用QR分解解方程

vector<double> QR_least_square(vector<vector<double>>, vector<double>);//求最小二乘解

//第四章古典迭代解法所要求的函数

void jacobi_itr(vector<vector<double>>, vector<double>, vector<double>&, double tol = 1e-6);//jacobi 迭代

void GS_itr(vector<vector<double>> A, vector<double> b, vector<double>& x, double tol = 1e-6);//GS迭代

void SOR_itr(vector<vector<double>> A, vector<double> b, vector<double>& x, double w, double tol = 1e-6);//SOR迭代



//第五章古典迭代法要求的函数

void CG_method(vector<vector<double>> A, vector<double> b, vector<double>& x, double tol = 1e-7); //共轭梯度法

//第六章非对称特征值问题的计算方法所要求的函数

double find_largest_root(const vector<double> a, int times = 1000);  //幂法找最大根

tuple<vector<vector<double>>, vector<double>> hessenberg_decomp(vector<vector<double>>& A);
//上Hessenberg分解，Hessenberg阵存任A中,Hessenberg阵存在A中，Householder变换的v和beta存在V和b中

void two_step_displacement_qr(vector<vector<double>>& H);//双重步位移的QR代,存在H中，并直接对H12和H23作变换

vector<complex<double>> implicit_qr(vector<vector<double>>& A);//隐式QR算法

//第七章对称特征值问题的计算方法所要求的函数


void Jacobi_Eig(vector<vector<double>>& A, vector<vector<double>>& Q); //jocabi方法计算特征值

double bisection_eig(vector<vector<double>>A, int m);//二分法求指定特征值

vector<double> inverse_power_method(vector<vector<double>> A, double mu);//反幂法计算特征向量

//SVD算法所要求的函数

tuple<vector<vector<double>>, vector<double>, vector<vector<double>>, vector<double>> bidiagnize(vector<vector<double>>& A);//二对角化

tuple<vector<vector<double>>, vector<vector<double>>> SVD_itr(vector<double>& D, vector<double>& S);//SVD迭代


vector<vector<double>> modify(vector<double>& D, vector<double>& S, int ind); //将对角元为0的情况调整成整行为0


tuple<vector<double>, vector<vector<double>>, vector<vector<double>>> SVD(vector<vector<double>> A); //SVD算法最终版