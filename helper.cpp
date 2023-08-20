#include<vector>
#include<iostream>
#include<stdlib.h>
#include<math.h>
#include <algorithm>
using namespace std;
vector<double> get_col(vector<vector<double>> A, int n)
{
	//求矩阵指定列的列向量
	//输入：矩阵A 所求的n
	//输出：指定列的列向量 vector形式
	if (A[0].size() < n) exit(100);
	vector<double> result;
	for (int i = 0; i < A.size(); i++)
		result.push_back(A[i][n - 1]);
	return result;
}

vector<double> get_col(vector<vector<double>>A, int n, int p, int q)
{
	//输入：列标n，行p-q
	if (A[0].size() < n) exit(100);
	vector<double> result;
	for (int i = p; i <= q; i++)
		result.push_back(A[i][n]);
	return result;
}//截取矩阵指定列的部分列向量

vector<double> get_row(vector<vector<double>>A, int n, int p, int q)
{
	//输入：行标n，列p-q
	if (A.size() < n) exit(100);
	vector<double> result;
	for (int i = p; i <= q; i++)
		result.push_back(A[n][i]);
	return result;
}

double get_innerproduct(vector<double> x, vector<double> y)
{
	//求两个向量的内积
	//输入：向量x，y
	//输出：x和y的内积
	if (x.size() != y.size()) exit(100);
	double result = 0;
	for (int i = 0; i < x.size(); i++)
		result = result + x[i] * y[i];
	return result;
}

vector<double> vector_minus(vector<double>x, vector<double>y)
{
	//求向量差（前-后）
	int n = x.size();
	vector <double> result(n);
	for (int i = 0; i < x.size(); i++)
		result[i] = x[i] - y[i];
	return result;
}

vector<double> vector_plus(vector<double>x, vector<double>y)
{
	//求向量差（前-后）
	int n = x.size();
	vector <double> result(n);
	for (int i = 0; i < x.size(); i++)
		result[i] = x[i] + y[i];
	return result;
}

void scalar_mul(vector<double>& vec, double x)
{
	//向量数乘
	for (int i = 0; i < vec.size(); i++)
		vec[i] = x * vec[i];
}

vector<double> scalar_mul(double x, vector<double> vec)
{
	for (int i = 0; i < vec.size(); i++)
		vec[i] = x * vec[i];
	return vec;
}
vector<double> mat_vec_mult(vector<vector<double>> A, vector<double>& b)
{
	//求矩阵向量乘积
	vector<double> x(A.size(), 0);
	for (int i = 0; i < A.size(); i++)
		x[i] = get_innerproduct(A[i], b);
	return x;
}
vector<vector<double>> martrix_mult(vector<vector<double>> A, vector<vector<double>> B)
{
	//求两个矩阵乘积
	//输入：矩阵A，B
	//输出：A*B
	if (A[0].size() != B.size()) exit(100);
	int m = A.size(), n = B.size(), q = B[0].size();
	vector<vector<double>> result(m, vector<double>(q, 0));
	for (int i = 0; i < m; i++)
		for (int j = 0; j < q; j++)
			for (int k = 0; k < n; k++)
				result[i][j] = result[i][j] + A[i][k] * B[k][j];
	return result;
}

vector<vector<double>> martrix_minus(vector<vector<double>>A, vector<vector<double>>B)
{
	//求矩阵减法
	int m = A.size(), n =A[0].size();
	vector<vector<double>> result(m, vector<double>(n, 0));
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			result[i][j] = A[i][j] - B[i][j];
	return result;

}


vector<double> sub_vector(vector <double> vec, int m, int n)
{
	//求子向量
	//输入：原向量vec,截取片段开始下标m结束下标n
	//输出：子向量 vec[m,m+1,...n]
	vector<double> sub(vec.begin() + m, vec.begin() + n + 1);
	return sub;
}


vector<vector<double>> sub_matrix(vector<vector<double>>A, int m, int n, int p, int q)
{
	//求子矩阵
	//输入：原矩阵A，截取m到n行，p到q列的矩阵
	//输出：截取的子矩阵
	vector<vector<double>> sub(n - m + 1, vector<double>(q - p + 1));
	for (int i = m; i <= n; i++)
		for (int j = p; j <= q; j++)
			sub[i - m][j - p] = A[i][j];
	return sub;
}
void swap_col(vector<vector<double>>& A, int m, int n)
{
	//交换矩阵的两列
	double temp;
	for (int i = 0; i < A.size(); i++)
	{
		temp = A[i][m];
		A[i][m] = A[i][n];
		A[i][n] = temp;
	}

}

double vector_1norm(vector<double>x)
{
	//求向量的1范数
	double sum = 0;
	for (int i = 0; i < x.size(); i++)
		sum = sum + fabs(x[i]);
	return sum;
}

double vector_infnorm(vector<double>y)
{
	//求向量的无穷范数
	double max = 0;
	for (int i = 0; i < y.size(); i++)
		if (max < abs(y[i]))
		{
			max = abs(y[i]);
		}
	return max;
}

double vector_infnorm(vector<double>y, int& index)
{
	//求向量的无穷范数
	//返回下标及无穷范数
	double max = 0;
	for (int i = 0; i < y.size(); i++)
		if (max < abs(y[i]))
		{
			max = abs(y[i]); index = i;
		}
	return max;
}

double matrix_infnorm(vector<vector<double>> A)
{
	//求向量的无穷范数
	double max = 0;
	for (int i = 0; i < A.size(); i++)
		if (max < vector_1norm(A[i]))
			max= vector_1norm(A[i]);
	return max;

}

void matrix_dispaly(vector<vector<double>>A)
{
	//打印矩阵
	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A[0].size(); j++)
			cout << A[i][j] << ' ';
		cout << endl;
	}
	cout << endl;
}

void vector_display(vector<double>b)
{
	for (int i = 0; i < b.size(); i++)
		cout << b[i] << ' ';
	cout << '\n' << endl;
}

void square_transpose(vector<vector<double>>& A)

{
	//方阵转置
	double temp;
	for (int i = 0; i < A.size(); i++)
		for (int j = i + 1; j < A.size(); j++)
		{
			temp = A[i][j];
			A[i][j] = A[j][i];
			A[j][i] = temp;
		}

}

double compute_error(vector<vector<double>>A, vector<double> b, vector<double> x)
{
	//计算方程求解误差
	//输入：矩阵A,方程右端项b，方程的解x
	//输出： Ax-b的无穷范数
	vector<double> y = mat_vec_mult(A, x);
	for (int i = 0; i < b.size(); i++)
		y[i] = abs(y[i] - b[i]);
	vector<double>::iterator biggest = max_element(begin(y), end(y));
	return *biggest;
}

double sgn(double x)
{
	if (x > 0) return 1;
	if (x < 0) return -1;
	if (x == 0) return 0;
}

double bidiag_inf_norm(vector<double> D, vector<double> S)
{
	//二对角阵的无穷范数（m>=n)
	int n = D.size();
	double mymax = 0;
	for (int i = 0; i < n; i++)
	{
		if (i < n - 1)
		{
			double mysum = abs(D[i]) + abs(S[i]);
			if (mysum > mymax)
				mymax = mysum;
		}
		else
			if (abs(D[n - 1]) > mymax)
				mymax = abs(D[n - 1]);
	}
	return mymax;
}

double biggest_elim(vector<vector<double>> A)
{
	//找矩阵的极大元
	int m = A.size(), n = A[0].size();
	double mymax = 0;
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			if (abs(A[i][j]) > mymax)
				mymax = A[i][j];
	return mymax;
}

void diagnize(vector<double> d,vector<vector<double>> &A)
{
	//求diag（d）
	int n = d.size();
	for (int i = 0; i < n; i++)
		A[i][i] = d[i];
}

vector<vector<double>> eye(int n)
{	
	//n阶单位阵

	vector<vector<double>> A(n, vector<double>(n, 0));
	for (int i = 0; i < n; i++)
		A[i][i] = 1;
	return A;

}


