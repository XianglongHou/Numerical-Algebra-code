#include<vector>
#include<iostream>
#include<stdlib.h>
#include<math.h>
#include <algorithm>
using namespace std;
vector<double> get_col(vector<vector<double>> A, int n)
{
	//�����ָ���е�������
	//���룺����A �����n
	//�����ָ���е������� vector��ʽ
	if (A[0].size() < n) exit(100);
	vector<double> result;
	for (int i = 0; i < A.size(); i++)
		result.push_back(A[i][n - 1]);
	return result;
}

vector<double> get_col(vector<vector<double>>A, int n, int p, int q)
{
	//���룺�б�n����p-q
	if (A[0].size() < n) exit(100);
	vector<double> result;
	for (int i = p; i <= q; i++)
		result.push_back(A[i][n]);
	return result;
}//��ȡ����ָ���еĲ���������

vector<double> get_row(vector<vector<double>>A, int n, int p, int q)
{
	//���룺�б�n����p-q
	if (A.size() < n) exit(100);
	vector<double> result;
	for (int i = p; i <= q; i++)
		result.push_back(A[n][i]);
	return result;
}

double get_innerproduct(vector<double> x, vector<double> y)
{
	//�������������ڻ�
	//���룺����x��y
	//�����x��y���ڻ�
	if (x.size() != y.size()) exit(100);
	double result = 0;
	for (int i = 0; i < x.size(); i++)
		result = result + x[i] * y[i];
	return result;
}

vector<double> vector_minus(vector<double>x, vector<double>y)
{
	//�������ǰ-��
	int n = x.size();
	vector <double> result(n);
	for (int i = 0; i < x.size(); i++)
		result[i] = x[i] - y[i];
	return result;
}

vector<double> vector_plus(vector<double>x, vector<double>y)
{
	//�������ǰ-��
	int n = x.size();
	vector <double> result(n);
	for (int i = 0; i < x.size(); i++)
		result[i] = x[i] + y[i];
	return result;
}

void scalar_mul(vector<double>& vec, double x)
{
	//��������
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
	//����������˻�
	vector<double> x(A.size(), 0);
	for (int i = 0; i < A.size(); i++)
		x[i] = get_innerproduct(A[i], b);
	return x;
}
vector<vector<double>> martrix_mult(vector<vector<double>> A, vector<vector<double>> B)
{
	//����������˻�
	//���룺����A��B
	//�����A*B
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
	//��������
	int m = A.size(), n =A[0].size();
	vector<vector<double>> result(m, vector<double>(n, 0));
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			result[i][j] = A[i][j] - B[i][j];
	return result;

}


vector<double> sub_vector(vector <double> vec, int m, int n)
{
	//��������
	//���룺ԭ����vec,��ȡƬ�ο�ʼ�±�m�����±�n
	//����������� vec[m,m+1,...n]
	vector<double> sub(vec.begin() + m, vec.begin() + n + 1);
	return sub;
}


vector<vector<double>> sub_matrix(vector<vector<double>>A, int m, int n, int p, int q)
{
	//���Ӿ���
	//���룺ԭ����A����ȡm��n�У�p��q�еľ���
	//�������ȡ���Ӿ���
	vector<vector<double>> sub(n - m + 1, vector<double>(q - p + 1));
	for (int i = m; i <= n; i++)
		for (int j = p; j <= q; j++)
			sub[i - m][j - p] = A[i][j];
	return sub;
}
void swap_col(vector<vector<double>>& A, int m, int n)
{
	//�������������
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
	//��������1����
	double sum = 0;
	for (int i = 0; i < x.size(); i++)
		sum = sum + fabs(x[i]);
	return sum;
}

double vector_infnorm(vector<double>y)
{
	//�������������
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
	//�������������
	//�����±꼰�����
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
	//�������������
	double max = 0;
	for (int i = 0; i < A.size(); i++)
		if (max < vector_1norm(A[i]))
			max= vector_1norm(A[i]);
	return max;

}

void matrix_dispaly(vector<vector<double>>A)
{
	//��ӡ����
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
	//����ת��
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
	//���㷽��������
	//���룺����A,�����Ҷ���b�����̵Ľ�x
	//����� Ax-b�������
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
	//���Խ�����������m>=n)
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
	//�Ҿ���ļ���Ԫ
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
	//��diag��d��
	int n = d.size();
	for (int i = 0; i < n; i++)
		A[i][i] = d[i];
}

vector<vector<double>> eye(int n)
{	
	//n�׵�λ��

	vector<vector<double>> A(n, vector<double>(n, 0));
	for (int i = 0; i < n; i++)
		A[i][i] = 1;
	return A;

}


