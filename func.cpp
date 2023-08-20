#include<iostream>
#include<vector>
#include"helper.h"
#include<math.h>
#include<tuple>
#include<complex>

using namespace std;
void forward_subs(vector<vector<double>>& L, vector<double>& b)
{
	//前代法解下三角方程组
	//输入：下三角矩阵L，方程右端项b
	//输出：方程的解b
	int i, j, n = L.size();
	for (i = 0; i < n - 1; i++)
	{
		b[i] = b[i] / L[i][i];
		for (j = i + 1; j < n; j++)
			b[j] = b[j] - b[i] * L[j][i];
	}
	b[n - 1] = b[n - 1] / L[n - 1][n - 1];
}

void forward_subs1(vector<vector<double>>& L, vector<double>& b)
{
	//对角元为1的前代法
	//输入：单位下三角矩阵L，方程右端项b
	//输出：方程的解b
	int i, j, n = L.size();
	for (i = 0; i < n; i++)
		for (j = i + 1; j < n; j++)
			b[j] = b[j] - b[i] * L[j][i];
	//cnm!!!
}

void back_subs(vector<vector<double>>& U, vector<double>& b)
{
	//回代法
	//输入：上三角矩阵U，方程右端项b
	//输出：方程的解b
	int n = U[0].size();
	for (int i = n - 1; i >= 1; i--) //搞到0就错了
	{
		b[i] = b[i] / U[i][i];
		for (int j = i - 1; j >= 0; j--)
			b[j] = b[j] - b[i] * U[j][i];
	}
	b[0] = b[0] / U[0][0];
}

void back_subs1(vector<vector<double>>& U, vector<double>& b)
{
	//对角元为1的回代法
	//输入：单位下三角矩阵U，方程右端项b
	//输出：方程的解b
	int n = U.size();
	for (int i = n - 1; i >= 0; i--)
		for (int j = i - 1; j >= 0; j--)
			b[j] = b[j] - b[i] * U[j][i];
}

void gauss_elim(vector<vector<double>>& A)
{
	//Gauss消去法（LU分解）
	//输入：1到n-1阶顺序主子式非0的方阵A
	//输出：LU分解中L矩阵的下三角部分（不包括对角线）和U矩阵的上三角部分，均存储在A中
	int n = A.size();
	for (int k = 0; k < n - 1; k++)
		for (int i = k + 1; i < n; i++)
		{
			A[i][k] = A[i][k] / A[k][k];
			for (int j = k + 1; j < n; j++)
				A[i][j] = A[i][j] - A[i][k] * A[k][j];
		}
}

void gauss_elim_solvequa(vector<vector<double>> A, vector<double>& b)
{
	//利用gauss消去法解方程
	//输入：可进行LU分解的矩阵A，方程右端项b
	//输出：方程的解b
	gauss_elim(A);
	forward_subs1(A, b);
	back_subs(A, b);
}

void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v)
{
	//全主元Gauss消去法
	//输入：矩阵A，已初始化的向量u，v
	//输出：：LU分解中L矩阵的下三角部分（不包括对角线）和U矩阵的上三角部分，均存储在A中
	int n = A.size(), max_row, max_col;
	for (int k = 0; k < n - 1; k++)
	{
		max_row = k;
		max_col = k;
		double mymax = abs(A[k][k]);
		for (int i = k; i < n; i++)
			for (int j = k; j < n; j++)
				if (mymax < abs(A[i][j])) { max_row = i; max_col = j; mymax = abs(A[i][j]); }
		A[k].swap(A[max_row]);
		swap_col(A, k, max_col);
		u[k] = max_row;
		v[k] = max_col;
		if (A[k][k] != 0)
		{
			for (int i = k + 1; i < n; i++)
			{
				A[i][k] = A[i][k] / A[k][k];
				for (int j = k + 1; j < n; j++)
					A[i][j] = A[i][j] - A[i][k] * A[k][j];
			}//matrix_dispaly(A);
		}

		else exit(100);
	}
}

void gauss_elim_full_pivoting_solvequa(vector<vector<double>> A, vector<double>& b)
{
	////全主元Gauss消去法解方程
	//输入：矩阵A,方程右端项b
	//输出：方程的解b
	vector<int> u(A.size() - 1, 0), v(A.size() - 1, 0);
	gauss_elim_full_pivoting(A, u, v);
	for (int i = 0; i < u.size(); i++)
		swap(b[i], b[u[i]]);
	forward_subs1(A, b);  //subs1!!!!
	back_subs(A, b);
	for (int i = v.size() - 1; i >= 0; i--)
		swap(b[i], b[v[i]]);
}

void gauss_elim_col_pivoting(vector<vector<double>>& A, vector<int>& u)
{
	//列主元Gauss消去法
	//输入：矩阵A，已初始化的向量u
	//输出：：LU分解中L矩阵的下三角部分（不包括对角线）和U矩阵的上三角部分，均存储在A中
	int n = A.size(), max_row;
	for (int k = 0; k < n - 1; k++)
	{
		max_row = k;
		double mymax = abs(A[k][k]);
		for (int i = k + 1; i < n; i++)
			if (mymax < abs(A[i][k])) { max_row = i; mymax = abs(A[i][k]); }
		A[k].swap(A[max_row]);
		u[k] = max_row;
		if (A[k][k] != 0)
		{
			for (int i = k + 1; i < n; i++)
			{
				A[i][k] = A[i][k] / A[k][k];
				for (int j = k + 1; j < n; j++)
					A[i][j] = A[i][j] - A[i][k] * A[k][j];
			}
		}

		else exit(100);
	}
}

void gauss_elim_col_pivoting_solvequa(vector<vector<double>> A, vector<double>& b)
{
	//列主元法解方程
	//输入：矩阵A,方程右端项b
	//输出：方程的解b
	int n = A.size();
	vector<int> u(A.size() - 1, 0);
	gauss_elim_col_pivoting(A, u);
	for (int i = 0; i < u.size(); i++)
		swap(b[i], b[u[i]]);
	forward_subs1(A, b);  //subs1!!!!
	back_subs(A, b);
}

void gauss_col_piv_solvequa(vector<vector<double>> A, vector<int> u, vector<double>& b)
{
	//在已知A的列主元LU分解的情况下，解方程
	for (int i = 0; i < u.size(); i++)
		swap(b[i], b[u[i]]);
	forward_subs1(A, b);  //subs1!!!!
	back_subs(A, b);
}


void cholesky_decomp(vector<vector<double>>& A)
{
	//对称正定阵标准Cholesky分解
	//输入：正定对称阵A
	//输出：Cholesky因子L，存储在A的下三角中
	int n = A.size();
	for (int k = 0; k < n; k++)
	{
		A[k][k] = sqrt(A[k][k]);
		for (int j = k + 1; j < n; j++)
			A[j][k] = A[j][k] / A[k][k]; //得到第k列元素的最终值
		for (int j = k + 1; j < n; j++)
			for (int i = j; i < n; i++)
				A[i][j] = A[i][j] - A[i][k] * A[j][k];
	}
}

void cholesky_decomp_solvequa(vector<vector<double>> A, vector<double>& b)
{
	//利用Cholesky分解解方程
	//输入：正定阵A，方程右端项b
	//输出：方程的解b
	cholesky_decomp(A);
	forward_subs(A, b);
	square_transpose(A);
	back_subs(A, b);
}

void modified_cholesky_decomp(vector<vector<double>>& A)
{
	//改进的平方根法A=LDL^T
	//输入：正定对称阵A
	//输出：单位下三角阵L,对角阵D (L除对角元外存储在A的下三角，D储存在A的对角线上）
	int n = A.size();
	vector<double> v(n, 0);
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i <= j - 1; i++)   //这里控制条件要取等
			v[i] = A[j][i] * A[i][i];
		for (int i = 0; i <= j - 1; i++)
			A[j][j] = A[j][j] - A[j][i] * v[i];
		for (int k = j + 1; k < n; k++)
		{
			for (int i = 0; i <= j - 1; i++)
				A[k][j] = A[k][j] - A[k][i] * v[i];
			A[k][j] = A[k][j] / A[j][j];
		}
	}
}

void modified_cholesky_decomp_solvequa(vector<vector<double>> A, vector<double>& b)
{
	//用改进的平方根法解方程
	//输入：正定阵A，方程右端项b
	//输出：方程的解b
	modified_cholesky_decomp(A);
	forward_subs1(A, b);
	for (int i = 0; i < b.size(); i++)
		b[i] = b[i] / A[i][i];
	square_transpose(A);
	back_subs1(A, b);
}

double matrix_1norm(vector<vector<double>> B)
{
	//估计矩阵的1范数（优化法）
	//输入：矩阵B
	//输出：矩阵1范数的估计值
	int n = B.size();
	int index = 0;
	vector<vector<double>> Bt(B);
	square_transpose(Bt);//B的转置
	vector<double>x(n, 1 / double(n)); //这个地方float
	vector<double>w(n), v(n), z(n);
	while (1)
	{
		w = mat_vec_mult(B, x);
		for (int i = 0; i < n; i++)
			v[i] = sgn(w[i]);
		z = mat_vec_mult(Bt, v);
		double z_norm = vector_infnorm(z, index);
		if (z_norm <= get_innerproduct(z, x))
		{
			return(vector_1norm(w));
		}
		else
		{
			for (int j = 0; j < n; j++)
				x[j] = 0;
			x[index] = 1;
		}
	}
}

double matrix_inverse_infnorm(vector<vector<double>>A)
{
	//估计矩阵逆的1范数（优化法）
	//输入：矩阵A
	//输出：矩阵A^(-1)的1范数
	int n = A.size();
	int index = 0;
	vector<vector<double>> At(A);
	square_transpose(At);//A的转置
	vector<double>x(n, 1 / double(n)); //这个地方double
	vector<double>w(n), v(n), z(n);
	vector<int> u1(A.size() - 1, 0);
	gauss_elim_col_pivoting(A, u1);
	vector<int> u2(At.size() - 1, 0);
	gauss_elim_col_pivoting(At, u2);//列主元的Gauss分解


	while (1)
	{
		w.assign(x.begin(), x.end());
		for (int i = 0; i < u2.size(); i++)
			swap(w[i], w[u2[i]]);
		forward_subs1(At, w);
		back_subs(At, w);//利用分解解方程At*w=x
		for (int i = 0; i < n; i++)
			v[i] = sgn(w[i]);
		z.assign(v.begin(), v.end());
		for (int i = 0; i < u1.size(); i++)
			swap(z[i], z[u1[i]]);
		forward_subs1(A, z);
		back_subs(A, z);//利用分解解方程A*z=v
		double z_norm = vector_infnorm(z, index);
		if (z_norm <= get_innerproduct(z, x))
			return(vector_1norm(w));
		else
		{
			for (int j = 0; j < n; j++)
				x[j] = 0;
			x[index] = 1;
		}
	}
}

double inf_condition_num(vector<vector<double>>A)
{
	//估计矩阵的条件数
	//输入：矩阵A
	//输出：矩阵A的条件数
	vector<vector<double>>At(A);
	square_transpose(At);
	return(matrix_1norm(At) * matrix_inverse_infnorm(A));
}

double equasolv_error(vector<vector<double>>A, vector<double>b, vector<double>x)
{
	//计算解的精度估计
	//输入：矩阵A，方程右端项b，方程的计算解x
	//输出：计算解的精度估计
	int n = A.size();
	double b_norm = vector_infnorm(b);
	vector<double> r(n);
	r = vector_minus(b, mat_vec_mult(A, x));
	double r_norm = vector_infnorm(r);
	return inf_condition_num(A) * r_norm / b_norm;
}

//第三章主要用到的函数

tuple<vector<double>, double> house(vector<double>x)
{
	//求使向量除第一个元素全变为0的household变换
	int n = x.size();
	double coe = 1, norm_square;
	double inf_norm = vector_infnorm(x);
	scalar_mul(x, 1 / inf_norm);
	vector<double> v(x);
	double s = get_innerproduct(sub_vector(x, 1, n - 1), sub_vector(x, 1, n - 1));
	if (s == 0)
		coe = 0;
	else
	{
		norm_square = pow(x[0] * x[0] + s, 0.5);
		if (x[0] <= 0) v[0] = x[0] - norm_square;
		else v[0] = -s / (x[0] + norm_square);
		coe = (2 * v[0] * v[0]) / (s + v[0] * v[0]);
		scalar_mul(v, 1 / v[0]);
	}
	return make_tuple(v, coe);
}

void QR(vector<vector<double>>& A, vector<double>& d)
{
	//QR分解
	//输入：原矩阵A,已初始化用于储存household矩阵系数的向量d
	//输出：A下三角部分存储Household，上三角部分存放R，d存放household矩阵系数
	int m = A.size(), n = A[0].size();
	vector<double> v;
	for (int j = 0; j < n; j++)
	{
		if (j < m - 1)
		{
			tuple<vector<double>, double> household = house(get_col(A, j, j, m - 1));
			for (int i = j; i < n; i++)
			{
				v = get<0>(household);
				double s = get<1>(household) * get_innerproduct(v, get_col(A, i, j, m - 1));
				for (int k = j; k < m; k++)
					A[k][i] = A[k][i] - s * v[k - j];
			}
			d[j] = get<1>(household);
			for (int k = j + 1; k < m; k++)
				A[k][j] = v[k - j];
		}
	}
}

void QR_equsolve(vector<vector<double>>A, vector<double>& b)
{
	//利用QR分解解方程
	int n = A[0].size();
	vector<double> d(n);
	QR(A, d);
	for (int j = 0; j < n - 1; j++)
	{
		vector<double> v = get_col(A, j, j + 1, n - 1);
		v.insert(v.begin(), 1.0);
		double s = d[j] * get_innerproduct(v, sub_vector(b, j, n - 1));
		for (int k = j; k < n; k++)
			b[k] = b[k] - s * v[k - j];
	}//对b进行household变换
	back_subs(A, b);
}

vector<double> QR_least_square(vector<vector<double>> A, vector<double> b)
{
	//求最小二乘解
	int m = A.size(), n = A[0].size();
	vector<double> d(n);
	QR(A, d);

	for (int j = 0; j < n; j++)
	{
		if (j < m - 1)
		{
			vector<double> v = get_col(A, j, j + 1, m - 1);
			v.insert(v.begin(), 1.0);
			double s = d[j] * get_innerproduct(v, sub_vector(b, j, m - 1));
			for (int k = j; k < m; k++)
				b[k] = b[k] - s * v[k - j];
		}
	}//对b进行household变换
	vector<double> r = sub_vector(b, 0, n - 1);
	back_subs(A, r);
	return r;
}

void jacobi_itr(vector<vector<double>> A, vector<double> b, vector<double>& x, double tol = 1e-6)
{
	//jacobi 迭代
	//输入：系数矩阵A,迭代初始值b,迭代终止条件tol（以x_{k+1}-x_k的无穷范数）
	// 输出：方程的解b，并打印迭代次数
	double d, err = 10;
	int count = 0;
	for (int i = 0; i < A.size(); i++)
	{
		d = A[i][i];
		A[i][i] = 0;
		b[i] = b[i] / d;
		for (int j = 0; j < A.size(); j++)
			A[i][j] = -A[i][j] / d;
	}//计算迭代矩阵及偏置项(搞错了差个负号
	while (err > tol)
	{
		vector<double> temp(x);
		x = vector_plus(mat_vec_mult(A, x), b);
		err = vector_infnorm(vector_minus(x, temp));
		count++;
	}
	cout << "jacobi迭代次数：" << count << endl;
}

void GS_itr(vector<vector<double>> A, vector<double> b, vector<double>& x, double tol = 1e-6)
{
	double d, err;
	int count = 0;
	for (int i = 0; i < A.size(); i++)
	{
		d = A[i][i];
		A[i][i] = 0;
		b[i] = b[i] / d;
		for (int j = 0; j < A.size(); j++)
			A[i][j] = -A[i][j] / d;
	}//计算迭代矩阵及偏置项
	do
	{
		vector<double> temp(x);
		for (int i = 0; i < A.size(); i++)
			x[i] = get_innerproduct(A[i], x) + b[i];
		err = vector_infnorm(vector_minus(x, temp));
		count++;
	} while (err > tol);
	cout << "G-S迭代次数：" << count << endl;
}

void SOR_itr(vector<vector<double>> A, vector<double> b, vector<double>& x, double w, double tol = 1e-6)
{
	double d, err = 10;
	int count = 0;
	for (int i = 0; i < A.size(); i++)
	{
		d = A[i][i];
		A[i][i] = 0;
		b[i] = b[i] / d;
		for (int j = 0; j < A.size(); j++)
			A[i][j] = -A[i][j] / d;
	}//计算迭代矩阵及偏置项
	do
	{
		vector<double> temp(x);
		for (int i = 0; i < A.size(); i++)
			x[i] = (1 - w) * x[i] + w * (get_innerproduct(A[i], x) + b[i]);
		err = vector_infnorm(vector_minus(x, temp));
		count++;
	} while (err > tol);
	cout << "SOR迭代次数：" << count << "  ，松弛因子" << w << endl;
}

void CG_method(vector<vector<double>> A, vector<double> b, vector<double>& x, double tol = 1e-7)
{
	//共轭梯度法
	//输入：系数矩阵（正定）A，右端项b，初始迭代值x，容忍度tol
	//输出：解x
	int count = 1;
	vector<double> r = vector_minus(b, mat_vec_mult(A, x)), p, w, x_temp;
	double rou = get_innerproduct(r, r), rou_temp = 0, alpha = 0, beta = 0;
	while (1)
	{
		x_temp = x;
		if (count != 1)
		{
			beta = rou / rou_temp;
			p = vector_plus(r, scalar_mul(beta, p));
		}
		else
			p = r;
		w = mat_vec_mult(A, p);
		alpha = rou / get_innerproduct(p, w);
		x = vector_plus(x, scalar_mul(alpha, p));
		r = vector_minus(r, scalar_mul(alpha, w));
		rou_temp = rou;
		rou = get_innerproduct(r, r);
		if (vector_infnorm(vector_minus(x_temp, x)) < tol)
			break;
		count++;
	}
	cout << "共轭梯度法迭代次数：" << count << endl;
}

double find_largest_root(const vector<double> a, int times = 1000)
{
	//幂法找多项式最大根（不写成矩阵乘法形式，算出成友方阵乘迭代向量的表达式计算直接迭代）
	//输入：多项式各项系数（从0次到n-1次）
	//输出：多项式的模最大根
	int n = a.size(), ind=0;
	vector<double> x(n, 1);
	double u = 1, temp = 0, temp1;
	for (int i = 0; i < times; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (j != 0)
			{
				temp1 = x[j];
				x[j] = temp - a[j] * x[n - 1];
				temp = temp1;
			}
			else
			{
				temp = x[0];
				x[0] = -a[0] * x[n - 1];
			}
		}
		u = vector_infnorm(x, ind);
		scalar_mul(x, 1 / u);
	}
	if (x[ind] < 0)
		return -u;
	return u;
}

tuple<vector<vector<double>>, vector<double>> hessenberg_decomp(vector<vector<double>>& A)
{
	//上Hessenberg分解，Hessenberg阵存任A中,Hessenberg阵存在A中，Householder变换的v和beta存在V和b中
	//输入：待上Hessenberg化的矩阵A
	//输出：Hessenberg矩阵A，Householder变换的v和beta以元组形式返回
	vector<vector<double>> V;
	vector<double> b, v;
	double beta;
	int n = A.size();
	for (int j = 0; j < n - 2; j++)
	{
		tuple<vector<double>, double> household = house(get_col(A, j, j + 1, n - 1));
		v = get<0>(household);
		beta = get<1>(household);
		for (int i = j; i < n; i++)
		{

			double s = beta * get_innerproduct(v, get_col(A, i, j + 1, n - 1));
			for (int k = j + 1; k < n; k++)
				A[k][i] = A[k][i] - s * v[k - j - 1];
		}
		for (int i = 0; i < n; i++)
		{
			double s = beta * get_innerproduct(v, get_row(A, i, j + 1, n - 1));
			for (int k = j + 1; k < n; k++)
				A[i][k] = A[i][k] - s * v[k - j - 1];
		}
		b.push_back(beta);
		V.push_back(v);
	}
	return{ V, b };
}

void two_step_displacement_qr(vector<vector<double>>& H)
{
	//双重步位移的QR代,存在H中，并直接对H12和H23作变换
	int n = H.size(), m = n - 1;
	double s = H[m - 1][m - 1] + H[n - 1][n - 1], t = H[m - 1][m - 1] * H[n - 1][n - 1] - H[m - 1][n - 1] * H[n - 1][m - 1]
		, x = H[0][0] * H[0][0] + H[0][1] * H[1][0] - s * H[0][0] + t, y = H[1][0] * (H[0][0] + H[1][1] - s),
		z = H[1][0] * H[2][1];
	tuple<vector<double>, double> household;
	vector<double> to_house, v;
	double beta;
	for (int k = -1; k < n - 3; k++)
	{
		to_house = { x,y,z };
		household = house(to_house);
		v = get<0>(household);
		beta = get<1>(household);
		int q = max(0, k);
		for (int i = q; i < n; i++)
		{

			double s = beta * get_innerproduct(v, get_col(H, i, k + 1, k + 3));
			for (int j = k + 1; j <= k + 3; j++)
				H[j][i] = H[j][i] - s * v[j - k - 1];

		}

		int r = min(k + 4, n - 1);
		for (int i = 0; i <= r; i++)
		{

			double s = beta * get_innerproduct(v, get_row(H, i, k + 1, k + 3));
			for (int j = k + 1; j <= k + 3; j++)
				H[i][j] = H[i][j] - s * v[j - k - 1];
		}
		x = H[k + 2][k + 1];
		y = H[k + 3][k + 1];
		if (k < n - 4)
			z = H[k + 4][k + 1];

	}
	to_house = { x,y };
	household = house(to_house);
	v = get<0>(household);
	beta = get<1>(household);
	for (int i = n - 3; i < n; i++)
	{
		double s = beta * get_innerproduct(v, get_col(H, i, n - 2, n - 1));
		for (int j = n - 2; j < n; j++)
			H[j][i] = H[j][i] - s * v[j - n + 2];
	}
	for (int i = 0; i < n; i++)
	{
		double s = beta * get_innerproduct(v, get_row(H, i, n - 2, n - 1));
		for (int j = n - 2; j < n; j++)
			H[i][j] = H[i][j] - s * v[j - n + 2];
	}
}

vector<complex<double>> get_root(vector<vector<double>> A);

vector<complex<double>> implicit_qr(vector<vector<double>>& A)
{
	int n = A.size();
	double acc = 1e-6;
	vector<vector<double>> A_copy;
	hessenberg_decomp(A);
	while (1)
	{
		for (int i = 1; i < n; i++)
		{
			if (abs(A[i][i - 1]) < (abs(A[i][i]) + abs(A[i - 1][i - 1])) * acc)
				A[i][i - 1] = 0;
		}

		int m = n - 1;
		for (m = n - 1; m >= 1; )
		{
			if (abs(A[m][m - 1]) < acc)
				m--;
			else if (m == 1)
				break;
			else if (abs(A[m - 1][m - 2]) < acc)
				m = m - 2;
			else
				break;
		}
		int l;
		for (l = m; l > 0; l--)
			if (abs(A[l][l - 1]) < acc)
				break;
		if (m == 0 || m == 1)
			break;
		else

		{
			A_copy = sub_matrix(A, l, m, l, m);
			two_step_displacement_qr(A_copy);
			for (int i = l; i <= m; i++)
				for (int j = l; j <= m; j++)
					A[i][j] = A_copy[i - l][j - l];
		}

	}
	vector<complex<double>> z;
	z = get_root(A);
	return z;
}//隐式QR算法

vector<complex<double>> get_root(vector<vector<double>> A)
{
	vector<complex<double>> z;
	int n = A.size();
	for (int m = 0; m < n; )
	{
		if (m == n - 1)
		{
			z.push_back(A[m][m]);
			m++;
		}

		else if (A[m + 1][m] == 0)
		{
			z.push_back(A[m][m]);
			m++;
		}
		else
		{
			double delta = (A[m][m] - A[m + 1][m + 1]) * (A[m][m] - A[m + 1][m + 1]) + 4 * A[m][m + 1] * A[m + 1][m];
			double real = (A[m][m] + A[m + 1][m + 1]) / 2;
			double img = pow(-delta, 0.5) / 2;
			z.push_back({ real,img });
			z.push_back({ real,-img });
			m = m + 2;
		}
	}
	return z;
}


void Givens(vector<vector<double>>& A, vector<vector<double>>& Q, int i, int j)
{
	//Givens矩阵左乘A
	//输入：A,Givens矩阵Q，作用列i,j
	//输出：变换后的矩阵A
	int Matrix_Size;
	Matrix_Size = A.size();

	double ta;
	ta = (A[j][j] - A[i][i]) / (2 * A[i][j]);

	double t;
	t = 1.0 / (fabs(ta) + sqrt(1 + ta * ta));
	if (ta < 0) { t = -t; }

	double c = 1.0/ sqrt(1 + t * t);
	double s = t * c;

	for (int i = 0; i < Matrix_Size; i++)
		Q[i][i] = 1;

	Q[i][i] = c;
	Q[i][j] = s;
	Q[j][i] = -s;
	Q[j][j] = c;
}


// givens 变换右乘
void Givens_matrix_times_right(vector<vector<double>>& A, vector<vector<double>>& Q, int i, int j)
{
	// givens 变换右乘
	int Matrix_Size;
	Matrix_Size = A.size();
	double c = Q[i][i];
	double s = -Q[i][j];

	for (int k = 0; k < Matrix_Size; k++)
	{
		double temp = A[k][i];
		A[k][i] = A[k][i] + A[k][j] * (s / c);
		A[k][j] = A[k][j] + temp * (-s / c);
		A[k][i] = A[k][i] * c;
		A[k][j] = A[k][j] * c;
	}
}

void Givens_matrix_times_right(vector<vector<double>>& A, vector<double> Q, int i, int j)
{
	// givens 变换右乘
	int n = A.size();
	double c = Q[0],s = Q[1];

	for (int k = 0; k < n; k++)
	{
		double temp = A[k][i];
		A[k][i] = A[k][i] + A[k][j] * (-s / c);
		A[k][j] = A[k][j] + temp * (s / c);
		A[k][i] = A[k][i] * c;
		A[k][j] = A[k][j] * c;
	}
}

void Givens_matrix_times_left(vector<vector<double>>& A, vector<vector<double>>& Q, int i, int j)
{
	//givens变换左乘
	int n;
	n = A.size();
	double c = Q[i][i],s = -Q[i][j];

	for (int k = 0; k < n; k++)
	{
		double temp = A[i][k];
		A[i][k] = A[i][k] + A[j][k] * (-s / c);
		A[j][k] = A[j][k] + temp * (s / c);
		A[i][k] = A[i][k] * c;
		A[j][k] = A[j][k] * c;
	}
}

void Givens_matrix_times_left(vector<vector<double>>& A, vector<double>& Q, int i, int j)
{
	//givens变换左乘
	int n;
	n = A.size();
	double c = Q[0],s = Q[1];
	for (int k = 0; k < n; k++)
	{
		double temp = A[i][k];
		A[i][k] = A[i][k] + A[j][k] * (s / c);
		A[j][k] = A[j][k] + temp * (-s / c);
		A[i][k] = A[i][k] * c;
		A[j][k] = A[j][k] * c;
	}
}

void Jacobi_Eig(vector<vector<double>>& A, vector<vector<double>>& Q)
{
	int n = A.size();

	vector<vector<double>> A1((n), vector<double>(n));
	
	double delta=0.0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (i != j)
				delta = delta + pow(A[i][j], 2);

	delta = sqrt(delta);
	int times = 0;
	while (delta >= 1e-7)
	{
		int flag = 0;
		for (int i = 0; i < n; i++)
		{
			for (int j = i + 1; j < n; j++)
			{
				if (fabs(A[i][j]) > delta)
				{
					vector<vector<double>> Q1((n), vector<double>(n));
					Givens(A, Q1, i, j);
					Givens_matrix_times_right(Q, Q1, i, j);
					Givens_matrix_times_right(A, Q1, i, j);
					square_transpose(Q1);
					Givens_matrix_times_left(A, Q1, i, j);
					flag = 1;
					times++;
				}
			}
		}
		if (flag == 0)
			delta = delta / n;
	}
	cout << "迭代次数为：" << times << endl;
}





int number_changes(vector<double>& A, vector<double>& B, double mu)
{
	//求mu的变号数
	//输入：对角元A，次对角元B，mu
	//输出：变号数s
	int n = A.size();
	int s = 0;
	double q = A[0] - mu;
	for (int i = 0; i < n; i++)
	{
		if (q < 0)
			s++;
		if (i < n - 1)
		{
			if (q == 0)
				q = fabs(B[i + 1]) * (1e-10);
			q = A[i + 1] - mu - B[i + 1] * B[i + 1] / q;
		}
	}
	return s;
}


double bisection_eig(vector<vector<double>>A, int m)
{
	//二分法求指定特征值
	//输入：矩阵A，特征值大小序m
	//输出：第m个特征值
	int n = A.size(),times=0;
	vector<double> diag(n), sub_diag(n);
	sub_diag[0] = 0;
	for (int i = 0; i < n; i++)
	{
		if (i != n - 1)
			sub_diag[i+1] = A[i][i + 1];
		diag[i] = A[i][i];
	}
	double r = matrix_infnorm(A),l=-r, mid = (l + r) / 2,length=2*r;
	double tol = 1e-7;
	
	while (length>tol)
	{
		mid = (l + r) / 2;
		if (number_changes(diag, sub_diag, mid) >= m)
			r = mid;
		else
			l = mid;
		length = length / 2;
		times++;
	}
	cout << "二分法迭代次数：" << times;
	return mid;
}




vector<double>inverse_power_method(vector<vector<double>> A, double mu)
{
	// 利用反幂法计算特征向量
	// 输入：矩阵A，特征值估计mu
	// 输出：特征向量 b
	int n = A.size();
	vector<double> b(n);
	for (int i = 0; i < n; i++)
		b[i] = 1.0;
	for (int i = 0; i < n; i++)
		A[i][i] = A[i][i] - mu;

	vector<int> u(n);
	gauss_elim_col_pivoting(A, u);

	// 反幂法
	int times = 0,k=0;
	double eig1 = 0, eig2 = 0.1;
	while (fabs(eig1 - eig2) > 1e-3)
	{
		times++;
		eig1 = eig2;
		gauss_col_piv_solvequa(A, u, b);
		vector_infnorm(b, k);
		eig2 = b[k];
		for (int j = 0; j < n; j++)
			b[j] = b[j] / eig2;
	}
	/*cout << "迭代次数为：" << times << endl;*/
	return b;
}



//SVD算法所要求的函数

tuple<vector<vector<double>>, vector<double>, vector<vector<double>>,vector<double>> bidiagnize(vector<vector<double>>& A)
{
	//二对角化A
	//输入：矩阵A
	//输出：左乘的household变换集合P,b1，右乘的Household变换集合Q,b2
	int m=A.size(),n=A[0].size(),mymin=min(m,n);
	vector<vector<double>> P, Q;
	vector<double> b1, b2;
	tuple<vector<double>, double> household;
	for (int j = 0; j < mymin; j++)
	{
		vector<double>  v;
		double beta;
		if (m - 1 - j > 0)
		{
			household = house(get_col(A, j, j , m - 1));
			v = get<0>(household);
			beta = get<1>(household);
			for (int i = j; i < n; i++)
			{
				double s = beta * get_innerproduct(v, get_col(A, i, j , m - 1));
				for (int k = j ; k < m; k++)
					A[k][i] = A[k][i] - s * v[k - j];
			}
			P.push_back(v);
			b1.push_back(beta);
		}

		if (n - 1 - j > 1)
		{
			household = house(get_row(A, j, j+1, n - 1));
			v = get<0>(household);
			beta = get<1>(household);
			for (int i = j; i < m; i++)
			{
				double s = beta * get_innerproduct(v, get_row(A, i, j + 1, n - 1));
				for (int k = j + 1; k < n; k++)
					A[i][k] = A[i][k] - s * v[k - j - 1];
			}
			Q.push_back(v);
			b2.push_back(beta);
		}
	}
	return{ P,b1,Q,b2 };


}

vector<double> givens_l(double a, double b)
{
	//计算将（a，b）^T中b化0的givens变换
	vector<double> bina(2, 0);
	if (b == 0)
		bina[0] = 1;
	else
	{
		if (abs(b) > abs(a))
		{
			if (b < 0)
			{
				double tao = a / b;
				bina[1] = -1 / sqrt(1 + tao * tao);
				bina[0] = bina[1] * tao;
			}
			else
			{
				double tao = a / b;
				bina[1] = 1 / sqrt(1 + tao * tao);
				bina[0] = bina[1] * tao;
			}
		}
		else
		{
			if (a < 0)
			{
				double tao = b / a;
				bina[0] = -1 / sqrt(1 + tao * tao);
				bina[1] = bina[0] * tao;
			}
			else
			{
				double tao = b / a;
				bina[0] = 1 / sqrt(1 + tao * tao);
				bina[1] = bina[0] * tao;
			}
			
		}
	}
	return bina;
}

vector<double> givens_r(double a, double b)
{
	//计算将（a，b）中b化0的givens变换
	vector<double> bina = givens_l(a, b);
	bina[1] = -bina[1];
	return bina;
}

void givens_times_l(vector<double> givens, double& a, double& b)
{
	//求Givens作用后的二维向量的值（左乘）
	double temp = a;
	a = givens[0] * a + givens[1] * b;
	b = -givens[1] * temp + givens[0] * b;
}

void givens_times_r(vector<double> givens, double& a, double& b)
{
	//求Givens作用后的二维向量的值（右乘）
	double temp = a;
	a = givens[0] * a - givens[1] * b;
	b = givens[1] * temp + givens[0] * b;
}

tuple<vector<vector<double>>, vector<vector<double>>> SVD_itr(vector<double>& D, vector<double>& S)
{
	//带wilkinson位移的SVD迭代
	//输入：主对角元D，次对角元S
	//输出：一次迭代后的D,S；左乘的Givens变换集合P(不考虑其转置），右乘的Givens变换集合Q
	int n = D.size();
	double alpha, beta, delta, mu, y, z;
	
	if (n != 2)
	{
		  alpha = D[n - 1] * D[n - 1] + S[n - 2] * S[n - 2], delta = (D[n - 2] * D[n - 2] + S[n - 3] * S[n - 3] - alpha) / 2,
			beta = D[n - 2] * S[n - 2], mu = alpha - beta * beta / (delta + sgn(delta) * sqrt(pow(delta, 2) + pow(beta, 2))),
			y = D[0] * D[0] - mu, z = D[0] * S[0];
	}
	if (n == 2)
	{
		alpha = D[n - 1] * D[n - 1] + S[n - 2] * S[n - 2], delta = (D[n - 2] * D[n - 2]  - alpha) / 2,
			beta = D[n - 2] * S[n - 2], mu = alpha - beta * beta / (delta + sgn(delta) * sqrt(pow(delta, 2) + pow(beta, 2))),
			y = D[0] * D[0] - mu, z = D[0] * S[0];
	}
	vector<vector<double>> P,Q;
	vector<double> given;
	for (int k = 0; k < n - 1; k++)
	{
		given = givens_r(y, z);
		if (k != 0)
			S[k - 1] = given[0] * y - given[1] * z;
		Q.push_back(given);
		givens_times_r(given, D[k], S[k]);
		y = D[k]; z = 0;
		givens_times_r(given, z, D[k+1]);
		given=givens_l(y, z);
		D[k] = given[0] * y + given[1] * z;
		P.push_back(given);
		if (k < n - 2)
		{
			givens_times_l(given, S[k], D[k + 1]);
			y = S[k]; z = 0;
			givens_times_l(given, z, S[k + 1]);
		}
		else
			givens_times_l(given, S[k], D[k + 1]);
	}
	return{ P,Q };
}

vector<vector<double>> modify(vector<double>& D, vector<double>& S,int ind)
{
	//将对角元为0的情况调整成整行为0
	//输入：对角元D，次对角元S，0元素的位置ind
	//输出：调整后的元素及左乘的givens阵
	int n = D.size();
	double z=S[ind];
	S[ind] = 0;
	vector<double> given;
	vector<vector<double>> P;
	for (int i = ind+1; i < n; i++)
	{
		given = givens_l(D[i], z);
		given[1] = -given[1];
		givens_times_l(given, z, D[i]);
		
		z = 0;
		if (i < n - 1)
			givens_times_l(given, z, D[i + 1]);
		P.push_back(given);
	}
	return P;
}

tuple<vector<double>,vector<vector<double>>, vector<vector<double>>> SVD(vector<vector<double>> A)
{
	//SVD算法最终版(长>=宽）
	//输入：奇异值分解矩阵A
	//输出：奇异值集合D，左正交阵P 和 右正交阵Q。

	//初始化
	int m = A.size(), n = A[0].size();
	vector<vector<double>> P(m, vector<double>(m, 0)), Q(n, vector<double>(n, 0));
	vector<double> D, S;
	for (int i = 0; i < m; i++)
		P[i][i] = 1;
	for (int i = 0; i < n; i++)
		Q[i][i] = 1;
	
	
	//二对角化
	tuple<vector<vector<double>>, vector<double>, vector<vector<double>>, vector<double>> hs = bidiagnize(A);
	vector<vector<double>> H1 = get<0>(hs), H2 = get<2>(hs);

	vector<double> b1 = get<1>(hs), b2 = get<3>(hs);  //household阵对应向量及系数

	for (int j = 0; j < b1.size(); j++)
	{
		vector<double> v = H1[j];
		double beta = b1[j];
		for (int i = 0; i < P.size(); i++)
		{
			double s = beta * get_innerproduct(v, get_col(P, i, j, m - 1));
			for (int k = j; k < m; k++)
				P[k][i] = P[k][i] - s * v[k - j];
		}
	}//左乘积累

	for (int j = 0; j < b2.size(); j++)
	{
		vector<double> v = H2[j];
		double beta = b2[j];
		for (int i = 0; i < Q.size(); i++)
		{
			double s = beta * get_innerproduct(v, get_row(Q, i, j + 1, n - 1));
			for (int k = j + 1; k < n; k++)
				Q[i][k] = Q[i][k] - s * v[k - j - 1];
		}
	}//右乘积累

	for (int j = 0; j < n; j++)
	{
		if (j != n - 1)
			S.push_back(A[j][j + 1]);
		D.push_back(A[j][j]);
	} //分离出对角元D和次对角元S


	//收敛性判定+SVD迭代
	double acc = 1e-15;
	int p, q,times=0;
	
	while (1)
	{
		double inf_norm = bidiag_inf_norm(D, S);
		for (int i = 0; i < n-1; i++)
		{
			if(abs(S[i])<=acc*(abs(D[i])+abs(D[i+1])))
				S[i]=0.0;
			if (abs(D[i]) <= acc * inf_norm)
				D[i] = 0;
		}
		if (abs(D[n - 1]) <= acc * inf_norm) D[n - 1] = 0;

		for (q = n - 1; q > 0; q--)
			if (S[q-1] != 0) break;
		if (q == 0)
		{
			cout << "SVD迭代次数： " << times << endl;
			cout << "机器精度" << acc << endl;
			square_transpose(P);
			square_transpose(Q);
			return { D,P,Q }; //终止条件

		}
		for (p = q; p > 0; p--)
			if (S[p - 1] == 0) break;

		int flag = 0;
		vector<double>D_sub(D.begin() + p, D.begin() + q + 1), S_sub(S.begin() + p, S.begin() + q);
		//检查是否有对角元为0，有则调整
		for(int i=D_sub.size()-2;i>=0;i-- )
			if (D_sub[i] == 0)
			{
				flag = 1;
				vector<vector<double>> givens = modify(D_sub, S_sub,i);
				for (int j = i + 1; j <= q - p; j++)
					Givens_matrix_times_left(P, givens[j - i - 1], p + i, p + j);
				for (int j = 0; j <= q - p - 1; j++)
				{
					D[j + p] = D_sub[j];
					S[j + p] = S_sub[j];
				}
				D[q] = D_sub[q - p];
				break;
			}
		if (flag == 1) continue;

		//SVD迭代
		tuple<vector<vector<double>>, vector<vector<double>>> two_givens=SVD_itr(D_sub, S_sub);
		times++;
		vector<vector<double>> givens_left = get<0>(two_givens), givens_right = get<1>(two_givens);

		for (int j = 0; j <= q - p - 1; j++)
		{
			Givens_matrix_times_left(P, givens_left[j], p + j, p + j + 1);
			Givens_matrix_times_right(Q, givens_right[j], p + j, p + j + 1);
			//积累正交阵
			D[j + p] = D_sub[j];
			S[j + p] = S_sub[j];
		}
		D[q] = D_sub[q - p];

	}
}
