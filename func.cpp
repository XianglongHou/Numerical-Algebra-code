#include<iostream>
#include<vector>
#include"helper.h"
#include<math.h>
#include<tuple>
#include<complex>

using namespace std;
void forward_subs(vector<vector<double>>& L, vector<double>& b)
{
	//ǰ�����������Ƿ�����
	//���룺�����Ǿ���L�������Ҷ���b
	//��������̵Ľ�b
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
	//�Խ�ԪΪ1��ǰ����
	//���룺��λ�����Ǿ���L�������Ҷ���b
	//��������̵Ľ�b
	int i, j, n = L.size();
	for (i = 0; i < n; i++)
		for (j = i + 1; j < n; j++)
			b[j] = b[j] - b[i] * L[j][i];
	//cnm!!!
}

void back_subs(vector<vector<double>>& U, vector<double>& b)
{
	//�ش���
	//���룺�����Ǿ���U�������Ҷ���b
	//��������̵Ľ�b
	int n = U[0].size();
	for (int i = n - 1; i >= 1; i--) //�㵽0�ʹ���
	{
		b[i] = b[i] / U[i][i];
		for (int j = i - 1; j >= 0; j--)
			b[j] = b[j] - b[i] * U[j][i];
	}
	b[0] = b[0] / U[0][0];
}

void back_subs1(vector<vector<double>>& U, vector<double>& b)
{
	//�Խ�ԪΪ1�Ļش���
	//���룺��λ�����Ǿ���U�������Ҷ���b
	//��������̵Ľ�b
	int n = U.size();
	for (int i = n - 1; i >= 0; i--)
		for (int j = i - 1; j >= 0; j--)
			b[j] = b[j] - b[i] * U[j][i];
}

void gauss_elim(vector<vector<double>>& A)
{
	//Gauss��ȥ����LU�ֽ⣩
	//���룺1��n-1��˳������ʽ��0�ķ���A
	//�����LU�ֽ���L����������ǲ��֣��������Խ��ߣ���U����������ǲ��֣����洢��A��
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
	//����gauss��ȥ���ⷽ��
	//���룺�ɽ���LU�ֽ�ľ���A�������Ҷ���b
	//��������̵Ľ�b
	gauss_elim(A);
	forward_subs1(A, b);
	back_subs(A, b);
}

void gauss_elim_full_pivoting(vector<vector<double>>& A, vector<int>& u, vector<int>& v)
{
	//ȫ��ԪGauss��ȥ��
	//���룺����A���ѳ�ʼ��������u��v
	//�������LU�ֽ���L����������ǲ��֣��������Խ��ߣ���U����������ǲ��֣����洢��A��
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
	////ȫ��ԪGauss��ȥ���ⷽ��
	//���룺����A,�����Ҷ���b
	//��������̵Ľ�b
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
	//����ԪGauss��ȥ��
	//���룺����A���ѳ�ʼ��������u
	//�������LU�ֽ���L����������ǲ��֣��������Խ��ߣ���U����������ǲ��֣����洢��A��
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
	//����Ԫ���ⷽ��
	//���룺����A,�����Ҷ���b
	//��������̵Ľ�b
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
	//����֪A������ԪLU�ֽ������£��ⷽ��
	for (int i = 0; i < u.size(); i++)
		swap(b[i], b[u[i]]);
	forward_subs1(A, b);  //subs1!!!!
	back_subs(A, b);
}


void cholesky_decomp(vector<vector<double>>& A)
{
	//�Գ��������׼Cholesky�ֽ�
	//���룺�����Գ���A
	//�����Cholesky����L���洢��A����������
	int n = A.size();
	for (int k = 0; k < n; k++)
	{
		A[k][k] = sqrt(A[k][k]);
		for (int j = k + 1; j < n; j++)
			A[j][k] = A[j][k] / A[k][k]; //�õ���k��Ԫ�ص�����ֵ
		for (int j = k + 1; j < n; j++)
			for (int i = j; i < n; i++)
				A[i][j] = A[i][j] - A[i][k] * A[j][k];
	}
}

void cholesky_decomp_solvequa(vector<vector<double>> A, vector<double>& b)
{
	//����Cholesky�ֽ�ⷽ��
	//���룺������A�������Ҷ���b
	//��������̵Ľ�b
	cholesky_decomp(A);
	forward_subs(A, b);
	square_transpose(A);
	back_subs(A, b);
}

void modified_cholesky_decomp(vector<vector<double>>& A)
{
	//�Ľ���ƽ������A=LDL^T
	//���룺�����Գ���A
	//�������λ��������L,�Խ���D (L���Խ�Ԫ��洢��A�������ǣ�D������A�ĶԽ����ϣ�
	int n = A.size();
	vector<double> v(n, 0);
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i <= j - 1; i++)   //�����������Ҫȡ��
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
	//�øĽ���ƽ�������ⷽ��
	//���룺������A�������Ҷ���b
	//��������̵Ľ�b
	modified_cholesky_decomp(A);
	forward_subs1(A, b);
	for (int i = 0; i < b.size(); i++)
		b[i] = b[i] / A[i][i];
	square_transpose(A);
	back_subs1(A, b);
}

double matrix_1norm(vector<vector<double>> B)
{
	//���ƾ����1�������Ż�����
	//���룺����B
	//���������1�����Ĺ���ֵ
	int n = B.size();
	int index = 0;
	vector<vector<double>> Bt(B);
	square_transpose(Bt);//B��ת��
	vector<double>x(n, 1 / double(n)); //����ط�float
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
	//���ƾ������1�������Ż�����
	//���룺����A
	//���������A^(-1)��1����
	int n = A.size();
	int index = 0;
	vector<vector<double>> At(A);
	square_transpose(At);//A��ת��
	vector<double>x(n, 1 / double(n)); //����ط�double
	vector<double>w(n), v(n), z(n);
	vector<int> u1(A.size() - 1, 0);
	gauss_elim_col_pivoting(A, u1);
	vector<int> u2(At.size() - 1, 0);
	gauss_elim_col_pivoting(At, u2);//����Ԫ��Gauss�ֽ�


	while (1)
	{
		w.assign(x.begin(), x.end());
		for (int i = 0; i < u2.size(); i++)
			swap(w[i], w[u2[i]]);
		forward_subs1(At, w);
		back_subs(At, w);//���÷ֽ�ⷽ��At*w=x
		for (int i = 0; i < n; i++)
			v[i] = sgn(w[i]);
		z.assign(v.begin(), v.end());
		for (int i = 0; i < u1.size(); i++)
			swap(z[i], z[u1[i]]);
		forward_subs1(A, z);
		back_subs(A, z);//���÷ֽ�ⷽ��A*z=v
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
	//���ƾ����������
	//���룺����A
	//���������A��������
	vector<vector<double>>At(A);
	square_transpose(At);
	return(matrix_1norm(At) * matrix_inverse_infnorm(A));
}

double equasolv_error(vector<vector<double>>A, vector<double>b, vector<double>x)
{
	//�����ľ��ȹ���
	//���룺����A�������Ҷ���b�����̵ļ����x
	//����������ľ��ȹ���
	int n = A.size();
	double b_norm = vector_infnorm(b);
	vector<double> r(n);
	r = vector_minus(b, mat_vec_mult(A, x));
	double r_norm = vector_infnorm(r);
	return inf_condition_num(A) * r_norm / b_norm;
}

//��������Ҫ�õ��ĺ���

tuple<vector<double>, double> house(vector<double>x)
{
	//��ʹ��������һ��Ԫ��ȫ��Ϊ0��household�任
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
	//QR�ֽ�
	//���룺ԭ����A,�ѳ�ʼ�����ڴ���household����ϵ��������d
	//�����A�����ǲ��ִ洢Household�������ǲ��ִ��R��d���household����ϵ��
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
	//����QR�ֽ�ⷽ��
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
	}//��b����household�任
	back_subs(A, b);
}

vector<double> QR_least_square(vector<vector<double>> A, vector<double> b)
{
	//����С���˽�
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
	}//��b����household�任
	vector<double> r = sub_vector(b, 0, n - 1);
	back_subs(A, r);
	return r;
}

void jacobi_itr(vector<vector<double>> A, vector<double> b, vector<double>& x, double tol = 1e-6)
{
	//jacobi ����
	//���룺ϵ������A,������ʼֵb,������ֹ����tol����x_{k+1}-x_k���������
	// ��������̵Ľ�b������ӡ��������
	double d, err = 10;
	int count = 0;
	for (int i = 0; i < A.size(); i++)
	{
		d = A[i][i];
		A[i][i] = 0;
		b[i] = b[i] / d;
		for (int j = 0; j < A.size(); j++)
			A[i][j] = -A[i][j] / d;
	}//�����������ƫ����(����˲������
	while (err > tol)
	{
		vector<double> temp(x);
		x = vector_plus(mat_vec_mult(A, x), b);
		err = vector_infnorm(vector_minus(x, temp));
		count++;
	}
	cout << "jacobi����������" << count << endl;
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
	}//�����������ƫ����
	do
	{
		vector<double> temp(x);
		for (int i = 0; i < A.size(); i++)
			x[i] = get_innerproduct(A[i], x) + b[i];
		err = vector_infnorm(vector_minus(x, temp));
		count++;
	} while (err > tol);
	cout << "G-S����������" << count << endl;
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
	}//�����������ƫ����
	do
	{
		vector<double> temp(x);
		for (int i = 0; i < A.size(); i++)
			x[i] = (1 - w) * x[i] + w * (get_innerproduct(A[i], x) + b[i]);
		err = vector_infnorm(vector_minus(x, temp));
		count++;
	} while (err > tol);
	cout << "SOR����������" << count << "  ���ɳ�����" << w << endl;
}

void CG_method(vector<vector<double>> A, vector<double> b, vector<double>& x, double tol = 1e-7)
{
	//�����ݶȷ�
	//���룺ϵ������������A���Ҷ���b����ʼ����ֵx�����̶�tol
	//�������x
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
	cout << "�����ݶȷ�����������" << count << endl;
}

double find_largest_root(const vector<double> a, int times = 1000)
{
	//�ݷ��Ҷ���ʽ��������д�ɾ���˷���ʽ��������ѷ���˵��������ı��ʽ����ֱ�ӵ�����
	//���룺����ʽ����ϵ������0�ε�n-1�Σ�
	//���������ʽ��ģ����
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
	//��Hessenberg�ֽ⣬Hessenberg�����A��,Hessenberg�����A�У�Householder�任��v��beta����V��b��
	//���룺����Hessenberg���ľ���A
	//�����Hessenberg����A��Householder�任��v��beta��Ԫ����ʽ����
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
	//˫�ز�λ�Ƶ�QR��,����H�У���ֱ�Ӷ�H12��H23���任
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
}//��ʽQR�㷨

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
	//Givens�������A
	//���룺A,Givens����Q��������i,j
	//������任��ľ���A
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


// givens �任�ҳ�
void Givens_matrix_times_right(vector<vector<double>>& A, vector<vector<double>>& Q, int i, int j)
{
	// givens �任�ҳ�
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
	// givens �任�ҳ�
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
	//givens�任���
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
	//givens�任���
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
	cout << "��������Ϊ��" << times << endl;
}





int number_changes(vector<double>& A, vector<double>& B, double mu)
{
	//��mu�ı����
	//���룺�Խ�ԪA���ζԽ�ԪB��mu
	//����������s
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
	//���ַ���ָ������ֵ
	//���룺����A������ֵ��С��m
	//�������m������ֵ
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
	cout << "���ַ�����������" << times;
	return mid;
}




vector<double>inverse_power_method(vector<vector<double>> A, double mu)
{
	// ���÷��ݷ�������������
	// ���룺����A������ֵ����mu
	// ������������� b
	int n = A.size();
	vector<double> b(n);
	for (int i = 0; i < n; i++)
		b[i] = 1.0;
	for (int i = 0; i < n; i++)
		A[i][i] = A[i][i] - mu;

	vector<int> u(n);
	gauss_elim_col_pivoting(A, u);

	// ���ݷ�
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
	/*cout << "��������Ϊ��" << times << endl;*/
	return b;
}



//SVD�㷨��Ҫ��ĺ���

tuple<vector<vector<double>>, vector<double>, vector<vector<double>>,vector<double>> bidiagnize(vector<vector<double>>& A)
{
	//���Խǻ�A
	//���룺����A
	//�������˵�household�任����P,b1���ҳ˵�Household�任����Q,b2
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
	//���㽫��a��b��^T��b��0��givens�任
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
	//���㽫��a��b����b��0��givens�任
	vector<double> bina = givens_l(a, b);
	bina[1] = -bina[1];
	return bina;
}

void givens_times_l(vector<double> givens, double& a, double& b)
{
	//��Givens���ú�Ķ�ά������ֵ����ˣ�
	double temp = a;
	a = givens[0] * a + givens[1] * b;
	b = -givens[1] * temp + givens[0] * b;
}

void givens_times_r(vector<double> givens, double& a, double& b)
{
	//��Givens���ú�Ķ�ά������ֵ���ҳˣ�
	double temp = a;
	a = givens[0] * a - givens[1] * b;
	b = givens[1] * temp + givens[0] * b;
}

tuple<vector<vector<double>>, vector<vector<double>>> SVD_itr(vector<double>& D, vector<double>& S)
{
	//��wilkinsonλ�Ƶ�SVD����
	//���룺���Խ�ԪD���ζԽ�ԪS
	//�����һ�ε������D,S����˵�Givens�任����P(��������ת�ã����ҳ˵�Givens�任����Q
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
	//���Խ�ԪΪ0���������������Ϊ0
	//���룺�Խ�ԪD���ζԽ�ԪS��0Ԫ�ص�λ��ind
	//������������Ԫ�ؼ���˵�givens��
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
	//SVD�㷨���հ�(��>=��
	//���룺����ֵ�ֽ����A
	//���������ֵ����D����������P �� ��������Q��

	//��ʼ��
	int m = A.size(), n = A[0].size();
	vector<vector<double>> P(m, vector<double>(m, 0)), Q(n, vector<double>(n, 0));
	vector<double> D, S;
	for (int i = 0; i < m; i++)
		P[i][i] = 1;
	for (int i = 0; i < n; i++)
		Q[i][i] = 1;
	
	
	//���Խǻ�
	tuple<vector<vector<double>>, vector<double>, vector<vector<double>>, vector<double>> hs = bidiagnize(A);
	vector<vector<double>> H1 = get<0>(hs), H2 = get<2>(hs);

	vector<double> b1 = get<1>(hs), b2 = get<3>(hs);  //household���Ӧ������ϵ��

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
	}//��˻���

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
	}//�ҳ˻���

	for (int j = 0; j < n; j++)
	{
		if (j != n - 1)
			S.push_back(A[j][j + 1]);
		D.push_back(A[j][j]);
	} //������Խ�ԪD�ʹζԽ�ԪS


	//�������ж�+SVD����
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
			cout << "SVD���������� " << times << endl;
			cout << "��������" << acc << endl;
			square_transpose(P);
			square_transpose(Q);
			return { D,P,Q }; //��ֹ����

		}
		for (p = q; p > 0; p--)
			if (S[p - 1] == 0) break;

		int flag = 0;
		vector<double>D_sub(D.begin() + p, D.begin() + q + 1), S_sub(S.begin() + p, S.begin() + q);
		//����Ƿ��жԽ�ԪΪ0���������
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

		//SVD����
		tuple<vector<vector<double>>, vector<vector<double>>> two_givens=SVD_itr(D_sub, S_sub);
		times++;
		vector<vector<double>> givens_left = get<0>(two_givens), givens_right = get<1>(two_givens);

		for (int j = 0; j <= q - p - 1; j++)
		{
			Givens_matrix_times_left(P, givens_left[j], p + j, p + j + 1);
			Givens_matrix_times_right(Q, givens_right[j], p + j, p + j + 1);
			//����������
			D[j + p] = D_sub[j];
			S[j + p] = S_sub[j];
		}
		D[q] = D_sub[q - p];

	}
}
