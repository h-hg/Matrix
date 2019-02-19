#include <iostream>
#include <stdexcept>
#include <cstring>
#include <algorithm>
#include <sstream>
#ifndef Matrix_H_
#define Matrix_H_
template <typename T>
class Matrix
{
	template <typename U>
	friend std::ostream& operator<<(std::ostream &os,const Matrix<U> &matrix);
	template <typename U>
	friend std::istream& operator>>(std::istream &is,Matrix<U> &matrix);
	template <typename U>
	friend const Matrix<U> operator*(const U &lhs,const Matrix<U> &rhs);
private:
	//allocate the memory according to the row and column
	void build();
	//free the memory
	void clear();
	void copy(const Matrix<T> &matrix);
	void movecopy(Matrix<T> &matrix);
	
	//This function convert the matrix to Row-Echlon Form, and return the number of pivot
	//pivot_pos is a pointer to data of the column of the pivot, it should not be less than sizeof(T)*column
	//该函数将矩阵转化为阶梯型，返回主元的个数，limitCol表示高斯消元最多消到limitCol列，pivot_pos[i]储存第i行主元的列位置
	size_t ref(size_t limitCol,size_t *pivot_pos);//Gauss-Jordan Elimination, 高斯若当消元
	
	size_t ref(size_t *pivot_pos){return ref(column,pivot_pos);};

	//convert Row-Echlon Form to Reduced Row-Echlon Form
	//将阶梯型矩阵转化为标准型
	void ref_to_rref(size_t *pivot_pos,size_t pivot_num);

	//获取逆序数
	template <typename U>
	static size_t getInverseNumber(U const arr[],size_t n);
	//将多项式转为字符串
	static std::string poly_to_string(T const poly[],size_t n);

	static void multi(T ans[],T tmp[],size_t n,const T &multiplier);
	
	static void extend(T const elems[],size_t elems_num,T ans[]);

	size_t row;
	size_t column;
	T **data;
public:
	Matrix(size_t row,size_t column):row(row),column(column){build();}
	Matrix(const Matrix<T> &matrix){copy(matrix);}
	Matrix(Matrix<T> &&matrix){movecopy(matrix);}
	~Matrix(){clear();}
	Matrix& operator=(const Matrix<T> &matrix);
	Matrix& operator=(Matrix<T> &&matrix);

	const Matrix<T> operator+(const Matrix<T> &rhs)const;
	Matrix<T>& operator+=(const Matrix<T> &rhs);
	const Matrix<T> operator-(const Matrix<T> &rhs)const;
	Matrix<T>& operator-=(const Matrix<T> &rhs);
	const Matrix<T> operator*(const Matrix<T> &rhs)const;

	Matrix<T>& operator*=(const T &rhs);
	const Matrix<T> operator*(const T&rhs)const;

	bool operator!=(const Matrix &rhs)const;
	bool operator==(const Matrix &rhs)const;

	T& operator()(size_t row,size_t column){return data[row][column];}
	T operator()(size_t row,size_t column)const{return data[row][column];}
	size_t getCol()const{return column;}
	size_t getRow()const{return row;}

	void rowInterchange(const size_t &row1,const size_t &row2);
	void swap(Matrix<T> &matrix2);

	void rref(size_t limitCol);
	void rref(){rref(column);};

	T det()const;

	//转置矩阵
	Matrix<T> transpose()const;
	
	//inverse, 逆矩阵
	Matrix<T> inv()const;
	
	//求秩
	size_t rank()const;

	//伴随矩阵
	Matrix<T> adjugate()const;/* 未测试 */

	//characteristic polynomial,特征多项式
	const std::string poly()const;/* 未测试 */
};

template <typename T>
template <typename U>
size_t Matrix<T>::getInverseNumber(U const arr[],size_t n)
{
	size_t ans = 0;
	for(size_t i = 1;i < n;++i)
		for(size_t j = 0;j < i;++j)
			if(arr[j] > arr[i])
				++ans;
	return ans;
}
template <typename U>
std::ostream& operator<<(std::ostream &os,const Matrix<U> &matrix)
{
	for(size_t i = 0;i < matrix.row;++i)
	{
		for(size_t j = 0; j < matrix.column - 1;++j)
			os << matrix.data[i][j] << '\t';
		os << matrix.data[i][matrix.column-1] << '\n';
	}
	return os;
}
template <typename U>
std::istream& operator>>(std::istream &is,Matrix<U> &matrix)
{
	for(size_t i = 0;i < matrix.row;++i)
		for(size_t j = 0;j < matrix.column;++j)
			is >> matrix.data[i][j];
	return is;
}
template<typename T>
void Matrix<T>::build()
{
	data = new T*[row];
	for(size_t sub = 0;sub < row;++sub)
		data[sub] = new T[column];
}
template<typename T>
void Matrix<T>::clear()
{
	for(size_t sub = 0;sub < row;++sub)
		delete[] data[sub];
	delete[] data;
}
template<typename T>
void Matrix<T>::copy(const Matrix<T> &matrix)
{
	row = matrix.row;
	column = matrix.column;
	build();
	for(auto sub = 0;sub < row;++sub)
		memcpy(data[sub],matrix.data[sub],sizeof(T)*column);
}
template<typename T>
void Matrix<T>::movecopy(Matrix<T> &matrix)
{
	row = matrix.row;
	column = matrix.column;
	data = matrix.data;
	matrix.data = nullptr;
}
template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &matrix)
{
	if(this == &matrix)
		return *this;
	clear();
	copy(matrix);
	return *this;
}
template <typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T> &&matrix)
{
	if(this == &matrix)
		return *this;
	clear();
	movecopy(matrix);
	return *this;
}
template <typename T>
bool Matrix<T>::operator!=(const Matrix<T> &rhs)const
{
	if(row != rhs.row || column != rhs.column)
		return true;
	for(size_t i = 0;i < row;++i)
		for(size_t j = 0;j < column;++j)
			if(data[i][j] != rhs.data[i][j])
				return true;
	return false;
}
template <typename T>
bool Matrix<T>::operator==(const Matrix<T> &rhs)const
{
	return !(*this != rhs);
}
template<typename T>
void Matrix<T>::rowInterchange(const size_t &row1,const size_t &row2)
{
	std::swap(data[row1],data[row2]);
}
template<typename T>
void Matrix<T>::swap(Matrix<T> &matrix2)
{
	std::swap(row,matrix2.row);
	std::swap(column,matrix2.column);
	std::swap(data,matrix2.data);
}
namespace std{
	template <typename T>
	void swap(Matrix<T> &matrix1,Matrix<T> &matrix2)
	{
		matrix1.swap(matrix2);
	}
}
template<typename T>
size_t Matrix<T>::ref(size_t limitCol,size_t *pivot_pos)
{
	size_t currentCol = 0,currentRow = 0;//currentRow <=> pivot_num;
	//the followed steps is converting the matrix to Row-Echelon Form
	for(;currentRow < row && currentCol < limitCol;++currentCol)
	{
		size_t non_zero_row = currentRow;
		//find first row which has non-zero entry at column currentCol
		for(;non_zero_row < row;++non_zero_row)
			if(data[non_zero_row][currentCol])
				break;
		if(non_zero_row == row)//can't find the non-zero iterm, enter the next column
			continue;
		if(non_zero_row != currentRow)
			rowInterchange(non_zero_row,currentRow);
		//minus
		for(auto i = currentRow + 1;i < row;++i)
		{
			auto foo = data[i][currentCol] / data[currentRow][currentCol];
			for(auto j = currentCol;j < column;++j)
				data[i][j] -=  foo * data[currentRow][j];
		}
		pivot_pos[currentRow++] = currentCol;//recod tht column of pivot and enter next row
	}
	return currentRow;
}
template <typename T>
void Matrix<T>::ref_to_rref(size_t *pivot_pos,size_t pivot_num)
{
	//pivot_num <=> currentRow
  	while(pivot_num--)
	{
		size_t currentCol = pivot_pos[pivot_num];//first non-zero entry 
		//minus
		for(size_t i = 0; i < pivot_num;++i)
		{
			auto foo = data[i][currentCol] / data[pivot_num][currentCol];
			if(data[i][currentCol])
				for(size_t j = currentCol;j < column;++j)
					data[i][j] -= foo * data[pivot_num][j];
		}
		// unit 1
		auto foo = data[pivot_num][currentCol];
		for(size_t j = currentCol;j < column;++j)
			data[pivot_num][j] /= foo;
	}
}
template <typename T>
void Matrix<T>::rref(size_t limitCol)
{
	size_t *pivot_pos = new size_t[row];
	size_t pivot_num = ref(limitCol,pivot_pos);
	ref_to_rref(pivot_pos,pivot_num);
	delete[] pivot_pos;
}
template <typename T>
T Matrix<T>::det()const
{
	if(column != row)
		throw std::logic_error("Only the square matrix can find the determinant");
	auto m_copy = *this;
	bool negative = false;
	T ans(1);
	for(size_t currentCol = 0,currentRow = 0;currentRow < m_copy.row && currentCol < m_copy.column;++currentCol)
	{
		size_t non_zero_row = currentRow;
		//find first row which has non-zero entry at column currentCol
		for(;non_zero_row < m_copy.row;++non_zero_row)
			if(m_copy.data[non_zero_row][currentCol])
				break;
		if(non_zero_row == row)//can't find the non-zero iterm, enter the next column
			return T(0);
		if(non_zero_row != currentRow){
			m_copy.rowInterchange(non_zero_row,currentRow);
			negative = !negative;
		}
		//minus
		for(auto i = currentRow + 1;i < row;++i)
		{
			auto foo = m_copy.data[i][currentCol] / m_copy.data[currentRow][currentCol];
			for(auto j = currentCol;j < column;++j)
				m_copy.data[i][j] -=  foo * m_copy.data[currentRow][j];
		}
		ans *=  m_copy.data[currentRow++][currentCol];//recod tht column of pivot and enter next row
	}
	return negative ? -ans : ans;
}
template<typename T>
Matrix<T> Matrix<T>::transpose()const
{
	Matrix ans(column,row);
	for(size_t i = 0;i < row;++i)
		for(size_t j = 0;j < column;++j)
			ans.data[j][i] = data[i][j];
	return ans;
}
template<typename T>
Matrix<T> Matrix<T>::inv()const
{
	if(row != column)
		throw std::logic_error("The operator of inverse must be square matrix.");
	Matrix<T> foo(row,column+column);
	for(size_t i = 0;i < row;++i)
	{
		memcpy(foo.data[i],data[i],sizeof(T)*column);
		foo.data[i][column + i] = T(1);
	}
	size_t *pivot_pos = new size_t[column];
	size_t pivot_num = foo.ref(column,pivot_pos);
	if(pivot_num != row)
	{
		delete[] pivot_pos;
		throw std::logic_error("The matrix can't be inversed.");
	}
	foo.ref_to_rref(pivot_pos,pivot_num);
	delete[] pivot_pos;
	Matrix<T> ans(row,column);
	for(size_t i = 0;i < row;++i)
		memcpy(ans.data[i],foo.data[i] + column,sizeof(T)*column);
	return ans;
}
//伴随矩阵
template <typename T>
Matrix<T> Matrix<T>::adjugate()const
{
	if(row != column)
		throw std::logic_error("The operator of adjugate must be square matrix.");
	Matrix<T> temp(row,column+column);
	for(size_t i = 0;i < row;++i)
	{
		memcpy(temp.data[i],data[i],sizeof(T)*column);
		temp.data[i][column + i] = T(1);
	}
	size_t *pivot_pos = new size_t[column];
	size_t pivot_num = temp.ref(column,pivot_pos);

	Matrix<T> P(row,column);//zero matrix
	if(pivot_num < row -1)//秩小于(n-1)，任何一个(n-1)阶的子式都是零
	{
		delete[] pivot_pos;
		return P;//return zero-matrix
	}
	temp.ref_to_rref(pivot_pos,pivot_num);
	delete[] pivot_pos;
	for(size_t i = 0;i < row;++i)
		memcpy(P.data[i],temp.data[i] + column,sizeof(T)*column);
	//P == (*this)^(-1)
	if(pivot_num == row)//row == n
	{
		T detA = (*this).det();
		return P * detA;
	}
	else//rank == row -1
	{
		Matrix<T> J(row,column);
		for(int i = 0;i < row;++i)
			memcpy(J.data[i],temp.data[i],sizeof(T)*column);
		Matrix<T> JT = J.transpose();
		//Matrix<T> temp2(row,column + column);
		for(size_t i = 0;i < row;++i)
		{
			memcpy(temp.data[i],JT.data[i],sizeof(T)*column);
			temp.data[i][column + i] = T(1);
		}
		temp.rref(column);
		Matrix<T> QT(row,column);
		for(size_t i = 0;i < row;++i)
			memcpy(QT.data[i],temp.data[i]+column,sizeof(T)*column);
		Matrix<T> Q = QT.transpose();
		Matrix _Ek(row,column);
		_Ek.data[row-1][column-1] = T(1);
		T t = T(1)/P.det()/Q.det();
		return t * Q * _Ek * P;
	}
}

// convert poly to the string
//poly[n-1] * x ^(n-1) + ... + poly[0]
template <typename T>
std::string Matrix<T>::poly_to_string(T const poly[],size_t n)
{
	std::ostringstream os;
	//找到第一个非零的项
	for(--n;n && !poly[n];--n)
	//避免 0 0 0 0没有输出
	if(!n) {
		os << poly[0];
		return os.str();
	}
	os << '(' << poly[n] << ")x^" << n;
	while(--n)
		if(poly[n])
			os << " + (" << poly[n] << ")x^" << n;
	if(poly[0])
		os << " + (" << poly[0] << ')';
	return os.str();
}

//calc ( ans[0] + ans[1] * λ^1 + ... + ans[ans_num-1] * λ^(ans_num-1) ) * ( multiplier - λ ), make sure ans and tmp has enough space
template <typename T>
void Matrix<T>::multi(T ans[],T tmp[],size_t n,const T &multiplier)
{
    //* multiplier
    for(size_t i = 0;i < n;++i)
        tmp[i] = ans[i] * multiplier;
    //* (-λ)
    for(size_t i = n;i > 0;--i)
        ans[i] = -ans[i-1];
	ans[0] = T(0);
    //add
    for(size_t i = 0;i < n;++i)
        ans[i] += tmp[i]; 
}
/* 
elems[0] * (elems[1] - λ) * (elems[2] - λ) * ... * (elems[elems_num-1] - λ)
->
ans[0] + ans[1] * λ^1 + ... + ans[ans_num-1] * λ^(ans_num-1)
 */
template <typename T>
void Matrix<T>::extend(T const elems[],size_t elems_num,T ans[])
{
	ans[0] = T(1);
	T *tmp = new T[elems_num];
	for(size_t i = 1;i < elems_num;++i)
		multi(ans,tmp,i,elems[i]);
	for(size_t i = 0;i < elems_num;++i)
		ans[i] *= elems[0];
	delete[] tmp;
}

//characteristic polynomial,特征多项式
template <typename T>
const std::string Matrix<T>::poly()const
{
	if(row != column)
		throw std::logic_error("The operator of poly must be square matrix.");
	//ans表示的式子：ans[0] + ans[1] * λ^1 + ... + ans[ans_num-1] * λ^(ans_num-1)
	T *ans = new T[row+1];//0-row 阶, 必须默认申请的时候全为零，不然要自己初始化
	//在求行列式，该数组的值由列指标的序列构成，ColNum[i]表示选取第i行第ColNum[i]列
	size_t *ColNum = new size_t[row];
	//初始化列指标，123...n
	for(size_t j = 0;j < column;++j)
		ColNum[j] = j;
	//elems表示式子：elems[0] * (elems[1] - λ) * (elems[2] - λ) * ... * (elems[elems_num-1] - λ)
	T *elems = new T[row+1];
	do{
		size_t elems_num = 0;
		//生成elems所表示的式子
		elems[0] = T(1);
		for(size_t i = 0; i < column;++i)
			if(i == ColNum[i])
				elems[++elems_num] = data[i][i];
			else
				elems[0] *= data[i][ColNum[i]];
		++elems_num;//正真的数组长度
		//计算逆序数
		if(getInverseNumber(ColNum,column) & 1)
			elems[0] = -elems[0];
		//expr_extend表示的式子同ans
		T *expr_extend = new T[elems_num];
		//展开并计算elems所表示的式子
		extend(elems,elems_num,expr_extend);
		//将展开结果加到最终结果
		for(size_t i = 0;i < elems_num;++i)
			ans[i] += expr_extend[i];
		delete[] expr_extend;
		
	}while(std::next_permutation(ColNum,ColNum+row));
	delete[] ColNum;
	delete[] elems;
	//将多项式转化为string
	auto retVal = poly_to_string(ans,row + 1);
	delete[] ans;
	return retVal;
}
template <typename T>
size_t Matrix<T>::rank()const
{
	Matrix<T> foo = *this;
	size_t *pivot_pos = new size_t[column];
	size_t pivot_num = foo.ref(column,pivot_pos); 
	delete [] pivot_pos;
	return pivot_num;
}
template <typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T> &rhs)
{
	if(row != rhs.row || column != rhs.column)
		throw std::logic_error("The operator+= of Matrix must have the same row and column.");
	for(size_t i = 0;i < row;++i)
		for(size_t j = 0;j < column;++j)
			data[i][j] += rhs.data[i][j];
	return *this;
}
template <typename T>
const Matrix<T> Matrix<T>::operator+(const Matrix<T> &rhs)const
{
	auto ans = *this;
	return ans += rhs;
}
template <typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T> &rhs)
{
	if(row != rhs.row || column != rhs.column)
		throw std::logic_error("The operator- of Matrix must have the same row and column.");
	for(size_t i = 0;i < row;++i)
		for(size_t j = 0;j < column;++j)
			data[i][j] -= rhs.data[i][j];
	return *this;
}
template <typename T>
const Matrix<T> Matrix<T>::operator-(const Matrix<T> &rhs)const
{
	auto ans = *this;
	return ans -= rhs;
}
template <typename T>
const Matrix<T> Matrix<T>::operator*(const Matrix<T> &rhs)const
{
	if(column != rhs.row)
		throw std::logic_error("The operator* of Matrix requires that lhs.row and rhs.column are equal.");
	Matrix<T> product(row,rhs.column);
	for(size_t i = 0;i < row;++i)
		for(size_t j = 0;j < rhs.column;++j)
		{
			product.data[i][j] = T(0);
			for(size_t k = 0;k < column;++k)
				product.data[i][j] += data[i][k] * rhs.data[k][j];
		}
	return product;
}
template <typename T>
Matrix<T>& Matrix<T>::operator*=(const T &rhs)
{
	for(size_t i = 0;i < row;++i)
		for(size_t j = 0;j < column;++j)
			data[i][j] *= rhs;
	return *this;
}
template <typename T>
const Matrix<T> Matrix<T>::operator*(const T&rhs)const
{
	auto ans = *this;
	return ans *= rhs;
}
template <typename U>
const Matrix<U> operator*(const U &lhs,const Matrix<U>& rhs)
{
	return rhs * lhs;
}

#endif