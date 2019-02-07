#include "Matrix.hpp"
#include "Fraction.h"
int main()
{
	std::cout << "本计算器支持分数，分数输入格式：1/3 表示三分之一，数值(分子、分母、整数)范围-2^31到 +(2^31-1)\n";
	int model,m,n,p;
    Matrix<Fraction> *pA = nullptr,*pB = nullptr;
	while(1)
	{
		std::cout << "请输入计算类型编号：\n"
					"1.行列式\n"
					"2.矩阵乘法\n"
					"3.矩阵的秩\n"
					"4.伴随矩阵\n"
					"5.逆矩阵\n"
					"6.线性方程组化标准型\n"
					"7.特征多项式\n"
					"8.退出\n";
		while(std::cin >> model,model != 0 && model != 8)
		{
			switch(model)
			{
			case 1://行列式
				std::cout << "请输入矩阵A的阶数：";
				std::cin >> n;
                pA = new Matrix<Fraction>(n,n);
				std::cout << "请输入矩阵 A" << n <<"×" << n << '\n';
				std::cin >> (*pA);
				std::cout << "行列式结果为：" << pA->det() << std::endl;
				break;
			case 2://矩阵乘法
				std::cout << "请输入矩阵 Am×n 和Bn×p 中的m n p：";
				std::cin >> m >> n >> p;
                pA = new Matrix<Fraction>(m,n);
                pB = new Matrix<Fraction>(n,p);
				std::cout << "请输入矩阵 A" << m << "×" << n << '\n';
				std::cin >> (*pA);
				std::cout << "请输入矩阵 B" << n << "×" << p << '\n';
				std::cin >> (*pB);
				std::cout << "乘法结果如下：\n" << (*pA) * (*pB);
				break;
			case 3://矩阵的秩
				std::cout << "请输入矩阵 Am×n 中的m n：";
				std::cin >> m >> n;
                pA = new Matrix<Fraction>(m,n);
				std::cout << "请输入矩阵 A" << m << "×" << n << '\n';
				std::cin >> (*pA);
				std::cout << "矩阵A的秩为:" << pA->rank() << std::endl;
				break;
			case 4://伴随矩阵
				std::cout << "请输入矩阵A的阶数：";
				std::cin >> n;
                pA = new Matrix<Fraction>(n,n);
				std::cout << "请输入矩阵 A" << n <<"×" << n << '\n';
				std::cin >> (*pA);
				std::cout << "伴随矩阵结果为：\n" << pA->adjugate() << std::endl;
				break;
			case 5://逆矩阵
				std::cout << "请输入矩阵的阶数：";
				std::cin >> n;
                pA = new Matrix<Fraction>(n,n);
				std::cout << "请输入矩阵 A" << n << "×" << n << '\n';
				std::cin >> (*pA);
				if(!pA->det())
					std::cout << "矩阵A 的行列式 = 0，不可逆\n";
				else
					std::cout << "矩阵A的逆矩阵如下：\n" << pA->inv();
				break;
			case 6://线性方程组
				std::cout << "请输入线性方程组的增广矩阵 Am×n 中的m n：";
				std::cin >> m >> n;
                pA = new Matrix<Fraction>(m,n);
				std::cout << "请输入矩阵 A" << m << "×" << n << '\n';
				std::cin >> (*pA);
				pA->rref(n-1);
				std::cout << "增广矩阵的标准矩阵如下：\n" << (*pA);
				break;
			case 7://特征多项式
				std::cout << "请输入矩阵A的阶数：";
				std::cin >> n;
                pA = new Matrix<Fraction>(n,n);
				std::cout <<"请输入矩阵 A" << n << "×" << n << '\n';
				std::cin >> (*pA);
				std::cout << "特征多项式为：" << pA->poly() << std::endl;
				break;
			default:
				std::cout << "不存在此选项\n";
				break;
			}
			std::cout << "请输入计算类型编号，当前类型编号为 " << model <<" ，输入 0 表示 回到功能选择,8 表示退出:";
		}
        delete pA;
        delete pB;
        pA = pB = nullptr;
		if(model == 8)
			break;
	}
	std::cout << "感谢使用，注意保存结果" << std::endl;
	system("pause");
	return 0;
}