#include "Matrix.hpp"
#include "Fraction.h"
int main()
{
	std::cout << "��������֧�ַ��������������ʽ��1/3 ��ʾ����֮һ����ֵ(���ӡ���ĸ������)��Χ-2^31�� +(2^31-1)\n";
	int model,m,n,p;
    Matrix<Fraction> *pA = nullptr,*pB = nullptr;
	while(1)
	{
		std::cout << "������������ͱ�ţ�\n"
					"1.����ʽ\n"
					"2.����˷�\n"
					"3.�������\n"
					"4.�������\n"
					"5.�����\n"
					"6.���Է����黯��׼��\n"
					"7.��������ʽ\n"
					"8.�˳�\n";
		while(std::cin >> model,model != 0 && model != 8)
		{
			switch(model)
			{
			case 1://����ʽ
				std::cout << "���������A�Ľ�����";
				std::cin >> n;
                pA = new Matrix<Fraction>(n,n);
				std::cout << "��������� A" << n <<"��" << n << '\n';
				std::cin >> (*pA);
				std::cout << "����ʽ���Ϊ��" << pA->det() << std::endl;
				break;
			case 2://����˷�
				std::cout << "��������� Am��n ��Bn��p �е�m n p��";
				std::cin >> m >> n >> p;
                pA = new Matrix<Fraction>(m,n);
                pB = new Matrix<Fraction>(n,p);
				std::cout << "��������� A" << m << "��" << n << '\n';
				std::cin >> (*pA);
				std::cout << "��������� B" << n << "��" << p << '\n';
				std::cin >> (*pB);
				std::cout << "�˷�������£�\n" << (*pA) * (*pB);
				break;
			case 3://�������
				std::cout << "��������� Am��n �е�m n��";
				std::cin >> m >> n;
                pA = new Matrix<Fraction>(m,n);
				std::cout << "��������� A" << m << "��" << n << '\n';
				std::cin >> (*pA);
				std::cout << "����A����Ϊ:" << pA->rank() << std::endl;
				break;
			case 4://�������
				std::cout << "���������A�Ľ�����";
				std::cin >> n;
                pA = new Matrix<Fraction>(n,n);
				std::cout << "��������� A" << n <<"��" << n << '\n';
				std::cin >> (*pA);
				std::cout << "���������Ϊ��\n" << pA->adjugate() << std::endl;
				break;
			case 5://�����
				std::cout << "���������Ľ�����";
				std::cin >> n;
                pA = new Matrix<Fraction>(n,n);
				std::cout << "��������� A" << n << "��" << n << '\n';
				std::cin >> (*pA);
				if(!pA->det())
					std::cout << "����A ������ʽ = 0��������\n";
				else
					std::cout << "����A����������£�\n" << pA->inv();
				break;
			case 6://���Է�����
				std::cout << "���������Է������������� Am��n �е�m n��";
				std::cin >> m >> n;
                pA = new Matrix<Fraction>(m,n);
				std::cout << "��������� A" << m << "��" << n << '\n';
				std::cin >> (*pA);
				pA->rref(n-1);
				std::cout << "�������ı�׼�������£�\n" << (*pA);
				break;
			case 7://��������ʽ
				std::cout << "���������A�Ľ�����";
				std::cin >> n;
                pA = new Matrix<Fraction>(n,n);
				std::cout <<"��������� A" << n << "��" << n << '\n';
				std::cin >> (*pA);
				std::cout << "��������ʽΪ��" << pA->poly() << std::endl;
				break;
			default:
				std::cout << "�����ڴ�ѡ��\n";
				break;
			}
			std::cout << "������������ͱ�ţ���ǰ���ͱ��Ϊ " << model <<" ������ 0 ��ʾ �ص�����ѡ��,8 ��ʾ�˳�:";
		}
        delete pA;
        delete pB;
        pA = pB = nullptr;
		if(model == 8)
			break;
	}
	std::cout << "��лʹ�ã�ע�Ᵽ����" << std::endl;
	system("pause");
	return 0;
}