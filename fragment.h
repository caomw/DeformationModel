#ifndef __FRAGMENT_H 
#define __FRAGMENT_H

#include "MathCore.h"
/*
*�������� ������� "������", "������", "������� ����������"
*� ������� ������ � ����
*/
namespace model
{

	class Fragment
	{
	public:
		Tensor d;						//������ ���������� ��������
		Tensor w;						//������ �����
		Tensor d_in;					//������ ��������� ������������ ����������
		Tensor o;						//�������������� ������
		Tensor om;						//������ ����� �������
		Tensor sgm;						//������ ����������
		Tensor dsgm;					//������ ��������� ����������
		Tensor e;						//������ ����������
		Tensor4 p;					//������ ������� �������
		
		int SS_count;					//���-�� ������ ����������
		SlipSystem *SS;					//������� ����������

		int material;					//��� ���������
		double size;					//�������� ������ ���������
		double stress;					//������������� ����������
		double strain;					//������������� ����������
		
		double norm;//����� ��� ������
	
		Fragment *surrounds;			//������ �� ���������� ���������
		Vector *normals;				//������� ������� � ���������� ����������
		Vector *moments;				//������������� ������� �� ������
		int *contact;					
		/*
		contact - ������ � ����������� � ���, ��� ������������� ���������
		��������� ��������:
		|-1	| ������� ��� �� �����						|
		| 0	| ��� ��������								|
		| 1	| ������� �� ������ ������� ������� (100 %)	|
		| 2	| ������� �� ����� (10 %)					|
		| 3	| ������� �� ������� (5 %)					|
		*/

		bool isRotate;					//��������� �� ������� �� ������� ����
		double rot_speed;				//�������� �������� �������
		double sum_angle;				//����������� ���� �������� �������
		double mc;						//��������� ����������� ������
		double rot_energy;				//������� ������� ���������

		void setMaterialParams(int);	//������� ������������ ����������
		void Orientate(double, double,
			double, double);			//������� ��������� ����������
		void NDScalc();					//���������� ��� ���������
		void Rotate(double, Vector);	//������� �������

		Fragment();
		~Fragment();


	private:

	};


}

#endif