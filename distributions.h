#ifndef __DISTRIB_H 
#define __DISTRIB_H

namespace model
{
	/*
	��������� ������������� ��������� ��������
	��� ������������������� - ������ �������� - 
	���.��������, ������ - ���������
	*/
	double UniformDistrib(double, double);		//����������� �������������
	double NormalDistrib(double m, double d);	//���������� �������������
	double LogNormalDistrib(double m, double d);//������������� �������������
}

#endif __DISTRIB_H