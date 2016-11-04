#ifndef __MATERIALS_H 
#define __MATERIALS_H
/*
*� ���� ����� ���������� ��� ������������ ���������-���������
*/
namespace model
{
	const int SS_COUNT_BCC = 96;		//���-�� �� � ���
	const int SS_COUNT_FCC = 24;		//���-�� �� � ���
	const int SS_COUNT_HCP = 36;		//���-�� �� � ���

	/********************************************
	*************    ����� 45     ***************
	********************************************/

	/******       ������� ���������        *****/
	const double STEEL_P1 = 2.2e11;
	const double STEEL_P2 = 1.66e11;
	const double STEEL_P3 = 8.7e10;

	/***** ��������� ����������� ���������� *****/
	const double STEEL_TC1 = 0.1e9;
	const double STEEL_TC2 = 0.47e9;
	const double STEEL_TC3 = 2.47e9;



	/********************************************
	***********    ������ ����     **************
	********************************************/

	/******       ������� ���������        *****/
	const double CUPR_P1 = 1.684e11;
	const double CUPR_P2 = 1.214e11;
	const double CUPR_P3 = 7.54e10;

	/***** ��������� ����������� ���������� *****/
	const double CUPR_TC = 1.75e7;


	/********************************************
	***********    �����-�����     **************
	********************************************/
	const double TITAN_TC1 = 15e7;
	const double TITAN_TC2 = 3e7;
	const double TITAN_TC3 = 12e7;

	const double TITAN_P1 = 18.07e10;
	const double TITAN_P2 = 16.24e10;
	const double TITAN_P3 = 6.9e10;
	const double TITAN_P4 = 9.2e10;
	const double TITAN_P5 = 4.67e10;

}

#endif __MATERIALS_H