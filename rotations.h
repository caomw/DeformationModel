#ifndef __ROTATIONS_H 
#define __ROTATIONS_H

#include "fragment.h"

/*
*����������� ���������
*/

namespace model
{
	
	void Taylor_rotations(Fragment*);		//������ ���������� �������� �� �������
	void Trusov_rotations(Fragment*);		//������, ��������� � ��������������� �������
	void Rotation_hardening(Fragment*);		//������ ������������ ����������

	/*********************************************************
	********	    ��������� �������� ����� � ���	   *******
	*********************************************************/

	void SavePoints(Tensor,	char*,
		int, int, int);			//������ �������� �� � ����
	void GetPoleFig(Fragment*);				//���������� �������� ������
	
	void SaveSSTPoints(Tensor&, float,		//������ �������� ��� � ����
		char*, int, int, int);
	void GetSST(Fragment&);					//���������� ���
}

#endif __ROTATIONS_H