#ifndef __TENSION_H 
#define __TENSION_H

#include "MathCore.h"

/*
*������� ��� ���������� �������� ���������� � ���������� ��� ��������� ����������
*/

namespace model
{
	Tensor TensionStrainCalc(Tensor4 P, 
		Tensor D_in, double tens_comp);					//���������� ������� ���������� ��� ���������
	Tensor TensionStressCalc(Tensor4 &P,
		Tensor &D_in, Tensor &D);						//���������� ������� ���������� ��� ����������
	Tensor UnloadingStrainCalc(Tensor4 &P, 
		Tensor &D_in, Tensor &Sgm, double lam);			//��������� ������ ���������� ��� ���������
}
#endif __TENSION_H