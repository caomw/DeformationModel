#include "stdafx.h"
#include <cmath>

#include "params.h"
#include "materials.h"
#include "hardening.h"
namespace model
{
	void Base_hardening(Fragment *f)
	{
		double nsd = 0;	//����������� �����
		for (int k = 0; k < f->SS_count; k++)
		{
			nsd += f->SS[k].gmm;		//��������� �� ���� ��
		}
		for (int k = 0; k < f->SS_count; k++)
		{
			double osn = 0;
			for (int j = 0; j < f->SS_count; j++)
			{
				double a = (j == k ? HARD_BASE_A : 1.25 * HARD_BASE_A);
				if ((f->SS[j].dgm > EPS) && (fabs(nsd) > EPS))
				{
					osn += a*pow(f->SS[j].gmm, HARD_BASE_PSI)*pow((f->SS[j].dgm / dgm0), HARD_BASE_DELTA) / nsd;
				}
			}
			double napr = 0;
			switch (f->material)
			{
			case 0://����� 45 (���)
			{
				if (k < 24) napr = STEEL_TC1;
				else if (k < 48) napr = STEEL_TC2;
				else napr = STEEL_TC3;
				break;
			}
			case 1://���� (���)
			{
				napr = CUPR_TC;
				break;
			}
			
			}
			f->SS[k].tc += napr*osn*dt;
		}
	}
	
	void Boundary_hardening(Fragment *f)
	{
		for (int k = 0; k < f->SS_count; k++)	//���� �� �� �������� ���������
		{
			Vector b1 = ScalMult(f->o, f->SS[k].b);//�������� ������ b ������� �� ������� ����� � ���
			double zgu = 0;
			double zguk;
			
			for (int h = 0; h < surround_count; h++)	//���� �� ��������			
			{
				if (f->contact[h] == 0) continue;//���� ��� �������� - ����������
				if (f->SS[k].b.ScalMult(f->normals[h]) < 0) continue; //���������� �� ������� - ����������
				zguk = HARD_BOUND_K * f->SS[k].dgm * f->SS[k].gmm / f->size;
				double min = 1.0;//�������
				for (int p = 0; p < f->surrounds[h].SS_count; p++)	//���� �� �������� ��������� �����
				{
					Vector b2 = ScalMult(f->surrounds[h].o, f->surrounds[h].SS[p].b);//�������� ������ b p-�� �� ��������� ����� � ���
					Vector diff = b1 - b2;
					diff.Normalize();
					double M = fabs(diff.ScalMult(f->normals[h]));

					if (M < min)
					{
						min = M;
					}
				}
				zgu += zguk*min;
			}

			/*if (!isnan(zgu))*/ f->SS[k].tc += zgu;
		}
	}

}