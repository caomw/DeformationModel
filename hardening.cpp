#include "stdafx.h"
#include <cmath>

#include "params.h"
#include "materials.h"
#include "hardening.h"
namespace model
{
	void Base_hardening(Fragment *f)
	{
		double nsd = 0;	//Накопленный сдвиг
		for (int k = 0; k < f->SS_count; k++)
		{
			nsd += f->SS[k].gmm;		//Суммируем по всем СС
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
			case 0://Сталь 45 (ОЦК)
			{
				if (k < 24) napr = STEEL_TC1;
				else if (k < 48) napr = STEEL_TC2;
				else napr = STEEL_TC3;
				break;
			}
			case 1://Медь (ГЦК)
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
		Vector **dbb = new Vector*[f->SS_count];
		for (int i = 0; i < f->SS_count; i++)
		{
			dbb[i] = new Vector[f->SS_count];
		}
		Vector *s1 = new Vector[f->SS_count];
		Vector *s2 = new Vector[f->SS_count];
		double M = 0;		//мера разориентации
		double zgu = 0;		//Приращение критического напряжения

		for (int h = 0; h < surround_count; h++)		//цикл по нормалям			
		{
			if (f->contact[h] == 0) continue;//Если нет контакта - пропускаем

			for (int k = 0; k < f->SS_count; k++)
			{
				double min = 1;
				int from;
				for (int p = 0; p < f->surrounds[h].SS_count; p++)//по системам соседнего зерна
				{
					s1[p] = ScalMult(f->o, f->SS[p].b);
					s2[p] = ScalMult(f->surrounds[h].o, f->SS[p].b);
					dbb[k][p] = s1[k] - s2[p];
					M += f->normals[h].ScalMult(dbb[k][p]);
					M = fabs(M);
					if (M < min)
					{
						min = fabs(M);
						from = k;
					}
					M = 0;
				}
				if ((min >= 0) && (min <= 1))
				{
					zgu += HARD_BOUND_K * f->SS[from].dgm * f->SS[from].gmm * min / f->size;
					f->SS[from].tc += zgu*dt;
				}
			}
		}
		for (int i = 0; i < f->SS_count; i++)
		{
			delete[] dbb[i];
		}
		delete[] dbb;
		delete[] s1;
		delete[] s2;
	}
}