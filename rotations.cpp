#include "stdafx.h"
#include <fstream>

#include "rotations.h"
#include "params.h"
#include "fragment.h"


namespace model
{
	void Taylor_rotations(Fragment *f)
	{
		f->om.setZero();
		//f->om = f->w - f->d_in.getAntiSymmetryPart();
		for (int i = 0; i < DIM; i++)//���� �������
		{
			for (int j = 0; j < DIM; j++)
			{
				for (int k = 0; k < f->SS_count; k++)
				{
					f->om.C[i][j] -= f->SS[k].dgm * (f->SS[k].n.C[i] * f->SS[k].b.C[j] - f->SS[k].b.C[i] * f->SS[k].n.C[j]);
				}
			}
		}
		f->om /= 2.0;
		f->om += f->w;
		
		Vector e;				//��� �������� �������
		e.set(f->om.C[1][2], f->om.C[2][0], f->om.C[0][1]);
		
		
		double dFi = e.getNorm();
		f->isRotate = dFi > EPS;
		if (f->isRotate)
		{
			f->rot_speed = dFi;
			dFi *= dt;			//�������� ���� ��������
			f->sum_angle += dFi;
			e.Normalize();
			f->Rotate(dFi, e);
	
		}
		else
		{
			f->rot_speed = 0;
		
		}
	}

	void Rotation_hardening(Fragment *f)
	{
		if (f->sum_angle > 100*EPS)
		{
			double HardRotK1 = 1e-8;//����������� ����� �����������
			double HardRotK2 = 1e3;//����������� ������ ����������
			double vol = pow(f->size, 3);//����� ��������
			double dmc = HardRotK1 / vol * exp( - HardRotK2 * f->sum_angle);//�������� ����������
			f->dmc = dmc;
			f->mc += dmc*dt;//���������� ������������ �������
		}
	}

	void Trusov_rotations(Fragment *f)
	{
	
		Vector *dm = new Vector[surround_count];
		Vector M;
		Vector dM;
		double S = f->size*f->size;		//������� ������� (������)

		for (int h = 0; h < surround_count; h++)//��������� �� ���� ������� ���������
		{
			if (f->contact[h] == 0) continue;//���� ��� �������� - ����������
			
			Tensor Lp = f->d_in - f->surrounds[h].d_in;//������ ������������ ����������
			Lp.Transp();
					
			Tensor buf = VectMult(f->normals[h], Lp);
			Vector m = ScalMult(buf, f->normals[h]);//������������� ������-������ 
			Vector b1 = ScalMult(f->om, m);//(������������� �����������)
			Vector b2 = ScalMult(m, f->om);
			dm[h] = m - b1 + b2;
			dm[h] *= ROT_L;
	/*	}
	
		for (int h = 0; h < surround_count; h++)	
		{*/
			dM += dm[h];
			dm[h] *= dt;
			f->moments[h] += dm[h];
			double c;		//���������� ������� �������� (� ����� �� ������ ������� ������� ����)
			if (h < 6) c = 1;
			else if (h < 14) c = 0.1;
			else c = 0.05;
			f->moments[h] *= S*c;
			M += f->moments[h];
		}
		double volume = pow(f->size, 3);		//����� ���������
		dM /= f->size;
		M /= volume;
		
		double pr = M.ScalMult(dM);
		
		M /= f->mc;
		dM /= f->mc;
		double norm = M.getNorm();
		double dMnorm = dM.getNorm();
		f->norm = norm;
		double dFi = 0;
		if (norm >= f->mc && pr >= 0)	//������������ � ������� ���������
		{
			dFi = ROT_A * dMnorm + ROT_H * norm ;
		}
		else
		{
			dFi = ROT_A * dMnorm;		//������ ������� ���������
		}
		
		f->isRotate = dFi > EPS * 1e5;
		if (f->isRotate)
		{
			Vector e;					//��� �������� �������
			e = M;						//������������ � �������� �������
			e.Normalize();		
			f->rot_speed = dFi;
			dFi *= dt;
			f->sum_angle += dFi;		//����������� ���� �������� �������������
			f->Rotate(dFi, e);			//�������� �������
			f->rot_energy = norm*dFi;	//������� ������� �������

			f->om.setZero();			//���� �������
			for (int i = 0; i < DIM; i++)//���� ������� (����� �� �������)
			{
				for (int j = 0; j < DIM; j++)
				{
					
					for (int k = 0; k < f->SS_count; k++)
					{
						f->om.C[i][j] -= f->SS[k].dgm * (f->SS[k].n.C[i] * f->SS[k].b.C[j] - f->SS[k].b.C[i] * f->SS[k].n.C[j]);
					}
					f->om.C[i][j] /= 2.0;
				
				}
			}
			f->om += f->w;

			for (int i = 0; i < DIM; i++)    //���� ������� (������ �����)
			{
				for (int j = 0; j < DIM; j++)
				{
					for (int k = 0; k < DIM; k++)
						f->om.C[i][j] -= LeviCivit(i, j, k) * e.C[k] * f->rot_speed;
				}
			}
			
		}
		else
		{
			f->rot_speed = 0;		//������� �� ���������
			f->rot_energy = 0;		//������� �������� ����� ����
		}

		delete[] dm;
	}

	void inline SavePoints(Tensor O, const char *file, const int i, const int j, const int k)
	{
		Vector e, e1;
		e.set(i, j, k);
		e1 = ScalMult(e, O);//������� � ���
		e1.Normalize();
		std::ofstream of;
		of.open(file, std::ios::out | std::ios_base::app | std::ios::binary);
		//������ � ���� x,y,z
		of.write((char *)&e1.C[0], sizeof(double));
		of.write((char *)&e1.C[1], sizeof(double));
		of.write((char *)&e1.C[2], sizeof(double));
		of.close();
	}

	void GetPoleFig(Fragment *f)
	{
		/*---��������� ����������� [001]---*/
		SavePoints(f->o, "Polus\\S001.dat", 0, 0, 1);
		SavePoints(f->o, "Polus\\S010.dat", 0, 1, 0);
		SavePoints(f->o, "Polus\\S100.dat", 1, 0, 0);

		/*---��������� ����������� [011]---*/
		SavePoints(f->o, "Polus\\S10-1.dat", 1, 0, -1);
		SavePoints(f->o, "Polus\\S01-1.dat", 0, 1, -1);
		SavePoints(f->o, "Polus\\S1-10.dat", 1, -1, 0);
		SavePoints(f->o, "Polus\\S011.dat", 0, 1, 1);
		SavePoints(f->o, "Polus\\S110.dat", 1, 1, 0);
		SavePoints(f->o, "Polus\\S101.dat", 1, 0, 1);

		/*---��������� ����������� [111]---*/
		SavePoints(f->o, "Polus\\S11-1.dat", 1, 1, -1);
		SavePoints(f->o, "Polus\\S1-11.dat", 1, -1, 1);
		SavePoints(f->o, "Polus\\S-111.dat", -1, 1, 1);
		SavePoints(f->o, "Polus\\S111.dat", 1, 1, 1);
	}

	void inline GetSSTPoint(Tensor O, float dFi, const char *FileName, const int i, const int j, const int k)
	{
		Vector a, v;
		a.set(i, j, k);
		a.Normalize();

		double CosFi = cos(dFi);
		double SinFiX = sin(dFi)*a.C[0];
		double SinFiY = sin(dFi)*a.C[1];
		double SinFiZ = sin(dFi)*a.C[2];
		double COS = (1.0 - CosFi);

		Tensor dO;			//������ �������� �� ����

		dO.C[0][0] = CosFi + COS * a.C[0] * a.C[0];
		dO.C[0][1] = COS * a.C[0] * a.C[1] - SinFiZ;
		dO.C[0][2] = COS * a.C[0] * a.C[2] + SinFiY;

		dO.C[1][0] = COS * a.C[0] * a.C[1] + SinFiZ;
		dO.C[1][1] = CosFi + COS * a.C[1] * a.C[1];
		dO.C[1][2] = COS * a.C[1] * a.C[2] - SinFiX;

		dO.C[2][0] = COS * a.C[0] * a.C[2] - SinFiY;
		dO.C[2][1] = COS * a.C[1] * a.C[2] + SinFiX;
		dO.C[2][2] = CosFi + COS * a.C[2] * a.C[2];
		
		Tensor R = dO*O;
		v = ScalMult(a, R);
		v.Normalize();
	
		std::ofstream Of;
		Of.open(FileName, std::ios::out | std::ios_base::app | std::ios::binary);
		//������ ��������� � ������� [xyz] ������ ��� ���� ���� 
		Of.write((char *)&v.C[0], sizeof(double));
		Of.write((char *)&v.C[1], sizeof(double));
		Of.write((char *)&v.C[2], sizeof(double));

		Of.close();
	}

	void GetSST(Fragment *f)
	{
		for (float fi = 0; fi < 2 * PI; fi += PI_2)
		{
			GetSSTPoint(f->o, fi, "Polus\\SST001.dat", 0, 0, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST001.dat", 0, 1, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST001.dat", 1, 0, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST001.dat", 0, 0, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST001.dat", 0, -1, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST001.dat", -1, 0, 0);
		}
		for (float fi = 0; fi < 2 * PI; fi += PI)
		{
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 0, 1, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 0, 1, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 1, 0, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 1, 1, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 0, -1, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", -1, 1, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", -1, 0, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", -1, 0, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 1, 0, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 1, -1, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", -1, -1, 0);
			GetSSTPoint(f->o, fi, "Polus\\SST011.dat", 0, -1, -1);
		}
		for (float fi = 0; fi < 2 * PI; fi += 2.0 * PI / 3.0)
		{
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", 1, 1, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", 1, 1, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", 1, -1, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", -1, 1, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", 1, -1, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", -1, 1, -1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", -1, -1, 1);
			GetSSTPoint(f->o, fi, "Polus\\SST111.dat", -1, -1, -1);
		}

	}
}