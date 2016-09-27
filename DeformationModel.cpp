#include "stdafx.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <Windows.h>
#include <omp.h>

#include "MathCore.h"
#include "params.h"
#include "fragment.h"
#include "rotations.h"
#include "hardening.h"
#include "functions.h"
#include "tension.h"
#include "distributions.h"

using namespace model;

int _tmain(int argc, _TCHAR* argv[])
{
	if (argc == 1) return 1;	//��������� ���������, ���� ������� ��� ����������
	
	char* param_file = new char[256];
	wcstombs(param_file, argv[1], 256);//�������� ��� ����� � �����������
	ReadParams(param_file);		//������� ��������� �� �����
	
		
	if (!isDirectoryExists(L"Plot"))
	{
		CreateDirectory(L"Plot", NULL);
	}
	if (!isDirectoryExists(L"Polus"))
	{
		CreateDirectory(L"Polus", NULL);
	}
	if (!isDirectoryExists(L"DBG"))
	{
		CreateDirectory(L"DBG", NULL);
	}
	const int total_fragm_count = (int)pow(fragm_count, 3);	//����� ���-�� ����������

	/*****************************************************************
	********     ��������� �����/������ ���������� ������     ********
	*****************************************************************/
	std::cout << " Build on " << __DATE__ << " " << __TIME__ << std::endl;
	std::cout << " Parameters file: " << param_file << std::endl;
	delete param_file;//������ �� �����
	std::cout << " __________________________________________" << std::endl;
	std::cout << " Fragments count: " << total_fragm_count << std::endl;
	std::cout << " Max. strain: " << strain_max << std::endl;
	std::cout << " Integration step: " << dt << std::endl;
	if (ROTATIONS_HARDENING)
	{
		std::cout << " Using rotation hardening model" << std::endl;
	}
	if (ROTATIONS_TAYLOR)
	{
		std::cout << " Using Taylor's rotation model" << std::endl;
	}
	if (ROTATIONS_TRUSOV)
	{
		std::cout << " Using Trusov's rotation model" << std::endl;
		std::cout << "       A: " << ROT_A << std::endl;
		std::cout << "       H: " << ROT_H << std::endl;
		std::cout << "       L: " << ROT_L << std::endl;
		std::cout << "       MC: " << ROT_MC << std::endl;
	}
	if (HARDENING_BASE)
	{
		std::cout << " Using basic hardening" << std::endl;
		std::cout << "       A: " << HARD_BASE_A << std::endl;
		std::cout << "       Delta: " << HARD_BASE_DELTA << std::endl;
		std::cout << "       Psi: " << HARD_BASE_PSI << std::endl;
	}
	if (HARDENING_BOUND)
	{
		std::cout << " Using boundary hardening" << std::endl;
		std::cout << "       K: " << HARD_BOUND_K << std::endl;
	}
	if (debug_period > 0)
	{
		std::cout << " ++++++++++++++DEBUG MODE++++++++++++++" << std::endl;
		std::cout << "       Period: " << debug_period << std::endl;
	}
	if (fix_orient == 1)
	{
		std::cout << " Saving current orientations and normals" << std::endl;
	}
	if (fix_orient == 2)
	{
		std::cout << " Reading saved orientations and normals" << std::endl;
	}

	Fragment ***PC = new Fragment**[fragm_count];		//������������	
	for (int i = 0; i < fragm_count; i++)
	{
		PC[i] = new Fragment*[fragm_count];
		for (int j = 0; j < fragm_count; j++)
		{
			PC[i][j] = new Fragment[fragm_count];
		}
	}
	unsigned long t1, t2;			//������� �������

	Tensor macro_D;					//������ ���������� �������� �����������
	Tensor macro_W;					//������ ����� �����������
	Tensor macro_D_in;				//������ ��������� ����� ���������� �����������
	Tensor macro_E;					//������ ���������������
	Tensor macro_dSgm;				//������ ��������� ���������������
	Tensor macro_Sgm;				//������ ���������������
	Tensor4 macro_P;				//���������� ������ ������� ��������
	
	std::cout << " Initializing all fragments... ";
	t1 = clock();

	macro_D = gradV.getSymmetryPart();
	macro_W = gradV.getAntiSymmetryPart();

	std::ofstream TestStream[4];
	TestStream[0].open("Test0.txt", std::ios_base::out | std::ios_base::trunc);
	TestStream[1].open("Test1.txt", std::ios_base::out | std::ios_base::trunc);
	TestStream[2].open("Test2.txt", std::ios_base::out | std::ios_base::trunc);
	TestStream[3].open("Test3.txt", std::ios_base::out | std::ios_base::trunc);
	
	std::srand(time(NULL));

	switch (SurroundsGrade)			//������� ����� �������� ���������
	{
	case 0:
	{
		surround_count = 6;			//������� �������
		break;
	}
	case 1:
	{
		surround_count = 18;		//���������� �������
		break;
	}
	case 2:
	{
		surround_count = 26;		//����� ������� �������
		break;
	}
	}
	
	for (int q1 = 0; q1 < fragm_count; q1++)
	{
		for (int q2 = 0; q2 < fragm_count; q2++)
		{
			for (int q3 = 0; q3 < fragm_count; q3++)
			{
				//������� ��������� 
				int another_material;//��������� ����
				another_material = (material == 1) ? 0 : 1;

				int a = (int)(((double)rand() / RAND_MAX) * 100);//�� �� ���� �����
				if (a <= material_purity)
				{
					PC[q1][q2][q3].setMaterialParams(material);
				}
				else
				{
					PC[q1][q2][q3].setMaterialParams(another_material);
				}
				PC[q1][q2][q3].mc = ROT_MC;//������� ��������� ����������� ��������

				if (RAND_ORIENT)//��������� ����������� �����
				{

					double a = ((double)rand() / RAND_MAX) * (PI);
					double g = ((double)rand() / RAND_MAX) * (PI);
					double y1 = ((double)rand() / RAND_MAX);
					double y2 = ((double)rand() / RAND_MAX);
					PC[q1][q2][q3].Orientate(a, g, y1, y2);
				}
				else//���=���
				{
					PC[q1][q2][q3].o.setUnit();
				}
				//������� �������� ����������
				switch (fragm_size_law)
				{
				case 0://�����������
				{
					PC[q1][q2][q3].size = UniformDistrib(fragm_size_m,fragm_size_dsp);
					break;
				}
				case 1://����������
				{
					PC[q1][q2][q3].size = NormalDistrib(fragm_size_m, fragm_size_dsp);
					break;
				}
				case 2://�������������
				{
					PC[q1][q2][q3].size = LogNormalDistrib(fragm_size_m, fragm_size_dsp);
					break;
				}
				}
				
				PC[q1][q2][q3].surrounds = new Fragment[surround_count];
				PC[q1][q2][q3].normals = new Vector[surround_count];
				PC[q1][q2][q3].moments = new Vector[surround_count];
				PC[q1][q2][q3].contact = new int[surround_count];

				for (int h = 0; h < surround_count; h++)
				{
					PC[q1][q2][q3].contact[h] = -1;		//���������� ������� �� �����
				}
			}
		}
	}
	
	for (int q1 = 0; q1 < fragm_count; q1++)
	{
		for (int q2 = 0; q2 < fragm_count; q2++)
		{
			for (int q3 = 0; q3 < fragm_count; q3++)
			{

				for (int h = 0; h < surround_count; h++)
				{
					//���� ������� ��� ��� ����� - ����������
					if (PC[q1][q2][q3].contact[h] != -1) continue;
					//����������, �������� �� ���������
					//������ 6, �.�. ������� �����, �������� ������
					double a = h < 6 ? 1 : ((double)rand() / RAND_MAX);//�� �� ���� �����
					if (a < 0.5)
					{
						//�������� ��� - ���� ����������
						PC[q1][q2][q3].contact[h] = 0;
						continue;
					}

					int qq1 = q1, qq2 = q2, qq3 = q3, y;
					//qq1, qq2, qq3 - ���������� ����� ������
					//y - ����� ������� � �������� ����� � ����������� ������� �����
					double fi = ((double)rand() / RAND_MAX) * (M_PI / 12);
					switch (h)
					{
					case 0://�����
					{
						PC[q1][q2][q3].normals[h].set(-sin(fi), sin(fi) / cos(fi), 1 / cos(fi));
						qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
						y = 5;
						break;
					}
					case 1://�� ���
					{
						PC[q1][q2][q3].normals[h].set(-1 / cos(fi), sin(fi), sin(fi) / cos(fi));
						qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
						y = 3;
						break;
					}
					case 2://������
					{
						PC[q1][q2][q3].normals[h].set(sin(fi) / cos(fi), 1 / cos(fi), -sin(fi));
						qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
						y = 4;
						break;
					}
					case 3://�� ���
					{
						PC[q1][q2][q3].normals[h].set(1 / cos(fi), -sin(fi), sin(fi) / cos(fi));
						qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
						y = 1;
						break;
					}
					case 4://�����
					{
						PC[q1][q2][q3].normals[h].set(sin(fi), -1 / cos(fi), -sin(fi) / cos(fi));
						qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
						y = 2;
						break;
					}
					case 5://����
					{
						PC[q1][q2][q3].normals[h].set(sin(fi) / cos(fi), sin(fi), -1 / cos(fi));
						qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
						y = 0;
						break;
					}
					//����� ���� ��� �������������� ������
					/**************           и��� ����          ***************************/
					case 6://���� �� ���
					{
						PC[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * cos(fi)*cos(PI_2), -cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
						qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
						qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
						y = 9;
						break;
					}
					case 7://���� �� ���
					{
						PC[q1][q2][q3].normals[h].set(cos(fi + PI_2) * cos(fi)*cos(PI_2), -cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
						qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
						qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
						y = 8;
						break;
					}
					case 8://����� �� ���
					{
						PC[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
						qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
						qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
						y = 7;
						break;
					}
					case 9://����� �� ���
					{
						PC[q1][q2][q3].normals[h].set(cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi + PI_2) * cos(fi)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2));
						qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
						qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
						y = 6;
						break;
					}
					case 10://���� ����
					{
						PC[q1][q2][q3].normals[h].set(cos(fi + PI_2) * sin(fi), -cos(fi + PI_2) * cos(fi), sin(fi + PI_2)*cos(fi));
						qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
						qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
						y = 15;
						break;
					}
					case 11://���� �����
					{
						PC[q1][q2][q3].normals[h].set(cos(fi + PI_2) * sin(fi), cos(fi + PI_2) * cos(fi), sin(fi + PI_2)*cos(fi));
						qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
						qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
						y = 14;
						break;
					}
					case 12://���� �� ���
					{
						PC[q1][q2][q3].normals[h].set(cos(fi + PI_2) * cos(fi), cos(fi + PI_2) * sin(fi), sin(fi + PI_2)*cos(fi));
						qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
						qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
						y = 17;
						break;
					}
					case 13://���� �� ���
					{
						PC[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * cos(fi), cos(fi + PI_2) * sin(fi), sin(fi + PI_2)*cos(fi));
						qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
						qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
						y = 16;
						break;
					}
					case 14://��� ����
					{
						PC[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * sin(fi), -cos(fi + PI_2) * cos(fi), -sin(fi + PI_2)*cos(fi));
						qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
						qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
						y = 11;
						break;
					}
					case 15://��� �����
					{
						PC[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * sin(fi), cos(fi + PI_2) * cos(fi), -sin(fi + PI_2)*cos(fi));
						qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
						qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
						y = 10;
						break;
					}
					case 16://��� �� ���
					{
						PC[q1][q2][q3].normals[h].set(cos(fi + PI_2) * cos(fi), -cos(fi + PI_2) * sin(fi), -sin(fi + PI_2)*cos(fi));
						qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
						qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
						y = 13;
						break;
					}
					case 17://��� �� ���
					{
						PC[q1][q2][q3].normals[h].set(-cos(fi + PI_2) * cos(fi), -cos(fi + PI_2) * sin(fi), -sin(fi + PI_2)*cos(fi));
						qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
						qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
						y = 12;
						break;
					}
					/**************      �������     *****************/
					case 18://���� ���� �� ���
					{
						PC[q1][q2][q3].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
						qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
						qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
						qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
						y = 25;
						break;
					}
					case 19://���� ���� �� ���
					{
						PC[q1][q2][q3].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
						qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
						qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
						qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
						y = 24;
						break;
					}
					case 20://���� ����� �� ���
					{
						PC[q1][q2][q3].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
						qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
						qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
						qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
						y = 23;
						break;
					}
					case 21://���� ����� �� ���
					{
						PC[q1][q2][q3].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
						qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
						qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
						qq3 = q3 == fragm_count - 1 ? 0 : q3 + 1;
						y = 22;
						break;
					}
					case 22://��� ���� �� ���
					{
						PC[q1][q2][q3].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
						qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
						qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
						qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
						y = 21;
						break;
					}
					case 23://��� ���� �� ���
					{
						PC[q1][q2][q3].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
						qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
						qq2 = q2 == 0 ? fragm_count - 1 : q2 - 1;
						qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
						y = 20;
						break;
					}
					case 24://��� ����� �� ���
					{
						PC[q1][q2][q3].normals[h].set(-cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));

						qq1 = q1 == 0 ? fragm_count - 1 : q1 - 1;
						qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
						qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
						y = 19;
						break;
					}
					case 25://��� ����� �� ���
					{
						PC[q1][q2][q3].normals[h].set(cos(fi) * cos(fi)*cos(PI_2)*cos(PI_2), cos(fi)*cos(fi) * cos(PI_2)*cos(PI_2), -cos(fi)*cos(fi)*cos(PI_2)*cos(PI_2));
						qq1 = q1 == fragm_count - 1 ? 0 : q1 + 1;
						qq2 = q2 == fragm_count - 1 ? 0 : q2 + 1;
						qq3 = q3 == 0 ? fragm_count - 1 : q3 - 1;
						y = 18;
						break;
					}
					}

					PC[q1][q2][q3].surrounds[h] = PC[qq1][qq2][qq3];//����������, �����!
					PC[qq1][qq2][qq3].surrounds[y] = PC[q1][q2][q3];//������� �������������!
					PC[q1][q2][q3].normals[h].Normalize();

					for (int i = 0; i < DIM; i++)
					{
						PC[qq1][qq2][qq3].normals[y].C[i] = -PC[q1][q2][q3].normals[h].C[i];//�������� ��������
					}

					if (h < 6) PC[q1][q2][q3].contact[h] = 1;//100 % ���������������
					else if (h < 14) PC[q1][q2][q3].contact[h] = 3;//5 % ���������������
					else PC[q1][q2][q3].contact[h] = 2;//10 % ���������������
				}

			}
		}
	}

	if (fix_orient == 2)	//���������� ���������� ����������
	{
		std::ifstream StreamO("DBG\\o.txt", std::ios_base::in);
		std::ifstream StreamNorm("DBG\\Norm.txt", std::ios_base::in);
		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{
					for (int i = 0; i < DIM; i++)	//��������� �������� �������������� ��������
					{
						for (int j = 0; j < DIM; j++)
						{
							StreamO >> PC[q1][q2][q3].o.C[i][j];
						}
					}
					for (int h = 0; h < surround_count; h++)//��������� �������� ��������
					{
						for (int i = 0; i < DIM; i++)
						{
							StreamNorm >> PC[q1][q2][q3].normals[h].C[i];
						}
					}
				}
			}
		}
		StreamO.close();
		StreamNorm.close();
	}

	if (fix_orient == 1)//����������� ��������� ����������
	{
		std::ofstream StreamO("DBG\\o.txt", std::ios_base::out | std::ios_base::trunc);
		std::ofstream StreamNorm("DBG\\Norm.txt", std::ios_base::out | std::ios_base::trunc);
		for (int q1 = 0; q1 < fragm_count; q1++)
		{
			for (int q2 = 0; q2 < fragm_count; q2++)
			{
				for (int q3 = 0; q3 < fragm_count; q3++)
				{
					WriteDebugInfo(StreamO, PC[q1][q2][q3].o.C);//���������� �������� �������� ����������
					for (int h = 0; h < surround_count; h++)//���������� �������� ��������
					{
						for (int i = 0; i < DIM; i++)
						{
							StreamNorm << PC[q1][q2][q3].normals[h].C[i] << " ";
						}
						StreamNorm << std::endl;
					}
				}
			}
		}
		StreamO.close();
		StreamNorm.close();
	}
	t2 = clock();
	std::cout << (t2 - t1) / 1000.0 << " sec" << std::endl;

	int CURR_STEP = 0;				//������� ��� ��������������
	int PLOT_STEP;					//��� ���������� ��������
	int POLUS_STEP;					//��� ���������� ��
	int PROC_STEP = 0;				//��� ����������� ���������
	int DEBUG_STEP = 0;				//��� ������ ���������� ������
	int proc_period = 400;			//������ ���������� �������� ����������
	double macro_stress = 0;		//������������� ���������������
	double macro_strain = 0;		//������������� ���������������

	/*������ � ������� ����������*/
	double tension_component = 0;	//������������ ���������� � ��������
	double final_stress = 1e3;		//��������, �� �������� ����������
	double lam = 2.5;				//����������� � ���������
	double addition_strain = 1e-4;	//���������� ��������� ��� ����������� �������
	if (REAL_UNIAX)
	{
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				if (macro_D.C[i][j] != 0)	//����������� ������������� ����������
				{
					tension_component = macro_D.C[i][j];
					break;
				}
			}
		}
	}
	

	/*******************************************************
	**********       ������ � ������� ������       *********
	*******************************************************/

	TruncPoleFiles();				//������� ���� ������ �������� �����
	TruncSSTFiles();
	std::ofstream StrainStream("Plot\\X.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
	std::ofstream StressStream("Plot\\Y.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
	std::ofstream StrainStreamAll("Plot\\Xall.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
	std::ofstream StressStreamAll("Plot\\Yall.dat", std::ios::out | std::ios_base::trunc | std::ios::binary);
	std::ofstream ActiveSysStream("Plot\\ActiveSS.dat", std::ios_base::out | std::ios_base::trunc | std::ios::binary);

	std::ofstream dbgstream[15];
	if (debug_period > 0)				//�������� ������ ��� ���������� ������
	{
		dbgstream[0].open("DBG\\o.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[1].open("DBG\\e.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[2].open("DBG\\d.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[3].open("DBG\\sgm.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[4].open("DBG\\om.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[5].open("DBG\\dsgm.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[6].open("DBG\\din.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[7].open("DBG\\w.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[8].open("DBG\\dgamma.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[9].open("DBG\\t.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[10].open("DBG\\Macro_D.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[11].open("DBG\\Macro_Din.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[12].open("DBG\\Macro_Sgm.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[13].open("DBG\\Macro_dSgm.txt", std::ios_base::out | std::ios_base::trunc);
		dbgstream[14].open("DBG\\Macro_E.txt", std::ios_base::out | std::ios_base::trunc);
	}
	std::cout << " Saving pole figures... ";
	t1 = clock();
	//���������� ��������� �������� ����� � ���
	for (int q1 = 0; q1 < fragm_count; q1++)
	{
		for (int q2 = 0; q2 < fragm_count; q2++)
		{
			for (int q3 = 0; q3 < fragm_count; q3++)
			{
				GetPoleFig(&PC[q1][q2][q3]);
				//GetSST(&PC[q1][q2][q3]);
			}
		}
	}
	t2 = clock();
	std::cout << (t2 - t1) / 1000.0 << " sec" << std::endl;
	/*******************************************************
	**********      ���� �� ������ ����������      *********
	*******************************************************/
	t1 = clock();
	
	for (int cycle = 0; cycle < cycle_count; cycle++)
	{
		PLOT_STEP = 0;
		POLUS_STEP = 0;
		if (cycle_count > 1) std::cout << std::endl << " Cycle # " << cycle+1;
		std::cout << std::endl << "        0.00%";

		/***************************************
		********    ���� �� �������    *********
		***************************************/
		while (macro_strain < strain_max)
		{
			//����������������� ������� ������ �����
			if (REAL_UNIAX)	//��������� ����������
			{
				//����������
				macro_P.setZero();
				macro_D_in.setZero();
				for (int q1 = 0; q1 < fragm_count; q1++)
				{
					for (int q2 = 0; q2 < fragm_count; q2++)
					{
						for (int q3 = 0; q3 < fragm_count; q3++)
						{
							macro_D_in += PC[q1][q2][q3].d_in;
							macro_P += PC[q1][q2][q3].p.ToLSK(PC[q1][q2][q3].o);
						}
					}
				}

				macro_D_in /= (total_fragm_count);
				macro_P /= (total_fragm_count);

				//�������������
				macro_P.Symmetrize();

				macro_D = TensionStrainCalc(macro_P, macro_D_in, macro_D.C[0][0]);

				Tensor b = macro_D;
				b *= dt;
				macro_E += b;

				macro_strain = macro_E.doubleScalMult(macro_E);
				macro_strain = SQRT2_3*sqrt(macro_strain);

				macro_dSgm = TensionStressCalc(macro_P, macro_D_in, macro_D);
				macro_dSgm *= dt;				//���������� ���������� �� ����
				macro_Sgm += macro_dSgm;
				macro_stress = macro_Sgm.doubleScalMult(macro_Sgm);
				macro_stress = SQRT3_2*sqrt(macro_stress);
			}
			omp_set_num_threads(thread_count);
			#pragma omp parallel for
			//�����, ������� ����� ����������
			//����� ���������� ������������� ������ ������ ������� ���������
			//�� ���������� ������� �����������
			for (int q1 = 0; q1 < fragm_count; q1++)
			{
				for (int q2 = 0; q2 < fragm_count; q2++)
				{
					for (int q3 = 0; q3 < fragm_count; q3++)
					{

						/**************************************************
						************       ��������� � ���       **********
						**************************************************/

						Tensor O = PC[q1][q2][q3].o;
						Tensor OT = O;
						OT.Transp();
						PC[q1][q2][q3].d = O*macro_D*OT;//�������� ������
						PC[q1][q2][q3].w = O*macro_W*OT;//�����������

						PC[q1][q2][q3].sgm = O*PC[q1][q2][q3].sgm*OT;
						PC[q1][q2][q3].d_in = O*PC[q1][q2][q3].d_in*OT;


						/***************************************************
						***********       ������������� ���      ***********
						***************************************************/

						PC[q1][q2][q3].NDScalc();

						if (HARDENING_BASE)			//������� ����������
						{
							Base_hardening(&PC[q1][q2][q3]);
						}
						if (HARDENING_BOUND)
						{
							Boundary_hardening(&PC[q1][q2][q3]);
						}
						if (ROTATIONS_TAYLOR)		//������� �� �������
						{
							Taylor_rotations(&PC[q1][q2][q3]);
						}
						if (ROTATIONS_TRUSOV)		//������� �� �������
						{
							Trusov_rotations(&PC[q1][q2][q3]);
						}
						if (ROTATIONS_TRUSOV && ROTATIONS_HARDENING)	//����������� ����������
						{
							Rotation_hardening(&PC[q1][q2][q3]);
						}

						/**************************************************
						************       ��������� � ���       **********
						**************************************************/

						PC[q1][q2][q3].sgm = OT*PC[q1][q2][q3].sgm*O;
						PC[q1][q2][q3].d_in = OT*PC[q1][q2][q3].d_in*O;

					}
				}
			}

			if (!REAL_UNIAX)	//���������� ��������������
			{
				macro_stress = 0;
				macro_strain = 0;
				for (int q1 = 0; q1 < fragm_count; q1++)
				{
					for (int q2 = 0; q2 < fragm_count; q2++)
					{
						for (int q3 = 0; q3 < fragm_count; q3++)
						{
							macro_strain += PC[q1][q2][q3].strain;
							macro_stress += PC[q1][q2][q3].stress;
						}
					}
				}
				macro_strain /= total_fragm_count;
				macro_stress /= total_fragm_count;
			}

			if (!REAL_UNIAX)		//���� ���� ����� ������������� ��� ������ � ��������!
			{
				macro_Sgm.setZero();
				for (int q1 = 0; q1 < fragm_count; q1++)
				{
					for (int q2 = 0; q2 < fragm_count; q2++)
					{
						for (int q3 = 0; q3 < fragm_count; q3++)
						{
							macro_Sgm += PC[q1][q2][q3].sgm;
						}
					}
				}
				//macro_Sgm = macro_dSgm;
				//macro_Sgm *= dt;
			}


			/************************************************************
			***********	        �������� ���������� 	      ***********
			************************************************************/
						
			double progress = (macro_strain / strain_max) * 100;
			
			if (!(cycle_count == 1 || cycle == 0))	//������������� ����������
			{
				progress /= 2;
				if (macro_E.C[0][0] > 0)			//���� ��� ��������� ���������
				{									//���������� �������� ����������,
					if (macro_Sgm.C[0][0] > 0)		//����� �� ������� ����� �����������
					{
						progress += 50.0;
					}
					else
					{
						progress = 50.0 - progress;
					}
				}
				else
				{
					if (macro_Sgm.C[0][0] > 0)
					{
						progress = 50.0 - progress;
					}
					else
					{
						progress += 50.0;
					}
				}
			}

			if (PROC_STEP == proc_period)
			{
				PROC_STEP = 0;
				std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
				printf("        %2.2f", progress);
				std::cout << "%";
			}

			/************************************************************
			***********	    ������ ������ ��� �������� ���    ***********
			************************************************************/
			
			if (progress - PLOT_STEP > plot_period && plot_period > 0)
			{
				if (!REAL_UNIAX)
				{
					StrainStream.write((char *)&macro_strain, sizeof macro_strain);
					StressStream.write((char *)&macro_stress, sizeof macro_stress);
				}
				else
				{
					StrainStream.write((char *)&macro_E.C[0][0], sizeof macro_E.C[0][0]);
					StressStream.write((char *)&macro_Sgm.C[0][0], sizeof macro_Sgm.C[0][0]);
				}
				double ActiveSysCount = 0;			//������� ���-�� �������� ������ ���������� �� ����
				double RotEnergy = 0;				//������� ������� �� ����
				double RotSpeed = 0;				//������� �������� �������� �� ����
				int RotCount = 0;					//���-�� ����������� ����������
				double norma = 0;
				for (int q1 = 0; q1 < fragm_count; q1++)
				{
					for (int q2 = 0; q2 < fragm_count; q2++)
					{
						for (int q3 = 0; q3 < fragm_count; q3++)
						{
							if (!REAL_UNIAX)
							{
								StrainStreamAll.write((char *)&PC[q1][q2][q3].strain, sizeof PC[q1][q2][q3].strain);
								StressStreamAll.write((char *)&PC[q1][q2][q3].stress, sizeof PC[q1][q2][q3].stress);
							}
							else
							{
								StrainStreamAll.write((char *)&PC[q1][q2][q3].e.C[0][0], sizeof PC[q1][q2][q3].e.C[0][0]);
								StressStreamAll.write((char *)&PC[q1][q2][q3].sgm.C[0][0], sizeof PC[q1][q2][q3].sgm.C[0][0]);
							}
							
							for (int i = 0; i < PC[q1][q2][q3].SS_count; i++)
							{
								if (PC[q1][q2][q3].SS[i].dgm > EPS) ActiveSysCount++;//������� �������� ��
							}
							
							if (PC[q1][q2][q3].isRotate) RotCount++;		//������� ����������� �������
							RotEnergy += PC[q1][q2][q3].rot_energy;				//������������ ������� ��������
							RotSpeed += PC[q1][q2][q3].rot_speed;				//������������ ��������� ��������
							norma += PC[q1][q2][q3].norm;
						}
					}
				}
				norma /= total_fragm_count;
				ActiveSysCount /= total_fragm_count;
				ActiveSysStream.write((char *)&ActiveSysCount, sizeof ActiveSysCount);//������ ���-�� �������� ��
				if (RotCount != 0)
				{
					RotSpeed /= RotCount;
				}
				else RotSpeed = 0;
				
				/*******************************************************
				********* 		   ������ � ��������          **********
				*******************************************************/

				//������ ������� �������������� - ����� ������������ ������� �� ������ ����
				//������������ ������� - ������ ���������� � ����������� ����������
				//������� ������� - ������*���������� ����

				Tensor dE = macro_D;
				dE *= dt;			//���������� ���������� �� ����

				double StepEnergy = macro_Sgm.doubleScalMult(dE);	//������ ������� �� ����

				TestStream[0] << RotCount << std::endl;	
				TestStream[1] << /*RotSpeed*/norma << std::endl;
				TestStream[2] << RotEnergy << std::endl;
				TestStream[3] << StepEnergy << std::endl;

				PLOT_STEP =  progress;
			}
			
			/************************************************************
			***********	      ���������� �������� �����	      ***********
			************************************************************/
			if (progress - POLUS_STEP > polus_period && polus_period > 0)
			{
				for (int q1 = 0; q1 < fragm_count; q1++)
				{
					for (int q2 = 0; q2 < fragm_count; q2++)
					{
						for (int q3 = 0; q3 < fragm_count; q3++)
						{
							GetPoleFig(&PC[q1][q2][q3]);
						//	GetSST(&PC[q1][q2][q3]); //�����������������, ���� ����� ��� (BETA)
						}
					}
				}
				POLUS_STEP = progress;
			}

			/************************************************************
			***********	       ������ ��������� ������	      ***********
			************************************************************/
			if (DEBUG_STEP == debug_period)
			{
				DEBUG_STEP = 0;

				for (int i = 0; i < 15; i++)
				{
					dbgstream[i] << "********      STEP " << CURR_STEP << "      ********" << std::endl << std::endl;
				}

				for (int q1 = 0; q1 < fragm_count; q1++)
				{
					for (int q2 = 0; q2 < fragm_count; q2++)
					{
						for (int q3 = 0; q3 < fragm_count; q3++)
						{
							WriteDebugInfo(dbgstream[0], PC[q1][q2][q3].o.C);
							WriteDebugInfo(dbgstream[1], PC[q1][q2][q3].e.C);
							WriteDebugInfo(dbgstream[2], PC[q1][q2][q3].d.C);
							WriteDebugInfo(dbgstream[3], PC[q1][q2][q3].sgm.C);
							WriteDebugInfo(dbgstream[4], PC[q1][q2][q3].om.C);
							WriteDebugInfo(dbgstream[5], PC[q1][q2][q3].dsgm.C);
							WriteDebugInfo(dbgstream[6], PC[q1][q2][q3].d_in.C);
							WriteDebugInfo(dbgstream[7], PC[q1][q2][q3].w.C);
							for (int f = 0; f < PC[q1][q2][q3].SS_count; f++)
							{
								dbgstream[8] << PC[q1][q2][q3].SS[f].dgm << " ";
							}
							dbgstream[8] << std::endl << std::endl;
							for (int f = 0; f < PC[q1][q2][q3].SS_count; f++)
							{
								dbgstream[9] << PC[q1][q2][q3].SS[f].t << " ";
							}
							dbgstream[9] << std::endl << std::endl;
						}
					}
				}
				WriteDebugInfo(dbgstream[10], macro_D.C);
				WriteDebugInfo(dbgstream[11], macro_D_in.C);
				WriteDebugInfo(dbgstream[12], macro_Sgm.C);
				WriteDebugInfo(dbgstream[13], macro_dSgm.C);
				WriteDebugInfo(dbgstream[14], macro_E.C);
			}
			CURR_STEP++;
			PROC_STEP++;
			DEBUG_STEP++;
		}

		if (UNLOADING)	//������� ��������� ����������������� ������
		{
			std::cout << std::endl << " Unloading # " << cycle;
			std::cout << std::endl << "        0.00%";
			while (abs(macro_stress) > final_stress)
			{
				//����������
				macro_P.setZero();
				macro_D_in.setZero();
				for (int q1 = 0; q1 < fragm_count; q1++)
				{
					for (int q2 = 0; q2 < fragm_count; q2++)
					{
						for (int q3 = 0; q3 < fragm_count; q3++)
						{
							macro_D_in += PC[q1][q2][q3].d_in;
							macro_P += PC[q1][q2][q3].p.ToLSK(PC[q1][q2][q3].o);
						}
					}
				}
				macro_D_in /= total_fragm_count;
				macro_P /= total_fragm_count;

				//�������������
				macro_P.Symmetrize();
				macro_D = UnloadingStrainCalc(macro_P, macro_D_in, macro_Sgm, lam);
				Tensor b = macro_D;
				b *= dt;
				macro_E += b;
				macro_strain = macro_E.doubleScalMult(macro_E);
				macro_strain = SQRT2_3*sqrt(macro_strain);
				
				macro_dSgm = TensionStressCalc(macro_P, macro_D_in, macro_D);
				macro_dSgm *= dt;
				macro_Sgm += macro_dSgm;
				macro_stress = macro_Sgm.doubleScalMult(macro_Sgm);
				macro_stress = SQRT3_2*sqrt(macro_stress);

				for (int q1 = 0; q1 < fragm_count; q1++)
				{
					for (int q2 = 0; q2 < fragm_count; q2++)
					{
						for (int q3 = 0; q3 < fragm_count; q3++)
						{

							/**************************************************
							************       ��������� � ���       **********
							**************************************************/

							Tensor O = PC[q1][q2][q3].o;	//�������� ������
							Tensor OT = O;
							OT.Transp();
							PC[q1][q2][q3].d = O*macro_D*OT;
							PC[q1][q2][q3].w = O*macro_W*OT;
							PC[q1][q2][q3].sgm = O*PC[q1][q2][q3].sgm*OT;
							PC[q1][q2][q3].d_in = O*PC[q1][q2][q3].d_in*OT;


							/***************************************************
							***********       ������������� ���      ***********
							***************************************************/

							PC[q1][q2][q3].NDScalc();

							if (HARDENING_BASE)			//������� ����������
							{
								Base_hardening(&PC[q1][q2][q3]);
							}
							if (HARDENING_BOUND)
							{
								Boundary_hardening(&PC[q1][q2][q3]);
							}
							if (ROTATIONS_TAYLOR)		//������� �� �������
							{
								Taylor_rotations(&PC[q1][q2][q3]);
							}
							if (ROTATIONS_TRUSOV)		//������� �� �������
							{
								Trusov_rotations(&PC[q1][q2][q3]);
							}
							if (ROTATIONS_TRUSOV && ROTATIONS_HARDENING)	//����������� ����������
							{
								Rotation_hardening(&PC[q1][q2][q3]);
							}

							/**************************************************
							************       ��������� � ���       **********
							**************************************************/
							PC[q1][q2][q3].sgm = OT*PC[q1][q2][q3].sgm*O;
							PC[q1][q2][q3].d_in = OT*PC[q1][q2][q3].d_in*O;

						}
					}
				}

				/************************************************************
				***********	        �������� ���������� 	      ***********
				************************************************************/

				double progress = (final_stress / macro_stress) * 100;

				if (PROC_STEP == proc_period/4)
				{
					PROC_STEP = 0;
					std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
					printf("        %2.2f", progress);
					std::cout << "%";
				}

				/************************************************************
				***********	    ������ ������ ��� �������� ���    ***********
				************************************************************/
				if (progress - PLOT_STEP > plot_period && plot_period > 0)
				{
					//PLOT_STEP = 0;
					StrainStream.write((char *)&macro_E.C[0][0], sizeof macro_E.C[0][0]);
					StressStream.write((char *)&macro_Sgm.C[0][0], sizeof macro_Sgm.C[0][0]);
					double ActiveSysCount = 0;		//������� ���-�� �������� ������ ����������
					for (int q1 = 0; q1 < fragm_count; q1++)
					{
						for (int q2 = 0; q2 < fragm_count; q2++)
						{
							for (int q3 = 0; q3 < fragm_count; q3++)
							{
								StrainStreamAll.write((char *)&PC[q1][q2][q3].e.C[0][0], sizeof PC[q1][q2][q3].e.C[0][0]);
								StressStreamAll.write((char *)&PC[q1][q2][q3].sgm.C[0][0], sizeof PC[q1][q2][q3].sgm.C[0][0]);
								for (int i = 0; i < PC[q1][q2][q3].SS_count; i++)
								{
									if (PC[q1][q2][q3].SS[i].dgm != 0) ActiveSysCount++;//������� �������� ��
								}
							
							}
						}
					}
					ActiveSysCount /= total_fragm_count;
					ActiveSysStream.write((char *)&ActiveSysCount, sizeof ActiveSysCount);
					PLOT_STEP = progress;
				}

				/************************************************************
				***********	      ���������� �������� �����	      ***********
				************************************************************/
				if (progress - POLUS_STEP > polus_period && polus_period > 0)
				{
					for (int q1 = 0; q1 < fragm_count; q1++)
					{
						for (int q2 = 0; q2 < fragm_count; q2++)
						{
							for (int q3 = 0; q3 < fragm_count; q3++)
							{
								GetPoleFig(&PC[q1][q2][q3]);
						//		GetSST(&PC[q1][q2][q3]); //�����������������, ���� ����� ��� (BETA)
							}
						}
					}
					POLUS_STEP = progress;
				}
				/************************************************************
				***********	       ������ ��������� ������	      ***********
				************************************************************/
				if (DEBUG_STEP == debug_period)
				{
					DEBUG_STEP = 0;

					for (int i = 0; i < 15; i++)
					{
						dbgstream[i] << "********      STEP " << CURR_STEP << "      ********" << std::endl << std::endl;
					}

					for (int q1 = 0; q1 < fragm_count; q1++)
					{
						for (int q2 = 0; q2 < fragm_count; q2++)
						{
							for (int q3 = 0; q3 < fragm_count; q3++)
							{
								WriteDebugInfo(dbgstream[0], PC[q1][q2][q3].o.C);
								WriteDebugInfo(dbgstream[1], PC[q1][q2][q3].e.C);
								WriteDebugInfo(dbgstream[2], PC[q1][q2][q3].d.C);
								WriteDebugInfo(dbgstream[3], PC[q1][q2][q3].sgm.C);
								WriteDebugInfo(dbgstream[4], PC[q1][q2][q3].om.C);
								WriteDebugInfo(dbgstream[5], PC[q1][q2][q3].dsgm.C);
								WriteDebugInfo(dbgstream[6], PC[q1][q2][q3].d_in.C);
								WriteDebugInfo(dbgstream[7], PC[q1][q2][q3].w.C);
								for (int f = 0; f < PC[q1][q2][q3].SS_count; f++)
								{
									dbgstream[8] << PC[q1][q2][q3].SS[f].dgm << " ";
								}
								dbgstream[8] << std::endl << std::endl;
								for (int f = 0; f < PC[q1][q2][q3].SS_count; f++)
								{
									dbgstream[9] << PC[q1][q2][q3].SS[f].t << " ";
								}
								dbgstream[9] << std::endl << std::endl;
							}
						}
					}
					WriteDebugInfo(dbgstream[10], macro_D.C);
					WriteDebugInfo(dbgstream[11], macro_D_in.C);
					WriteDebugInfo(dbgstream[12], macro_Sgm.C);
					WriteDebugInfo(dbgstream[13], macro_dSgm.C);
					WriteDebugInfo(dbgstream[14], macro_E.C);
				}
				
				CURR_STEP++;
				PROC_STEP++;
				DEBUG_STEP++;
			}
		}

		if (cycle_count > 1 && REAL_UNIAX)
		{
			macro_D.C[0][0] = pow(-1, cycle + 1) * tension_component;	//������ ���� ������������� ����������
			strain_max += strain_max * addition_strain;					//�������� ������ �������������
		}
	}

	//���������� �������� �������� �����
/*	for (int q1 = 0; q1 < fragm_count; q1++)
	{
		for (int q2 = 0; q2 < fragm_count; q2++)
		{
			for (int q3 = 0; q3 < fragm_count; q3++)
			{
				GetPoleFig(PC[q1][q2][q3]);
				GetSST(PC[q1][q2][q3]);
			}
		}
	}
	*/
	t2 = clock();//��������� ������� �������

	/***********************************************************
	*********      �������� ���� �������� �������      *********
	***********************************************************/

	StrainStream.close();
	StressStream.close();
	StrainStreamAll.close();
	StressStreamAll.close();
	ActiveSysStream.close();
	for (int i = 0; i < 4; i++)
	{
		TestStream[i].close();
	}


	if (debug_period > 0)
	{
		for (int i = 0; i < 15; i++)
		{
			dbgstream[i].close();
		}
	}
	/***********************************************************
	*********       ���������� � ������� � �����		********
	*********	(������������� ������������ � ������	********
	*********	���������� � ��������� �����������,		********
	*********		���� ��� ���� �������� �� ��		********
	***********************************************************/

	if (isnan(macro_strain)) std::cout << std::endl << " Calculation ERROR!" << std::endl;
	else std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << " Done    " << std::endl;
	std::cout << " __________________________________________________" << std::endl;
	std::cout << " Processing time: " << (t2 - t1) / 1000.0 << " sec" << std::endl;
	std::cout << " Number of steps: " << CURR_STEP << std::endl;
	if (isnan(macro_strain))//���� �� ������������� ������ - ��������
	{
		std::cout << " __________________________________________________" << std::endl;
		std::cout << " Press any key or STOP button to exit...";
		std::system("title Done");
		std::cin.get();
	}

	/************************************************************
	***********      ������������ ������� ������      ***********
	************************************************************/

	for (int i = 0; i < fragm_count; i++)
	{
		for (int j = 0; j < fragm_count; j++)
		{
			delete[] PC[i][j];
		}
		delete[] PC[i];
	}
	delete[] PC;


	/************************************************************
	*******      �������� ������ � ������ ����������      *******
	************************************************************/
	HWND hwnd;
	hwnd = ::FindWindow(NULL, L"������ ���������� �������");
	if (hwnd != NULL)
	{
		::SendMessage(hwnd, WM_USER + 1, CURR_STEP, (t2 - t1));
		//�������� 1 - ���������� �����
		//�������� 2 - ����������� ����� (��)
	}
	return 0;
}

