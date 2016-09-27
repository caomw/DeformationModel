#ifndef __MATHCORE_H 
#define __MATHCORE_H

namespace model
{
	/**********************************************************************
	***********				�������������� ���������			***********
	**********************************************************************/
	const double SQRT3 = 1.732050807568;		//������ �� 3-�
	const double SQRT2 = 1.414213562373;		//������ �� 2-�
	const double SQRT2_3 = 0.816496580927;		//������ �� 2/3
	const double SQRT3_2 = 1.224744871391;		//������ �� 3/2

	const double PI = 3.141592653589;			//����� ��
	const double PI_2 = 1.570796326794;			//��/2
	const double PIx2 = 6.283185307179;			//2��

	const double EPS = 1e-10;					//����� ��������

	const int DIM = 3;	//����������� ������������


	/***********************************************************************
	****************                  ������               *****************
	***********************************************************************/

	class Vector
	{
	public:
		double C[DIM];					//���������� �������

		double getNorm();				//���������� ����� �������
		int Normalize();				//����������� ������
		void setZero();					//�������� ���������� �������
		void set(double,
			double, double);			//����� �������� ���������
		double ScalMult(Vector v);		//��������� ������������ ��������
		Vector VectMult(Vector v);		//��������� ������������ ��������
			
		Vector operator + (Vector);		//�������� ��������
		Vector operator - (Vector);		//�������� ���������
		void operator += (Vector);		//�������� ����������� �������
		void operator -= (Vector);		//�������� ��������� �������
		void operator *= (double);		//�������� ��������� ������� �� �����
		int operator /= (double);		//�������� ������� ������� �� �����

		Vector();
		~Vector();

	private:

	};


	/***********************************************************************
	****************                  ������               *****************
	***********************************************************************/

	class Tensor
	{
	public:
		double C[DIM][DIM];				//���������� �������

		double getDet();				//���������� ������������ ������� ���������
		void setZero();					//�������� ���������� �������
		void setUnit();					//������ ������� ��������� ������� ���������
		void Transp();					//������������� ������� ��������� �������
		double doubleScalMult(Tensor);	//������ (������� ��������� ������������ ��������)
		Tensor getSymmetryPart();		//���������� ������������ ����� �������
		Tensor getAntiSymmetryPart();	//���������� ���������������� ����� �������
		Vector getRow(int);				//������ �� ��������� �������� ������
		Vector getCol(int);				//������ �� ��������� ��������� ������
		
		Tensor operator + (Tensor);		//�������� �������� ��������
		Tensor operator - (Tensor);		//�������� ��������� ��������
		Tensor operator * (Tensor);		//�������� ��������� ��������
		void operator += (Tensor);		//�������� ����������� �������
		void operator -= (Tensor);		//�������� ��������� �������
		void operator *= (Tensor);		//�������� ���������� �� ������
		void operator *= (double);		//�������� ��������� ������� �� �����
		int operator /= (double);		//�������� ������� ������� �� �����
		
		Tensor();
		~Tensor();

	private:

	};

	/***********************************************************************
	****************         ������� ����������            *****************
	***********************************************************************/

	class SlipSystem
	{
	public:
		Vector n;						//������ �������
		Vector b;						//������ ��������
		double t;						//����������� ����������� ����������
		double tc;						//����������� ����������� ����������
		double dgm;						//�������� ������
		double gmm;						//����������� �����

		void Initialize(
			double, double,	double,
			double,	double, double);	//������������� ��������
		
		SlipSystem();
		~SlipSystem();
	};

	/***********************************************************************
	**************         ������ ��������� �����            **************
	***********************************************************************/

	class Tensor4
	{
	public:
		double C[DIM][DIM][DIM][DIM];	//���������� �������

		void setZero();					//��������� ��������� �������
		void Symmetrize();				//������������� ��������� �������

		Tensor4 ToLSK(Tensor O);		//������� ��������� ������� � ���

		void operator += (Tensor4);		//�������� ����������� �������
		void operator -= (Tensor4);		//�������� ��������� �������
		void operator *= (double);		//�������� ��������� ������� �� �����
		int operator /= (double);		//�������� ������� ������� �� �����

		Tensor4();
		~Tensor4();
	private:
	};


	int LeviCivit(int i, int j, int k);	//������-������ ����-������


	Tensor VectMult(Vector, Tensor);	//��������� ������������ ������� �� ������
	Tensor VectMult(Tensor, Vector);	//��������� ������������ ������� �� ������
	Vector ScalMult(Vector, Tensor);	//��������� ������������ ������� �� ������
	Vector ScalMult(Tensor, Vector);	//��������� ������������ ������� �� ������
}

#endif __MATHCORE_H