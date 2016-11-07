﻿#include "stdafx.h"

#include "tension.h"
#include "Eigen/Eigen"

namespace model
{
	Tensor TensionStrainCalc(const Tensor4 P,const Tensor D_in, const double tens_comp)
	{
		Eigen::MatrixXd A(5, 5);
		Eigen::VectorXd B(5), X(5);

		A(0, 0) = P.C[1][1][0][1] + P.C[1][1][1][0];
		A(0, 1) = P.C[1][1][0][2] + P.C[1][1][2][0];
		A(0, 2) = P.C[1][1][1][2] + P.C[1][1][2][1];
		A(0, 3) = P.C[1][1][1][1];
		A(0, 4) = P.C[1][1][2][2];

		A(1, 0) = P.C[2][2][0][1] + P.C[2][2][1][0];
		A(1, 1) = P.C[2][2][0][2] + P.C[2][2][2][0];
		A(1, 2) = P.C[2][2][1][2] + P.C[2][2][2][1];
		A(1, 3) = P.C[2][2][1][1];
		A(1, 4) = P.C[2][2][2][2];

		A(2, 0) = P.C[0][1][0][1] + P.C[0][1][1][0];
		A(2, 1) = P.C[0][1][0][2] + P.C[0][1][2][0];
		A(2, 2) = P.C[0][1][1][2] + P.C[0][1][2][1];
		A(2, 3) = P.C[0][1][1][1];
		A(2, 4) = P.C[0][1][2][2];

		A(3, 0) = P.C[0][2][0][1] + P.C[0][2][1][0];
		A(3, 1) = P.C[0][2][0][2] + P.C[0][2][2][0];
		A(3, 2) = P.C[0][2][1][2] + P.C[0][2][2][1];
		A(3, 3) = P.C[0][2][1][1];
		A(3, 4) = P.C[0][2][2][2];

		A(4, 0) = P.C[1][2][0][1] + P.C[1][2][1][0];
		A(4, 1) = P.C[1][2][0][2] + P.C[1][2][2][0];
		A(4, 2) = P.C[1][2][1][2] + P.C[1][2][2][1];
		A(4, 3) = P.C[1][2][1][1];
		A(4, 4) = P.C[1][2][2][2];

		B(0) = P.C[1][1][0][0] * (-tens_comp + D_in.C[0][0]) + (P.C[1][1][0][1] + P.C[1][1][1][0]) * D_in.C[0][1] +
			(P.C[1][1][0][2] + P.C[1][1][2][0]) * D_in.C[0][2] + (P.C[1][1][1][2] + P.C[1][1][2][1]) * D_in.C[1][2] +
			P.C[1][1][1][1] * D_in.C[1][1] + P.C[1][1][2][2] * D_in.C[2][2];

		B(1) = P.C[2][2][0][0] * (-tens_comp + D_in.C[0][0]) + (P.C[2][2][0][1] + P.C[2][2][1][0])*D_in.C[0][1] +
			(P.C[2][2][0][2] + P.C[2][2][2][0])*D_in.C[0][2] + (P.C[2][2][1][2] + P.C[2][2][2][1])*D_in.C[1][2] +
			P.C[2][2][1][1] * D_in.C[1][1] + P.C[2][2][2][2] * D_in.C[2][2];

		B(2) = P.C[0][1][0][0] * (-tens_comp + D_in.C[0][0]) + (P.C[0][1][0][1] + P.C[0][1][1][0])*D_in.C[0][1] +
			(P.C[0][1][0][2] + P.C[0][1][2][0])*D_in.C[0][2] + (P.C[0][1][1][2] + P.C[0][1][2][1])*D_in.C[1][2] +
			P.C[0][1][1][1] * D_in.C[1][1] + P.C[0][1][2][2] * D_in.C[2][2];

		B(3) = P.C[0][2][0][0] * (-tens_comp + D_in.C[0][0]) + (P.C[0][2][0][1] + P.C[0][2][1][0])*D_in.C[0][1] +
			(P.C[0][2][0][2] + P.C[0][2][2][0])*D_in.C[0][2] + (P.C[0][2][1][2] + P.C[0][2][2][1])*D_in.C[1][2] +
			P.C[0][2][1][1] * D_in.C[1][1] + P.C[0][2][2][2] * D_in.C[2][2];

		B(4) = P.C[1][2][0][0] * (-tens_comp + D_in.C[0][0]) + (P.C[1][2][0][1] + P.C[1][2][1][0])*D_in.C[0][1] +
			(P.C[1][2][0][2] + P.C[1][2][2][0])*D_in.C[0][2] + (P.C[1][2][1][2] + P.C[1][2][2][1])*D_in.C[1][2] +
			P.C[1][2][1][1] * D_in.C[1][1] + P.C[1][2][2][2] * D_in.C[2][2];

		X = A.fullPivHouseholderQr().solve(B);

		Tensor res;
		res.C[0][1] = res.C[1][0] = X(0);
		res.C[0][2] = res.C[2][0] = X(1);
		res.C[1][2] = res.C[2][1] = X(2);
		res.C[1][1] = X(3);
		res.C[2][2] = X(4);
		res.C[0][0] = tens_comp;
		return res;
	}
	Tensor TensionStressCalc(Tensor4 &P, Tensor &D_in, Tensor &D)
	{
		Tensor res;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						res.C[i][j] += P.C[i][j][k][l] * (D.C[l][k] - D_in.C[l][k]);
					}
				}
			}
		}
		return res;
	}
	Tensor UnloadingStrainCalc(Tensor4 &P, Tensor &D_in, Tensor &Sgm, double lam)
	{
		Eigen::MatrixXd Ar(6, 6);
		Eigen::VectorXd Br(6), Xr(6);

		Ar(0, 0) = P.C[1][1][0][1] + P.C[1][1][1][0];
		Ar(0, 1) = P.C[1][1][0][2] + P.C[1][1][2][0];
		Ar(0, 2) = P.C[1][1][1][2] + P.C[1][1][2][1];
		Ar(0, 3) = P.C[1][1][1][1];
		Ar(0, 4) = P.C[1][1][2][2];
		Ar(0, 5) = P.C[1][1][0][0];

		Ar(1, 0) = P.C[2][2][0][1] + P.C[2][2][1][0];
		Ar(1, 1) = P.C[2][2][0][2] + P.C[2][2][2][0];
		Ar(1, 2) = P.C[2][2][1][2] + P.C[2][2][2][1];
		Ar(1, 3) = P.C[2][2][1][1];
		Ar(1, 4) = P.C[2][2][2][2];
		Ar(1, 5) = P.C[2][2][0][0];

		Ar(2, 0) = P.C[0][1][0][1] + P.C[0][1][1][0];
		Ar(2, 1) = P.C[0][1][0][2] + P.C[0][1][2][0];
		Ar(2, 2) = P.C[0][1][1][2] + P.C[0][1][2][1];
		Ar(2, 3) = P.C[0][1][1][1];
		Ar(2, 4) = P.C[0][1][2][2];
		Ar(2, 5) = P.C[0][1][0][0];

		Ar(3, 0) = P.C[0][2][0][1] + P.C[0][2][1][0];
		Ar(3, 1) = P.C[0][2][0][2] + P.C[0][2][2][0];
		Ar(3, 2) = P.C[0][2][1][2] + P.C[0][2][2][1];
		Ar(3, 3) = P.C[0][2][1][1];
		Ar(3, 4) = P.C[0][2][2][2];
		Ar(3, 5) = P.C[0][2][0][0];

		Ar(4, 0) = P.C[1][2][0][1] + P.C[1][2][1][0];
		Ar(4, 1) = P.C[1][2][0][2] + P.C[1][2][2][0];
		Ar(4, 2) = P.C[1][2][1][2] + P.C[1][2][2][1];
		Ar(4, 3) = P.C[1][2][1][1];
		Ar(4, 4) = P.C[1][2][2][2];
		Ar(4, 5) = P.C[1][2][0][0];

		Ar(5, 0) = P.C[0][0][0][1] + P.C[0][0][1][0];
		Ar(5, 1) = P.C[0][0][0][2] + P.C[0][0][2][0];
		Ar(5, 2) = P.C[0][0][1][2] + P.C[0][0][2][1];
		Ar(5, 3) = P.C[0][0][1][1];
		Ar(5, 4) = P.C[0][0][2][2];
		Ar(5, 5) = P.C[0][0][0][0];

		Br(0) = -lam*Sgm.C[1][1] + (P.C[1][1][0][0] * D_in.C[0][0] + (P.C[1][1][0][1] + P.C[1][1][1][0])*D_in.C[0][1] +
			(P.C[1][1][0][2] + P.C[1][1][2][0])*D_in.C[0][2] + (P.C[1][1][1][2] + P.C[1][1][2][1])*D_in.C[1][2] +
			P.C[1][1][1][1] * D_in.C[1][1] + P.C[1][1][2][2] * D_in.C[2][2]);
		Br(1) = -lam*Sgm.C[2][2] + (P.C[2][2][0][0] * D_in.C[0][0] + (P.C[2][2][0][1] + P.C[2][2][1][0])*D_in.C[0][1] +
			(P.C[2][2][0][2] + P.C[2][2][2][0])*D_in.C[0][2] + (P.C[2][2][1][2] + P.C[2][2][2][1])*D_in.C[1][2] +
			P.C[2][2][1][1] * D_in.C[1][1] + P.C[2][2][2][2] * D_in.C[2][2]);
		Br(2) = -lam*Sgm.C[0][1] + (P.C[0][1][0][0] * D_in.C[0][0] + (P.C[0][1][0][1] + P.C[0][1][1][0])*D_in.C[0][1] +
			(P.C[0][1][0][2] + P.C[0][1][2][0])*D_in.C[0][2] + (P.C[0][1][1][2] + P.C[0][1][2][1])*D_in.C[1][2] +
			P.C[0][1][1][1] * D_in.C[1][1] + P.C[0][1][2][2] * D_in.C[2][2]);
		Br(3) = -lam*Sgm.C[0][2] + (P.C[0][2][0][0] * D_in.C[0][0] + (P.C[0][2][0][1] + P.C[0][2][1][0])*D_in.C[0][1] +
			(P.C[0][2][0][2] + P.C[0][2][2][0])*D_in.C[0][2] + (P.C[0][2][1][2] + P.C[0][2][2][1])*D_in.C[1][2] +
			P.C[0][2][1][1] * D_in.C[1][1] + P.C[0][2][2][2] * D_in.C[2][2]);
		Br(4) = -lam*Sgm.C[1][2] + (P.C[1][2][0][0] * D_in.C[0][0] + (P.C[1][2][0][1] + P.C[1][2][1][0])*D_in.C[0][1] +
			(P.C[1][2][0][2] + P.C[1][2][2][0])*D_in.C[0][2] + (P.C[1][2][1][2] + P.C[1][2][2][1])*D_in.C[1][2] +
			P.C[1][2][1][1] * D_in.C[1][1] + P.C[1][2][2][2] * D_in.C[2][2]);
		Br(5) = -lam*Sgm.C[0][0] + (P.C[0][0][0][0] * D_in.C[0][0] + (P.C[0][0][0][1] + P.C[0][0][1][0])*D_in.C[0][1] +
			(P.C[0][0][0][2] + P.C[0][0][2][0])*D_in.C[0][2] + (P.C[0][0][1][2] + P.C[0][0][2][1])*D_in.C[1][2] +
			P.C[0][0][1][1] * D_in.C[1][1] + P.C[0][0][2][2] * D_in.C[2][2]);

		Xr = Ar.fullPivHouseholderQr().solve(Br);
		Tensor res;
		res.C[0][1] = res.C[1][0] = Xr(0);
		res.C[0][2] = res.C[2][0] = Xr(1);
		res.C[1][2] = res.C[2][1] = Xr(2);
		res.C[1][1] = Xr(3);
		res.C[2][2] = Xr(4);
		res.C[0][0] = Xr(5);
		return res;
	}
}