/*
 * tensor.cpp
 *
 *  Created on: Oct 6, 2018
 *      Author: fabricio fernandez <fabricio.hmf@gmail.com>
 */

#include <iostream>
#include <float.h>
#include <math.h>

#include "matrix.h"


static void JacobiCyclicMethod(double eigenvalues[], double *eigenvectors, double *A, int n)
{
	int i, j, k, m;
	double *pAk, *pAm, *p_r, *p_e;
	double threshold_norm;
	double threshold;
	double tan_phi, sin_phi, cos_phi, tan2_phi, sin2_phi, cos2_phi;
	double sin_2phi, cot_2phi;
	double dum1;
	double dum2;
	double dum3;
	//double r;
	double max;

	// Take care of trivial cases

	if ( n < 1) return;
	if ( n == 1) {
	  eigenvalues[0] = *A;
	  *eigenvectors = 1.0;
	  return;
	}

	// Initialize the eigenvalues to the identity matrix.

	for (p_e = eigenvectors, i = 0; i < n; i++)
	  for (j = 0; j < n; p_e++, j++)
		 if (i == j) *p_e = 1.0; else *p_e = 0.0;

	// Calculate the threshold and threshold_norm.

	for (threshold = 0.0, pAk = A, i = 0; i < ( n - 1 ); pAk += n, i++)
	  for (j = i + 1; j < n; j++) threshold += *(pAk + j) * *(pAk + j);
	threshold = sqrt(threshold + threshold);
	threshold_norm = threshold * DBL_EPSILON;
	max = threshold + 1.0;
	while (threshold > threshold_norm) {
	  threshold /= 10.0;
	  if (max < threshold) continue;
	  max = 0.0;
	  for (pAk = A, k = 0; k < (n-1); pAk += n, k++) {
		 for (pAm = pAk + n, m = k + 1; m < n; pAm += n, m++) {
			if ( fabs(*(pAk + m)) < threshold ) continue;

				 // Calculate the sin and cos of the rotation angle which
				 // annihilates A[k][m].

				cot_2phi = 0.5 * ( *(pAk + k) - *(pAm + m) ) / *(pAk + m);
				dum1 = sqrt( cot_2phi * cot_2phi + 1.0);
				if (cot_2phi < 0.0) dum1 = -dum1;
				tan_phi = -cot_2phi + dum1;
				tan2_phi = tan_phi * tan_phi;
				sin2_phi = tan2_phi / (1.0 + tan2_phi);
				cos2_phi = 1.0 - sin2_phi;
				sin_phi = sqrt(sin2_phi);
				if (tan_phi < 0.0) sin_phi = - sin_phi;
				cos_phi = sqrt(cos2_phi);
				sin_2phi = 2.0 * sin_phi * cos_phi;
				// cos_2phi = cos2_phi - sin2_phi;

				// Rotate columns k and m for both the matrix A
				//     and the matrix of eigenvectors.

				p_r = A;
				dum1 = *(pAk + k);
				dum2 = *(pAm + m);
				dum3 = *(pAk + m);
				*(pAk + k) = dum1 * cos2_phi + dum2 * sin2_phi + dum3 * sin_2phi;
				*(pAm + m) = dum1 * sin2_phi + dum2 * cos2_phi - dum3 * sin_2phi;
				*(pAk + m) = 0.0;
				*(pAm + k) = 0.0;
				for (i = 0; i < n; p_r += n, i++) {
				   if ( (i == k) || (i == m) ) continue;
				   if ( i < k ) dum1 = *(p_r + k); else dum1 = *(pAk + i);
				   if ( i < m ) dum2 = *(p_r + m); else dum2 = *(pAm + i);
				   dum3 = dum1 * cos_phi + dum2 * sin_phi;
				   if ( i < k ) *(p_r + k) = dum3; else *(pAk + i) = dum3;
				   dum3 = - dum1 * sin_phi + dum2 * cos_phi;
				   if ( i < m ) *(p_r + m) = dum3; else *(pAm + i) = dum3;
			}
			for (p_e = eigenvectors, i = 0; i < n; p_e += n, i++) {
			   dum1 = *(p_e + k);
			   dum2 = *(p_e + m);
			   *(p_e + k) = dum1 * cos_phi + dum2 * sin_phi;
			   *(p_e + m) = - dum1 * sin_phi + dum2 * cos_phi;
			}
		 }
		 for (i = 0; i < n; i++)
			if ( i == k ) continue;
			else if ( max < fabs(*(pAk + i))) max = fabs(*(pAk + i));
	  }
	}
	for (pAk = A, k = 0; k < n; pAk += n, k++) eigenvalues[k] = *(pAk + k);
}

void eigen(Matrix3 s, Vector3& eival, Matrix3& eivec)
{
	double A[3][3];
	double eigenvalues[3];
	double eigenvectors[3][3];

	A[0][0]= s.xx;
	A[0][1]= s.xy;
	A[0][2]= s.xz;

	A[1][0]= s.yx;
	A[1][1]= s.yy;
	A[1][2]= s.yz;

	A[2][0]= s.zx;
	A[2][1]= s.zy;
	A[2][2]= s.zz;

	JacobiCyclicMethod(eigenvalues, (double*)eigenvectors, (double*)A,3);

	eival.x = eigenvalues[0];
	eival.y = eigenvalues[1];
	eival.z = eigenvalues[2];

	eivec.xx = eigenvectors[0][0];
	eivec.xy = eigenvectors[0][1];
	eivec.xz = eigenvectors[0][2];

	eivec.yx = eigenvectors[1][0];
	eivec.yy = eigenvectors[1][1];
	eivec.yz = eigenvectors[1][2];

	eivec.zx = eigenvectors[2][0];
	eivec.zy = eigenvectors[2][1];
	eivec.zz = eigenvectors[2][2];
}

void eigen(Matrix3 s, Vector3& eival,Vector3& n1,Vector3& n2,Vector3& n3)
{
	double A[3][3];
	double eigenvalues[3];
	double eigenvectors[3][3];

	A[0][0]= s.xx;
	A[0][1]= s.xy;
	A[0][2]= s.xz;

	A[1][0]= s.yx;
	A[1][1]= s.yy;
	A[1][2]= s.yz;

	A[2][0]= s.zx;
	A[2][1]= s.zy;
	A[2][2]= s.zz;

	JacobiCyclicMethod(eigenvalues, (double*)eigenvectors, (double*)A,3);

	double val1,val2,val3;

	if (eigenvalues[0]<=eigenvalues[1]&&eigenvalues[0]<=eigenvalues[2])
	{
	 	val1 = eigenvalues[0];
		n1.x = eigenvectors[0][0];
		n1.y = eigenvectors[1][0];
		n1.z = eigenvectors[2][0];

		if (eigenvalues[1]<=eigenvalues[2])
		{
			val2 = eigenvalues[1];
			n2.x = eigenvectors[0][1];
			n2.y = eigenvectors[1][1];
			n2.z = eigenvectors[2][1];

			val3 = eigenvalues[2];
			n3.x = eigenvectors[0][2];
			n3.y = eigenvectors[1][2];
			n3.z = eigenvectors[2][2];
		}
		if (eigenvalues[2]<=eigenvalues[1])
		{
			val2 = eigenvalues[2];
			n2.x = eigenvectors[0][2];
			n2.y = eigenvectors[1][2];
			n2.z = eigenvectors[2][2];

			val3 = eigenvalues[1];
		    n3.x = eigenvectors[0][1];
			n3.y = eigenvectors[1][1];
			n3.z = eigenvectors[2][1];
		}
	}
    if (eigenvalues[1]<=eigenvalues[0]&&eigenvalues[1]<=eigenvalues[2])
	{
	 	val1 = eigenvalues[1];
		n1.x = eigenvectors[0][1];
		n1.y = eigenvectors[1][1];
		n1.z = eigenvectors[2][1];

		if (eigenvalues[0]<=eigenvalues[2])
		{
			val2 = eigenvalues[0];
			n2.x = eigenvectors[0][0];
			n2.y = eigenvectors[1][0];
			n2.z = eigenvectors[2][0];

			val3 = eigenvalues[2];
			n3.x = eigenvectors[0][2];
			n3.y = eigenvectors[1][2];
			n3.z = eigenvectors[2][2];
		}
		if (eigenvalues[2]<=eigenvalues[0])
		{
			val2 = eigenvalues[2];
			n2.x = eigenvectors[0][2];
			n2.y = eigenvectors[1][2];
			n2.z = eigenvectors[2][2];

			val3 = eigenvalues[0];
			n3.x = eigenvectors[0][0];
			n3.y = eigenvectors[1][0];
			n3.z = eigenvectors[2][0];
		}
	}
	if (eigenvalues[2]<=eigenvalues[0]&&eigenvalues[2]<=eigenvalues[1])
	{
	 	val1 = eigenvalues[2];
		n1.x = eigenvectors[0][2];
		n1.y = eigenvectors[1][2];
		n1.z = eigenvectors[2][2];

	 	if (eigenvalues[0]<=eigenvalues[1])
		{
			val2 = eigenvalues[0];
			n2.x = eigenvectors[0][0];
			n2.y = eigenvectors[1][0];
			n2.z = eigenvectors[2][0];

			val3 = eigenvalues[1];
			n3.x = eigenvectors[0][1];
			n3.y = eigenvectors[1][1];
			n3.z = eigenvectors[2][1];
		}
		if (eigenvalues[1]<=eigenvalues[0])
		{
			val2 = eigenvalues[1];
			n2.x = eigenvectors[0][1];
			n2.y = eigenvectors[1][1];
			n2.z = eigenvectors[2][1];

			val3 = eigenvalues[0];
			n3.x = eigenvectors[0][0];
			n3.y = eigenvectors[1][0];
			n3.z = eigenvectors[2][0];
		}
	}

	eival.x = val1;
	eival.y = val2;
	eival.z = val3;

	// eival.x = eigenvalues[0];
	// eival.y = eigenvalues[1];
	// eival.z = eigenvalues[2];

	// eivec.xx = eigenvectors[0][0];
	// eivec.xy = eigenvectors[0][1];
	// eivec.xz = eigenvectors[0][2];

	// eivec.yx = eigenvectors[1][0];
	// eivec.yy = eigenvectors[1][1];
	// eivec.yz = eigenvectors[1][2];

	// eivec.zx = eigenvectors[2][0];
	// eivec.zy = eigenvectors[2][1];
	// eivec.zz = eigenvectors[2][2];

	// std::cout<<"1="<<eival.x<<", n1=("<<eivec.xx<<","<<eivec.yx<<","<<eivec.zx<<")\n";
	// std::cout<<"2="<<eival.y<<", n1=("<<eivec.xy<<","<<eivec.yy<<","<<eivec.zy<<")\n";
	// std::cout<<"3="<<eival.z<<", n1=("<<eivec.xz<<","<<eivec.yz<<","<<eivec.zz<<")\n";
}
