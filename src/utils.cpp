/*
 * utils.cpp
 *
 *  Created on: 23 de nov de 2018
 *      Author: fabricio
 */

#include "tensor.h"

namespace Utils
{
	double volumeTetrahedonGet(Vector3 p1,Vector3 p2,Vector3 p3,Vector3 p4)
	{
		double x1 = p1.x;
		double x2 = p2.x;
		double x3 = p3.x;
		double x4 = p4.x;
	 
		double y1 = p1.y;
		double y2 = p2.y;
		double y3 = p3.y;
		double y4 = p4.y;
	 
		double z1 = p1.z;
		double z2 = p2.z;
		double z3 = p3.z;
		double z4 = p4.z;
	 
		double volume = (x1*y3*z2)/6.0 - (x1*y2*z3)/6.0 +
						(x2*y1*z3)/6.0 - (x2*y3*z1)/6.0 - 
						(x3*y1*z2)/6.0 + (x3*y2*z1)/6.0 + 
						(x1*y2*z4)/6.0 - (x1*y4*z2)/6.0 - 
						(x2*y1*z4)/6.0 + (x2*y4*z1)/6.0 + 
						(x4*y1*z2)/6.0 - (x4*y2*z1)/6.0 - 
						(x1*y3*z4)/6.0 + (x1*y4*z3)/6.0 + 
						(x3*y1*z4)/6.0 - (x3*y4*z1)/6.0 - 
						(x4*y1*z3)/6.0 + (x4*y3*z1)/6.0 + 
						(x2*y3*z4)/6.0 - (x2*y4*z3)/6.0 - 
						(x3*y2*z4)/6.0 + (x3*y4*z2)/6.0 + 
						(x4*y2*z3)/6.0 - (x4*y3*z2)/6.0;
						
		return std::abs(volume);
	}
}
