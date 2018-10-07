/*
 * model.h
 *
 *  Created on: Oct 6, 2018
 *      Author: fabricio
 */

#ifndef INC_MODEL_H_
#define INC_MODEL_H_

#include "tensor.h"

namespace Model
{
	struct DemBox
	{
		Vector3 pos;
		double zbase;
		double dx,dy,dz;
		int matid;
		DemBox():pos(0.0),zbase(0.0),dx(0.0),dy(0.0),dz(0.0),matid(0){}
	};

	struct Particle
	{
		Vector3 pos;
		double vol;
		double lp;
		int matid;
		Particle():pos(0.0),vol(0.0),lp(0.0),matid(0){}
	};
	
	std::vector<DemBox>& GetDemVector();
	std::vector<Particle>& GetParticleVector();
	
	void CreateMPMmodel();
}

#endif /* INC_MODEL_H_ */
