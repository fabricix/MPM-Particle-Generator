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
	// domain partition for fracture generation
	struct DemBox
	{
	    Vector3 pos;
	    double zbase;
	    double dx;
	    double dy;
	    int matid;
	    DemBox():pos(0.0),zbase(0.0),dx(0.0),dy(0.0),matid(0){}
	};
	std::vector<DemBox>& GetDemVector();
}

#endif /* INC_MODEL_H_ */
