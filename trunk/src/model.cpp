/*
 * model.cpp
 *
 *  Created on: Oct 6, 2018
 *      Author: fabricio
 */


#include <iostream>
#include <vector>

#include <model.h>

namespace Model
{
	static std::vector<DemBox> DemBoxVector;
	std::vector<DemBox>& GetDemVector(){return DemBoxVector;}
}
