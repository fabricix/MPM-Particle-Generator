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
	using namespace std;

	static vector<DemBox> DemBoxVector;
	vector<DemBox>& GetDemVector(){return DemBoxVector;}

	static vector<Particle> ParticleVector;
	vector<Particle>& GetParticleVector(){return ParticleVector;}

	static void CreateMpmModelByDEM()
	{ 
		if(DemBoxVector.size()==0) { return; }

		ParticleVector.clear();
		int pid = 0;

		for (size_t i = 0; i < DemBoxVector.size(); ++i)
		{
			DemBox ibox = DemBoxVector.at(i);
			double dx = ibox.dx; 
			double dy = ibox.dy;
			double dz = ibox.dz;
			double zbase = ibox.zbase;
			double ztop  = ibox.pos.z;
			int matid = ibox.matid;
			double meanL = (dx+dy+dz)/3.0;
            double ppc = 2.0;
            double iz = zbase + (dz*0.25);

			Particle prt1,prt2,prt3,prt4;
			do
			{	
				prt1.pos.x = ibox.pos.x-(dx*0.25);
				prt1.pos.y = ibox.pos.y-(dy*0.25);
				prt1.pos.z = iz;
				prt1.vol = (dx*dy*dz)/8.0;
				prt1.lp  = (meanL/ppc)*0.5;
				prt1.matid = matid;
				pid++;
				prt1.id = pid; 
				ParticleVector.push_back(prt1);

				prt2.pos.x = ibox.pos.x+(dx*0.25);
				prt2.pos.y = ibox.pos.y-(dy*0.25);
				prt2.pos.z = iz;
			    prt2.vol = (dx*dy*dz)/8.0;
				prt2.lp  = (meanL/ppc)*0.5;
				prt2.matid = matid;
				pid++;
				prt2.id = pid; 
				ParticleVector.push_back(prt2);

				prt3.pos.x = ibox.pos.x+(dx*0.25);
				prt3.pos.y = ibox.pos.y+(dy*0.25);
				prt3.pos.z = iz;
			    prt3.vol = (dx*dy*dz)/8.0;
				prt3.lp  = (meanL/ppc)*0.5;
				prt3.matid = matid;
				pid++;
				prt3.id = pid; 
				ParticleVector.push_back(prt3);

				prt4.pos.x = ibox.pos.x-(dx*0.25);
				prt4.pos.y = ibox.pos.y+(dy*0.25);
				prt4.pos.z = iz;
				prt4.vol = (dx*dy*dz)/8.0;
				prt4.lp  = (meanL/ppc)*0.5;
				prt4.matid = matid;
				pid++;
				prt4.id = pid; 
				ParticleVector.push_back(prt4);
				
				iz += (dz*0.5);

			}
			while(iz<ztop);
		}

		cout<<"it were read "<<DemBoxVector.size()<<" DEM boxes...\n";
		DemBoxVector.clear();
	}

	void CreateMPMmodel()
	{
		CreateMpmModelByDEM();
	}
}
