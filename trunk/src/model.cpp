/*
 * model.cpp
 *
 *  Created on: Oct 6, 2018
 *      Author: fabricio
 */


#include <iostream>
#include <vector>

#include <model.h>
#include <utils.h>

namespace Model
{
	using namespace std;

	static vector<DemBox> DemBoxVector;
	vector<DemBox>& GetDemVector(){return DemBoxVector;}

	static vector<Particle> ParticleVector;
	vector<Particle>& GetParticleVector(){return ParticleVector;}

	static vector<mshPoint> meshPointsVector;
	vector<mshPoint>& GetmshPointsVector(){return meshPointsVector;}

	static vector<mshElement> mshElementsVector;
	vector<mshElement>& GetmshElementsVector(){return mshElementsVector;}

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

	static void CreateMpmModelBySMESH()
	{ 
		if(mshElementsVector.size()==0) { return; }

		ParticleVector.clear();
		int pid = 0;
		double ppc = 2.0;
		double nparticles = 5;

		std::vector<Vector3> coord;
		coord.push_back(Vector3(0.0));
		coord.push_back(Vector3(0.0));
		coord.push_back(Vector3(0.0));
		coord.push_back(Vector3(0.0));
		
		double min_x(0.0),max_x(0.0),min_y(0.0),max_y(0.0),min_z(0.0),max_z(0.0);


		for( size_t i=0; i < mshElementsVector.size(); i++ )
		{
			// get global cords of nodes
			for( size_t j=0; j < mshElementsVector[i].points.size(); j++ )
			{
				coord[j].x = meshPointsVector[mshElementsVector[i].points[j]-1].pos.x;
				coord[j].y = meshPointsVector[mshElementsVector[i].points[j]-1].pos.y;
				coord[j].z = meshPointsVector[mshElementsVector[i].points[j]-1].pos.z;
			
				if(j==0){
					min_x = max_x = coord[j].x;
					min_y = max_y = coord[j].y;
					min_z = max_z = coord[j].z;
				}

				if(coord[j].x < min_x) min_x = coord[j].x;
				if(coord[j].x > max_x) max_x = coord[j].x;
				if(coord[j].y < min_y) min_y = coord[j].y;
				if(coord[j].y > max_y) max_y = coord[j].y;
				if(coord[j].z < min_z) min_z = coord[j].z;
				if(coord[j].z > max_z) max_z = coord[j].z;

				//std::cout<<"elem:"<<i<<"\n";
			}

			// volume of the element and particles
			double vol  = Utils::volumeTetrahedonGet(coord[0],coord[1],coord[2],coord[3]);
			double pvol = vol/nparticles;

			// particle lenght for load and pressure conditions
			Vector3 L;
			L.x = 0.5 * (max_x - min_x);
			L.y = 0.5 * (max_y - min_y);
			L.z = 0.5 * (max_z - min_z);

			double meanL = (L.x+L.y+L.z)/3.0;

			// area coords
			double L1[5] = {0.1,0.1,0.1,0.7,0.25};
			double L2[5] = {0.1,0.1,0.7,0.1,0.25};
			double L3[5] = {0.1,0.7,0.1,0.1,0.25};
			double L4[5] = {0.7,0.1,0.1,0.1,0.25};

			// element node coords
			double x21 = coord[0].x;
			double x22 = coord[1].x;
			double x23 = coord[2].x;
			double x24 = coord[3].x;

			double x31 = coord[0].y;
			double x32 = coord[1].y;
			double x33 = coord[2].y;
			double x34 = coord[3].y;    

			double x41 = coord[0].z;
			double x42 = coord[1].z;
			double x43 = coord[2].z;
			double x44 = coord[3].z;
			
			// mat id - TODO
			int matid = 1;

			Particle prti;
			
			for (int ip = 0;ip<nparticles;++ip)
			{
				prti.pos.x = (L1[ip]*x21 + L2[ip]*x22 + L3[ip]*x23 + L4[ip]*x24);
				prti.pos.y = (L1[ip]*x31 + L2[ip]*x32 + L3[ip]*x33 + L4[ip]*x34);
				prti.pos.z = (L1[ip]*x41 + L2[ip]*x42 + L3[ip]*x43 + L4[ip]*x44);
				prti.vol   = pvol;
				prti.lp    = (meanL/ppc)*0.5;
				prti.matid = matid;
				pid++;
				prti.id = pid; 
				ParticleVector.push_back(prti);
			}
		}
	}

	void CreateMPMmodel()
	{
		CreateMpmModelByDEM();
		CreateMpmModelBySMESH();
	}
}
