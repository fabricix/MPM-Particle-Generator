/*
 * model.cpp
 *
 *  Created on: Oct 6, 2018
 *      Author: fabricio fernandez <fabricio.hmf@gmail.com>
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

	static vector<vector<HorizontPoint> > HorizontPointVector;
	vector<vector<HorizontPoint> >& GetHorizontVector(){return HorizontPointVector;}

	static mshGrid Grid;
	mshGrid& GetmshGrid(){return Grid;}

	static std::vector<bool> cellid;
	static std::vector<int> materialidcell;

	static int nHorizonts;
	void SetHorizonNumber(int nhor){nHorizonts = nhor;}
	int  GetHorizonNumber( ){return nHorizonts;}

	static BoundingBoxDem dembndbox;
	void BoundingBoxDemSet(Vector3 pmin,Vector3 pmax){dembndbox.min=pmin; dembndbox.max=pmax;}
	void BoundingBoxDemGet(Vector3& lmin,Vector3& lmax){lmin=dembndbox.min; lmax=dembndbox.max; }

	bool SEPARATE_FILES = false;

	void setWriteParticlesSeparateFiles(bool value){ SEPARATE_FILES = value; }
	bool getWriteParticlesSeparateFiles(){return SEPARATE_FILES;}

	static void mapHorizonPoints2DemPoints()
	{	
		for (size_t ih = 0; ih < HorizontPointVector.size(); ++ih)
		{
			for (size_t i = 0; i < HorizontPointVector.at(ih).size(); ++i)
			{
				int dem_index = -1;
				double mindist = 1e300;	
				Vector3 ihpos = HorizontPointVector.at(ih).at(i).pos;

				for (size_t j = 0; j < DemBoxVector.size(); ++j)
				{
					Vector3 pos = DemBoxVector.at(j).pos;
					double dist = std::sqrt(std::pow((pos.x-ihpos.x),2.0)+std::pow((pos.y-ihpos.y),2.0));

					if (dist<=mindist)
					{
						mindist = dist;
						dem_index = j;
					}
				}

				DemBoxVector.at(dem_index).hpnt.push_back(HorizontPointVector.at(ih).at(i));
			}
		}
	}

	static int inCell(Vector3 pos)
	{
		int i = floor((pos.x-Grid.regionBegin.x)/Grid.dx+Grid.Ng);
      	int j = floor((pos.y-Grid.regionBegin.y)/Grid.dy+Grid.Ng);
      	int k = floor((pos.z-Grid.regionBegin.z)/Grid.dz+Grid.Ng);

      	if(i>=0 && i<Grid.I && j>=0 && j<Grid.J && k>=0 && k<Grid.K)
      	{
        	return ((j*Grid.I+i)+(Grid.I*Grid.J*k));
      	}
      	else
      	{
			std::cout<<"particle in cell not found... line "<<__LINE__<<", in file"<< __FILE__<<"\n";

		    std::cout<<"xp = "<<pos.x<<"\n";
			std::cout<<"yp = "<<pos.y<<"\n";
			std::cout<<"zp = "<<pos.z<<"\n";
			
         	return(-1); 
      	} 
	}

	static void CreateMpmModelByDEM()
	{ 
		if(DemBoxVector.size()==0) { return; }
		if(HorizontPointVector.size()!=0) { return; }
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
				prt1.lp  = (meanL/ppc);
				prt1.matid = matid;
				pid++;
				prt1.id = pid; 
				ParticleVector.push_back(prt1);

				prt2.pos.x = ibox.pos.x+(dx*0.25);
				prt2.pos.y = ibox.pos.y-(dy*0.25);
				prt2.pos.z = iz;
			    prt2.vol = (dx*dy*dz)/8.0;
				prt2.lp  = (meanL/ppc);
				prt2.matid = matid;
				pid++;
				prt2.id = pid; 
				ParticleVector.push_back(prt2);

				prt3.pos.x = ibox.pos.x+(dx*0.25);
				prt3.pos.y = ibox.pos.y+(dy*0.25);
				prt3.pos.z = iz;
			    prt3.vol = (dx*dy*dz)/8.0;
				prt3.lp  = (meanL/ppc);
				prt3.matid = matid;
				pid++;
				prt3.id = pid; 
				ParticleVector.push_back(prt3);

				prt4.pos.x = ibox.pos.x-(dx*0.25);
				prt4.pos.y = ibox.pos.y+(dy*0.25);
				prt4.pos.z = iz;
				prt4.vol = (dx*dy*dz)/8.0;
				prt4.lp  = (meanL/ppc);
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

	static void CreateMpmModelByDEMHorizonts()
	{ 
		if(DemBoxVector.size()==0) { return; }
		if(HorizontPointVector.size()==0) { return; }

		// map the horizont points to dem points
		mapHorizonPoints2DemPoints();

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
				for (size_t ih = 0; ih < DemBoxVector.at(i).hpnt.size(); ++ih)
				{
					if (iz > DemBoxVector.at(i).hpnt.at(ih).pos.z)
					{
						matid = DemBoxVector.at(i).hpnt.at(ih).matid;
					}
				}

				prt1.pos.x = ibox.pos.x-(dx*0.25);
				prt1.pos.y = ibox.pos.y-(dy*0.25);
				prt1.pos.z = iz;
				prt1.vol = (dx*dy*dz)/8.0;
				prt1.lp  = (meanL/ppc);
				prt1.matid = matid;
				pid++;
				prt1.id = pid; 
				ParticleVector.push_back(prt1);

				prt2.pos.x = ibox.pos.x+(dx*0.25);
				prt2.pos.y = ibox.pos.y-(dy*0.25);
				prt2.pos.z = iz;
			    prt2.vol = (dx*dy*dz)/8.0;
				prt2.lp  = (meanL/ppc);
				prt2.matid = matid;
				pid++;
				prt2.id = pid; 
				ParticleVector.push_back(prt2);

				prt3.pos.x = ibox.pos.x+(dx*0.25);
				prt3.pos.y = ibox.pos.y+(dy*0.25);
				prt3.pos.z = iz;
			    prt3.vol = (dx*dy*dz)/8.0;
				prt3.lp  = (meanL/ppc);
				prt3.matid = matid;
				pid++;
				prt3.id = pid; 
				ParticleVector.push_back(prt3);

				prt4.pos.x = ibox.pos.x-(dx*0.25);
				prt4.pos.y = ibox.pos.y+(dy*0.25);
				prt4.pos.z = iz;
				prt4.vol = (dx*dy*dz)/8.0;
				prt4.lp  = (meanL/ppc);
				prt4.matid = matid;
				pid++;
				prt4.id = pid; 
				ParticleVector.push_back(prt4);
				
				iz += (dz*0.5);

			}
			while(iz<ztop);
		}

		cout<<"it were read "<<DemBoxVector.size()<<" DEM boxes...\n";
		cout<<"it were read "<<HorizontPointVector.size()<<" Horizons...\n";
		DemBoxVector.clear();
	}
	
	static void MpmModel_SMESH_Element_Mapping()
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

			// area cords
			double L1[5] = {0.1,0.1,0.1,0.7,0.25};
			double L2[5] = {0.1,0.1,0.7,0.1,0.25};
			double L3[5] = {0.1,0.7,0.1,0.1,0.25};
			double L4[5] = {0.7,0.1,0.1,0.1,0.25};

			// element node cords
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
				prti.lp    = (meanL/ppc);
				prti.matid = matid;
				pid++;
				prti.id = pid; 
				ParticleVector.push_back(prti);
			}
		}
	}

	static void AddParticlesInActiveCellGrid(std::vector<bool> cellid,std::vector<int>materialcell)
	{
		ParticleVector.clear();
		int pid = 0;
		double ppc = 2.0;
		
		// cell dimension
		double dx = Grid.dx;
		double dy = Grid.dy;
		double dz = Grid.dz;

		double meanL = (dx+dy+dz)/3.0;

		// particle volume
		double pvol = (dx*dy*dz)/8.0;

		// length for pressure calculations
		Vector3 L;
		L.x = 0.5*dx;
		L.y = 0.5*dy;
		L.z = 0.5*dz;

		// add 8 particle in each active cell grid
		for (size_t i = 0; i < cellid.size(); ++i)
		{
			if (cellid.at(i)==true)
			{
				int idn1 = i;
				int idn2 = idn1+1;
				int idn3 = idn1+Grid.I;
				int idn4 = idn3+1;
				int idn5 = idn1 + (Grid.I*Grid.J);
				int idn6 = idn2 + (Grid.I*Grid.J);
				int idn7 = idn3 + (Grid.I*Grid.J);
				int idn8 = idn4 + (Grid.I*Grid.J);

				int matid = materialcell.at(i);
				
				Vector3 pt1 = Grid.gnodes[idn1];
				Vector3 pt2 = Grid.gnodes[idn2];
				Vector3 pt3 = Grid.gnodes[idn3];
				Vector3 pt4 = Grid.gnodes[idn4];

				Vector3 pt5 = Grid.gnodes[idn5];
				Vector3 pt6 = Grid.gnodes[idn6];
				Vector3 pt7 = Grid.gnodes[idn7];
				Vector3 pt8 = Grid.gnodes[idn8];

				pt1.x += (dx*0.25);
				pt1.y += (dy*0.25);
				pt1.z += (dz*0.25);

				pt2.x -= (dx*0.25);
				pt2.y += (dy*0.25);
				pt2.z += (dz*0.25);

				pt3.x += (dx*0.25);
				pt3.y -= (dy*0.25);
				pt3.z += (dz*0.25);

				pt4.x -= (dx*0.25);
				pt4.y -= (dy*0.25);
				pt4.z += (dz*0.25);

				pt5.x += (dx*0.25);
				pt5.y += (dy*0.25);
				pt5.z -= (dz*0.25);

				pt6.x -= (dx*0.25);
				pt6.y += (dy*0.25);
				pt6.z -= (dz*0.25);

				pt7.x += (dx*0.25);
				pt7.y -= (dy*0.25);
				pt7.z -= (dz*0.25);

				pt8.x -= (dx*0.25);
				pt8.y -= (dy*0.25);
				pt8.z -= (dz*0.25);
				
				Particle prt1;
				prt1.pos   = pt1;
				prt1.vol   = pvol;
				prt1.lp    = (meanL/ppc);
				prt1.matid = matid;
				pid++;
				prt1.id = pid; 
				ParticleVector.push_back(prt1);

				Particle prt2;
				prt2.pos   = pt2;
				prt2.vol   = pvol;
				prt2.lp    = (meanL/ppc);
				prt2.matid = matid;
				pid++;
				prt2.id = pid; 
				ParticleVector.push_back(prt2);

				Particle prt3;
				prt3.pos   = pt3;
				prt3.vol   = pvol;
				prt3.lp    = (meanL/ppc);
				prt3.matid = matid;
				pid++;
				prt3.id = pid; 
				ParticleVector.push_back(prt3);

				Particle prt4;
				prt4.pos   = pt4;
				prt4.vol   = pvol;
				prt4.lp    = (meanL/ppc);
				prt4.matid = matid;
				pid++;
				prt4.id = pid; 
				ParticleVector.push_back(prt4);

				Particle prt5;
				prt5.pos   = pt5;
				prt5.vol   = pvol;
				prt5.lp    = (meanL/ppc);
				prt5.matid = matid;
				pid++;
				prt5.id = pid; 
				ParticleVector.push_back(prt5);

				Particle prt6;
				prt6.pos   = pt6;
				prt6.vol   = pvol;
				prt6.lp    = (meanL/ppc);
				prt6.matid = matid;
				pid++;
				prt6.id = pid; 
				ParticleVector.push_back(prt6);

				Particle prt7;
				prt7.pos   = pt7;
				prt7.vol   = pvol;
				prt7.lp    = (meanL/ppc);
				prt7.matid = matid;
				pid++;
				prt7.id = pid; 
				ParticleVector.push_back(prt7);

				Particle prt8;
				prt8.pos   = pt8;
				prt8.vol   = pvol;
				prt8.lp    = (meanL/ppc);
				prt8.matid = matid;
				pid++;
				prt8.id = pid; 
				ParticleVector.push_back(prt8);
			}
		}

		std::cout<<"was created "<<ParticleVector.size()<<" particles by grid mapping..\n";
	}

	static void MpmModel_SMESH_Grid_Mapping()
	{ 	
		cellid.clear();
		cellid.resize(Grid.Nnodes,false);
		
		materialidcell.clear();
		materialidcell.resize(Grid.Nnodes,0);

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
			}

			// volume of the element and particles
			double vol  = Utils::volumeTetrahedonGet(coord[0],coord[1],coord[2],coord[3]);
			double pvol = vol/nparticles;

			// particle lenght for load and pressure conditions
			Vector3 L;
			L.x = 0.5 * (max_x - min_x);
			L.y = 0.5 * (max_y - min_y);
			L.z = 0.5 * (max_z - min_z);

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
				Vector3 pt(0,0,0);
				pt.x = (L1[ip]*x21 + L2[ip]*x22 + L3[ip]*x23 + L4[ip]*x24);
				pt.y = (L1[ip]*x31 + L2[ip]*x32 + L3[ip]*x33 + L4[ip]*x34);
				pt.z = (L1[ip]*x41 + L2[ip]*x42 + L3[ip]*x43 + L4[ip]*x44);

				cellid[inCell(pt)] = true;
				materialidcell[inCell(pt)] = mshElementsVector[i].mat;	
			}
		}

		AddParticlesInActiveCellGrid(cellid,materialidcell);
	}

	static void CreateMpmModelBySMESH()
	{ 
		//MpmModel_SMESH_Element_Mapping();
		MpmModel_SMESH_Grid_Mapping();
	}

	void initGrid( Vector3 limmin,Vector3 limmax, Vector3 celldim)
	{
		Grid.regionBegin = limmin;
		Grid.regionEnd   = limmax;
		Grid.dx = celldim.x; //(limmax.x-limmin.x)/double(Ncell.x);
		Grid.dy = celldim.y; //(limmax.y-limmin.y)/double(Ncell.y);
		Grid.dz = celldim.z; //(limmax.z-limmin.z)/double(Ncell.z);
		Grid.Ng = 2.0;

		Vector3 Ncell;

		Ncell.x = (limmax.x-limmin.x)/Grid.dx;
		Ncell.y = (limmax.y-limmin.y)/Grid.dy;
		Ncell.z = (limmax.z-limmin.z)/Grid.dz;

		Grid.I = (Ncell.x+2*Grid.Ng+1);
        Grid.J = (Ncell.y+2*Grid.Ng+1);
        Grid.K = (Ncell.z+2*Grid.Ng+1);

        Grid.gnodes.clear();
        Grid.gnodes.resize((Grid.I*Grid.J*Grid.K),Vector3(0.0));
        Grid.Nnodes = (int) Grid.gnodes.size();

        // Grid3D
		for (int k=0; k < Grid.K; ++k)
		{
		    double  z = (k-Grid.Ng)*Grid.dz+limmin.z;
		    for(int j=0; j<Grid.J; ++j)
		    {
		        double y = (j-Grid.Ng)*Grid.dy+limmin.y;
		        for(int i=0; i<Grid.I; ++i)
		        {
		            double x = (i-Grid.Ng)*Grid.dx+limmin.x;
		            int n=(j*Grid.I+i)+(Grid.J*Grid.I*k);

		            Grid.gnodes[n]=Vector3(x,y,z);
		        }
		    }
		}
		std::cout<<"was created a grid, with "<< Grid.gnodes.size()<<" nodes...\n";
	}

	int NnodesGet()
	{
		return Grid.Nnodes;
	}

	std::vector<Vector3> GridNodesGet()
	{
		return Grid.gnodes;
	}

	Vector3 SpacingbyDirectionGet()
	{
		return Vector3(Grid.dx,Grid.dy,Grid.dz);
	}
	
	Vector3 NnodesbyDirectionGet()
	{
		return Vector3(Grid.I,Grid.J,Grid.K);
	}

	Vector3 RegionBeginGet()
	{
		return Grid.regionBegin;
	}

	Vector3 MinGridCoordsGet()
	{
		double minx(0.0), miny (0.0), minz (0.0);
		Vector3 pos(0.0);

		for (size_t i = 0; i < Grid.gnodes.size(); ++i)
		{
			pos = Grid.gnodes.at(i);

			if (pos.x<minx){ minx = pos.x ;}
			if (pos.y<miny){ miny = pos.y ;}
			if (pos.z<minz){ minz = pos.z ;}
		}

		return Vector3(minx,miny,minz);
	}

	void CreateMPMmodel()
	{
		CreateMpmModelByDEM();
		CreateMpmModelByDEMHorizonts();
		CreateMpmModelBySMESH();
	}
}
