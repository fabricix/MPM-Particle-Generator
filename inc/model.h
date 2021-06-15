/*
 * model.h
 *
 *  Created on: Oct 6, 2018
 *      Author: fabricio
 */

#ifndef INC_MODEL_H_
#define INC_MODEL_H_

#include "matrix.h"

namespace Model
{
	enum elemtypes { none, tet4 };

	struct HorizontPoint
	{
		Vector3 pos;
		int matid;
		HorizontPoint():pos(0.0),matid(0){}
	};
	
	struct DemBox
	{
		Vector3 pos;
		double zbase;
		double dx,dy,dz;
		int matid;
		std::vector<HorizontPoint> hpnt;
		DemBox():pos(0.0),zbase(0.0),dx(0.0),dy(0.0),dz(0.0),matid(0){}
	};

	struct Particle
	{
		int id;
		Vector3 pos;
		double vol;
		double lp;
		int matid;
		Particle():id(0),pos(0.0),vol(0.0),lp(0.0),matid(0){}
	};
	
	struct mshPoint
	{
		int id;
		Vector3 pos;
		int bndry;
	};

	struct mshElement
	{
		int id;
		int bndry;
		int mat;
		std::vector<int> points;
	};

	struct mshGrid
	{
		Vector3 regionBegin;       		// minimum x,y,z of region
		Vector3 regionEnd;         		// maximum x,y,z of region
		double dx;                 		// x-dir cell size
		double dy;                 		// y-dir cell size
		double dz;                 		// z-dir cell size
		int Ng;  	              		// number of rings of nodes outside region
		int I;                     		// number of grid nodes in x direction
		int J;                     		// number of grid nodes in y direction
		int K;                     		// number of grid nodes in z direction
		int Nnodes;
		std::vector<Vector3> gnodes;	// grid nodes
	};

	struct BoundingBoxDem
	{
		Vector3 min;
		Vector3 max;
		BoundingBoxDem():min(0.0),max(0){}
	};

	std::vector<DemBox>& GetDemVector();
	std::vector<Particle>& GetParticleVector();
	std::vector<mshPoint>& GetmshPointsVector();
	std::vector<mshElement>& GetmshElementsVector();
	std::vector<std::vector<HorizontPoint> >& GetHorizontVector();
	
	mshGrid& GetmshGrid();
	int NnodesGet();
	std::vector<Vector3> GridNodesGet();
	Vector3 SpacingbyDirectionGet();
	Vector3 NnodesbyDirectionGet();
	Vector3 RegionBeginGet();
	Vector3 MinGridCoordsGet();

	void CreateMPMmodel();
	void initGrid(Vector3 limmin,Vector3 limmax, Vector3 celldim);
	
	// horizonts	
	void SetHorizonNumber(int nhor);
	int  GetHorizonNumber( );

	void BoundingBoxDemSet(Vector3 pmin,Vector3 pmax);
	void BoundingBoxDemGet(Vector3& lmin,Vector3& lmax);

	void setWriteParticlesSeparateFiles(bool value);
	bool getWriteParticlesSeparateFiles();
}

#endif /* INC_MODEL_H_ */
