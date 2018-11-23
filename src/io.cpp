/*
 * io.cpp
 *
 *  Created on: Oct 6, 2018
 *      Author: fabricio
 */

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>

#include "io.h"
#include "model.h"

namespace IO {

	static std::string input_fname;
	static std::string out_fname;
	static std::string path_out_fname;
	static std::string path;
	static std::string line;

	std::ifstream inputfile;
	std::ofstream outfile;

	//** statics

	static void read_dem2mpm()
	{
		std::string auxline;
		std::getline(inputfile,auxline);
		
		// get lines number
		int nbox = 0;
		std::sscanf(auxline.c_str(),"%d",&nbox);
		
		// read lines
		std::getline(inputfile,auxline);

		std::vector<Model::DemBox>& demvector = Model::GetDemVector();

		if(nbox==0){
			std::cout<<"ERROR: number of DEM box is = "<<nbox<<"\n";
			return;
		}

		for( int i = 0; i < nbox; i++ )
		{
			Model::DemBox ibox;
			std::sscanf(auxline.c_str(),"%lf%lf%lf%lf%lf%lf%lf%d", &ibox.pos.x,&ibox.pos.y,&ibox.pos.z,&ibox.zbase,&ibox.dx,&ibox.dy,&ibox.dz,&ibox.matid);
			demvector.push_back(ibox);
			std::getline(inputfile,auxline);
		}
	}


	static void read_smesh_points()
	{
		std::getline(inputfile,line);
		
		// get lines number
		int npnts;	// # of points
		int dim;	// dimension
		int attr;	// attributes
		int bndry;	// # of boundary markers

		std::sscanf(line.c_str(),"%d %d %d %d",&npnts,&dim,&attr,&bndry);

		std::vector<Model::mshPoint>& pntsvector = Model::GetmshPointsVector();

		if(npnts==0){
			std::cout<<"ERROR: number of points is = "<<npnts<<"\n";
			return;
		}

		for( int i = 0; i < npnts; i++ )
		{
			std::getline(inputfile,line);
			Model::mshPoint ipnt;
			std::sscanf(line.c_str(),"%d %lf %lf %lf %d", &ipnt.id,&ipnt.pos.x,&ipnt.pos.y,&ipnt.pos.z,&ipnt.bndry);
			pntsvector.push_back(ipnt);
		}

		std::cout<<"reading "<<pntsvector.size()<<" points...\n";
	}

	static void read_smesh_elements()
	{
		std::getline(inputfile,line);
		
		// get lines number
		int nelem;		// # of elements
		int nodes;	    // # nodes or corners
		int bndry;	    // # of boundary markers

		std::sscanf(line.c_str(),"%d %d",&nelem, &bndry);

		std::vector<Model::mshElement>& elemsvector = Model::GetmshElementsVector();

		if(nelem==0){
			std::cout<<"ERROR: number of elements is = "<<nelem<<"\n";
			return;
		}

		for( int i = 0; i < nelem; i++ )
		{
			std::getline(inputfile,line);
			Model::mshElement iele;
			int nnodes,n1,n2,n3,n4,bndry;
			
			std::sscanf(line.c_str(),"%d %d %d %d %d %d",&nnodes,&n1,&n2,&n3,&n4,&bndry);
			
			iele.points.clear();
			iele.points.push_back(n1);
			iele.points.push_back(n2);
			iele.points.push_back(n3);
			iele.points.push_back(n4);
			
			iele.bndry = bndry;

			elemsvector.push_back(iele);
		}

		std::cout<<"reading "<<elemsvector.size()<<" elements...\n";
	}

	static void read_smesh2mpm()
	{	
		std::cout<<"reading SMESH file...\n";
		
		while (std::getline(inputfile,line))
		{
			// points
			std::size_t found1 = line.find("SMESH.POINTS");
			if (found1!=std::string::npos)
			{
				std::cout<<"reading points in SMESH file...\n";
				read_smesh_points();
			}

			// elements
			std::size_t found2 = line.find("SMESH.ELEMENTS");
			if (found2!=std::string::npos)
			{	
				std::cout<<"reading elements in SMESH file...\n";
				read_smesh_elements();
			}
		}
	}

	static void init_file_name(int argc, char **argv)
	{
		// path set
		if(argc>0){
			std::string argv_str(argv[1]);
			std::string base = argv_str.substr(0, argv_str.find_last_of("/"));
			path = base;
		}

		// name of file to process
		if (argc > 1) {
			input_fname = argv[1];
		}
		else {
			std::cout<<"ERROR: please insert the file to process...\n";
			return;
	    }

		// name and path of output file
		out_fname = "/mpm.part";
		path_out_fname = path + out_fname;

	}


	//** globals

	void ReadInputFile(int argc, char **argv)
	{
		init_file_name(argc, argv);

		inputfile.open (input_fname.c_str(), std::ifstream::in);
		//std::string line;

		while (std::getline(inputfile,line))
		{
			// DEM to MPM
			std::size_t found1 = line.find("DEM.TO.MPM");
			if (found1!=std::string::npos)
			{
				read_dem2mpm();
			}

			// SMESH to MPM
			std::size_t found2 = line.find("SMESH");
			if (found2!=std::string::npos)
			{
				read_smesh2mpm();
			}

		}
		// close the file
		if(inputfile.is_open()) inputfile.close();
	}

	void WriteOutputFile()
	{
		using namespace std;
		outfile.open(path_out_fname.c_str(), fstream::out);
		
		// test
		vector<Model::Particle>& parvec = Model::GetParticleVector();

		// header for plot
		outfile<<"id x y z vol lp matid\n";
		outfile<<"%PARTICLES\n";
		outfile<<parvec.size()<<"\n";

		for (size_t i = 0; i < parvec.size(); ++i)
		{
			outfile<<parvec.at(i).id<<" ";
			outfile<<parvec.at(i).pos.x<<" ";
			outfile<<parvec.at(i).pos.y<<" ";
			outfile<<parvec.at(i).pos.z<<" ";
			outfile<<parvec.at(i).vol<<" ";
			outfile<<parvec.at(i).lp<<" ";
			outfile<<parvec.at(i).matid<<"\n";
		}

		cout<<"writing "<<parvec.size()<<" particles...\n";
		cout<<"please see the output file, "<<out_fname.c_str()<<"...\n";

		// close the file
		if(outfile.is_open()) outfile.close();
	}
}
