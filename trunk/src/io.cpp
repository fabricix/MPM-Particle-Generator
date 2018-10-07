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

	//** globals

	void InitFileName(int argc, char **argv)
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
		out_fname = "mpm.part";
		path_out_fname = path + out_fname;

	}

	void ReadInputFile()
	{
		inputfile.open (input_fname.c_str(), std::ifstream::in);
		std::string line;

		while (std::getline(inputfile,line))
		{
			// DEM to MPM
			std::size_t found = line.find("DEM.TO.MPM");
			if (found!=std::string::npos)
			{
				read_dem2mpm();
			}

		}
		// close the file
		if(inputfile.is_open()) inputfile.close();
	}

	void WriteOutputFile()
	{
		using namespace std;
		outfile.open(out_fname.c_str(), fstream::out);
		
		// test
		vector<Model::Particle>& parvec = Model::GetParticleVector();

		// header for plot
		outfile<<"id x y z vol lp matid\n";

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

		cout<<"it was writing "<<parvec.size()<<" particles...\n";
		cout<<"please see the output file, "<<out_fname.c_str()<<"...\n";

		// close the file
		if(outfile.is_open()) outfile.close();
	}
}
