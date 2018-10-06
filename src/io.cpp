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
	static std::string output_fname;
	static std::string path_out_fname;
	static std::string path_input_fname;
	static std::string path;

	std::ifstream inputfile;
	std::ifstream outputfile;

	static std::string TimeDataGet(){

   		time_t now = time(0);
   		char* dt = ctime(&now);
   		return std::string(dt);
	}

	void InitFileName(int argc, char **argv){

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
		output_fname = "//mpm.part";
		path_out_fname = path + output_fname;

	}

	 void ReadInputFile(){

		 inputfile.open (input_fname.c_str(), std::ifstream::in);
		 std::string line;

		 while (std::getline(inputfile,line))
		 {
			 std::size_t found = line.find("DEM.TO.MPM");
			 if (found!=std::string::npos)
			 {
				 std::string auxline;
				 std::getline(inputfile,auxline);
				 int nbox = 0;
				 std::sscanf(auxline.c_str(), "%d", &nbox);

				 std::getline(inputfile,auxline);

				 std::vector<Model::DemBox> demvector = Model::GetDemVector();
				 demvector.clear();
				  if(nbox==0){
					  std::cout<<"ERROR: number of DEM box is = "<<nbox<<"\n";
					  return;
				  }

				  for( int i = 0; i < nbox; i++ )
				  {
					  Model::DemBox ibox;
				      std::sscanf( auxline.c_str(), "%lf%lf%lf%lf%lf%lf%d", &ibox.pos.x,&ibox.pos.y,&ibox.pos.z,&ibox.zbase,&ibox.dx,&ibox.dy,&ibox.matid);
				      demvector.push_back(ibox);
				  }
			  }

		 }
#if 0
		 char str[256];
		 char c;
		   while (inputfile.get(c))          // loop getting single characters
		     std::cout << c;
#endif

		 inputfile.close();
	 }

	 void WriteOutputFile(){



	 }
}
