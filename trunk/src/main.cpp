/*
 * main.cpp
 *
 *  Created on: Oct 6, 2018
 *      Author: fabricio
 */

#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv)
{

	// name of file to process
	std::string file_name;
	if (argc > 1) {
	  file_name = argv[1];
	  std::cout<<file_name<<"\n";
	}
	else {
		std::cout<<"ERROR: please insert the file to process...\n";
		return 0;
    }

	// set the current path


	// process the file


	// write the particle file

	//exit
  return 0;
}
