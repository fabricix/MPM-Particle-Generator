/*
 * main.cpp
 *
 *  Created on: Oct 6, 2018
 *      Author: fabricio
 */

// C/C++
#include <iostream>
#include <string>
#include <vector>

// local
#include "io.h"
#include "model.h"

int main(int argc, char **argv)
{	
    IO::ReadInputFile(argc, argv);

    Model::CreateMPMmodel();

    IO::WriteOutputFile();

    return 0;
}
