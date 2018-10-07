/*
 * main.cpp
 *
 *  Created on: Oct 6, 2018
 *      Author: fabricio
 */

#include <iostream>
#include <string>
#include <vector>
#include "io.h"
#include "model.h"

int main(int argc, char **argv)
{
    IO::InitFileName(argc, argv);

    IO::ReadInputFile();

    Model::CreateMPMmodel();

    IO::WriteOutputFile();

    return 0;
}
