/*
 * paraview.cpp
 *
 *  Created on: 23 de nov de 2018
 *      Author: fabricio
 */

#include <string>
#include <fstream>
#include <vector>
#include <stdint.h>

#include "paraview.h"
#include "io.h"
#include "model.h"

namespace ParaView
{
	static std::string Edian = "";

	static void defineEdian()
	{
	   int16_t i = 1;
	   int8_t *p = (int8_t *) &i;
	   if (p[0] == 1) { Edian = "LittleEndian"; }
	   else { Edian = "BigEndian"; }
	}

	void writegrid( )
	{
		if(Edian==""){ defineEdian(); }

		// open the file
		std::ofstream gridfile;
		std::string fname = IO::pathGet() + "grid" + ".vtu";
		gridfile.open(fname.c_str());
		gridfile.precision(4);

		std::cout<<"writing grid file name: "+fname+"\n";
		// write results
		gridfile <<"<?xml version=\"1.0\"?>\n";
		gridfile <<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\""<<Edian.c_str()<<"\">\n";
		gridfile <<"<UnstructuredGrid>\n";

		// points
		int ngpoints = Model::NnodesGet();
		gridfile <<std::fixed<<"\t<Piece NumberOfPoints=\""<<ngpoints<<"\" NumberOfCells=\""<<ngpoints<<"\">\n";
		gridfile <<"\t\t<Points>\n";

		// position
		std::vector<Vector3> gnodes = Model::GridNodesGet();
		gridfile <<"\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
		for (int i = 0; i < ngpoints; ++i){
		  gridfile<<std::scientific<<"\t\t\t\t"<<gnodes[i].x<<"\t"<<gnodes[i].y<<"\t"<<gnodes[i].z<<"\n";
		}
		gridfile <<"\t\t\t</DataArray>\n";

		gridfile <<"\t\t</Points>\n";
	
		// scalar values 
		gridfile <<"\t\t\t<PointData Scalars=\"scalars\">\n";
		
		// user fixed points
		gridfile <<"\t\t\t\t<DataArray type=\"Float32\" Name=\"Fixed\" Format=\"ascii\">\n";
		for (int i = 0; i < ngpoints; ++i){
			gridfile<<std::scientific<<"\t\t\t\t\t"<<0.0<<"\n";
		}
		gridfile <<"\t\t\t\t</DataArray>\n";

		// id
		gridfile <<"\t\t\t\t<DataArray type=\"Float32\" Name=\"Id\" Format=\"ascii\">\n";
		for (int i = 0; i < ngpoints; ++i){
			gridfile<<std::scientific<<"\t\t\t\t\t"<<i<<"\n";
		}
		gridfile <<"\t\t\t\t</DataArray>\n";

		gridfile <<"\t\t\t</PointData>\n";

		// cells -->
		gridfile <<"\t\t\t<Cells>\n";
		gridfile <<"\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";
		for (int i = 0; i < ngpoints; ++i){
			gridfile<<std::scientific<<"\t\t\t\t\t"<<i<<"\n";
		}
		gridfile <<"\t\t\t\t</DataArray>\n";

		gridfile <<"\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n";
		for (int i = 0; i < ngpoints; ++i){
			gridfile <<"\t\t\t\t\t"<<i+1<<"\n";
		}
		gridfile <<"\t\t\t\t</DataArray>\n";

		gridfile <<"\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" Format=\"ascii\">\n";
		for (int i = 0; i < ngpoints; ++i){
			gridfile <<"\t\t\t\t\t"<<1<<"\n";
		}
		gridfile <<"\t\t\t\t</DataArray>\n";
		gridfile <<"\t\t\t</Cells>\n";
		//<-- cells

		gridfile <<"\t</Piece>\n";
		gridfile <<"</UnstructuredGrid>\n";
		gridfile <<"</VTKFile>\n";

		// close the file
		gridfile.close();
	}

	void writegridcell( )
	{
		if(Edian==""){ defineEdian(); }

		// open the file
		std::ofstream gridfile;
		std::string fname = IO::pathGet() + "grid-cells" + ".vtk";
		gridfile.open(fname.c_str());
		gridfile.precision(4);

		Vector3 spacing = Model::SpacingbyDirectionGet();
		Vector3 ncelldir = Model::NnodesbyDirectionGet();
		Vector3 origin = Model::MinGridCoordsGet();

		Vector3 nbydir = Model::NnodesbyDirectionGet();
		int ngpoints = (nbydir.x-1)*(nbydir.y-1)*(nbydir.z-1);

		gridfile<<"# vtk DataFile Version 4.0\n";
		gridfile<<"vtk output\n";
		gridfile<<"ASCII\n";
		gridfile<<"DATASET STRUCTURED_POINTS\n";
		gridfile<<"DIMENSIONS "<<ncelldir.x<<" "<<ncelldir.y<<" "<<ncelldir.z<<"\n";
		gridfile<<"SPACING "<<spacing.x<<" "<<spacing.y<<" "<<spacing.z<<"\n";
		gridfile<<"ORIGIN  "<<origin.x<<" "<<origin.y<<" "<<origin.z<<" \n";
		gridfile<<"CELL_DATA  "<<ngpoints<<"\n";
		gridfile<<"FIELD FieldData "<<1<<"\n";
		gridfile<<"Id "<<1<<" "<<ngpoints<<" int\n";

		std::cout<<"writing grid file name: "+fname+"\n";

		// id
		for (int i = 0; i < ngpoints; ++i){
			gridfile <<" "<<i;
		}
		gridfile <<"\n";

		// close the file
		gridfile.close();
	}

	void writeparticles()
	{
		if(Edian==""){ defineEdian(); }

		// open the file
		std::ofstream partfile;
		std::string fname = IO::pathGet() + "particles" + ".vtu";
		partfile.open(fname.c_str());
		partfile.precision(4);

		std::cout<<"writing vtk particles file: "+fname+"\n";
		
		// write results
		partfile <<"<?xml version=\"1.0\"?>\n";
		partfile <<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\""<<Edian.c_str()<<"\">\n";
		partfile <<"<UnstructuredGrid>\n";

		// points
		int ngpoints = Model::GetParticleVector().size();
		partfile <<std::fixed<<"\t<Piece NumberOfPoints=\""<<ngpoints<<"\" NumberOfCells=\""<<ngpoints<<"\">\n";
		partfile <<"\t\t<Points>\n";

		// position
		std::vector<Model::Particle> particles = Model::GetParticleVector();
		partfile <<"\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
		for (int i = 0; i < ngpoints; ++i){
		  partfile<<std::scientific<<"\t\t\t\t"<<particles[i].pos.x<<"\t"<<particles[i].pos.y<<"\t"<<particles[i].pos.z<<"\n";
		}
		partfile <<"\t\t\t</DataArray>\n";

		partfile <<"\t\t</Points>\n";
	
		// scalar values 
		partfile <<"\t\t\t<PointData Scalars=\"scalars\">\n";
		
		// user fixed points
		partfile <<"\t\t\t\t<DataArray type=\"Float32\" Name=\"Mat\" Format=\"ascii\">\n";
		for (int i = 0; i < ngpoints; ++i){
			partfile<<std::scientific<<"\t\t\t\t\t"<<particles[i].matid<<"\n";
		}
		partfile <<"\t\t\t\t</DataArray>\n";

		// id
		partfile <<"\t\t\t\t<DataArray type=\"Float32\" Name=\"Id\" Format=\"ascii\">\n";
		for (int i = 0; i < ngpoints; ++i){
			partfile<<std::scientific<<"\t\t\t\t\t"<<i<<"\n";
		}
		partfile <<"\t\t\t\t</DataArray>\n";

		// Zo
		partfile <<"\t\t\t\t<DataArray type=\"Float32\" Name=\"Zo\" Format=\"ascii\">\n";
		for (int i = 0; i < ngpoints; ++i){
			partfile<<std::scientific<<"\t\t\t\t\t"<<particles[i].pos.z<<"\n";
		}
		partfile <<"\t\t\t\t</DataArray>\n";


		partfile <<"\t\t\t</PointData>\n";

		// cells -->
		partfile <<"\t\t\t<Cells>\n";
		partfile <<"\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";
		for (int i = 0; i < ngpoints; ++i){
			partfile<<std::scientific<<"\t\t\t\t\t"<<i<<"\n";
		}
		partfile <<"\t\t\t\t</DataArray>\n";

		partfile <<"\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n";
		for (int i = 0; i < ngpoints; ++i){
			partfile <<"\t\t\t\t\t"<<i+1<<"\n";
		}
		partfile <<"\t\t\t\t</DataArray>\n";

		partfile <<"\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" Format=\"ascii\">\n";
		for (int i = 0; i < ngpoints; ++i){
			partfile <<"\t\t\t\t\t"<<1<<"\n";
		}
		partfile <<"\t\t\t\t</DataArray>\n";
		partfile <<"\t\t\t</Cells>\n";
		//<-- cells

		partfile <<"\t</Piece>\n";
		partfile <<"</UnstructuredGrid>\n";
		partfile <<"</VTKFile>\n";

		// close the file
		partfile.close();
	}
}
