#include "ismrmrd/dataset.h"
#include "ismrmrd/xml.h"
#include "unistd.h"
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char *argv[]){//take the filename as an argument

	int pos;
	unsigned long long int num_echos,namelen;
	struct ISMRMRD::SequenceParameters seq_par;

	if(argc!=2){
	std::cout<<"This file takes the name of an ismrmrd file as its one argument. Program exiting."<<std::endl; //since this is automated what's the point
	return -1;
	}
	int childPID;

	
	ISMRMRD::Dataset d(argv[1], "/images", false);

	std::string xml;
	d.readHeader(xml);
	ISMRMRD::IsmrmrdHeader hdr;
	ISMRMRD::deserialize(xml.c_str(),hdr);
	std::stringstream attributes;
	hdr.measurementInformation.get().seriesDescription.set("[testing]");

	ISMRMRD::serialize(hdr, attributes);
      
	d.writeHeader(attributes.str());

}
