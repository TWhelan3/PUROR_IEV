#include "ismrmrd/dataset.h"
#include "ismrmrd/xml.h"
#include "unistd.h"
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <complex>

int main(int argc, char *argv[]){//take the filename as an argument

	int pos;
	unsigned long long int num_echos,namelen;
	struct ISMRMRD::SequenceParameters seq_par;

	std::stringstream ss;
	std::string id;

	if(argc!=2){
	std::cout<<"This file takes the name of an ismrmrd file as its one argument. Program exiting."<<std::endl; //since this is automated what's the point
	return -1;
	}
	ISMRMRD::Dataset datt("newset.ismrmrd", "/images", true);
	ISMRMRD::Dataset d(argv[1], "/images", false);
	ISMRMRD::Image<std::complex<float> > image;
	ISMRMRD::ImageHeader hdr;
	int numImages = d.getNumberOfImages("complex");
	for (uint32_t i=0; i<numImages; i++) {
	 try {
            d.readImage("complex", i, image);
         }
         catch (std::exception &ex) {
	std::cout<<"FAIL\n";
                  return -1;
         }
	ss.str("");
	ss<<i;
	id.clear();
	id+="<?xml version=\"1.0\"?> <ismrmrdMeta><meta><name>GADGETRON_ImageNumber</name><value>";
	id+=ss.str();
	id+="</value></meta></ismrmrdMeta>";
	image.setAttributeString(id);
	datt.appendImage("complex", image);
	

	}
	std::string xml;
	d.readHeader(xml);
   	datt.writeHeader(xml);

}














