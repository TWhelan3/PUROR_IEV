#include "ismrmrd/dataset.h"
#include "ismrmrd/xml.h"
#include <iostream>
#include <string>
int main(int argc, char *argv[]){//take the filename as an argument

	ISMRMRD::IsmrmrdHeader hdr;
	int numImages_data, numImages_header;	//numbers to compare 
	int numSlices, numEchos;		//two parts of number of images in header	
	
	std::string xml;
	ISMRMRD::Dataset d(argv[1], "/images", false); //open dataset

	

	
	d.readHeader(xml);
	
	ISMRMRD::deserialize(xml.c_str(),hdr);

	numSlices = hdr.encoding[0].reconSpace.matrixSize.z;
	numEchos  = hdr.encoding[0].encodingLimits.contrast->maximum + 1; //0-based

	numImages_header = numSlices*numEchos;
	numImages_data	 = d.getNumberOfImages("complex");

	std::cout<<"Slices = "<< numSlices <<"\nEchos = "<<numEchos<<"\nHeader = "<<numImages_header<<"\nDataset = "<<numImages_data<<std::endl;	


	return (numImages_header!=numImages_data); //return 0 if they match

}
