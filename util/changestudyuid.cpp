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

	 ISMRMRD::StudyInformation study_info;
	    study_info.studyDate.set("19000101");
            study_info.studyTime.set("170000");
            study_info.studyID.set("smallteststudy");
            study_info.accessionNumber.set(0);
            study_info.referringPhysicianName.set("XXXXXXXX");
            study_info.studyDescription.set("Gadgetron^Gadgetron^SmallTest");
            study_info.studyInstanceUID.set("1.2.826.0.1.3680043.9.5687.1");//.1 for small //.2 for large (those two test sets. use unique numbers in real life)

	 hdr.studyInformation.set(study_info);


	ISMRMRD::serialize(hdr, attributes);
      
	d.writeHeader(attributes.str());

}
