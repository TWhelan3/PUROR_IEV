#include "ismrmrd/dataset.h"
#include "ismrmrd/xml.h"
#include "unistd.h"
#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
int main(int argc, char *argv[]){//take the filename as an argument

	int pos,posb;
	
	struct ISMRMRD::SequenceParameters seq_par;

	if(argc!=2){
	std::cout<<"This file takes the name of an ismrmrd file as its one argument. Program exiting."<<std::endl; 
	return -1;
	}


	ISMRMRD::Dataset d(argv[1], "/images", false);
	//ismrmrd_dataset_exists isn't in library yet but is coming soon I think (https://github.com/ismrmrd/ismrmrd/issues/29)
	/*if(d) 
	{
		std::cout<<"Input file not found. Check to make sure the path is correct."<<std::endl;
		return -1;
	}*/
	
	std::string xml;
	std::string new_xml;
	std::string starttag= "\n\t\t\t<value>";	//all of the things we're adding will be value tags, likely three tabs in
	std::string endtag= "</value>\n";		//end of value tag
	std::string value;				//what to add
	std::string reader;				//buffer for lines read in
	std::string tofollow; 				//trigger line
	std::string fileNameBase= argv[1];		//name for new file(s)

	pos=fileNameBase.find_first_of(".");
	posb=fileNameBase.find_last_of("/");
	if(posb==-1)
		posb=0;
	fileNameBase = fileNameBase.substr(posb,pos);

	new_xml=fileNameBase+".xml";		//cut off .ismrmrd and add .xml
	
	std::ifstream base_xml("base.xml");
	if(!base_xml.is_open())
	{
		std::cout<<"Failed to open base xml file\n"<<std::endl;
		return -1;
	}
	std::ofstream xml_output(new_xml.c_str());
	if(!xml_output.is_open())
	{
		std::cout<<"Failed to create new xml file\n"<<std::endl;
		return -1;
	}

	d.readHeader(xml);
	ISMRMRD::IsmrmrdHeader hdr;
	ISMRMRD::deserialize(xml.c_str(),hdr);
	seq_par=hdr.sequenceParameters.get();

	//Add input file to xml as LoadIsmrmrdDatasetImages' filename parameter

	tofollow = "<name>filename</name>";
	value = starttag + argv[1] + endtag;						//<value>xyz.ismrmrd</value>
	std::getline(base_xml, reader);									
	pos=reader.find_first_of("<"); 							//safe? should handle -1 (no < found)
	while(reader.substr(pos,std::string::npos).compare(tofollow)!=0){		//read line and compare to trigger line
		xml_output<<reader<<std::endl;							//copy
		std::getline(base_xml, reader);							//get next line
		pos=reader.find_first_of("<");							//find first <
	}
	xml_output<<reader<<value;							//if triggered, output the new value line

	//Add volume parameters for 3D unwrap/Autoscaling here
	

	//Add output file to xml as SaveDICOM's filename parameter

	tofollow = "<name>dicom_dir</name>";
	value = starttag + fileNameBase+ ".dicom" + endtag;		//<value>xyz.dicom</value>
	std::getline(base_xml, reader);
	pos=reader.find_first_of("<");
	
	while(reader.substr(pos,std::string::npos).compare(tofollow)!=0){		//read line and compare to trigger line
		xml_output<<reader<<std::endl;							//copy
		std::getline(base_xml, reader);							//get next line
		pos=reader.find_first_of("<");							//find first <
		
	}
	xml_output<<reader<<value;

	//look for end of file, and copy until it's found
	tofollow = "</gadgetronStreamConfiguration>";					
	std::getline(base_xml, reader);
	pos=reader.find_first_of("<");

	while(reader.substr(pos,std::string::npos).compare(tofollow)!=0){
	xml_output<<reader<<std::endl;
	std::getline(base_xml, reader);
	pos=reader.find_first_of("<");
	}
	xml_output<<reader;
 	
	xml_output.close();
	base_xml.close();
	///end of xml creation

	//figure out sqsub arguments from hdr info
	long int num_echos=hdr.encoding[0].encodingLimits.contrast().maximum+1;
	long int num_chan=hdr.acquisitionSystemInformation.get().receiverChannels.get();
	long int x = hdr.encoding[0].reconSpace.matrixSize.x;
	long int y = hdr.encoding[0].reconSpace.matrixSize.y;
	long int z = hdr.encoding[0].reconSpace.matrixSize.z;
	double run_minutes_d = x*y*z*num_echos*num_chan/120000000.0; 
	std::string run_minutes=std::to_string(ceil(run_minutes_d));
	char run_gigs[]="8";//how am I going to figure this out? 
	//std::cout<<x<<" "<<y<<" "<<z<<" "<<num_echos<<" "<<num_chan<<" "<<run_minutes_d <<" "<<run_minutes<<" "<<run_gigs<<std::endl;	
	if(execl("/work/twhelan5/local/bin/submitscript", "submitscript", fileNameBase.c_str(),run_minutes.c_str(),run_gigs,  0)==-1) //run the program
	//	std::cout<<"Whoops"<<std::endl;
	
		

	return 0;
}
