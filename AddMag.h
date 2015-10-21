//AddMag.h

#ifndef AddMag_H
#define AddMag_H

#include "PUROR.h"
#include <exception>

namespace Gadgetron
{

class EXPORTGADGETSPUROR AddMag : 
public Gadget1<ISMRMRD::ImageHeader>
{
 public:
  GADGET_DECLARE(AddMag)

    protected:

	virtual int process_config(ACE_Message_Block* mb);
  	virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1);

	GADGET_PROPERTY(filename, std::string, "Name of ISMRMRD Dataset file", "inputfilename");//@LoadIsmrmrdDatasetImages"); Gadgets can only see things after them for now
        GADGET_PROPERTY(groupname, std::string, "Name of group in ISMRMRD Dataset file", "groupname");//@LoadIsmrmrdDatasetImages");
	GADGET_PROPERTY(varname, std::string, "Name of variable within group in ISMRMRD Dataset file", "varname");//@LoadIsmrmrdDatasetImages");
	
	//std::string filename,groupname,varname;
	
	int xres;
	int yres;
	

	int numEchos;
	int numChan;
	int numSlices;

	ISMRMRD::Dataset *pDataset;
};

}  
#endif //AddMag_H
