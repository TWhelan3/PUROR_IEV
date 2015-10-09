//SWIGadget.h

#ifndef SWIGADGET_H
#define SWIGADGET_H

#include "PUROR.h"
#include <exception>

namespace Gadgetron
{

class EXPORTGADGETSPUROR SWIGadget : 
public Gadget1<ISMRMRD::ImageHeader>
{
 public:
  GADGET_DECLARE(SWIGadget)

    protected:

	virtual int process_config(ACE_Message_Block* mb);
  	virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1);

	GADGET_PROPERTY(filename, std::string, "Name of ISMRMRD Dataset file", "filename@LoadIsmrmrdDatasetImages");
        GADGET_PROPERTY(groupname, std::string, "Name of group in ISMRMRD Dataset file", "groupname@LoadIsmrmrdDatasetImages");
	GADGET_PROPERTY(varname, std::string, "Name of variable within group in ISMRMRD Dataset file", "varname@LoadIsmrmrdDatasetImages");
	GADGET_PROPERTY(timetimes, int, "How many times to multiply by phase mask", 4);
	
	int xres;
	int yres;
	int cres;

	int numEchos;
	int numChan;

	ISMRMRD::Dataset *pDataset;
};

}  
#endif //SWIGADGET_H
