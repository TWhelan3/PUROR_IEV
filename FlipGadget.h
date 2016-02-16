//FlipGadget.h

#ifndef FlipGadget_H
#define FlipGadget_H

#include "PUROR.h"
#include <exception>

namespace Gadgetron
{

class EXPORTGADGETSPUROR FlipGadget : 
public Gadget1<ISMRMRD::ImageHeader>
{
 public:
  GADGET_DECLARE(FlipGadget)

    protected:

	virtual int process_config(ACE_Message_Block* mb);
  	virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1);


	//std::string filename,groupname,varname;
	
	int xres;
	int yres;
	
	bool read_flip, phase_flip, swap;

	int numEchos;
	int numChan;
	int numSlices;

	ISMRMRD::Dataset *pDataset;
};

}  
#endif //FlipGadget_H
