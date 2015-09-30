//SWIGadget.h

#ifndef SWIGADGET_H
#define SWIGADGET_H

#include "PUROR.h"


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
	
	int xres;
	int yres;
	int cres;
};

}  
#endif //SWIGADGET_H
