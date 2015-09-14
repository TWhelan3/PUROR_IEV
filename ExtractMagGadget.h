//ExtractMagGadget.h

#ifndef ExtractMagGADGET_H
#define ExtractMagGADGET_H

#include "PUROR.h"


namespace Gadgetron
{

class EXPORTGADGETSPUROR ExtractMagGadget : 
public Gadget1<ISMRMRD::ImageHeader>
{
 public:
  GADGET_DECLARE(ExtractMagGadget)

    protected:

	virtual int process_config(ACE_Message_Block* mb);
  	virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1);
	
	int xres;
	int yres;
	int cres;
};

}  
#endif //ExtractMagGADGET_H
