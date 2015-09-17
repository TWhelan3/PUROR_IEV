#ifndef PhaseMaskGadget_H_
#define PhaseMaskGadget_H_

#include "PUROR.h"
namespace Gadgetron{

  class EXPORTGADGETSPUROR PhaseMaskGadget:
public Gadget1<ISMRMRD::ImageHeader>
    {

    public:
      	GADGET_DECLARE(PhaseMaskGadget);

    protected:
      
     	virtual int process_config(ACE_Message_Block*);
	virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1);
	
	GADGET_PROPERTY(min, float, "Completely suppresed", PI*-1);
	GADGET_PROPERTY(max, float, "Unity Value", 0);

	float FOVx;
	float FOVy;

     };
}

#endif /* PhaseMaskGadget_H_ */
