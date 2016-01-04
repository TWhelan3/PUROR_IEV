#ifndef ReorderGadget_H_
#define ReorderGadget_H_

#include "PUROR.h"
namespace Gadgetron{

  class EXPORTGADGETSPUROR ReorderGadget:
public Gadget1<ISMRMRD::ImageHeader>
    {

    public:
      	GADGET_DECLARE(ReorderGadget);

    protected:
      
     	virtual int process_config(ACE_Message_Block*);
	virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1);

	ISMRMRD::Dataset* temp_storage;
	int num_slices;
	int num_echos;
		
	std::vector<int> ordering;//how would I get this from inputs?
	std::vector<int> VolumeEnds;
	std::vector<std::vector<std::vector<float>>> slice_mean;
	GADGET_PROPERTY(limits, std::string, "Ends of Volumes", "");
	GADGET_PROPERTY(numVol, int, "Number of Volumes", 1);

     };
}

#endif /* ReorderGadget_H_ */
