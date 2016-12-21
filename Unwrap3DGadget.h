#ifndef Unwrap3DGadget_H_
#define Unwrap3DGadget_H_

#include "PUROR.h"
#include "hdf5.h"
namespace Gadgetron{

  class EXPORTGADGETSPUROR Unwrap3DGadget:
public Gadget1<ISMRMRD::ImageHeader>
    {

    public:
      	GADGET_DECLARE(Unwrap3DGadget);

    protected:
      
     	virtual int process_config(ACE_Message_Block*);
	virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1);
	GADGET_PROPERTY(filename, std::string, "Name of ISMRMRD Dataset file", "output");
	GADGET_PROPERTY(savephase, int, "Save Before Passing On", 0);
	ISMRMRD::Dataset* dsToWrite;


	hsize_t yres;
	hsize_t xres;
	hsize_t cres;
	ISMRMRD::Dataset* temp_storage;
	hsize_t num_slices;
	hsize_t num_echos;
	hsize_t num_ch;
	
	std::vector<int> ordering;//how would I get this from inputs?
	std::vector<int> VolumeEnds;
	std::vector<std::vector<std::vector<float>>> slice_mean;
	GADGET_PROPERTY(limits, std::string, "Ends of Volumes", "");
	GADGET_PROPERTY(numVol, int, "Number of Volumes", 1);
	GADGET_PROPERTY(do3D, int, "Unwrap Volumes Slicewise", 1); 

	float FOVx;
	float FOVy;

     };
}

#endif /* Unwrap3DGadget_H_ */
