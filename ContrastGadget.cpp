//ContrastGadget.cpp
//Written by Tim Whelan 2015
//Input ImageHeader->Float 3D Array (Phase Data) -> //
//Output ImageHeader->Float 3D Array (Phase Data) -> //
//Sums individual echos
#include "ContrastGadget.h"


namespace Gadgetron{

int ContrastGadget::process_config(ACE_Message_Block* mb)
{
	ISMRMRD::IsmrmrdHeader hdr;
        ISMRMRD::deserialize(mb->rd_ptr(),hdr);
	numEchos=hdr.encoding[0].encodingLimits.contrast().maximum +1; //number of echos is one more than highest numbers (0-based)
	numSlices=hdr.encoding[0].reconSpace.matrixSize.z; //number of slices (will this always work?)
	return GADGET_OK;
}
int ContrastGadget::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1)
{
	GadgetContainerMessage<hoNDArray< float > > *image_message =     AsContainerMessage<hoNDArray<float>>(m1->cont());
	int echo = m1->getObjectPtr()->contrast;	
	static hoNDArray< float >  *collapsed;
	static float* cptr;
	float* image_pixels = image_message->getObjectPtr()->get_data_ptr();
	int yres = image_message->getObjectPtr()->get_size(0);
	int xres = image_message->getObjectPtr()->get_size(1);
	if(echo==0)//If reordered properly, this will happen before others (within each slice)
	{
		collapsed = new hoNDArray< float >();
		try{collapsed->create(xres,yres);}
			catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to create collapsed space");
		return GADGET_FAIL;
			}
		cptr=collapsed->get_data_ptr();
		for (int i = 0; i < xres*yres; i++)
		{
	 	 cptr[i]=image_pixels[i];	//memcpy?
		}

	}
	else
	{
		for (int i = 0; i < xres*yres; i++)
		{
	  	cptr[i]+=image_pixels[i];		
		}
	}
	if(pass_on.value()==YES && (this->next()->putq(m1) == -1)) {//is this unclear?
		m1->release();
		GDEBUG("Unable to put collapsed images on next gadgets queue\n");
		return GADGET_FAIL;
	}

	if(echo==(numEchos-1))
	{	
		//divide by number of echos!
		GadgetContainerMessage<ISMRMRD::ImageHeader>* h1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
		
		*h1->getObjectPtr()= *m1->getObjectPtr();
		GadgetContainerMessage<hoNDArray< float > > *outimage = new GadgetContainerMessage< hoNDArray<float> >(collapsed);
	
		h1->cont(outimage);
		h1->getObjectPtr()->image_series_index=numEchos;//what should this be
		//h1->getObjectPtr()->image_index = (h1->getObjectPtr()->image_index-1) % numSlices;//+1;//v.s.
		h1->getObjectPtr()->image_index = (h1->getObjectPtr()->image_index);// % numSlices;//+1;//v.s.
		h1->getObjectPtr()->set =10;
		

		if (this->next()->putq(h1) == -1) {
		m1->release();
		GDEBUG("Unable to put collapsed images on next gadgets queue\n");
		return GADGET_FAIL;
		}
	}

}
GADGET_FACTORY_DECLARE(ContrastGadget)
}
