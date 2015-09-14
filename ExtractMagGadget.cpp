//ExtractMagGadget.cpp
//Written by Tim Whelan 2015
//Input ImageHeader->Complex Float 3D Array (Image Data) -> [MetaContainer]-> //
//Output ImageHeader->Float 3D Array (Phase Data) -> Int 3D Array (Support Mask) -> [Int 3D Array (Support Mask)] -> [MetaContainer] ->//
#include "ExtractMagGadget.h"
using namespace Gadgetron;


int ExtractMagGadget::process_config(ACE_Message_Block* mb)
{	
this->msg_queue()->high_water_mark(128);//This helps with memory. It's not a hard limit though. 

return GADGET_OK;
}



int ExtractMagGadget::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1)
{

    	//std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
      
	GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2 = AsContainerMessage< hoNDArray< std::complex<float> > > (m1->cont());
	GadgetContainerMessage<ISMRMRD::MetaContainer> *meta = AsContainerMessage<ISMRMRD::MetaContainer>(m2->cont());

	yres = m2->getObjectPtr()->get_size(0);
	xres = m2->getObjectPtr()->get_size(1);
	cres = m2->getObjectPtr()->get_size(3);

	int c;
	GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
	GadgetContainerMessage<hoNDArray< float > > *cm2 = new GadgetContainerMessage<hoNDArray< float > >();


	*cm1->getObjectPtr() = *m1->getObjectPtr();//correct way to copy? based on extract gadget, since does a similar task
	cm1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_FLOAT;//GADGET_IMAGE_REAL_FLOAT;
	cm1->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_PHASE;//BECAUSE I"M LAZY

	cm1->cont(cm2);



	hoNDArray< float > *mag = cm2->getObjectPtr();

	boost::shared_ptr< std::vector<size_t> > dims = m2->getObjectPtr()->get_dimensions();
	try{mag->create(dims.get());}
	catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to create itohx in Unwrap Gadget");
		return GADGET_FAIL;
	}


	std::complex<float>* src = m2->getObjectPtr()->get_data_ptr();
	
	float*  dst = cm2->getObjectPtr()->get_data_ptr();
	#pragma omp parallel for private(c)
	for(c = 0; c<cres; c++)
	{
		
		int ch_offset= xres*yres*c;

		for (int i = xres*yres*c; i < xres*yres*(c+1); i++) 
		{
		dst[i] = abs(src[i]);
		}
	}
	
	
	
	 if (this->next()->putq(cm1) < 0) { //pass to next gadget
		GERROR("Failed to initialize Mask\n");
		cm1->release();
	   return GADGET_FAIL;
	}
	m2->cont(NULL);//necessary so masks don't get deleted and stay connected. Release removes links further down chain.
	m1->release(); 
	return GADGET_OK;
	
}

GADGET_FACTORY_DECLARE(ExtractMagGadget)

