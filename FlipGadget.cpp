//FlipGadget.cpp
//Written by Tim Whelan 2015
//Throughput ImageHeader->Float 3D Array (Image Data) -> [MetaContainer]-> // 

//Output ImageHeader->Float 3D Array (Magnitude Data) -> [MetaContainer???] ->//

#include "FlipGadget.h"
using namespace Gadgetron;


int FlipGadget::process_config(ACE_Message_Block* mb)
{	
	this->msg_queue()->high_water_mark(128);//This helps with memory. It's not a hard limit though. 

return GADGET_OK;
}


int FlipGadget::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1)
{
	GadgetContainerMessage< hoNDArray< float > >* m2 = AsContainerMessage< hoNDArray< float > > (m1->cont());
	


	if(!m2){
		GERROR("Float image array expected and not found.\n");
		
		return GADGET_FAIL;
	}





	yres = m2->getObjectPtr()->get_size(1); //flipped from before
	xres = m2->getObjectPtr()->get_size(0);

	m1->getObjectPtr()->matrix_size[1]=yres;
	m1->getObjectPtr()->matrix_size[0]=xres;
	

	int c;
	ISMRMRD::Image<complex_float_t> image;

	GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
	GadgetContainerMessage<hoNDArray< float > > *cm2 = new GadgetContainerMessage<hoNDArray< float > >();


	*cm1->getObjectPtr() = *m1->getObjectPtr();
	cm1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_FLOAT;//GADGET_IMAGE_REAL_FLOAT;

		
	cm1->cont(cm2);
	cm2->cont(m2->cont());
	m2->cont(NULL);
	
	try{cm2->getObjectPtr()->create(yres*xres);}
	catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to allocate array in Flip Gadget.\n");
		return GADGET_FAIL;
	}
	
	float*  src;//= m2->getObjectPtr()->get_data_ptr();
	float*  dst;// = cm2->getObjectPtr()->get_data_ptr()
		
		
	dst=cm2->getObjectPtr()->get_data_ptr();
	src=m2->getObjectPtr()->get_data_ptr();
	for (int i = 0; i < xres; i++) 
	{
		for (int j =0; j< yres; j++)
		        dst[j*xres+i] = src[i*yres+j];
	}

	

	
	if (this->next()->putq(cm1) < 0) { //pass to next gadget
		GERROR("Failed to pass on magnitude\n");
		cm1->release();
	   return GADGET_FAIL;
	}
	return GADGET_OK;
	
}

GADGET_FACTORY_DECLARE(FlipGadget)

