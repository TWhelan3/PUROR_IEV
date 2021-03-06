//PhaseMaskGadget.cpp
//Written by Tim Whelan 2015
//Input ImageHeader->Float 3D Array (Phase Data) -> //
//Output ImageHeader->Float 3D Array (Phase Data) -> //

#include "PhaseMaskGadget.h"
namespace Gadgetron{

int PhaseMaskGadget::process_config(ACE_Message_Block* mb)
{
	//max.value(std::stof(dynamic_cast<Gadget*>(this->next())->find_property("numVol")->string_value()));
	GINFO("Max value =%f\n",max.value());
	//GINFO("Next gadget is %s\n", dynamic_cast<Gadget*>(this->next())->description());
	this->get_controller();//->find_gadget("Unwrap")->print_properties();
	
	return GADGET_OK;
}
int PhaseMaskGadget::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1)
{
	GadgetContainerMessage<hoNDArray<float > > *m2 =AsContainerMessage<hoNDArray<float >> (m1->cont());

	Gadget* g=this->get_controller()->find_gadget("Unwrap");

	if(g)
	g->print_properties();
	
	if(!m2){

		GERROR("Need phase image to mask!");
		return GADGET_FAIL;
	}

	int num_pixels=m2->getObjectPtr()->get_number_of_elements();
	float* pixel=m2->getObjectPtr()->get_data_ptr();
	float maxvalue,minvalue,range;

	maxvalue=max.value();
	minvalue=min.value();
	range=maxvalue-minvalue;
	for(int i=0;i<num_pixels;i++)
	{
	
		if(pixel[i]<minvalue)
			pixel[i]=0;	
		else if(pixel[i]>maxvalue)
			pixel[i]=1;
		else
			pixel[i]=(pixel[i]-minvalue)/range;


	}

	if (this->next()->putq(m1) == -1) {
		m1->release();
		GERROR("Unable to send image.\n");
	    	return GADGET_FAIL;
	}



	return GADGET_OK;
}
GADGET_FACTORY_DECLARE(PhaseMaskGadget)
}












