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

	float new_phase_dir[3];
	float new_read_dir[3];
	float*  src;
	float*  dst;

	int d;

	GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
	GadgetContainerMessage<hoNDArray< float > > *cm2 = new GadgetContainerMessage<hoNDArray< float > >();


	*cm1->getObjectPtr() = *m1->getObjectPtr();
	cm1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_FLOAT;//GADGET_IMAGE_REAL_FLOAT;

	cm1->cont(cm2);
	cm2->cont(m2->cont());
	m2->cont(NULL);

	//Want to make read_dir 1/0/0 and phase_dir 0/1/0
	//Assume we have it
	if(std::abs(m1->getObjectPtr()->phase_dir[0])>std::abs(m1->getObjectPtr()->read_dir[0]))
	{
		swap=true;
		for(d=0; d<3; d++)
		{
			new_read_dir[d]=m1->getObjectPtr()->phase_dir[d];
			new_phase_dir[d]=m1->getObjectPtr()->read_dir[d];
			
		}

		cm1->getObjectPtr()->field_of_view[0]=m1->getObjectPtr()->field_of_view[1];
		cm1->getObjectPtr()->field_of_view[1]=m1->getObjectPtr()->field_of_view[0];
		//pixel spacing?
		//slice location?
	}
	else
	{	
		swap=false;
		for(int d=0; d<3; d++)
		{
			new_read_dir[d]=m1->getObjectPtr()->read_dir[d];
			new_phase_dir[d]=m1->getObjectPtr()->phase_dir[d];
		}
	}

	if(new_read_dir[0]<0)//changed from <= Jun 28 2017
	{
		read_flip=true;
		for(d=0; d<3; d++)
		{
			cm1->getObjectPtr()->read_dir[d]=-new_read_dir[d];
		}
	}
	else
	{
		read_flip=false;
		for(d=0; d<3; d++)
		{
			cm1->getObjectPtr()->read_dir[d]=new_read_dir[d];
		}
	}

	if(new_phase_dir[1]<0)//changed from <= Jun 28 2017
	{
		phase_flip=true;
		for(d=0; d<3; d++)
		{
			cm1->getObjectPtr()->phase_dir[d]=-new_phase_dir[d];
		}
	}
	else
	{
		phase_flip=false;
		for(d=0; d<3; d++)
		{
			cm1->getObjectPtr()->phase_dir[d]=new_phase_dir[d];
		}
	}

	if(swap)
	{
		yres = m1->getObjectPtr()->matrix_size[0]; //flipped from before
		xres = m1->getObjectPtr()->matrix_size[1];
	}
	else
	{
		yres = m1->getObjectPtr()->matrix_size[1]; 
		xres = m1->getObjectPtr()->matrix_size[0];
	}

	cm1->getObjectPtr()->matrix_size[1]=yres;
	cm1->getObjectPtr()->matrix_size[0]=xres;

	try{cm2->getObjectPtr()->create(yres*xres);}
		catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to allocate array in Flip Gadget.\n");
		return GADGET_FAIL;
	}

	dst=cm2->getObjectPtr()->get_data_ptr();
	src=m2->getObjectPtr()->get_data_ptr();

	if(swap && !read_flip && !phase_flip)
	{
		for (int i = 0; i < xres; i++) 
		{
			for (int j =0; j< yres; j++)
				dst[j*xres+i] = src[i*yres+j];
		}
	}
	else if(!swap && read_flip && !phase_flip)
	{
		for (int i = 0; i < yres; i++) 
		{
			for (int j =0; j< xres; j++)
				dst[i*xres+j] = src[i*xres+xres-j-1];
		}
	}
	else //added catchall to make sure copy happens if swap or flip criteria are not met
	{
		for (int ij=0; ij<xres*yres; ij++)
			dst[ij]=src[ij];
	}



	//Acquisition Matrix might be screwed up but DicomFinishGadget might also screw it up

	if (this->next()->putq(cm1) < 0) { //pass to next gadget
		GERROR("Failed to pass on magnitude\n");
		cm1->release();
	   return GADGET_FAIL;
	}
	
	return GADGET_OK;
	
}

GADGET_FACTORY_DECLARE(FlipGadget)

