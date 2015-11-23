//ContrastGadget.cpp
//Written by Tim Whelan 2015
//Input ImageHeader->Float 3D Array (Phase Data) -> //
//Output ImageHeader->Float 3D Array (Phase Data) -> //
//Sums individual echos
#include "ContrastGadget.h"


namespace Gadgetron{

int ContrastGadget::process_config(ACE_Message_Block* mb)
{
	this->msg_queue()->high_water_mark(128);//This helps with memory. It's not a hard limit though. 
	ISMRMRD::IsmrmrdHeader hdr;
        ISMRMRD::deserialize(mb->rd_ptr(),hdr);
	numEchos=hdr.encoding[0].encodingLimits.contrast().maximum +1; //number of echos is one more than highest numbers (0-based)
	numSlices=hdr.encoding[0].reconSpace.matrixSize.z; //number of slices (will this always work?)

	/*phase_sum = new hoNDArray< float >();
	try{phase_sum->create(hdr.encoding[0].reconSpace.matrixSize.y,hdr.encoding[0].reconSpace.matrixSize.x);}
	catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to create phase_sum space.\n");
		return GADGET_FAIL;
	}
	phase_sum_ptr=phase_sum->get_data_ptr();

	mag_sum = new hoNDArray< float >();
	try{mag_sum->create(hdr.encoding[0].reconSpace.matrixSize.y,hdr.encoding[0].reconSpace.matrixSize.x);}
	catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to create mag_sum space.\n");
		return GADGET_FAIL;
	}
	mag_sum_ptr=mag_sum->get_data_ptr();*/


	numPixels=hdr.encoding[0].reconSpace.matrixSize.x *hdr.encoding[0].reconSpace.matrixSize.y;
	return GADGET_OK;
}
int ContrastGadget::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1)
{
	GadgetContainerMessage<hoNDArray< float > > *image_message =     AsContainerMessage<hoNDArray<float>>(m1->cont());
	GadgetContainerMessage<ISMRMRD::MetaContainer> *meta = AsContainerMessage<ISMRMRD::MetaContainer>(image_message->cont());
	
	if(!image_message){
		GERROR("Phase array (float/single) expected and not found.\n");
		
		return GADGET_FAIL;
	}
	int echo = m1->getObjectPtr()->contrast;	
	float* image_pixels = image_message->getObjectPtr()->get_data_ptr();
	
	std::string type = meta->getObjectPtr()->as_str(GADGETRON_DATA_ROLE);
	int typeIndex;
	bool found=false;	
	int ii;
	for(ii=0; ii<types.size(); ii++)
		if(strcmp(types[ii].c_str(),type.c_str())==0)
		{
			type=types[ii];
			found=true;
			break;
		}
	typeIndex=ii;
	if(!found)
	{
		types.push_back(type);
		sum_arrays.push_back(new hoNDArray<float>());
		
		try{sum_arrays[typeIndex]->create(numPixels);}
		catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to create mag_sum space.\n");
		return GADGET_FAIL;
		}
	}


	float* dst_ptr=sum_arrays[typeIndex]->get_data_ptr();	

	if(echo==0)
	{
		memcpy(dst_ptr, image_pixels, numPixels*sizeof(float));
	}
	else
	{
		#pragma omp parallel for	
		for (int i = 0; i < numPixels; i++)
		{
	  	dst_ptr[i]+=image_pixels[i];		
		}
			
	}
	

	/*
	if(echo==0)
	{
		memcpy(phase_sum_ptr, image_pixels, numPixels*sizeof(float));
	}
	else if(echo == 130)
	{
		memcpy(mag_sum_ptr, image_pixels, numPixels*sizeof(float));

	}
	else
	{
		
		if(m1->getObjectPtr()->image_type==ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE)
		{
			#pragma omp parallel for	
			for (int i = 0; i < numPixels; i++)
			{
		  	mag_sum_ptr[i]+=image_pixels[i];		
			}
		}
		else
		{	
			#pragma omp parallel for 
			for (int i = 0; i < numPixels; i++)
			{
		  	phase_sum_ptr[i]+=image_pixels[i];		
			}
		}
	}*/
 
	if(pass_on.value()==YES && (this->next()->putq(m1) == -1)) {//pass on individual echos if config says to, if pass fails, throw error, if sending to ismrmrd can throw off order
		m1->release();
		GERROR("Unable to put phase_sum images on next gadgets queue.\n");
		return GADGET_FAIL;
	}
	

	if(echo==(numEchos-1))//if(echo==(numEchos-1)||echo==(numEchos+129))//Phase or Magnitude
	{	
		
		//divide by number of echos!
		GadgetContainerMessage<ISMRMRD::ImageHeader>* h1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
		GadgetContainerMessage<hoNDArray< float > > *outimage;
		*h1->getObjectPtr()= *m1->getObjectPtr();

		/*if(m1->getObjectPtr()->image_type==ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE)
			outimage = new GadgetContainerMessage< hoNDArray<float> >(mag_sum);
		else
			outimage = new GadgetContainerMessage< hoNDArray<float> >(phase_sum);*/

		outimage= new GadgetContainerMessage<hoNDArray<float>>(sum_arrays[typeIndex]);
		
		float* out_ptr = outimage->getObjectPtr()->get_data_ptr();
		float meanterm=1/(float)numEchos;
	
		#pragma omp parallel for	
		for (int i = 0; i < numPixels; i++)
			out_ptr[i]*=meanterm;		
			
				


		h1->cont(outimage);
		h1->getObjectPtr()->slice =m1->getObjectPtr()->slice;
		h1->getObjectPtr()->image_index = m1->getObjectPtr()->image_index+1000;
		h1->getObjectPtr()->contrast=m1->getObjectPtr()->contrast+1;
	
		if(meta)
		{
			GadgetContainerMessage<ISMRMRD::MetaContainer> *newmeta = new GadgetContainerMessage<ISMRMRD::MetaContainer>();
			*(newmeta->getObjectPtr()) = *(meta->getObjectPtr());//actually copy instead of clone

			outimage->cont(newmeta);
		}		
		
		if (this->next()->putq(h1) == -1) {
		m1->release();
		GERROR("Unable to put phase_sum images on next gadget's queue.\n");
		return GADGET_FAIL;
		}
	}

}
GADGET_FACTORY_DECLARE(ContrastGadget)
}
