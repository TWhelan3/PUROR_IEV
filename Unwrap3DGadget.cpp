//Unwrap3DGadget.cpp
//Written by Tim Whelan 2015
//Input ImageHeader->Float 3D Array (Phase Data) -> // Int 3D Array (Support Mask) -> [Int 3D Array (Support Mask)] -> [MetaContainer]
//Output ImageHeader->Float 3D Array (Phase Data) -> // Int 3D Array (Support Mask) -> [Int 3D Array (Support Mask)] -> [MetaContainer]


#include "Unwrap3DGadget.h"

namespace Gadgetron{

int Unwrap3DGadget::process_config(ACE_Message_Block* mb)
{
	 boost::filesystem::remove("/tmp/temp.ismrmrd");
	temp_storage = new ISMRMRD::Dataset("/tmp/temp.ismrmrd","ooo_storage", true);

	ISMRMRD::IsmrmrdHeader hdr;
        ISMRMRD::deserialize(mb->rd_ptr(),hdr);
	//num_echos=6; num_slices=32; num_ch=32;
	num_echos=hdr.encoding[0].encodingLimits.contrast().maximum +1; //number of echos is one more than highest numbers (0-based)
	num_slices=hdr.encoding[0].reconSpace.matrixSize.z; //number of slices (will this always work?)
	//This is an optional part of the header and not sure if it will always be appropriate. # channels in acquisition may not be # channels in image gadget receives if they are already combined
	//for now leacing cres (actual image channels) separate
	num_ch=hdr.acquisitionSystemInformation.get().receiverChannels();


	//this->msg_queue()->high_water_mark(128);//This helps with memory. It's not a hard limit though. 
	ordering.resize(num_echos*num_slices);
	slice_mean.resize(num_echos);

	for(int e=0;e<num_echos; e++)
	{
		slice_mean[e].resize(num_ch);
		for(int ch=0;ch<num_ch; ch++)
			slice_mean[e][ch].resize(num_slices);
	}
	return GADGET_OK;
}

int Unwrap3DGadget::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1)
{
	GadgetContainerMessage<hoNDArray<float > > *m2 =AsContainerMessage<hoNDArray<float >> (m1->cont());
	GadgetContainerMessage<hoNDArray<int>> *supportmasks_block = AsContainerMessage<hoNDArray<int>>(m2->cont());
	GadgetContainerMessage<hoNDArray<int>> *masks_block;
	GadgetContainerMessage<ISMRMRD::MetaContainer> *meta;

	if((m2 && supportmasks_block)==0){

		GINFO("ERROR in u3dg\n");
		return GADGET_FAIL;
	}
	ISMRMRD::Image<float> image;
	ISMRMRD::ImageHeader *header=m1->getObjectPtr();
	int fullsignal =0;
	float tmp_mean=0;
	int k,ii,ch,e;

	int yres = m2->getObjectPtr()->get_size(0);	
	int xres = m2->getObjectPtr()->get_size(1);
	int cres = header->channels;
	float *unwrapped_x = m2->getObjectPtr()->get_data_ptr();
	int *supportmasks = supportmasks_block->getObjectPtr()->get_data_ptr();
	int *masks;

	float* im_ptr;
	//synthesize naming conventions between this and U2DG
	
	if(*supportmasks<0)//fullsignal masks
	{
		fullsignal=1;
		*supportmasks+=2;
		meta=  AsContainerMessage<ISMRMRD::MetaContainer>(supportmasks_block->cont());

	}
	else
	{
		masks_block = AsContainerMessage<hoNDArray<int>>(supportmasks_block->cont());
		if(!masks_block)
			GERROR("Soemthing's wrong");
		masks=masks_block->getObjectPtr()->get_data_ptr();
		meta = AsContainerMessage<ISMRMRD::MetaContainer>(supportmasks_block->cont()->cont());
	}



	image.setHead(*(header));
	/////
	im_ptr=image.getDataPtr();
	memcpy(im_ptr, unwrapped_x, image.getDataSize());
	



      std::stringstream attributes;
      
      if (meta) {
         ISMRMRD::serialize(*meta->getObjectPtr(), attributes);
	
      }
      image.setAttributeString(attributes.str());
	static int myplace=0;

	ordering[header->image_index-1]=myplace++;
	temp_storage->appendImage("here", image);


	/*if (this->next()->putq(m1) == -1) {
	m1->release();
	GDEBUG("Unable to put images on next gadgets queue\n");
	return GADGET_FAIL;
	}*/


	m1->release();
	

	if(myplace==num_slices*num_echos)
	{
		std::chrono::high_resolution_clock::time_point t1;
		std::chrono::high_resolution_clock::time_point t2;
		int duration;
		std::vector<float> slice_mean_original;
		std::vector<std::vector<std::vector<float>>> offset(slice_mean); //just need the size but in 3D
		float* im_ptr, *data_ptr;

		t1 = std::chrono::high_resolution_clock::now();
		//parallelize this
		//for outlerloop is echos or channels better?
		for(e=0; e<num_echos; e++)
		for(ch=0;ch<cres;ch++)
		{
			slice_mean_original=slice_mean[e][ch];
			unwrap(slice_mean[e][ch].data(), num_slices);
			tmp_mean=slice_mean[e][ch][num_slices/2-1];
			for(int sl=0; sl<num_slices; sl++)//number of slices
			{
				slice_mean[e][ch][sl]-=PI2*round(tmp_mean/PI2);
				if(fabs(slice_mean_original[sl]-slice_mean[e][ch][sl])>PI)
					offset[e][ch][sl]=slice_mean[e][ch][sl]-slice_mean_original[sl];
				else
					offset[e][ch][sl]=0;
			}
		}
 	
		///////////////////

	
	
			
		for(int i=0; i<num_slices; i++)//number of slices
		{
			for(int j=0; j<num_echos; j++)//number of echos/contrast
			{
				 // Read the image
				 try {
				    temp_storage->readImage("here", ordering[j*num_slices+i], image);///number of slices //ordering[num_slices*j+i]  <--necessary is thread >1
				 }
				 catch (std::exception &ex) {
				    GERROR("Error reading image %d\n",i);
				    return -1;
				 }
				 GadgetContainerMessage<ISMRMRD::ImageHeader>* mheader = new GadgetContainerMessage<ISMRMRD::ImageHeader >();
				 //GDEBUG("%d %d\n",  ordering[j*num_slices+i], j*num_slices+i);
				 header = mheader->getObjectPtr();
				 *header = image.getHead();
				 
				header->data_type = ISMRMRD::ISMRMRD_FLOAT;//GADGET_IMAGE_REAL_FLOAT;
				header->image_type = ISMRMRD::ISMRMRD_IMTYPE_PHASE;

				GadgetContainerMessage< hoNDArray<float> > *data = new GadgetContainerMessage< hoNDArray<float> >();
				
				data->getObjectPtr()->create(header->matrix_size[0], header->matrix_size[1],header->matrix_size[2], header->channels);
								
				//GINFO("%d %d %d %d\n",j*num_slices+i,header->matrix_size[0], header->matrix_size[1], header->channels);
				data_ptr=data->getObjectPtr()->get_data_ptr();
				
				memcpy(data_ptr, image.getDataPtr(), image.getDataSize()); //copy image 
				im_ptr=image.getDataPtr();
			
				std::string attributes2;
				image.getAttributeString(attributes2);
				 if (header->attribute_string_len > 0 && !attributes2.empty()) {
				    GadgetContainerMessage<ISMRMRD::MetaContainer> *meta = new GadgetContainerMessage<ISMRMRD::MetaContainer>();
				    
				    try {
				       ISMRMRD::deserialize(image.getAttributeString(), *meta->getObjectPtr());
				    }
				    catch (std::exception &ex) {
				       GERROR("Error parsing attribute string: %s\n", ex.what());
				       return -1;
				    }
				    data->cont(meta);
				 }
				 mheader->cont(data);
				 
				 // Send the image along the chain
				//while(this->msg_queue()->message_count()>1){}
				 if (this->next()->putq(mheader) == -1) {
				    mheader->release();
				    GERROR("Unable to send image.\n");
				    return -1;
				 }
				
			}
			if(i%10==0)
			GINFO("%d sent back out\n", i);
		}
		t2 = std::chrono::high_resolution_clock::now();
		 duration= std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
		GINFO("Total for 3D Unwrap:%d\n",duration/1000000);
	}
	
	//append (to particular set?)
	//possibly start figuring out order (how would I know???)
	//can save write order to make reading easier
	
//after some point
	//start loading back (possibly out of order)
	//either send pow pow pow or send image 'slab', latter may work better. 

return GADGET_OK;

}

GADGET_FACTORY_DECLARE(Unwrap3DGadget)

}