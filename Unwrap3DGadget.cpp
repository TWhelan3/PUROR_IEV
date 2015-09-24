//Unwrap3DGadget.cpp
//Written by Tim Whelan 2015
//Input ImageHeader->Float 3D Array (Phase Data) -> // Int 3D Array (Support Mask) -> [Int 3D Array (Support Mask)] -> [MetaContainer]
//Output ImageHeader->Float 3D Array (Phase Data) -> // Int 3D Array (Support Mask) -> [Int 3D Array (Support Mask)] -> [MetaContainer]


#include "Unwrap3DGadget.h"

namespace Gadgetron{

int Unwrap3DGadget::process_config(ACE_Message_Block* mb)
{
	 boost::filesystem::remove("/tmp/temp.ismrmrd");
	temp_storage = new ISMRMRD::Dataset("/tmp/temp.ismrmrd","storage", true);
	this->msg_queue()->high_water_mark(128);//This helps with memory. It's not a hard limit though. 
	ISMRMRD::IsmrmrdHeader hdr;
        ISMRMRD::deserialize(mb->rd_ptr(),hdr);

	std::string id = "MxDoe";

	
	num_echos=hdr.encoding[0].encodingLimits.contrast().maximum +1; //number of echos is one more than highest numbers (0-based)
	num_slices=hdr.encoding[0].reconSpace.matrixSize.z; //number of slices (will this always work?)


	yres=hdr.encoding[0].reconSpace.matrixSize.x; //match my (and MATLABs) unfortunate convention 
	xres=hdr.encoding[0].reconSpace.matrixSize.y;

	if(hdr.subjectInformation.is_present())
	{
		struct ISMRMRD::SubjectInformation sInfo=hdr.subjectInformation.get();
		if(sInfo.patientName.is_present())
			id= sInfo.patientName.get();//assuming this field is acutally a pseudoanonymous id
	}

	std::cout<<"PATIENT ID:" << id <<std::endl;
	if(savephase.value()==1)
	{
		dsToWrite = new ISMRMRD::Dataset((id+filename.value()+".ismrmrd").c_str(),"images", true);
	}




	//This is an optional part of the header and not sure if it will always be appropriate. # channels in acquisition may not be # channels in image gadget receives if they are already combined
	//for now leaving cres (actual image channels) separate
	num_ch=hdr.acquisitionSystemInformation.get().receiverChannels();

	ordering.resize(num_echos*num_slices);
	slice_mean.resize(num_echos);

	if(numVol.value()!=1)//I want to believe there's a better way. Gadget Properties can be vectors but can't make it work.
	{
		std::string buffer;
		std::stringstream ss(limits.value());
		while(getline(ss, buffer, ',')){
		VolumeEnds.push_back(stof(buffer));
		}
	}
	else
		VolumeEnds.push_back(num_slices);	



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
	GadgetContainerMessage<hoNDArray<int>> *supportmasks_msg = AsContainerMessage<hoNDArray<int>>(m2->cont());
	GadgetContainerMessage<hoNDArray<int>> *masks_msg;
	GadgetContainerMessage<ISMRMRD::MetaContainer> *meta;

	if(!(m2 && supportmasks_msg)){//may not NEED support mask. only if 3D unwrapping is requested

		GERROR("Image and/or support mask missing from message.\n");
		return GADGET_FAIL;
	}
	ISMRMRD::Image<float> image;
	ISMRMRD::ImageHeader *header=m1->getObjectPtr();
	int fullsignal =0;
	float tmp_mean=0;
	int ch, e;

	//yres = m2->getObjectPtr()->get_size(0);	
	//xres = m2->getObjectPtr()->get_size(1);
	cres = header->channels;
	float *unwrapped_x = m2->getObjectPtr()->get_data_ptr();
	int *supportmasks = supportmasks_msg->getObjectPtr()->get_data_ptr();
	int *masks;

	float* im_ptr;
	//synthesize naming conventions between this and U2DG
	
	if(*supportmasks<0)//fullsignal masks
	{
		fullsignal=1;
		*supportmasks+=2;
		meta=  AsContainerMessage<ISMRMRD::MetaContainer>(supportmasks_msg->cont());

	}
	else
	{
		masks_msg = AsContainerMessage<hoNDArray<int>>(supportmasks_msg->cont());
		if(!masks_msg){
			GERROR("Something's wrong with this message chain\n");
			return GADGET_FAIL;
		}
	
		masks=masks_msg->getObjectPtr()->get_data_ptr();
		meta = AsContainerMessage<ISMRMRD::MetaContainer>(supportmasks_msg->cont()->cont());
	}
	////////////////
	for(ch=0; ch<cres; ch++)
	{
		tmp_mean=0;
		int k=0;
		for(int ii=xres*yres*ch; ii < xres*yres*(ch+1); ii++)
		{
				if(supportmasks[ii]==1)
				{
					
					tmp_mean+=unwrapped_x[ii];
					k++;
				}

		
		}
		tmp_mean/=k;
		for(int ii=xres*yres*ch; ii < xres*yres*(ch+1); ii++)
		{
		
				if(fullsignal || masks[ii]==1)//should shortcircuit with fullsignal (i.e. not try to check mask)
					unwrapped_x[ii]-=PI2*round(tmp_mean/PI2);
		
		}

		slice_mean[header->contrast][ch][header->image_series_index%num_slices]=tmp_mean; //will this always work?
	}
	//////////////
	image.setHead(*(header));
	im_ptr=image.getDataPtr();
	memcpy(im_ptr, unwrapped_x, image.getDataSize());
	



      std::stringstream attributes;
      
      if (meta) {
         ISMRMRD::serialize(*meta->getObjectPtr(), attributes);
	
      }
      image.setAttributeString(attributes.str());
	static int myplace=0;

	ordering[header->image_index-1]=myplace++;//is this the appropriate thing to be ordering by? Would position be better?
	temp_storage->appendImage("here", image);

	m1->release();
	
	for(int v=0; v<numVol.value(); v++)//better way than for loop? maybe remove below to another function?
	{
	if(myplace==VolumeEnds[v]*num_echos)//num_slices*num_echos)
	{
		std::chrono::high_resolution_clock::time_point t1;
		std::chrono::high_resolution_clock::time_point t2;
		int duration;
		std::vector<float> slice_mean_original;
		std::vector<std::vector<std::vector<float>>> offset(slice_mean); //just need the size but in 3D
		float* im_pointer, *data_ptr;

		t1 = std::chrono::high_resolution_clock::now();

		int start=v==0?0:VolumeEnds[v-1];
		
		int end =VolumeEnds[v];

		//parallelize this
		//for outlerloop is echos or channels better?
		for(e=0; e<num_echos; e++)
		for(ch=0;ch<cres;ch++)
		{
			slice_mean_original=std::vector<float>(slice_mean[e][ch].begin()+start, slice_mean[e][ch].begin()+end);//end-1?
			unwrap(slice_mean[e][ch].data()+start, end-start);
			tmp_mean=slice_mean[e][ch][start+(end-start)/2-1];
			for(int sl=start; sl<end; sl++)//number of slices
			{
				slice_mean[e][ch][sl]-=PI2*round(tmp_mean/PI2);
				if(fabs(slice_mean_original[sl-start]-slice_mean[e][ch][sl])>PI)
					offset[e][ch][sl]=slice_mean[e][ch][sl]-slice_mean_original[sl-start];
				else
					offset[e][ch][sl]=0;
			}
		}
		/////////////////////////////
 		//if images are highpass filtered, don't think this helps
		///////////////////

	
	
	
			
		for(int i=start; i<end; i++)//number of slices
		{
			for(int j=0; j<num_echos; j++)//number of echos/contrast
			{
				 // Read the image
				 try {
				    temp_storage->readImage("here", ordering[j*num_slices+i], image);///number of slices //ordering[num_slices*j+i]  <--necessary if thread >1
				 }
				 catch (std::exception &ex) {
				    GERROR("Error reading image %d\n",i);
				    return GADGET_FAIL;
				 }
				 GadgetContainerMessage<ISMRMRD::ImageHeader>* mheader = new GadgetContainerMessage<ISMRMRD::ImageHeader >();
				  header = mheader->getObjectPtr();
				 *header = image.getHead();
				 
				GadgetContainerMessage< hoNDArray<float> > *data = new GadgetContainerMessage< hoNDArray<float> >();
				
				data->getObjectPtr()->create(header->matrix_size[0], header->matrix_size[1],header->matrix_size[2], header->channels);
								
				data_ptr=data->getObjectPtr()->get_data_ptr();
				
				im_pointer=image.getDataPtr();

				////
				for(ch=0; ch<cres; ch++)//if 3D unwrapping is unnecessary this is a waste (as are the rest of the u3D specific parts)
				{			//this gadget is mostly intended to reorder slices that are out of order from threading or interleaved scans(?)
							//maybe there is a better way to brek up the two ideas
					for (int p = xres*yres*ch; p < xres*yres*(ch+1); p++) 
					{
						//GDEBUG("%d %d %d %d\n",p,j, ch, i);
						im_pointer[p]+=offset[j][ch][i];
						
					}
				}
				///
				memcpy(data->getObjectPtr()->get_data_ptr(), image.getDataPtr(), image.getDataSize());
				std::string attributes2;


				image.getAttributeString(attributes2);
				 if (header->attribute_string_len > 0 && !attributes2.empty()) {
				    GadgetContainerMessage<ISMRMRD::MetaContainer> *meta = new GadgetContainerMessage<ISMRMRD::MetaContainer>();
				    
				    try {
				       ISMRMRD::deserialize(image.getAttributeString(), *meta->getObjectPtr());
				    }
				    catch (std::exception &ex) {
				       GERROR("Error parsing attribute string: %s\n", ex.what());
				       return GADGET_FAIL;
				    }
				    data->cont(meta);
				 }
				 mheader->cont(data);


				if(savephase.value()==1 && do3D.value()==1)//otherwise this is available in the temp file (or 2D unwrap)				
				{
					dsToWrite->appendImage("3DMultiChannelUnwrappedPhase", image);
				}


				 
				 // Send the image along the chain
				//while(this->msg_queue()->message_count()>1){}
				 if (this->next()->putq(mheader) == -1) {
				    mheader->release();
				    GERROR("Unable to send image.\n");
				    return GADGET_FAIL;
				 }
				
			}
			if(i%10==0)
			GINFO("%d sent back out\n", i);
		}
		t2 = std::chrono::high_resolution_clock::now();
		 duration= std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
		GINFO("Total time (s) for 3D Unwrap:%d\n",duration/1000000);
	}
	}
	
return GADGET_OK;

}

GADGET_FACTORY_DECLARE(Unwrap3DGadget)

}
