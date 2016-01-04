//ReorderGadget.cpp
//Written by Tim Whelan 2016
//Input ImageHeader->Float 3D Array (Phase Data) -> [MetaContainer]
//Output ImageHeader->Float 3D Array (Phase Data)-> [MetaContainer]


#include "ReorderGadget.h"

namespace Gadgetron{

int ReorderGadget::process_config(ACE_Message_Block* mb)
{
	 boost::filesystem::remove("/tmp/temp.ismrmrd");
	temp_storage = new ISMRMRD::Dataset("/tmp/temp.ismrmrd","storage", true);
	this->msg_queue()->high_water_mark(128);//This helps with memory. It's not a hard limit though. 
	ISMRMRD::IsmrmrdHeader hdr;
        ISMRMRD::deserialize(mb->rd_ptr(),hdr);

	std::string id = "MxDoe";

	
	num_echos=hdr.encoding[0].encodingLimits.contrast().maximum +1; //number of echos is one more than highest numbers (0-based)
	num_slices=hdr.encoding[0].reconSpace.matrixSize.z; //number of slices (will this always work?)

	if(hdr.subjectInformation.is_present())
	{
		struct ISMRMRD::SubjectInformation sInfo=hdr.subjectInformation.get();
		if(sInfo.patientName.is_present())
			id= sInfo.patientName.get();//assuming this field is acutally a pseudoanonymous id
	}

	ordering.resize(num_echos*num_slices);

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

	return GADGET_OK;
}

int ReorderGadget::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1)
{
	GadgetContainerMessage<hoNDArray<float > > *m2 =AsContainerMessage<hoNDArray<float >> (m1->cont());
	GadgetContainerMessage<ISMRMRD::MetaContainer> *meta=AsContainerMessage<ISMRMRD::MetaContainer>(m2->cont());

	ISMRMRD::Image<float> image;
	ISMRMRD::ImageHeader *header=m1->getObjectPtr();
	int fullsignal =0;
	float tmp_mean=0;
	int ch, e;

	int cres = header->channels;
	float *unwrapped_x = m2->getObjectPtr()->get_data_ptr();

	float* im_ptr;
	image.setHead(*(header));
	im_ptr=image.getDataPtr();
	memcpy(im_ptr, unwrapped_x, image.getDataSize());
	
        std::stringstream attributes;
      
        if (meta) {
          ISMRMRD::serialize(*meta->getObjectPtr(), attributes);
	
        }
        image.setAttributeString(attributes.str());
	static int myplace=0;

	ordering[header->image_index]=myplace++;//would be better to organize by position
	temp_storage->appendImage("here", image);

	m1->release();
	
	for(int v=0; v<numVol.value(); v++)//better way than for loop? maybe remove below to another function?
	{
	if(myplace==VolumeEnds[v]*num_echos)//num_slices*num_echos)
	{
		std::chrono::high_resolution_clock::time_point t1;
		std::chrono::high_resolution_clock::time_point t2;
		int duration;
		float* im_pointer, *data_ptr;

		t1 = std::chrono::high_resolution_clock::now();

		int start=v==0?0:VolumeEnds[v-1];
		
		int end =VolumeEnds[v];

		for(int i=start; i<end; i++)//number of slices
		{
			for(int j=0; j<num_echos; j++)//number of echos/contrast
			{
				 // Read the image
				 try {
				    temp_storage->readImage("here", ordering[j*num_slices+i], image);///number of slices 				 
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

				 
				 // Send the image along the chain
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

GADGET_FACTORY_DECLARE(ReorderGadget)

}
