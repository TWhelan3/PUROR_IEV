//UnwrapGadget.cpp
//Written by Tim Whelan 2015
//Input ImageHeader->Float 3D Array (Phase Data) -> Int 3D Array (Support Mask) -> [Int 3D Array (Support Mask)] -> [MetaContainer]
//Output ImageHeader->Float 3D Array (Phase Data) -> Int 3D Array (Support Mask) -> [MetaContainer]
#include "UnwrapGadget.h"

namespace Gadgetron{

int UnwrapGadget::process_config(ACE_Message_Block* mb)
{
	this->msg_queue()->high_water_mark(128);//This helps with memory. It's not a hard limit though.
	ISMRMRD::IsmrmrdHeader hdr;
	ISMRMRD::deserialize(mb->rd_ptr(),hdr);

	std::string id = "MxDoe";

	yres=hdr.encoding[0].reconSpace.matrixSize.x; //match my (and MATLABs) unfortunate convention
	xres=hdr.encoding[0].reconSpace.matrixSize.y;
	num_echos=hdr.encoding[0].encodingLimits.contrast().maximum +1;

	if(hdr.subjectInformation.is_present())
	{
		struct ISMRMRD::SubjectInformation sInfo=hdr.subjectInformation.get();
		if(sInfo.patientName.is_present())
			id= sInfo.patientName.get();//assuming this field is acutally a pseudoanonymous id
	}

	boost::filesystem::path fname((id+filename.value()+".ismrmrd"));

	boost::filesystem::remove((id+filename.value()+".ismrmrd"));
	if(savephase.value()==1)
	{
		dsToWrite = new ISMRMRD::Dataset((id+filename.value()+".ismrmrd").c_str(),"images", true);
	}

	return GADGET_OK;
}
int UnwrapGadget::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1)
{
	static int myid = 0;

	GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2 =AsContainerMessage<hoNDArray<std::complex<float> > > (m1->cont());
	GadgetContainerMessage<hoNDArray<int>> *supportmask_msg = AsContainerMessage<hoNDArray<int>>(m2->cont());
	GadgetContainerMessage<hoNDArray<int>> *mask_msg;
	GadgetContainerMessage<ISMRMRD::MetaContainer> *meta;

	if (!(m2 && supportmask_msg)) {
	GDEBUG("Wrong datatypes coming in! Gadget requires header, complex array and mask data\n");
	return GADGET_FAIL;
	}

	int channel_index;

	num_ch = m1->getObjectPtr()->channels;

	GadgetContainerMessage<ISMRMRD::ImageHeader>* new_header_msg = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
	GadgetContainerMessage<hoNDArray< float > > *new_image_msg = new GadgetContainerMessage<hoNDArray< float > >();

	bool fullsignal;

	fullsignal=0;

	//Copy the header
	*new_header_msg->getObjectPtr() = *m1->getObjectPtr();//correct way to copy? based on extract gadget, since does a similar task
	new_header_msg->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_FLOAT;//GADGET_IMAGE_REAL_FLOAT;
	new_header_msg->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_PHASE;

	boost::shared_ptr< std::vector<size_t> > dims = m2->getObjectPtr()->get_dimensions();

	try{new_image_msg->getObjectPtr()->create(dims.get());}
	catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to create itohx in Unwrap Gadget");
		return GADGET_FAIL;
	}

	float* unwrapped_phase = new_image_msg->getObjectPtr()->get_data_ptr();
	std::complex<float>* src = m2->getObjectPtr()->get_data_ptr();
	std::vector<float> original_phase(xres*yres*num_ch);
	new_header_msg->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_FLOAT;//GADGET_IMAGE_REAL_FLOAT;
	new_header_msg->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_PHASE;

	if(*(supportmask_msg->getObjectPtr()->get_data_ptr())<0)//fullsignal masks
	{
		fullsignal=1;
		*(supportmask_msg->getObjectPtr()->get_data_ptr())+=2;
		meta= AsContainerMessage<ISMRMRD::MetaContainer>(supportmask_msg->cont());
		supportmask_msg->cont(NULL);//allow support mask to be released without losing meta
	}
	else
	{
		mask_msg = AsContainerMessage<hoNDArray<int>>(supportmask_msg->cont());
		meta = AsContainerMessage<ISMRMRD::MetaContainer>(mask_msg->cont());
		mask_msg->cont(NULL);//allow mask to be released without losing meta
	}

	if(!meta)
		meta = new GadgetContainerMessage<ISMRMRD::MetaContainer>();

	#pragma omp parallel for private(channel_index)
	for(channel_index=0; channel_index<num_ch; channel_index++)
	{
		MaskData md(yres,xres);

		int ch_offset= xres*yres*channel_index;
		int next_offset= xres*yres*(channel_index+1);

		int *ch_mask;
		int *ch_sppt_mask=supportmask_msg->getObjectPtr()->get_data_ptr()+ch_offset;

		if(fullsignal)//fullsignal masks
		{
			md.fullsignal=1;
		}
		else
		{
			ch_mask=mask_msg->getObjectPtr()->get_data_ptr()+ch_offset;
			memcpy(md.MASK.data(),ch_mask, yres*xres*sizeof(int));
		}

		memcpy(md.support_MASK.data(),ch_sppt_mask, yres*xres*sizeof(int));

		md.iniMask();
		md.iniMask_y();

		for (int i = xres*yres*channel_index; i < xres*yres*(channel_index+1); i++)
		{
			original_phase[i] = arg(src[i]);
		}

		std::vector<float> phase_x(original_phase.data()+ch_offset,original_phase.data()+next_offset); 
		std::vector<float> phase_y(original_phase.data()+ch_offset,original_phase.data()+next_offset); 
		//Do 1D unwrapping
		unwrap_rows(phase_x,xres,yres); 
		unwrap_columns(phase_y,xres,yres);

		//Connect 'good' segments withing rows and columns
		//Align centers of rows and columns
		if(fullsignal)
		{
			shift_to_mean_brain(phase_x, md, xres, yres, 1);
			shift_to_mean_brain(phase_y, md, yres, xres, 2);
		}
		else
		{
			shift_to_mean(phase_x, md, xres, yres, 1);
			shift_to_mean(phase_y, md, yres, xres, 2);
		}

		std::vector<float> quality_x(xres);
		std::vector<float> quality_y(yres);

		int xy_start_L, xy_start_R, xy_start_up, xy_start_dw;

		//Find highest quality strips in the two images
		calc_quality_y(phase_x, md.support_MASK, quality_y,xy_start_dw,xy_start_up, xres, yres);
		calc_quality_x(phase_y, md.support_MASK, quality_x,xy_start_L,xy_start_R, xres, yres);

		if ((xy_start_R - xy_start_L)>=(xy_start_up - xy_start_dw))//need to change to a ratio if want to use a bias, in B0NICE it is >=0.5*
		{
			//Unwrapped columns have greater quality

			center_y(phase_y, phase_x, md.support_MASK_trans,xy_start_L,xy_start_R,xres, yres);//Align unwrapped rows with best part columns

			calc_quality_y(phase_x, md.support_MASK,quality_y,xy_start_dw,xy_start_up,xres, yres); //Recalculate best part of unwrapped rows

			center_x(phase_x, phase_y,md.support_MASK,xy_start_dw,xy_start_up, xres, yres);//Align unwrapped columns with best part of unwrapped rows

		}
		else
		{
			//Unwrapped rows of greater quality
			center_x(phase_x, phase_y, md.support_MASK,xy_start_dw,xy_start_up, xres, yres);//Align unwrapped columns with best part of unwrapped rows

			calc_quality_x(phase_y, md.support_MASK, quality_x,xy_start_L,xy_start_R, xres, yres);//Recalculate best part of unwrapped columns

			center_y(phase_y, phase_x, md.support_MASK_trans,xy_start_L,xy_start_R, xres, yres);//Align unwrapped rows with best part columns

		}

		if(fullsignal)
			final_compare_brain(phase_x, phase_y, 6, xres, yres);
		else
		{
			final_compare(phase_x, phase_y, 6,md, xres, yres);
		}

		diff_x(phase_x, fullsignal,md, xres, yres);

		/*double tmp_mean=0;   //*Changed over time, used for 3D?
		int k=0;
		for(int ii=0; ii < xres*yres; ii++)
		{
			if(md.support_MASK[ii]==1)
			{

				tmp_mean+=phase_x[ii];
				k++;
			}
		}
		tmp_mean/=k;
		for(int ii=0; ii < xres*yres; ii++)
		{

				if(fullsignal || md.MASK[ii]==1)//should shortcircuit with fullsignal (i.e. not try to check mask)
					phase_x[ii]-=PI2*round(tmp_mean/PI2);

		}*/

		for (int i = 0; i < xres*yres; i++)
		{
			unwrapped_phase[i+ch_offset] = phase_x[i];
		}
	}

	new_header_msg->cont(new_image_msg);

	meta->getObjectPtr()->append(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_PHASE);
	meta->getObjectPtr()->append(GADGETRON_IMAGENUMBER,(long) m1->getObjectPtr()->image_index);

	new_image_msg->cont(meta);

	if(savephase.value()==1)
	{
		ISMRMRD::Image<float> image;

		image.setHead(*(new_header_msg->getObjectPtr()));

		memcpy(image.getDataPtr(), unwrapped_phase, image.getDataSize());
		std::stringstream attributes;

		if (meta){
		ISMRMRD::serialize(*meta->getObjectPtr(), attributes);
		}
		image.setAttributeString(attributes.str());

		dsToWrite->appendImage("2DMultiChannelUnwrappedPhase", image);
	}

	if (this->next()->putq(new_header_msg) == -1) {
		new_header_msg->release();
		GERROR("Unable to put images on next gadgets queue\n");
		return GADGET_FAIL;
	}

	if(myid%(num_echos*4)==0)//logging value, shows approximately how long to do four sets of echos
		GINFO("Unwrapped %d \n",new_header_msg->getObjectPtr()->image_index);
	++myid;
	//m2->cont(NULL);//necessary so masks don't get deleted and stay connected. Release removes links further down chain.
	m1->release();

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(UnwrapGadget)
}


