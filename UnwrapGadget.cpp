//UnwrapGadget.cpp
//Written by Tim Whelan 2015
//Input ImageHeader->Float 3D Array (Phase Data) -> Int 3D Array (Support Mask) -> [Int 3D Array (Support Mask)] -> [MetaContainer]
//Output ImageHeader->Float 3D Array (Phase Data) -> Int 3D Array (Support Mask) -> [Int 3D Array (Support Mask)] -> [MetaContainer]
#include "UnwrapGadget.h"

namespace Gadgetron{

int UnwrapGadget::process_config(ACE_Message_Block* mb)
{

	this->msg_queue()->high_water_mark(128);//This helps with memory. It's not a hard limit though. 

	return GADGET_OK;
}
int UnwrapGadget::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1)
{
	//Gets the header, complex data, maskdata
	//Passes on the header, unwrapped x
	static int myid = 0;
	
	
	GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2 =AsContainerMessage<hoNDArray<std::complex<float> > > (m1->cont());
	 GadgetContainerMessage<hoNDArray<int>> *supportmasks = AsContainerMessage<hoNDArray<int>>(m2->cont());
	GadgetContainerMessage<hoNDArray<int>> *masks;
	GadgetContainerMessage<ISMRMRD::MetaContainer> *meta;
	if ((m1 && m2 && supportmasks)==0) {
	GDEBUG("Wrong datatypes coming in! Gadget requires header, complex array and mask data\n");
	return -1;
	}

	int channel_index;	
	yres = m2->getObjectPtr()->get_size(0);	
	xres = m2->getObjectPtr()->get_size(1);
	cres = m1->getObjectPtr()->channels;

	GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
	GadgetContainerMessage<hoNDArray< float > > *cm2 = new GadgetContainerMessage<hoNDArray< float > >();
	hoNDArray< float > *phase_x_block = cm2->getObjectPtr();	
	hoNDArray< float > *phase_y_block = new hoNDArray< float >;


	//std::chrono::high_resolution_clock::time_point t1;
	//std::chrono::high_resolution_clock::time_point t2;
	//int duration;
	
	bool fullsignal;


	fullsignal=0;

	//GINFO("m1 count %d m2 count %d support mask count %d\n",m1->reference_count(), m2->reference_count(), supportmasks->reference_count());
	//Copy the header
	*cm1->getObjectPtr() = *m1->getObjectPtr();//correct way to copy? based on extract gadget, since does a similar task
	cm1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_FLOAT;//GADGET_IMAGE_REAL_FLOAT;
	cm1->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_PHASE;

	boost::shared_ptr< std::vector<size_t> > dims = m2->getObjectPtr()->get_dimensions();
		//////

	//if(myid%106==0)
		//t1 = std::chrono::high_resolution_clock::now();
		//////
	try{phase_x_block->create(dims.get());}
	catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to create itohx in Unwrap Gadget");
		return GADGET_FAIL;
	}


	try{phase_y_block->create(yres*xres*cres);}
	catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to create itohy in Unwrap Gadget");
		return GADGET_FAIL;
	}

	std::complex<float>* src = m2->getObjectPtr()->get_data_ptr();
	std::vector<float>  dst(phase_x_block->get_number_of_elements());
	cm1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_FLOAT;//GADGET_IMAGE_REAL_FLOAT;
	cm1->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_PHASE;


	if(*(supportmasks->getObjectPtr()->get_data_ptr())<0)//fullsignal masks
	{
		fullsignal=1;
		*(supportmasks->getObjectPtr()->get_data_ptr())+=2;
		meta=  AsContainerMessage<ISMRMRD::MetaContainer>(supportmasks->cont());
	}
	else
	{
		masks = AsContainerMessage<hoNDArray<int>>(supportmasks->cont());
		meta = AsContainerMessage<ISMRMRD::MetaContainer>(masks->cont());
	}




	#pragma omp parallel for private(channel_index)	
	for(channel_index=0; channel_index<cres; channel_index++)
	{
		MaskData* md= new MaskData(yres,xres);

		int ch_offset= xres*yres*channel_index;
		float* phase_x, *phase_y;
		int *ch_mask;
		int *ch_sppt_mask=supportmasks->getObjectPtr()->get_data_ptr()+ch_offset;
		
		
	
		if(fullsignal)//fullsignal masks
		{
			md->fullsignal=1;
		}
		else
		{
			ch_mask=masks->getObjectPtr()->get_data_ptr()+ch_offset;
			memcpy(md->MASK,ch_mask, yres*xres*sizeof(int));
		}
	
			memcpy(md->support_MASK,ch_sppt_mask, yres*xres*sizeof(int));
		
		md->iniMask();
		md->iniMask_y();
	
		for (int i = xres*yres*channel_index; i < xres*yres*(channel_index+1); i++) 
		{
		dst[i] = arg(src[i]);
		}

		phase_x=phase_x_block->get_data_ptr()+ch_offset;
		phase_y=phase_y_block->get_data_ptr()+ch_offset;
		
		//Do 1D unwrapping
		unwrap_rows(dst.data()+ch_offset,phase_x); //was dst+ch_offset when dst was a 'c style array'
		unwrap_columns(dst.data()+ch_offset,phase_y); 
		
		//Connect 'good' segments withing rows and columns
		//Align centers of rows and columns 
		if(fullsignal)
		{
			shift_to_mean_brain(phase_x, md);
			shift_to_mean_y_brain(phase_y,md);
		}
		else
		{
			
			shift_to_mean(phase_x, md);
			shift_to_mean_y(phase_y,md);
		}
		
		std::vector<float> quality_x(xres);
		std::vector<float> quality_y(yres);;
		int xy_start_L, xy_start_R, xy_start_up, xy_start_dw;
		int* mask=md->support_MASK;
		int* t_mask=md->support_MASK_trans;

		//Find highest quality strips in the two images
		calc_quality_y(phase_x, mask, quality_y,xy_start_dw,xy_start_up);
		calc_quality_x(phase_y, mask, quality_x,xy_start_L,xy_start_R);
		//if(myid==0)
		//for (int i = 0; i < xres; i++) 
		//std::cout<<quality_x[i]<<" ";

		if ((xy_start_R - xy_start_L)>=(xy_start_up - xy_start_dw))//need to change to a ratio if want to use a bias
		{	
			//Unwrapped columns have greater quality

			
			center_y(phase_y, phase_x, t_mask,xy_start_L,xy_start_R);//Align unwrapped rows with best part columns
			
			
			calc_quality_y(phase_x, mask,quality_y,xy_start_dw,xy_start_up); //Recalculate best part of unwrapped rows
			
			center_x(phase_x, phase_y,mask,xy_start_dw,xy_start_up);//Align unwrapped columns with best part of unwrapped rows
			
		}
		else
		{
			//Unwrapped rows of greater quality
			center_x(phase_x, phase_y, mask,xy_start_dw,xy_start_up);//Align unwrapped columns with best part of unwrapped rows
			
			calc_quality_x(phase_y, mask, quality_x,xy_start_L,xy_start_R);//Recalculate best part of unwrapped columns
			
			center_y(phase_y, phase_x, t_mask,xy_start_L,xy_start_R);//Align unwrapped rows with best part columns
		
		}

		if(fullsignal)
			final_compare_brain(phase_x, phase_y, 6);
		else
		{
			final_compare(phase_x, phase_y, 6,md);
		}
		
		diff_x(phase_x, fullsignal,md);




		delete md;

	
	}
	//GINFO("m1 count %d m2 count %d support mask count %d\n",m1->reference_count(), m2->reference_count(), supportmasks->reference_count());

	delete phase_y_block;
	//need to deal with meta
	cm1->cont(cm2);
	//GadgetContainerMessage<hoNDArray<int>> *supportmasks_dup = supportmasks->duplicate();
	//GINFO("m1 count %d m2 count %d support mask count %d mask count %d\n",cm1->reference_count(), cm2->reference_count(), cm2->cont()->reference_count(),masks->reference_count());
	cm2->cont(supportmasks);
	//supportmasks->duplicate();
	//meta->duplicate();
	//masks->duplicate();
	if(fullsignal)
	{
		supportmasks->getObjectPtr()->get_data_ptr()[0]-=2;
	}

	//GINFO("m1 count %d m2 count %d support mask count %d mask count %d\n",cm1->reference_count(), cm2->reference_count(), supportmasks->reference_count(),masks->reference_count());
	if (this->next()->putq(cm1) == -1) {
		cm1->release();
		GDEBUG("Unable to put images on next gadgets queue\n");
		return GADGET_FAIL;
	}
		/*if(myid%106==0)
		{
		t2 = std::chrono::high_resolution_clock::now();
		 duration= std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
		//GINFO("~Total for %d: %d\n",myid, duration);
		//GDEBUG("%d Unwrapped Queue Length  %d HWM\n", this->msg_queue()->message_count(), this->msg_queue()->high_water_mark());
		}*/
	
	if(myid%24==0)//debugging value, shows approximately how long to do four sets of six echos
	GINFO("Unwrapped %d \n",cm1->getObjectPtr()->image_index);
	++myid;
	m2->cont(NULL);//necessary so masks don't get deleted and stay connected. Release removes links further down chain.
	m1->release(); 
	//GINFO("m1 count %d m2 count %d support mask count %d\n",m1->reference_count(), m2->reference_count(), supportmasks->reference_count());



	return GADGET_OK;
}

void  UnwrapGadget::unwrap_rows(float* toUnwrap, float* output)
{
	std::vector<float> phase_tmp(xres);
	int i,j;
	for(i=0; i<yres; i++)//Unwrap each row, alternating directions
	{

		if(i%2==0)
		{
			for(j=0; j<xres; j++)
				phase_tmp[xres-1-j]=toUnwrap[j*yres+i];

			unwrap(phase_tmp.data(), xres);

			for(j=0; j<xres; j++)
				output[j*yres+i]=phase_tmp[xres-1-j];

		}
		else
		{
			for(j=0; j<xres; j++)
				phase_tmp[j]=toUnwrap[j*yres+i];

			unwrap(phase_tmp.data(), xres);

			for(j=0; j<xres; j++)
				output[j*yres+i]=phase_tmp[j];
		}	

	}
}
void UnwrapGadget::unwrap_columns(float* toUnwrap, float* output)
{	
	std::vector<float> phase_tmp(yres);
	
	int i,j;
	for(i=0; i<xres; i++)//Unwrap each column alternating directions
	{
		if(i%2==0)
		{
			
			for(j=0; j<yres; j++)
			{
				phase_tmp[yres-j-1]=toUnwrap[i*yres+j];
	
			}
			unwrap(phase_tmp.data(),yres);
			for(j=0; j<yres; j++)
				output[i*yres+j]=phase_tmp[yres-j-1];
		
		}
		else
		{
			
			for(j=0; j<yres; j++)
				phase_tmp[j]=toUnwrap[i*yres+j];
			
			//memcpy(phase_tmp.data(), toUnwrap+i*yres, yres*sizeof(float));
	
			unwrap(phase_tmp.data(),yres);

			//memcpy(output+i*yres, phase_tmp.data(), yres*sizeof(float));
			for(j=0; j<yres; j++)
				output[i*yres+j]=phase_tmp[j];
		
		}	

	}
}

void UnwrapGadget::shift_to_mean(float* phase, MaskData* md)
{
	std::vector<float> phase_sample(xres);
	std::vector<float> phase_u(xres);
	std::vector<float> mean_connect(yres);
	std::vector<float> mean_unwrap(yres);
	int index_y, ii;
	std::vector<int> signal;
	std::vector<int> connect;
	float diff_test, correction, phi_good;
	int len_u, len_c;
	int midline=round(yres/2)-1;
	std::vector<int> good_segs;
	std::vector<int> index_ls;
	for(index_y=0; index_y<yres; index_y++)
	{
		
		good_segs=md->segX[index_y];

		len_u=md->signalX[index_y].size();
		signal = md->signalX[index_y];

		if(len_u>(xres/4))//if mask covers more than 1/4 of the row
		{
			for(ii=0; ii<len_u; ii++)
				phase_u[ii]=phase[signal[ii]*yres+index_y];

			phi_good=PI/2;
			index_ls = dePULM_1D_ls(phase_u,phi_good ,good_segs);//Find segments in phase_u with phase gradient < phi_good
		
			dePULM_1D(phase_u, signal,index_ls);//Reconnect these segments

			for(ii=0; ii<len_u; ii++)
				phase[signal[ii]*yres+index_y]=phase_u[ii];
		}

	}
	
	for(index_y=1; index_y<yres; index_y++)
	{
		//calculate connect mean
		len_u	=	md->signalX[index_y].size();
		len_c	=	md->connectXH[index_y].size();
		signal = 	md->signalX[index_y];
		connect =	md->connectXH[index_y];
		

		if(len_u==0)
		{
			mean_connect[index_y]=mean_connect[index_y-1];
		}
		else
		{
			//Find mean of each row
			if(len_c<=3)
			{
				for(ii=0; ii<len_u; ii++)
					phase_sample[ii]=phase[signal[ii]*yres+index_y];
				mean_connect[index_y]=mean(phase_sample.data(), len_u);
			}

			else
			{
				for(ii=0; ii<len_c; ii++)
					phase_sample[ii]=phase[connect[ii]*yres+index_y];
				mean_connect[index_y]=mean(phase_sample.data(), len_c);
			}

			//bring the average of the row to [-pi,pi] 
			correction=PI2*round(mean_connect[index_y]/PI2);

			for(ii=0; ii<len_u; ii++)
				phase[signal[ii]*yres+index_y]-=correction;

			if(len_c<=3)
			{
				for(ii=0;ii<len_u;ii++)
					phase_sample[ii]=phase[signal[ii]*yres+index_y];
				mean_connect[index_y]=mean(phase_sample.data(), len_u);

			}
			else
			{
				for(ii=0; ii<len_c; ii++)
					phase_sample[ii]=phase[connect[ii]*yres+index_y];
				mean_connect[index_y]=mean(phase_sample.data(), len_c);
			}	
		}

	}

	index_y=0;//This is special case of loop above
	len_u=md->signalX[0].size();
	len_c=md->connectXH[0].size();
	signal = md->signalX[0];
	connect = md->connectXH[0];
	if(len_u==0)
		mean_connect[index_y]=mean_connect[index_y+1];
	else
	{
		if(len_c<=3)
		{
			for(ii=0; ii<len_u; ii++)
				phase_sample[ii]=phase[signal[ii]*yres];
			mean_connect[index_y]=mean(phase_sample.data(), len_u);
		}
		else
		{
			for(ii=0; ii<len_c; ii++)
				phase_sample[ii]=phase[connect[ii]*yres];
			mean_connect[index_y]=mean(phase_sample.data(), len_c);
		}
		correction=PI2*round(mean_connect[index_y]/PI2);
		for(ii=0; ii<len_u; ii++)
			phase[signal[ii]*yres]-=correction;
		if(len_c<=3)
		{			for(ii=0; ii<len_u; ii++)
				phase_sample[ii]=phase[signal[ii]*yres];
			mean_connect[index_y]=mean(phase_sample.data(), len_u);
		}
		else
		{
			for(ii=0; ii<len_c; ii++)
				phase_sample[ii]=phase[connect[ii]*yres];
			mean_connect[index_y]=mean(phase_sample.data(), len_c);
		}	
	}

	mean_unwrap=mean_connect;
	//Actually want a copy here so mean_connect can be used as reference

	//unwrap the mean values before global shifting the data
	unwrap(mean_unwrap.data(), yres);
	correction=round(mean_unwrap[midline]/PI2)*PI2;
	for(ii=0; ii<yres; ii++)
		mean_unwrap[ii]-=correction;

	//shift phase data
	for(index_y = 0; index_y<yres; index_y++)
	{
				diff_test = mean_unwrap[index_y] - mean_connect[index_y];
		
		if (fabs(diff_test) > PI)
		{
			correction=PI2*round(diff_test/(PI2));

			len_u=md->signalX[index_y].size();
			signal = md->signalX[index_y];

			for(ii=0; ii<len_u; ii++)
				phase[signal[ii]*yres+index_y] +=correction;
		}
		
	}

}


void UnwrapGadget::shift_to_mean_y(float* phase,MaskData* md)
{
	int index_x, ii;
	float diff_test, correction, phi_good;
	int length2, len_u, len_c;

	std::vector<float> phase_sample(yres);
	std::vector<float> phase_u(yres);
	std::vector<float> mean_connect(xres);
	std::vector<float> mean_unwrap(xres);
	std::vector<int> good_segs;
	std::vector<int> index_ls;

	std::vector<int> signal;
	std::vector<int> connect;
	int midline=round(xres/2)-1;
	int columnOffset;
	for(index_x=0; index_x<xres; index_x++)
	{
		
		good_segs=md->segY[index_x];
		len_u=md->signalY[index_x].size();
		columnOffset=index_x*yres;
		signal = md->signalY[index_x];
		if(len_u>yres/4)//if mask covers more than 1/4 of the column		{
			for(ii=0; ii<len_u; ii++)
				phase_u[ii]=phase[columnOffset+signal[ii]];
			phi_good=PI/2;
			index_ls = dePULM_1D_ls(phase_u,phi_good,good_segs);//Find segments in phase_u with phase gradient < phi_good
		
			dePULM_1D(phase_u,signal,index_ls);//Reconnect these segments
			for(ii=0; ii<len_u; ii++)				phase[columnOffset+signal[ii]]=phase_u[ii];

		}		//took out else	}
	
	for(index_x=1; index_x<xres; index_x++)
	{		len_u=		md->signalY[index_x].size();
		len_c	=	md->connectYH[index_x].size();
		signal = 	md->signalY[index_x];
		connect =	md->connectYH[index_x];
		columnOffset=index_x*yres;
		//calculate the connect mean
		if(len_u==0)
		{
			mean_connect[index_x]=mean_connect[index_x-1];
		}
		else
		{			//Find mean of each column
			if(len_c<=3)
			{
				for(ii=0; ii<len_u; ii++)
					phase_sample[ii]=phase[columnOffset+signal[ii]];
				mean_connect[index_x]=mean(phase_sample.data(), len_u);
			}
			else
			{
				for(ii=0; ii<len_c; ii++)
					phase_sample[ii]=phase[columnOffset+connect[ii]];
				mean_connect[index_x]=mean(phase_sample.data(), len_c);
			}
			//bring the average of the row to [-pi,pi] 
			correction=round(mean_connect[index_x]/PI2)*PI2;
			for(ii=0; ii<len_u; ii++)				phase[columnOffset+signal[ii]]-=correction;
			//find new average of row			if(len_c<=3)
			{
				for(ii=0; ii<len_u; ii++)
					phase_sample[ii]=phase[columnOffset+signal[ii]];
				mean_connect[index_x]=mean(phase_sample.data(), len_u);			}
			else
			{
				for(ii=0; ii<len_c; ii++)
					phase_sample[ii]=phase[columnOffset+connect[ii]];
				mean_connect[index_x]=mean(phase_sample.data(), len_c);
			}
		}	}
	index_x=0;//This is special case of loop above
	len_u=		md->signalY[0].size();
	len_c	=	md->connectYH[0].size();
	signal = md->signalY[0];
	connect = md->connectYH[0];
	columnOffset=index_x*yres;
	if(len_u==0)
	{
		mean_connect[index_x]=mean_connect[index_x+1];
	}	else
	{
		//Find mean of each column
		if(len_c<=3)
		{
			for(ii=0; ii<len_u; ii++)
				phase_sample[ii]=phase[signal[ii]];			mean_connect[index_x]=mean(phase_sample.data(), len_u);
		}
		else
		{
			for(ii=0; ii<len_c; ii++)
				phase_sample[ii]=phase[connect[ii]];
			mean_connect[index_x]=mean(phase_sample.data(), len_c);
		}
		//bring the average of the row to [-pi,pi] 
		correction=round(mean_connect[index_x]/PI2)*PI2;
		for(ii=0; ii<len_u; ii++)
			phase[signal[ii]]-=correction;
		//find new average of row
		if(len_c<=3)
		{			for(ii=0; ii<len_u; ii++)
				phase_sample[ii]=phase[signal[ii]];
			mean_connect[index_x]=mean(phase_sample.data(), len_u);
		}
		else
		{
			for(ii=0; ii<len_c; ii++)
				phase_sample[ii]=phase[connect[ii]];
			mean_connect[index_x]=mean(phase_sample.data(), len_c);
		}
	}

	mean_unwrap=mean_connect;
	
	unwrap(mean_unwrap.data(), xres);
	
	correction=round(mean_unwrap[midline]/PI2)*PI2;
	for(ii=0; ii<xres; ii++)
	{
		mean_unwrap[ii]-=correction;
	}
	//shift phase data	for(index_x = 0; index_x<xres; index_x++)
	{
		diff_test = mean_unwrap[index_x] - mean_connect[index_x];
		signal = md->signalY[0];
		if (fabs(diff_test) > PI)
		{			correction=PI2*round(diff_test/(PI2));
			for(ii=0; ii<md->signalY[index_x].size(); ii++)
				phase[index_x*yres+signal[ii]] += correction;
		}
	}
}
void UnwrapGadget::shift_to_mean_brain(float* phase, MaskData* md)
{

	std::vector<float> phase_sample(xres);
	std::vector<float> mean_connect(yres);
	std::vector<float> mean_unwrap(yres);
	std::vector<int> index_ls;

	int index_y, ii;
	std::vector<int> signal;
	std::vector<int> connect;
	float diff_test, correction, phi_good;
	int len_c;
	int midline=round(yres/2)-1;
	
	for(index_y=0; index_y<yres; index_y++)
	{
		for(ii=0; ii<xres; ii++)
			phase_sample[ii]=phase[ii*yres+index_y];

		phi_good=PI/2;
		index_ls = dePULM_1D_ls_brain(phase_sample,phi_good);//Find segments in phase_sample with phase gradient < phi_good
	
		dePULM_1D_brain(phase_sample,index_ls);//Reconnect these segments

		for(ii=0; ii<xres; ii++)
				phase[ii*yres+index_y]=phase_sample[ii];
	

	}
	
	for(index_y=0; index_y<yres; index_y++)
	{
		//calculate connect mean
		len_c	=	md->connectXH[index_y].size();
		connect =	md->connectXH[index_y];

		//Find mean of each row
		if(len_c<=3)
		{
			for(ii=0; ii<xres; ii++)
				phase_sample[ii]=phase[ii*yres+index_y];
			mean_connect[index_y]=mean(phase_sample.data(), xres);
		}
		else
		{
			for(ii=0; ii<len_c; ii++)
				phase_sample[ii]=phase[connect[ii]*yres+index_y];

			mean_connect[index_y]=mean(phase_sample.data(), len_c);
		}

		//bring the average of the row to [-pi,pi] 
		correction=PI2*round(mean_connect[index_y]/PI2);

		for(ii=0; ii<xres; ii++)
			phase[ii*yres+index_y]-=correction;

		if(len_c<=3)
		{
			for(ii=0;ii<xres;ii++)
				phase_sample[ii]=phase[ii*yres+index_y];
			mean_connect[index_y]=mean(phase_sample.data(), xres);

		}
		else
		{
			for(ii=0; ii<len_c; ii++)
				phase_sample[ii]=phase[(int)(connect[ii])*yres+index_y];
			mean_connect[index_y]=mean(phase_sample.data(), len_c);
		}	
	

	}


	mean_unwrap=mean_connect;//Want a copy here so mean_connect can be used as reference

	//unwrap the mean values before global shifting the data
	unwrap(mean_unwrap.data(), yres);
	correction=round(mean_unwrap[midline]/PI2)*PI2;
	for(ii=0; ii<yres; ii++)
		mean_unwrap[ii]-=correction;

	//shift phase data
	for(index_y = 0; index_y<yres; index_y++)
	{
				diff_test = mean_unwrap[index_y] - mean_connect[index_y];

		if (fabs(diff_test) > PI)
		{
			correction=PI2*round(diff_test/(PI2));

			for(ii=0; ii<xres; ii++)
				phase[ii*yres+index_y] +=correction;
		}
		
	}

}
void UnwrapGadget::shift_to_mean_y_brain(float* phase,MaskData* md)
{


	int index_x, ii;
	std::vector<int> connect;
	float diff_test, correction, phi_good;
	int len_c;

	std::vector<float> phase_sample(yres);
	std::vector<float> mean_connect(xres);
	std::vector<float> mean_unwrap(xres);
	std::vector<int> index_ls;
	int midline=round(xres/2)-1;
	int columnOffset;
	for(index_x=0; index_x<xres; index_x++)
	{
		columnOffset=index_x*yres;
		for(ii=0; ii<yres; ii++)
			phase_sample[ii]=phase[ii+columnOffset];
		
		phi_good=PI/2;
		index_ls = dePULM_1D_ls_brain(phase_sample,phi_good);//Find segments in phase with phase gradient < phi_good
		dePULM_1D_brain(phase_sample,index_ls);//Reconnect these segments
		
		for(ii=0; ii<yres; ii++)
			phase[ii+columnOffset]=phase_sample[ii];



	}
	
	for(index_x=0; index_x<xres; index_x++)//should start at 0?
	{		len_c	=	md->connectYH[index_x].size();
		connect =	md->connectYH[index_x];
		columnOffset=index_x*yres;
		//calculate the connect mean
				//Find mean of each column
		if(len_c<=3)
		{
			mean_connect[index_x]=mean(phase+columnOffset, yres);
		}
		else
		{
			for(ii=0; ii<len_c; ii++)
				phase_sample[ii]=phase[columnOffset+connect[ii]];
			mean_connect[index_x]=mean(phase_sample.data(), len_c);
		}
		//bring the average of the row to [-pi,pi] 
		correction=round(mean_connect[index_x]/PI2)*PI2;
		for(ii=0; ii<yres; ii++)			phase[columnOffset+ii]-=correction;
		//find new average of row		if(len_c<=3)
		{
			mean_connect[index_x]=mean(phase+columnOffset, yres);		}
		else
		{
			for(ii=0; ii<len_c; ii++)
				phase_sample[ii]=phase[columnOffset+connect[ii]];
			mean_connect[index_x]=mean(phase_sample.data(), len_c);
		}
			}


	mean_unwrap=mean_connect;

	unwrap(mean_unwrap.data(), xres);
	correction=round(mean_unwrap[midline]/PI2)*PI2;
	for(ii=0; ii<xres; ii++)
	{
		mean_unwrap[ii]-=correction;
	}
	//shift phase data	for(index_x = 0; index_x<xres; index_x++)
	{
		diff_test = mean_unwrap[index_x] - mean_connect[index_x];
		if (fabs(diff_test) > PI)
		{			correction=PI2*round(diff_test/(PI2));
			for(ii=0; ii<yres; ii++)
				phase[index_x*yres+ii] += correction;
		}
	}
}
void UnwrapGadget::calc_quality_y(float* phase_x, int *mask, std::vector<float> &quality_y, int &xy_start_dw, int &xy_start_up)//eqv to matlabs quality_ch
{

	int row_index;
	int len_dw, len_up;

	int ii,jj;
	int flag_start, flag_end;
	int g_start, g_end;

	std::vector<int> pointsToUse;
	
	for(row_index=0; row_index<yres; row_index++)
	{
		len_dw=0; len_up=0;
		if(row_index!=0 && row_index!=(yres-1))
		{
			for(ii=0;ii<xres; ii++)//find spots in t_mask which has a point above and below to check quality
			{
				//if(t_mask[(row_index+1)*xres+ii] && t_mask[(row_index-1)*xres+ii] && t_mask[row_index*xres+ii])
				if(mask[row_index+1+yres*ii] && mask[row_index-1+yres*ii] && mask[row_index+yres*ii])
					pointsToUse.push_back(ii);
			}
			for(ii=0; ii<pointsToUse.size(); ii++)
			{
				//count number of jumps between rows
			
				if((fabs(phase_x[yres*pointsToUse[ii]+row_index]-phase_x[yres*pointsToUse[ii]+row_index+1]))>PI)
					len_dw++;
			}
			for(ii=0; ii<pointsToUse.size(); ii++)
			{
				//count number of jumps between rows
				if((fabs(phase_x[yres*pointsToUse[ii]+row_index]-phase_x[yres*pointsToUse[ii]+row_index-1]))>PI)
					len_up++;
			}
			if(!pointsToUse.empty())	//find quality of row			
			{
				quality_y[row_index]=(float)(pointsToUse.size()-len_up)*(float)(pointsToUse.size()-len_dw)/(float)(pointsToUse.size()*pointsToUse.size());
			}
			else
			{
				quality_y[row_index]=-1;

			}
		}	
		if(row_index==0)//special case
		{
			len_dw=0; 
			for(ii=0;ii<xres; ii++)
			{
				//if(t_mask[(row_index+1)*xres+ii] && t_mask[row_index*xres+ii])
				if(mask[row_index+1+yres*ii] && mask[row_index+yres*ii])
					pointsToUse.push_back(ii);
			}
		
			for(ii=0; ii<pointsToUse.size(); ii++)
			{
				if(fabs(phase_x[yres*pointsToUse[ii]] -  phase_x[yres*pointsToUse[ii]+1])>PI);
				len_dw++;
				//std::cout<<pointsToUse[ii]<<" ";
			}
			//std::cout<<pointsToUse.size()<<" "<<len_dw<<" "<<pointsToUse.size()<<"     ";
			//std::cout<<(pointsToUse.size()-len_dw)/(pointsToUse.size())<<"     ";
			if(!pointsToUse.empty())				
				quality_y[row_index]=(float)(pointsToUse.size()-len_dw)/((float)pointsToUse.size());
			else
				quality_y[row_index]=-1;
			//std::cout<<quality_y[row_index]<<" ";
		}
		if(row_index==(yres-1))//special case
		{
			for(ii=0;ii<xres; ii++)
			{
				if(mask[row_index-1+yres*ii] && mask[row_index+yres*ii])
					pointsToUse.push_back(ii);
				//if(t_mask[(row_index-1)*xres+ii] && t_mask[row_index*xres+ii])
					
			}
			for(ii=0; ii<pointsToUse.size(); ii++)
			{
				//count number of jumps between rows
				if((fabs(phase_x[yres*pointsToUse[ii]+row_index]-phase_x[yres*pointsToUse[ii]+row_index-1]))>PI)
					len_up++;
			}
			
			if(!pointsToUse.empty())			
				quality_y[row_index]=(float)(pointsToUse.size()-len_up)/(pointsToUse.size());
			else
				quality_y[row_index]=-1;

		}

		pointsToUse.clear();
		
	}



	xy_start_dw=0;
	xy_start_up=0;
	flag_start=1;
	flag_end=2;
	q_th=0.9;

	for(ii=0; ii<yres; ii++)//find highest quality strip
	{
		if(quality_y[ii]>q_th && flag_start==1)//quality above threshold, mark start
		{
			g_start=ii;
			flag_start=0;
			flag_end=1;
		}
		if(quality_y[ii]<q_th && flag_end ==1)//quality below threshold, if started, mark end
		{
			g_end=ii-1;
			flag_start=1;
			flag_end=0;
		}
		if(ii==(yres-1) && flag_end ==1)//last row, if started, save end
		{
			if (quality_y[ii]>=q_th)
				g_end=ii;
			else
				g_end=g_start;
			flag_end=0;
		}
		if(flag_end==0 && ii>0)//check if this strip is bigger than earlier strips
			if((g_end-g_start)>=(xy_start_up-xy_start_dw))
			{
				xy_start_dw=g_start;
				xy_start_up=g_end;
			}
	}





}
void UnwrapGadget::calc_quality_x(float* phase_y, int *mask, std::vector<float> &quality_x,int &xy_start_L, int &xy_start_R)//eqv to matlabs quality_ch_y
{
	int ii,jj;
	int flag_start, flag_end;
	int g_start, g_end;
	int len_dw, len_up;

	int col_index;
	std::vector<int> pointsToUse;
	for(col_index=0; col_index<xres; col_index++)
	{
		len_dw=0; len_up=0;
		if(col_index!=0 && col_index!=(xres-1))
		{
			for(ii=0;ii<yres; ii++)//find spots in mask which has a point left and right to check quality
			{
				if(mask[(col_index+1)*yres+ii] && mask[(col_index-1)*yres+ii] && mask[col_index*yres+ii])
					pointsToUse.push_back(ii);
			}
			for(ii=0; ii<pointsToUse.size(); ii++)
			{		
				//count jumps between columns
				if(fabs(phase_y[pointsToUse[ii]+yres*col_index] -  phase_y[pointsToUse[ii]+yres*(col_index+1)])>PI)
					len_dw++;
				if(fabs(phase_y[pointsToUse[ii]+yres*col_index] - phase_y[pointsToUse[ii]+yres*(col_index-1)])>PI)
					len_up++;
			}

			if(pointsToUse.size()>0)//get quality of column
				quality_x[col_index]=(float)(pointsToUse.size()-len_up)*(pointsToUse.size()-len_dw)/(pointsToUse.size()*pointsToUse.size());


			else

				quality_x[col_index]=-1;

		}	
		if(col_index==0)//special case of above loop
		{
			len_dw=0;
			for(ii=0;ii<yres; ii++)
			{
				if(mask[(col_index+1)*yres+ii] && mask[col_index*yres+ii])
					pointsToUse.push_back(ii);
			}
			for(ii=0; ii<pointsToUse.size(); ii++)
			{
				if(fabs(phase_y[pointsToUse[ii]+yres*col_index] -  phase_y[pointsToUse[ii]+yres*(col_index+1)])>PI)
				len_dw++;
			}

			if(!pointsToUse.empty())				
				quality_x[col_index]=(pointsToUse.size()-len_dw)/(float)(pointsToUse.size());
			else
				quality_x[col_index]=-1;

		}
		if(col_index==(xres-1))//special case of above loop
		{
			len_up=0;
			for(ii=0;ii<yres; ii++)
			{
				if(mask[(col_index-1)*yres+ii] && mask[col_index*yres+ii])
					pointsToUse.push_back(ii);
			}
			for(ii=0; ii<pointsToUse.size(); ii++)
			{
				for(ii=0; ii<pointsToUse.size(); ii++)
				{
				if(fabs(phase_y[pointsToUse[ii]+yres*col_index] -  phase_y[pointsToUse[ii]+yres*(col_index-1)])>PI)
				len_up++;
				}

				if(!pointsToUse.empty())				
					quality_x[col_index]=((float)pointsToUse.size()-len_up)/(pointsToUse.size());
				else
					quality_x[col_index]=-1;
			}
		}
		pointsToUse.clear();
	
	}
		xy_start_L=0;
		xy_start_R=0;
		flag_start=1;
		flag_end=2;
		q_th=0.9;
	for(ii=0; ii<xres; ii++)//find highest quality strip
	{
		if(quality_x[ii]>q_th && flag_start==1)//quality above threshold, mark start
		{
			g_start=ii;
			flag_start=0;
			flag_end=1;
		}
		if(quality_x[ii]<q_th && flag_end ==1)//quality below threshold, if started, mark end
		{
			g_end=ii-1;
			flag_start=1;
			flag_end=0;
		}
		if(ii==(xres-1) && flag_end ==1)//last column, if started, save end
		{
			if (quality_x[ii]>=q_th)
				g_end=ii;
			else
				g_end=g_start;
			flag_end=0;
		}
		if(flag_end==0 && ii>0)//check if this strip is bigger than earlier strips
			if((g_end-g_start)>=(xy_start_R-xy_start_L))
			{
				xy_start_L=g_start;
				xy_start_R=g_end;
			}
	}
	
}
void  UnwrapGadget::center_x(float* phase_x, float* phase_y, int* mask, int xy_start_dw, int xy_start_up){

	float diffmean;

	int col_index,ii;

	std::vector<int> index_l;
	std::vector<int> index_s;
	std::vector<float> diff;
	
	index_l.reserve(yres);
	index_s.reserve(yres);
	diff.reserve(yres);
	for(col_index=0; col_index<xres; col_index++)
	{
		for(ii=xy_start_dw;ii<=xy_start_up; ii++)
		{
			if(mask[col_index*yres+ii])
				index_s.push_back(ii);

		}
		for(ii=0;ii<yres; ii++)
		{
			if(phase_x[ii+yres*col_index]!=0)
				index_l.push_back(ii);
		}
		if(index_s.empty())//no good strip, use whole column
		{
			for(ii=0;ii<yres; ii++)
			{
				if(mask[col_index*yres+ii])
					index_s.push_back(ii);
			}

			for(ii=0;ii<index_s.size(); ii++)//find column mean difference
				diff.push_back(phase_x[index_s[ii]+yres*col_index]-phase_y[index_s[ii]+yres*col_index]);

			diffmean=mean(diff.data(), diff.size());
			diffmean=(PI*2)*round(diffmean/(PI*2));
			for(ii=0;ii<index_l.size(); ii++)//adjust column
				phase_y[index_l[ii]+yres*col_index]+=diffmean;
		}
		else//cross reference with good strip
		{
			//find column mean difference
			for(ii=0;ii<index_s.size(); ii++)
				diff.push_back(phase_x[index_s[ii]+yres*col_index]-phase_y[index_s[ii]+yres*col_index]);

			diffmean=mean(diff.data(), diff.size());
			if(fabs(diffmean)>PI)//adjust column
			{
				diffmean=(PI*2)*round(diffmean/(PI*2));
				for(ii=0;ii<index_l.size(); ii++)
					phase_y[index_l[ii]+yres*col_index]+=diffmean;
			}		//
		}

		index_l.clear();
		index_s.clear();
		diff.clear();

	}
}
void  UnwrapGadget::center_y(float* phase_y, float* phase_x, int* t_mask, int xy_start_L, int xy_start_R)
{	

	float diffmean;
	int row_index, ii;

	std::vector<int> index_l;
	std::vector<int> index_s;
	std::vector<float> diff;
	
	index_l.reserve(xres);
	index_s.reserve(xres);
	diff.reserve(xres);
	for(row_index=0; row_index<yres; row_index++)
	{

		for(ii=xy_start_L;ii<=xy_start_R; ii++)
		{
			if(t_mask[row_index*xres+ii])
				index_s.push_back(ii);
		}
		for(ii=0;ii<xres; ii++)
		{
			if(phase_y[row_index+yres*ii]!=0)
				index_l.push_back(ii);
		}
		
		if(index_s.empty())//no good strip use whole row
		{  
			for(ii=0;ii<xres; ii++)
			{
				if(t_mask[row_index*xres+ii])
					index_s.push_back(ii);
			}
			for(ii=0;ii<index_s.size(); ii++)//find mean row difference
				diff.push_back(phase_y[row_index+yres*index_s[ii]]-phase_x[row_index+yres*index_s[ii]]);
			
			diffmean=mean(diff.data(), diff.size());
			diffmean=(PI*2)*round(diffmean/(PI*2));
		
			for(ii=0;ii<index_l.size(); ii++)//adjust row
				phase_x[row_index+yres*index_l[ii]]+=diffmean;
			
		}
		else//cross reference with good strip
		{
			for(ii=0;ii<index_s.size(); ii++)//find mean row difference
				diff.push_back(phase_y[row_index+yres*index_s[ii]]-phase_x[row_index+yres*index_s[ii]]);
			
			diffmean=mean(diff.data(), diff.size());
			
			if(fabs(diffmean)>PI)
			{	diffmean=(PI*2)*round(diffmean/(PI*2));
				//adjust row
				
				for(ii=0;ii<index_l.size(); ii++)
					phase_x[row_index+yres*index_l[ii]]+=diffmean;
			}
		}
		index_l.clear();
		index_s.clear();
		diff.clear();
		

	}
}	
void UnwrapGadget::final_compare(float* phase_x, float* phase_y, int iterations, MaskData *md)
{

	float seg_phi=2*PI;
	float *line_tmp_phase,* diff_tmp, *seg_mean_shift;
	std::vector<int> g_seg;
	std::vector<int> index_s;
	int index_x,index_y, ii,length, loopcount;
	std::vector<float> ref_tmp;
	std::vector<float> trl_tmp;
	std::vector<int> index_ls_TR;
	if(xres>=yres)
	{

	line_tmp_phase= new float[xres];
	diff_tmp=	new float[xres];
	seg_mean_shift= new float[xres]; 
	}
	else
	{
  
	line_tmp_phase= new float[yres];
	diff_tmp=	new float[yres];
	seg_mean_shift= new float[yres];
	}
	for(loopcount=0;loopcount<iterations; loopcount ++)
	{
		ref_tmp.resize(yres);
		trl_tmp.resize(yres);
		seg_phi/=2;
		for(index_x=0; index_x<xres; index_x++)
		{
			if(!(md->segY[index_x].empty()) || (md->signalY[index_x].empty()))
			{
				g_seg=md->segY[index_x];
				index_s=md->signalY[index_x];

				for(ii=0; ii<index_s.size(); ii++)
				{
					ref_tmp[ii]=phase_x[index_s[ii]+yres*index_x];
					trl_tmp[ii]=phase_y[index_s[ii]+yres*index_x];
				}

				index_ls_TR=dePULM_1D_ls(trl_tmp, seg_phi, g_seg);//find phase segments with gradient < seg_phi
				dePULM_2D_merging(trl_tmp, index_ls_TR, ref_tmp, index_s.size(), line_tmp_phase, seg_mean_shift, diff_tmp);//Cross reference those segmentsmd->Ylengths[index_x][0] was yres
			
				for(ii=0; ii<index_s.size(); ii++)
					phase_y[index_s[ii]+index_x*yres]=trl_tmp[ii];

			}
			
		}
	//////////////////////
		ref_tmp.resize(xres);
		trl_tmp.resize(xres);
		for(index_y=0; index_y<yres; index_y++)
		{
			if(!(md->segX[index_y].empty()) || (md->signalX[index_y].empty()))
			{
					
				g_seg=md->segX[index_y];
				index_s=md->signalX[index_y];
	
				for(ii=0; ii<index_s.size(); ii++)
				{
					ref_tmp[ii]=phase_y[index_y+index_s[ii]*yres];
					trl_tmp[ii]=phase_x[index_y+index_s[ii]*yres];
				}

				index_ls_TR=dePULM_1D_ls(trl_tmp,seg_phi,g_seg);//Identify phase segments with gradient < seg_phi
			
				dePULM_2D_merging(trl_tmp,index_ls_TR,ref_tmp,index_s.size(), line_tmp_phase, seg_mean_shift, diff_tmp);//Cross reference those strips
			
				for(ii=0; ii<index_s.size(); ii++)
					phase_x[index_y+index_s[ii]*yres]=trl_tmp[ii];

				
			

			}
	
		}
	}
	delete[] line_tmp_phase;
	delete[] diff_tmp;
	delete[] seg_mean_shift;


}
void UnwrapGadget::final_compare_brain(float* phase_x, float* phase_y, int iterations)
{

	float seg_phi=2*PI;
	float *line_tmp_phase,* diff_tmp, *seg_mean_shift;

	int index_x,index_y, ii,length, loopcount;
	
	std::vector<float> ref_tmp;
	std::vector<float> trl_tmp;
 	std::vector<int> index_ls_TR;
	if(xres>=yres)
	{
 
	line_tmp_phase= new float[xres];
	diff_tmp=	new float[xres];
	seg_mean_shift= new float[xres]; 
	}
	else
	{
	line_tmp_phase= new float[yres];
	diff_tmp=	new float[yres];
	seg_mean_shift= new float[yres];
	}
	for(loopcount=0;loopcount<iterations; loopcount ++)
	{	
		ref_tmp.resize(yres);
		trl_tmp.resize(yres);
		seg_phi/=2;
		for(index_x=0; index_x<xres; index_x++)
		{
				
				
				memcpy(ref_tmp.data(), phase_x+index_x*yres, sizeof(float)*yres);
				memcpy(trl_tmp.data(), phase_y+index_x*yres, sizeof(float)*yres);
			
				index_ls_TR=dePULM_1D_ls_brain(trl_tmp, seg_phi);//find phase segments with gradient < seg_phi
				dePULM_2D_merging(trl_tmp, index_ls_TR, ref_tmp,yres,line_tmp_phase, seg_mean_shift, diff_tmp);//Cross reference those segmentsmd->Ylengths[index_x][0] was yres

				memcpy(phase_y+index_x*yres,trl_tmp.data(), sizeof(float)*yres);
			
		}
	//////////////////////
		ref_tmp.resize(xres);
		trl_tmp.resize(xres);
		for(index_y=0; index_y<yres; index_y++)
		{
				
				for(ii=0; ii<xres; ii++)
				{
					ref_tmp[ii]=phase_y[index_y+ii*yres];
					trl_tmp[ii]=phase_x[index_y+ii*yres];
				}

				index_ls_TR=dePULM_1D_ls_brain(trl_tmp,seg_phi);//Identify phase segments with gradient < seg_phi
			
				dePULM_2D_merging(trl_tmp, index_ls_TR,ref_tmp,xres, line_tmp_phase, seg_mean_shift, diff_tmp);//Cross reference those strips
			
				for(ii=0; ii<xres; ii++)
					phase_x[index_y+ii*yres]=trl_tmp[ii];
		}
	}
	delete[] line_tmp_phase;
	delete[] diff_tmp;
	delete[] seg_mean_shift;


}
void UnwrapGadget::diff_x(float* phase_x,bool fullsignal, MaskData *md)
{
	int index_y, ii,jj, len_dw,len_up, fixcount;
	float diff_up, diff_dw, quality_br_tmp;	

	float **trl_dw;
	
	bool* bad_up = new bool[xres];
	bool* bad_dw = new bool[xres];
	int* fix_index = new int[xres];
	trl_dw=new float*[yres];

	for(ii=0; ii<yres; ii++)
	{
		trl_dw[ii]= new float[xres];
		for(jj=0; jj<xres; jj++)
			trl_dw[ii][jj]=phase_x[ii+jj*yres];
	}

	index_y=0; len_dw=0;  

	for(ii=0; ii<xres; ii++)
	{
		if(fabs(phase_x[index_y+ii*yres] -  phase_x[index_y+1+ii*yres])>PI)
			len_dw++;
	}

	if(len_dw!=0)
	{
		quality_br_tmp=(float)((xres-len_dw)/(xres));

		if(quality_br_tmp<0.95)
		{
			for(ii=0; ii<xres; ii++)
			{
				diff_dw=phase_x[index_y+1+ii*yres] - phase_x[index_y+ii*yres];
				phase_x[index_y+ii*yres]+=PI2*(round(diff_dw/PI2));
			}
		}
	}

	for(index_y=1; index_y<(yres-1);index_y++)
	{	
		len_up=0; 
		len_dw=0; fixcount=0;

		for(ii=0;ii<xres; ii++)
		{
			bad_up[ii]=0;
			bad_dw[ii]=0;

			if(fabs(trl_dw[index_y][ii] -  trl_dw[index_y+1][ii])>PI)
			{
				len_dw++;
				bad_dw[ii]=1;
			}

			if(fabs(trl_dw[index_y][ii] -  trl_dw[index_y-1][ii])>PI)
			{
				len_up++;
				bad_up[ii]=1;
			}


		}

		if(len_dw || len_up)//if there are any jumps
		{
			quality_br_tmp=(float)((xres-len_up)*(xres-len_dw)/(xres*xres));
			if(quality_br_tmp<.95)
			{
				for(ii=0; ii<xres; ii++)
				{
					if(bad_up[ii] && bad_dw[ii])
						fix_index[fixcount++]=ii;
				}
				if(fixcount!=0)//for points >PI off from neighbours 
				{
					for(ii=0; ii<fixcount; ii++)
					{	
						jj=fix_index[ii];

						diff_up=trl_dw[index_y-1][jj] -  trl_dw[index_y][jj];
						diff_dw=trl_dw[index_y+1][jj] -  trl_dw[index_y][jj];
						phase_x[index_y+jj*yres]+=PI2*(round(0.5*(diff_up+diff_dw)/PI2));//split the difference

					}
				}
			}
		}
	}


	index_y=yres-1; len_dw=0;
	len_up=0;

	
		for(ii=0;ii<xres; ii++)
		{
			bad_up[ii]=0;
			bad_dw[ii]=0;
	
			if(fabs(trl_dw[index_y][ii] -  trl_dw[index_y-1][ii])>PI)
			{
				len_up++;
				bad_up[ii]=1;
			}
		}
	if(len_up!=0)
	{
		quality_br_tmp=(xres-len_up)/(xres);

		if(quality_br_tmp<0.95)
		{
			for(ii=0; ii<xres; ii++)
			{
				jj=ii;
				diff_up=phase_x[index_y-1+jj*yres] - phase_x[index_y+jj*yres];
				phase_x[index_y+jj*yres]+=PI2*(round(diff_up/PI2));
			}
		}
	}	
////////////////////
	float* filter_out_x=filter(phase_x, xres, yres);
	float tmp_mean;
	int k;
	for(ii=0; ii<yres*xres;ii++)
	{
			filter_out_x[ii]=phase_x[ii]-filter_out_x[ii];
			if(fabs(filter_out_x[ii])>=PI && md->support_MASK[ii]==1)
				phase_x[ii]-=round(filter_out_x[ii]/PI2)*PI2;//unwrapped_phase_x should be equivalent to out_x by now

	}

	
	tmp_mean=0;
	k=0;
	for(ii=0; ii<yres*xres; ii++)
	{
		
			if(md->support_MASK[ii]==1)
			{
				tmp_mean+=phase_x[ii];
				k++;
			}

		
	}
	tmp_mean/=k;
	for(ii=0; ii<yres*xres; ii++)
	{
		
			if(fullsignal || md->MASK[ii]==1)//should shortcircuit with fullsignal (i.e. not try to check mask)
				phase_x[ii]-=PI2*round(tmp_mean/PI2);
		
	}
	for(ii=0; ii<yres; ii++)
	{
		delete[] trl_dw[ii];
	}
	delete[] trl_dw;
	delete[] bad_up;
	delete[] bad_dw;
	delete[] fix_index;
	delete[] filter_out_x; //dyn allocated in filter

}

GADGET_FACTORY_DECLARE(UnwrapGadget)
}


