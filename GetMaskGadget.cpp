//GetMaskGadget.cpp
//Written by Tim Whelan 2015
//Input ImageHeader->Complex Float 3D Array (Image Data) -> [MetaContainer]-> //
//Output ImageHeader->Float 3D Array (Phase Data) -> Int 3D Array (Support Mask) -> [Int 3D Array (Support Mask)] -> [MetaContainer] ->//
#include "GetMaskGadget.h"
using namespace Gadgetron;


int GetMaskGadget::process_config(ACE_Message_Block* mb)
{	
this->msg_queue()->high_water_mark(128);//This helps with memory. It's not a hard limit though. 

return GADGET_OK;
}



int GetMaskGadget::process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1)
{

    	//std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
      
	GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2 = AsContainerMessage< hoNDArray< std::complex<float> > > (m1->cont());
	GadgetContainerMessage<ISMRMRD::MetaContainer> *meta = AsContainerMessage<ISMRMRD::MetaContainer>(m2->cont());

	yres = m2->getObjectPtr()->get_size(0);
	xres = m2->getObjectPtr()->get_size(1);
	cres = m2->getObjectPtr()->get_size(3);

	/////
	/*m1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_FLOAT;//GADGET_IMAGE_REAL_FLOAT;
	m1->getObjectPtr()->image_type = ISMRMRD::ISMRMRD_IMTYPE_PHASE;
	//std::complex<float>* src = m2->getObjectPtr()->get_data_ptr();
	//std::vector<float>  dst(m2->getObjectPtr()->get_number_of_elements());
	
	GadgetContainerMessage<hoNDArray< float > > *cm2 = new GadgetContainerMessage<hoNDArray< float > >();
	hoNDArray< float > *phase_y_block = cm2->getObjectPtr();
	try{phase_y_block->create(xres,yres,cres);}
	catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to create itohy in Unwrap Gadget");
		return GADGET_FAIL;
	}
		
	m1->cont(cm2);*/
	/////


	static int myid = 0;
	myid++;
	/*if(myid>6)
	{
		m1->release();
		return GADGET_OK;
	}*/
	
	int maskProcedure=maskflag.value();
	int c;	


	boost::shared_ptr< std::vector<size_t> > dims = m2->getObjectPtr()->get_dimensions();

	//GDEBUG("%d Saving Queue Length  %d My Queue Length\n", this->next()->msg_queue()->message_count(), this->msg_queue()->message_count());
	GadgetContainerMessage< hoNDArray<int>> *supportmasks = new GadgetContainerMessage<hoNDArray<int>>();
	
	try{supportmasks->getObjectPtr()->create(dims.get());} //create support mask
	catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to create support masks in GetMaskGadget");
		return GADGET_FAIL;
	}
		
		

	GadgetContainerMessage< hoNDArray<int>> *masks = new GadgetContainerMessage<hoNDArray<int>>();
	if(maskProcedure!=0)//0=default, 1=use thresholds provided, 2=load
	{
		masks = new GadgetContainerMessage<hoNDArray<int>>();
		try{masks->getObjectPtr()->create(dims.get());}
		catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to create masks in GetMaskGadget");
		return GADGET_FAIL;
		}

		
	}

	std::complex<float>* complex_data_ptr = m2->getObjectPtr()->get_data_ptr();
	
	
	#pragma omp parallel for private(c)
	for(c = 0; c<cres; c++)
	{
		float thr;
		int ch_offset= xres*yres*c;
		float** filter_mag;
		int sum=0;
		int *ch_sppt_mask=supportmasks->getObjectPtr()->get_data_ptr()+ch_offset;
		int *ch_mask;
		if(maskProcedure<2)
		{
			filter_mag=filter2(complex_data_ptr+ch_offset, xres,yres); //mean filter image
			
		
		

			if(maskProcedure==0)
			{
				thr=graythresh_more(filter_mag);		//find threshold to use
				mask_create(ch_sppt_mask,filter_mag,thr*1.5);
				//mask_create(ch_mask,filter_mag,thr*.05);
				
			}			
			else
			{
				ch_mask =masks->getObjectPtr()->get_data_ptr()+ch_offset;
				//mask_create(ch_sppt_mask,filter_mag,threshold.value());
				//mask_create(ch_mask,filter_mag,threshold2.value());

				thr=graythresh_more(filter_mag);
				mask_create(ch_sppt_mask,filter_mag,thr*1.5);
				mask_create(ch_mask,filter_mag,0);
			

			}

			/*for(int jj=0; jj<yres;jj++)
				for(int ii=0; ii<xres; ii++)
					phase_y_block->get_data_ptr()[ii*yres+jj] =filter_mag[ii][jj];//ch_sppt_mask[ii*yres+jj];*/
			for(int ii=0; ii<xres; ii++)
			{
				delete[] filter_mag[ii];
			}
			/*for(int jj=0; jj<yres;jj++)
			{
				for(int ii=0; ii<xres; ii++)
					sum+=ch_mask[ii*yres+jj];
			}
			if(sum!=(yres*xres))
				GINFO("sum error in %d %d\n", c, myid);*/
			delete[] filter_mag;
			
		}
		else
		{
			//read read read
			//load support mask from file
			//sppt_mask_ptr[c]=?
	
			if(maskProcedure==2)
			{
			//load signal mask from file
			//mask_ptr[c]=?
			}		

		}
	


	}
	//build message block chain to pass on
	m2->cont(supportmasks);
	if(maskProcedure==0 || threshold2.value()==0)
	{
		*(supportmasks->getObjectPtr()->get_data_ptr())-=2;
		if(meta)
		supportmasks->cont(meta);
	}
	else
	{
		supportmasks->cont(masks);
		if(meta)
		masks->cont(meta);
	}
//std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
//auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	
	
	 if (this->next()->putq(m1) < 0) { //pass to next gadget
		GERROR("Failed to initialize Mask\n");
		m1->release();
	   return GADGET_FAIL;
	}
	
	return GADGET_OK;
	
}

void GetMaskGadget::mask_create(int * tofill, float **filter_mag, float filter_th)
{
	std::vector<std::vector<int>> mag_tmp_x(yres, std::vector<int>(xres, 0));
	std::vector<std::vector<int>> mag_tmp_y(xres, std::vector<int>(yres, 0));
	std::vector<int> sum_tmp(xres>yres?xres:yres,0);
	
	int ii,jj;
	
	std::vector<int*> output(xres);
	std::vector<std::vector<int>> output_t(yres, std::vector<int>(xres, 0));

	for(ii=0; ii<xres; ii++)
	{
		output[ii]=tofill+ii*yres;
	}	

	for(ii=0; ii<xres;ii++)
		for(jj=0;jj<yres; jj++)
		{
			if(filter_mag[ii][jj]>=filter_th)//thresholding along rows
			{
				output[ii][jj]=1;
				output_t[jj][ii]=1;
			}
			else
				output[ii][jj]=0;
		}

	
	for(ii=0; ii<yres; ii++)//for mag_tmp_x mask, use only points which have enough neighbours above threshold
	{
		if(ii!=0 && ii!=(yres-1))
		{
			for(jj=0; jj<xres; jj++)
				sum_tmp[jj]=output_t[ii+1][jj]+output_t[ii][jj]+output_t[ii-1][jj];
			for(jj=1; jj<xres-1; jj++)
			{
				if(sum_tmp[jj]>=2 && (output_t[ii][jj-1]+output_t[ii][jj]+output_t[ii][jj+1])>=2)
					mag_tmp_x[ii][jj]=1;
			}	
		}
		else
		{

			for(jj=0; jj<xres; jj++)
				sum_tmp[jj]=output_t[ii][jj];
			for(jj=1; jj<xres-1; jj++)
			{
				if(sum_tmp[jj] && (output_t[ii][jj-1]+output_t[ii][jj]+output_t[ii][jj+1])>=2)
					mag_tmp_x[ii][jj]=1;
			}

		}

		if(mag_tmp_x[ii][1])
			mag_tmp_x[ii][0]=1;

		if(mag_tmp_x[ii][xres-2])
			mag_tmp_x[ii][xres-1]=1;

	}
	
	for(ii=0; ii<xres; ii++)//for mag_tmp_y mask, use only points which have enough neighbours above threshold
	{
		if(ii!=0 && ii!=(xres-1))
		{
			for(jj=0; jj<yres; jj++)
				sum_tmp[jj]=output[ii+1][jj]+output[ii][jj]+output[ii-1][jj];
			for(jj=1; jj<yres-1; jj++)
			{
				if(sum_tmp[jj]>=2 && (output[ii][jj-1]+output[ii][jj]+output[ii][jj+1])>=2)
					mag_tmp_y[ii][jj]=1;
			}	
		}
		else
		{

			for(jj=0; jj<yres; jj++)
				sum_tmp[jj]=output[ii][jj];
			for(jj=1; jj<yres-1; jj++)
			{
				if(sum_tmp[jj]>0 && (output[ii][jj-1]+output[ii][jj]+output[ii][jj+1])>=2)
					mag_tmp_y[ii][jj]=1;
			}

		}

		if(mag_tmp_y[ii][1])
			mag_tmp_y[ii][0]=1;
	
		if(mag_tmp_y[ii][yres-2])
			mag_tmp_y[ii][yres-1]=1;
	}
	
	for(ii=0; ii<xres;ii++){
	for(jj=0;jj<yres; jj++)
	{
		if(mag_tmp_y[ii][jj] || mag_tmp_x[jj][ii])//not sure about this one
			output[ii][jj]=1;
		else
			output[ii][jj]=0;
	}

	}



	

}
float  GetMaskGadget::graythresh_more(float** image)
{
	int ii,jj, numpix, max_index;
	float realmax=0;
	float maxmag=0;
	float maxsbs=0;
	int bins[256];
	float p[256];
	float omega[256];
	float mu[256];
	float sbs[256];	
	int count=0;
	float thr1;
	for(ii=0;ii<256;ii++)
	{
		bins[ii]=0;
		sbs[ii]=0;
	}
	numpix=xres*yres;
	std::vector<float> test_im(numpix);

	for(ii=0; ii<xres;ii++){
	for(jj=0; jj<yres; jj++)
		realmax=image[ii][jj]>realmax?image[ii][jj]:realmax;
	}
	
	for(ii=0; ii<xres;ii++){
	for(jj=0; jj<yres; jj++)
		test_im[jj+yres*ii]=image[ii][jj]/realmax;
	}

	for(ii=0; ii<numpix;ii++)
		bins[(int)(test_im[ii]*255)]++;
	
	for(ii=0; ii<256;ii++){
		p[ii]=(float)bins[ii]/numpix;
	
	}
	omega[0]=p[0];
	
	for(ii=1; ii<256;ii++)
	{
		omega[ii]=omega[ii-1]+p[ii];
			
	}
	mu[0]=p[0];
	for(ii=1; ii<256;ii++)
		mu[ii]=mu[ii-1]+p[ii]*(ii+1);
	
	for(ii=0; ii<255;ii++)//only to 254 to avoid divide by 0 (omega[255]=1) 
	{
		if(omega[ii]!=0 && omega[ii] != 1)
		sbs[ii]=(mu[255]*omega[ii]-mu[ii])*(mu[255]*omega[ii]-mu[ii])/(omega[ii]*(1-omega[ii]));
		
	}
	maxsbs=0;
	for(ii=0; ii<255;ii++)
	{
		if(sbs[ii]>maxsbs)
		{
		    maxsbs=sbs[ii];
		   max_index=ii;
		}
	}
	thr1 = (float)max_index/255;
	for(ii=0; ii<xres;ii++){
		for(jj=0;jj<yres;jj++)
		{
		if(image[ii][jj]<(thr1*realmax))
			test_im[jj+yres*ii]=image[ii][jj];
		else
			test_im[jj+yres*ii]=0;
		}
	}
	
	thr1=stdev(test_im.data(), numpix);
	return thr1;
}
GADGET_FACTORY_DECLARE(GetMaskGadget)

