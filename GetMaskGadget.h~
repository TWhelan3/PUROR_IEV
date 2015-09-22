//GetMaskGadget.h

#ifndef GETMASKGADGET_H
#define GETMASKGADGET_H

#include "PUROR.h"


namespace Gadgetron
{

class EXPORTGADGETSPUROR GetMaskGadget : 
public Gadget1<ISMRMRD::ImageHeader>
{
 public:
  GADGET_DECLARE(GetMaskGadget)

    protected:

	virtual int process_config(ACE_Message_Block* mb);
  	virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1);
	
	GADGET_PROPERTY(maskflag, int, "How to acquire initial masks", 0);//0=default, 1=use thresholds provided, 2=load, 3=load support only
	GADGET_PROPERTY(threshold, float, "Thresholding for mask", 0);
	GADGET_PROPERTY(threshold2, float, "Support thresholding for mask", 0);
	int xres;
	int yres;
	int cres;
	
	void mask_create(int* outputptr, float **filter_mag, float filter_th);
	float graythresh_more(float** image);

};

}  
#endif //GetMASKGADGET_H
