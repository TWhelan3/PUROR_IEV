//HPFGadget.h
//Written by Tim Whelan 2015

#ifndef HPFGadget_H_
#define HPFGadget_H_

#include "PUROR.h"
#include "hoNDFFT.h"

enum FilterProc 	{DERICHE, FFT};
enum OUTPUT {LFS, PHASE};
namespace Gadgetron{

  class EXPORTGADGETSPUROR HPFGadget:
public Gadget1<ISMRMRD::ImageHeader>
    {

    public:
      	GADGET_DECLARE(HPFGadget);

    protected:
      
     	virtual int process_config(ACE_Message_Block*);
	virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1);
	GADGET_PROPERTY(sigma, float, "Sigma for FFT in filter", .003); 
	GADGET_PROPERTY(output, int, "What to output?",0 ); //0=LFS 1=Phase
	int yres,xres,cres;
	int num_ch;

	std::vector<std::vector<double> > F;
	std::vector<float> x2;
	std::vector<float> y2;

	float FOVx;
	float FOVy;

     };
}

#endif /* HPFGadget_H_ */
