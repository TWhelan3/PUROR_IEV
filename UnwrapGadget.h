#ifndef UnwrapGadget_H_
#define UnwrapGadget_H_

#include "PUROR.h"
#include <vector>
namespace Gadgetron{

  class EXPORTGADGETSPUROR UnwrapGadget:
public Gadget1<ISMRMRD::ImageHeader>
    {

    public:
      	GADGET_DECLARE(UnwrapGadget);

    protected:
      
     	virtual int process_config(ACE_Message_Block*);
	virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader>* m1);
	int xres, yres;
	int  cres;
	
	float q_th; //might want to set this is config



	void unwrap_rows(float*,float*);
	void unwrap_columns(float*, float*);

	void shift_to_mean(float*, MaskData*);
	void shift_to_mean_y(float*, MaskData*);
	void shift_to_mean_brain(float*, MaskData*);
	void shift_to_mean_y_brain(float*, MaskData*);
	

	void calc_quality_x(float*, int*, std::vector<float>&,int&,int&);
	void calc_quality_y(float*, int*, std::vector<float>&,int&,int&);
	void center_x(float*, float*, int*, int,int);
	void center_y(float*, float*, int*, int,int);	

	void final_compare(float*, float*, int, MaskData *);
	void final_compare_brain(float*, float*, int);
	void diff_x(float*,bool,MaskData*);

     };
}

#endif /* UnwrapGadget_H_ */
