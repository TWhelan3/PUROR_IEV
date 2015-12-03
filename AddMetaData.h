#ifndef AddMetaData_H_
#define AddMetaData_H_

#include "PUROR.h"
#include <dcmtk/dcmdata/dcfilefo.h>
#include "meta_key_def.h"
#include <ismrmrd/xml.h>
namespace Gadgetron{

	class EXPORTGADGETSPUROR AddMetaData:
	public Gadget1<DcmFileFormat>
    	{
	
	public:
		GADGET_DECLARE(AddMetaData);

		virtual int process_config(ACE_Message_Block*);
		virtual int process(GadgetContainerMessage<DcmFileFormat> * m1);
		GADGET_PROPERTY(intercept, float, "Intercept for Dicom", 0);
		GADGET_PROPERTY(slope, float, "Slope for Dicom", 0);
	
		char generatedStudyUID[64];
	};

}

#endif
