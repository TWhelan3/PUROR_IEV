//
//  DicomSave.cpp
//  network
//
//  Created by Sahar Rabinoviz on 2015-07-17.
//  Copyright (c) 2015 CFMM. All rights reserved.
//

#include "DicomSave.h"

namespace Gadgetron {
    int DicomSave::process_config(ACE_Message_Block * mb)
    {       
        return GADGET_OK;
    }
    
    //save dicom file
    int DicomSave::process(GadgetContainerMessage<ISMRMRD::ImageHeader> * m1)
    {
        GadgetContainerMessage<DcmFileFormat>* dcm_file_message = AsContainerMessage<DcmFileFormat>(m1->cont());
        
        if (!dcm_file_message) {
            GERROR("Received invalid message object\n");
            m1->release();
            return GADGET_FAIL;
        }

        DcmFileFormat * dcm = dcm_file_message->getObjectPtr();

        std::string fname           = workingDirectory.value() + "/" +  dicom_dir.value();
        std::string mkdir_command   = "mkdir " + fname + + " > /dev/null 2>&1";
        
        system(mkdir_command.c_str());
        
        GadgetContainerMessage<std::string> * f = AsContainerMessage<std::string>(dcm_file_message->cont());
        
        std::string imageName = fname + "/" + *f->getObjectPtr() + ".dcm";
        
        OFCondition res = dcm->saveFile(imageName.c_str(), EXS_LittleEndianExplicit);
        
        if( res != EC_Normal){
            std::cout << "dicom file failed to save, code:" << res.code() << " reason:" << res.text() << std::endl;
            return GADGET_FAIL;
        }
        
        m1->release();
        
        return GADGET_OK;
    }
        
    GADGET_FACTORY_DECLARE(DicomSave);

}
