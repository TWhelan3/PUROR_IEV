#ifndef GADGETRON_PUROR_EXPORT_H_
#define GADGETRON_PUROR_EXPORT_H_

#if defined (WIN32)
#if defined (__BUILD_GADGETRON_PUROR__)
#define EXPORTGADGETSPUROR __declspec(dllexport)
#else
#define EXPORTGADGETSPUROR __declspec(dllimport)
#endif
#else
#define EXPORTGADGETSPUROR
#endif

#endif /* GADGETRON_PUROR_EXPORT_H_ */
