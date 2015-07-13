#ifndef SUBSAMPLEFEUT_H
#define SUBSAMPLEFEUT_H

#include "vtkSlicerRegistrationLogic.h"
#include "FeatureletRegistrationResult.h"

///	\fn			SubsampleVolume( regImage* Image, ImageType::SizeType FeatureletSize , ImageType::IndexType ImageIndex )
///	\brief		This function subsamples a volume
///	\param		Image The name of the #regImage object
///	\param		FeatureletSize The size of one featurelet
///	\param		ImageIndex The current index of the image
///	\return		This function returns an integer
///	This function performs subsampling of a volume using ITK classes. In order
///	to avoid aliasing artifacts, the volume must be processed by a low-pass
///	filter before resampling.

int SubsampleVolume( ImageType::Pointer Image,
                     ImageType::SizeType FeatureletSize,
                     ImageType::SizeType SearchRegionSize,
                     ImageType::IndexType ImageIndex) ;


///	\fn			CheckFeaturelet( regImage* Image )
///	\brief		This function checks a featurelet if it is necessary to be registered
///	\param		Image The name of the #regImage object
///	\return		This function returns an integer
///	This function checks a featurelet if it is necessary to be registered to
///	avoid registration errors or a bad result.

Featurelet::Status CheckFeaturelet( ImageType::Pointer Image ) ;

#endif // SUBSAMPLEFEUT_H
