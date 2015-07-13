#ifndef RESAMPLEVOLUME_H
#define RESAMPLEVOLUME_H
#include "RegisterVolumes.h"
//#include "vtkSlicerRegistrationLogic.h"
#include "regImage.h"



///	\fn			ResampleVolumesToBeIsotropic( regImage* Image )
///	\brief		This function resamples a voulme to be isotropic
///	\param		Image The name of the #regImage object
///	\return		This function returns an integer
///	It is unfortunate that it is still very common to find medical image
///	datasets that have been acquired with large inter-sclice spacings that
///	result in voxels with anisotropic shapes. In many cases these voxels have
///	ratios of [1:5] or even [1:10] between the resolution in the plane (x,y) and
///	the resolution along the z axis. Such datasets are close to useless for the
///	purpose of computer assisted image analysis.

int ResampleVolumesToBeIsotropic(regImage *Image);

int ResampleVolumes(regImage *Image) ;

#endif // RESAMPLEVOLUME_H
