#ifndef REGISTERVOLUMES_H
#define REGISTERVOLUMES_H

#include "vtkSlicerRegistrationLogic.h"

#include "CommandIterationUpdate.h"
#include "FeatureletRegistrationResult.h"
#include "regImage.h"

typedef FeatureletRegistrationResult FeatureletRegistrationResultType;
typedef FeatureletRegistrationResultType::Pointer FeatureletRegistrationResultPointer;

///	\fn			RegisterVolumes( ImageType* FixedImage, ImageType* MovingImage )
///	\brief		This function registers two volumes
///	\param		FixedImage The name of the fixed #ImageType object
///	\param		MovingImage The name of the moving #ImageType object
///	\return		This function returns an integer
///	This function uses of the VersorRigid3DTransform class for performing
///	registration of two 3D images. The class CenteredTransformInitializer is
///	used to initialize the center and translation of the transform. The case of
///	rigid registration of 3D images is probably one of the most commonly found
///	cases of image registration.

FeatureletRegistrationResult::Pointer RegisterVolumes(regImage* FixedImage, regImage* MovingImage);
FeatureletRegistrationResult::Pointer RegisterVolumes2(regImage* FixedImage, regImage* MovingImage);
FeatureletRegistrationResult::Pointer RegisterVolumesII2(regImage* FixedImage, regImage* MovingImage);
FeatureletRegistrationResult::Pointer RegisterVolumesII(regImage* FixedImage, regImage* MovingImage);

#endif // REGISTERVOLUMES_H
