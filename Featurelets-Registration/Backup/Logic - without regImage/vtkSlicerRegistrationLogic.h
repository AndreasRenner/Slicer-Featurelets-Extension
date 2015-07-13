/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// .NAME vtkSlicerRegistrationLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the volumes


#ifndef __vtkSlicerRegistrationLogic_h
#define __vtkSlicerRegistrationLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"
#include "vtkSlicerRegistrationModuleLogicExport.h"
class vtkSlicerVolumesLogic;
class vtkMRMLVolumeNode;
class vtkMRMLRegistrationNode;

// STD includes
#include <cstdlib>
#include "stdio.h"
#include "string"
#include "math.h"
#include "map"
#include "fstream"
#include "gdcmGlobal.h"

// ITK includes
#include "itkAmoebaOptimizer.h"
#include "itkBMPImageIO.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkCastImageFilter.h"
#include "itkCenteredTransformInitializer.h"
#include "itkCheckerBoardImageFilter.h"
#include "itkCommand.h"
#include "itkConstantPadImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkGradientDescentOptimizer.h"
#include "itkGradientDifferenceImageToImageMetric.h"
#include "itkGDCMImageIO.h"
#include "itkIdentityTransform.h"
#include <itkImage.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegistrationMethod.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkKappaStatisticImageToImageMetric.h"
#include "itkKullbackLeiblerCompareHistogramImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMatchCardinalityImageToImageMetric.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMeanReciprocalSquareDifferenceImageToImageMetric.h"
#include "itkMeanSquaresHistogramImageToImageMetric.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkMetaImageIO.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMutualInformationImageToImageMetric.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkNormalizeImageFilter.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkNormalVariateGenerator.h"
#include "itkOnePlusOneEvolutionaryOptimizer.h"
#include "itkPowellOptimizer.h"
#include "itkQuaternionRigidTransform.h"
#include "itkQuaternionRigidTransformGradientDescentOptimizer.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkSobelEdgeDetectionImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkVectorResampleImageFilter.h"
#include "itkVersorRigid3DTransform.h"
#include "itkVersorRigid3DTransformOptimizer.h"
#include "itkWarpImageFilter.h"
#include "itkWindowedSincInterpolateImageFunction.h"


typedef short PixelType;
const unsigned int Dimension = 3;

typedef itk::Image<PixelType, Dimension> ImageType;


/// \ingroup Slicer_QtModules_ExtensionTemplate
class VTK_SLICER_REGISTRATION_MODULE_LOGIC_EXPORT vtkSlicerRegistrationLogic :
  public vtkSlicerModuleLogic
{
public:

  static vtkSlicerRegistrationLogic *New();
  vtkTypeMacro(vtkSlicerRegistrationLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);
  void SetAndObserveRegistrationNode(vtkMRMLRegistrationNode* node);

  void SetVolumesLogic(vtkSlicerVolumesLogic* logic);
  vtkSlicerVolumesLogic* GetVolumesLogic();

  int RunClicked(vtkMRMLRegistrationNode*);
  int ShowVolume(vtkMRMLRegistrationNode*, bool);

  vtkGetObjectMacro(RegistrationNode, vtkMRMLRegistrationNode);

  bool ConvertVolumeNodeToItkImage(vtkMRMLVolumeNode*, itk::Image<PixelType, Dimension>::Pointer);


protected:
  vtkSlicerRegistrationLogic();
  virtual ~vtkSlicerRegistrationLogic();

  //static bool ComputeOrientationMatrixFromScanOrder(const char *order, vtkMatrix4x4 *outputMatrix);

  virtual void SetMRMLSceneInternal(vtkMRMLScene* newScene);

  /// Register MRML Node classes to Scene. Gets called automatically when the MRMLScene is attached to this logic class.
  virtual void RegisterNodes();

  virtual void UpdateFromMRMLScene();
  virtual void OnMRMLSceneNodeAdded(vtkMRMLNode* node);
  virtual void OnMRMLSceneNodeRemoved(vtkMRMLNode* node);
  virtual void OnMRMLSceneEndImport();
  virtual void OnMRMLSceneEndClose();

  vtkMRMLRegistrationNode* RegistrationNode;

private:

  vtkSlicerRegistrationLogic(const vtkSlicerRegistrationLogic&); // Not implemented
  void operator=(const vtkSlicerRegistrationLogic&); // Not implemented

  class vtkInternal;
  vtkInternal* Internal;
};

#endif
