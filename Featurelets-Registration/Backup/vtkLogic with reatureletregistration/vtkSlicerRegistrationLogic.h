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
//#include "FeatureletRegistrationResult.h"
#include "CommandIterationUpdate.h"
class vtkSlicerVolumesLogic;
class vtkMRMLVolumeNode;
class vtkMRMLRegistrationNode;
//class FeatureletRegistrationResult;

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
#include "itkObject.h"
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


namespace Featurelet{
enum Status
{
  OK=0,
  SizeError=1,
  sameColor=2,
  internalError=3
};
}


typedef short PixelType;
const unsigned int Dimension = 3;

typedef itk::Image<PixelType, Dimension> ImageType;

typedef FeatureletRegistrationResult FeatureletRegistrationResultType;
typedef FeatureletRegistrationResultType::Pointer FeatureletRegistrationResultPointer;


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
  int SubsampleVolume( ImageType::Pointer Image,
                       ImageType::SizeType FeatureletSize,
                       ImageType::SizeType SearchRegionSize,
                       ImageType::IndexType ImageIndex);
  Featurelet::Status CheckFeaturelet(ImageType::Pointer Image);
  FeatureletRegistrationResult::Pointer RegisterVolumes(ImageType* FixedImage, ImageType* MovingImage);
  FeatureletRegistrationResult::Pointer RegisterVolumes2(ImageType* FixedImage, ImageType* MovingImage);
  FeatureletRegistrationResult::Pointer RegisterVolumesII2(ImageType* FixedImage, ImageType* MovingImage);
  FeatureletRegistrationResult::Pointer RegisterVolumesII(ImageType* FixedImage, ImageType* MovingImage);


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


class FeatureletRegistrationResult :
  public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef FeatureletRegistrationResult          Self;
  typedef itk::Object                           Superclass;
  typedef itk::SmartPointer<Self>               Pointer;
  typedef itk::SmartPointer<const Self>         ConstPointer;

  typedef itk::TranslationTransform<double> TransformType;
  typedef TransformType::Pointer TransformPointer;
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef OptimizerType::Pointer OptimizerPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FeatureletRegistrationResult, itk::Object);

  virtual void PrintSelf(std::ostream& os,itk::Indent indent) const {
    Superclass::PrintSelf(os, indent);
    os << indent << "Offset: " << m_Transform->GetOffset() << "\n";
    os << indent << "Metric Value II: "<< m_Optimizer->GetValue()<<"\n";
    os << indent << "status Fixed featurelet: " << m_StatusFixed << "\n";
    os << indent << "status Moving featurelet: " << m_StatusMoving << "\n";
    os << indent << "Number of Featurelet: " << m_FeatureletNumber<< "\n";
  }

  itkGetObjectMacro(Transform, TransformType);
  itkSetObjectMacro(Transform, TransformType);
  itkSetObjectMacro(Optimizer, OptimizerType);
  itkGetObjectMacro(Optimizer, OptimizerType);

  Featurelet::Status GetStatusFixed() {
    return m_StatusFixed;
  }

  Featurelet::Status GetStatusMoving() {
    return m_StatusMoving;
  }

  void SetStatusFixed(Featurelet::Status status) {
    m_StatusFixed = status;
  }

  void SetStatusMoving(Featurelet::Status status) {
    m_StatusMoving = status;
  }

  int GetStatusRegistration() {
    return m_StatusRegistration;
  }

  void SetStatusRegistration(int error) {
    m_StatusRegistration = error;
  }

  int GetMetricRegistration() {
    return m_MetricRegistration;
  }

  void SetMetricRegistration(int error) {
    m_MetricRegistration = error;
  }

  int* GetFeatureletSize() {
    return m_FeatureletSize;
  }

  void SetFeatureletSize(int i, int j, int k) {
    m_FeatureletSize[0]=i;
    m_FeatureletSize[1]=j;
    m_FeatureletSize[2]=k;
  }

  int* GetFeatureletIndex() {
    return  m_FeatureletIndex;
  }

  void SetFeatureletIndex(int i, int j, int k) {
    m_FeatureletIndex[0]=i;
    m_FeatureletIndex[1]=j;
    m_FeatureletIndex[2]=k;
  }

  int* GetFeatureletCenter() {
    return  m_FeatureletCenter;
  }

  void SetFeatureletCenter(int i, int j, int k) {
    m_FeatureletCenter[0]=i;
    m_FeatureletCenter[1]=j;
    m_FeatureletCenter[2]=k;
  }

  int* GetFeatureletGridPosition() {
    return  m_FeatureletGridPosition;
  }

  void SetFeatureletGridPosition(int i, int j, int k) {
    m_FeatureletGridPosition[0]=i;
    m_FeatureletGridPosition[1]=j;
    m_FeatureletGridPosition[2]=k;
  }

  int GetFeatureletNumber() {
    return m_FeatureletNumber;
  }

  void SetFeatureletNumber(int numero) {
    m_FeatureletNumber= numero;
  }

  void SetImageSizeM(int i, int j, int k) {
    m_ImageSizeM[0]=i;
    m_ImageSizeM[1]=j;
    m_ImageSizeM[2]=k;
  }

protected:
  FeatureletRegistrationResult() {
    m_Transform = TransformType::New();
    m_Optimizer = OptimizerType::New();
  }
  virtual ~FeatureletRegistrationResult();

  int m_FeatureletNumber;
  TransformPointer m_Transform;
  OptimizerPointer m_Optimizer;

  Featurelet::Status m_StatusFixed;
  Featurelet::Status m_StatusMoving;
  int m_StatusRegistration;
  int m_MetricRegistration;
  int m_FeatureletSize[3];
  int m_FeatureletIndex[3];
  int m_FeatureletCenter[3];
  int m_FeatureletGridPosition[3];
  int m_ImageSizeM[3];
};

#endif
