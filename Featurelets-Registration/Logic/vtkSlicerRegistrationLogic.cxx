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

// Registration Logic includes
#include "vtkMRMLRegistrationNode.h"
#include "vtkSlicerRegistrationLogic.h"
#include "vtkSlicerVolumesLogic.h"

// MRML includes
//#include <vtkMRMLAnnotationROINode.h>    //Not yet implemented
#include <vtkMRMLDiffusionWeightedVolumeDisplayNode.h>
#include <vtkMRMLDiffusionWeightedVolumeNode.h>
#include <vtkMRMLDiffusionTensorVolumeNode.h>
#include <vtkMRMLLinearTransformNode.h>
#include <vtkMRMLMarkupsFiducialNode.h>
#include <vtkMRMLMarkupsNode.h>
#include <vtkMRMLRegistrationNode.h>
#include <vtkMRMLScalarVolumeDisplayNode.h>
#include <vtkMRMLScene.h>
#include <vtkMRMLSelectionNode.h>
#include <vtkMRMLVectorVolumeDisplayNode.h>
#include <vtkMRMLVectorVolumeNode.h>
#include <vtkMRMLVolumeNode.h>
#include <vtkMRMLTransformNode.h>

// VTK includes
#include <vtkAbstractTransform.h>
#include <vtkImageData.h>
#include <vtkImageExport.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>

// ITK includes
#include <itkImage.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkPoint.h>

// Featurelets includes
#include "CommandIterationUpdate.h"
#include "DeformationFieldGenerator.h"

// STD includes
#include <cassert>
#include <iostream>

//----------------------------------------------------------------------------
class vtkSlicerRegistrationLogic::vtkInternal {
public:
  vtkInternal();
  vtkSlicerVolumesLogic* VolumesLogic;
};

//----------------------------------------------------------------------------
vtkSlicerRegistrationLogic::vtkInternal::vtkInternal() {
  this->VolumesLogic = 0;
}

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerRegistrationLogic);

//----------------------------------------------------------------------------
vtkSlicerRegistrationLogic::vtkSlicerRegistrationLogic() {
  this->Internal = new vtkInternal;
  this->RegistrationNode = NULL;
}

//----------------------------------------------------------------------------
vtkSlicerRegistrationLogic::~vtkSlicerRegistrationLogic() {
  delete this->Internal;
  this->SetAndObserveRegistrationNode(NULL); //release the node object to avoid memory leaks
}

//----------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::SetVolumesLogic(vtkSlicerVolumesLogic* logic) {
  this->Internal->VolumesLogic = logic;
}

//----------------------------------------------------------------------------
vtkSlicerVolumesLogic* vtkSlicerRegistrationLogic::GetVolumesLogic() {
  return this->Internal->VolumesLogic;
}

//----------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::PrintSelf(ostream& os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::SetAndObserveRegistrationNode(vtkMRMLRegistrationNode *node) {
  vtkSetAndObserveMRMLNodeMacro(this->RegistrationNode, node);
}

//---------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene) {
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndImportEvent);
  events->InsertNextValue(vtkMRMLScene::EndCloseEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEvents(newScene, events.GetPointer());
}

//----------------------------------------------------------------------------
int vtkSlicerRegistrationLogic::RunClicked(vtkMRMLRegistrationNode* pnode) {
  vtkMRMLScene *scene = this->GetMRMLScene();

  bool debugMode = pnode->GetcheckBoxDebug();

  if(debugMode) {
    std::cerr << "Fixed Image Node:       " << pnode->GetFixedImageNodeID() << std::endl;
    std::cerr << "Moving Image Node:      " << pnode->GetMovingImageNodeID() << std::endl;
    std::cerr << "Deformed Image Node:    " << pnode->GetDeformedImageNodeID() << std::endl;
    std::cerr << "Deformation Field Node: " << pnode->GetDeformationFieldID() << std::endl;
  }

  //Get the needed data from the RegistrationNode
  vtkMRMLVolumeNode *fixedImage =
    vtkMRMLVolumeNode::SafeDownCast(scene->GetNodeByID(pnode->GetFixedImageNodeID()));
  vtkMRMLVolumeNode *movingImage =
    vtkMRMLVolumeNode::SafeDownCast(scene->GetNodeByID(pnode->GetMovingImageNodeID()));
  vtkMRMLVolumeNode *deformedImage =
    vtkMRMLVolumeNode::SafeDownCast(scene->GetNodeByID(pnode->GetDeformedImageNodeID()));
  vtkMRMLTransformNode *deformationField =
    vtkMRMLTransformNode::SafeDownCast(scene->GetNodeByID(pnode->GetDeformationFieldID()));
  vtkMRMLMarkupsFiducialNode *fiducialPoints = NULL;
  MaxStepLength = (double) pnode->GetMaxStepLength();
  MinStepLength = (double) pnode->GetMinStepLength();
  NumberIterations = (int) pnode->GetNumberIterations();
  bool correl = pnode->GetUseCorrelationForSimilarity();
  bool linear = pnode->GetUseLinearCorrelation();
  bool UseFiducialPoints = pnode->GetcheckBoxFiducial();
  bool rigid = pnode->GetcheckBoxRigid();
  bool zdifferent = pnode->GetcheckBoxZDifferent();
  if(UseFiducialPoints){
    fiducialPoints = vtkMRMLMarkupsFiducialNode::SafeDownCast(scene->GetNodeByID(pnode->GetFiducialPointsID()));
    if(debugMode) {
      std::cerr << "Fiducial Point Node:    " << pnode->GetFiducialPointsID() << std::endl;
    }
  }


  vtkImageData* imageData;
  if(debugMode) {
    imageData = fixedImage->GetImageData();
    std::cerr << "Fixed Pixeltype: " << imageData->GetScalarType() << " " << imageData->GetScalarTypeAsString() << std::endl;
    imageData = movingImage->GetImageData();
    std::cerr << "Moving Pixeltype: " << imageData->GetScalarType() << " " << imageData->GetScalarTypeAsString() << std::endl;
  }

  if(!fixedImage || !movingImage) {
    std::cerr << "Registration: Failed to look up fixed/moving image!" << std::endl;
    return -1;
  }
  if(!deformedImage || !deformationField) {
    std::cerr << "Registration: Deformed Image/Deformation field not initialized" << std::endl;
    return -1;
  }
  if(UseFiducialPoints && !fiducialPoints){
    std::cerr << "Registration: Failed to look up fiducial points!" << std::endl;
    return -1;
  }

  // check the volume type
  vtkMRMLDiffusionTensorVolumeNode *dtvnode = vtkMRMLDiffusionTensorVolumeNode::SafeDownCast(movingImage);
  vtkMRMLDiffusionWeightedVolumeNode *dwvnode = vtkMRMLDiffusionWeightedVolumeNode::SafeDownCast(movingImage);
  vtkMRMLVectorVolumeNode *vvnode = vtkMRMLVectorVolumeNode::SafeDownCast(movingImage);
  vtkMRMLScalarVolumeNode *svnode = vtkMRMLScalarVolumeNode::SafeDownCast(movingImage);

  if(!this->Internal->VolumesLogic) {
    std::cerr << "Registration: ERROR: failed to get hold of Volumes logic" << std::endl;
    return -2;
  }

  std::ostringstream outSS;
  outSS << movingImage->GetName() << "_deformed";
  deformedImage->SetName(outSS.str().c_str());

  if(dtvnode) {
    std::cerr << "Registration: ERROR: Diffusion tensor volume not supported by this module!" << std::endl;
    return -2;
  }
  else if(dwvnode) {
    std::cerr << "Registration: ERROR: Diffusion weighted volume not supported by this module!" << std::endl;
    return -2;
  }
  else if(vvnode) {
    std::cerr << "Registration: ERROR: Vector volume not supported by this module!" << std::endl;
    return -2;
  }
  else if(svnode)
    std::cerr << "Moving Image is a Scalar Volume - we can go on." << std::endl;
  else {
    std::cerr << "Moving Image not recognized!" << std::endl;
    return -1;
  }

  dtvnode = vtkMRMLDiffusionTensorVolumeNode::SafeDownCast(fixedImage);
  dwvnode = vtkMRMLDiffusionWeightedVolumeNode::SafeDownCast(fixedImage);
  vvnode = vtkMRMLVectorVolumeNode::SafeDownCast(fixedImage);
  svnode = vtkMRMLScalarVolumeNode::SafeDownCast(fixedImage);

  if(dtvnode) {
    std::cerr << "Registration: ERROR: Diffusion tensor volume not supported by this module!" << std::endl;
    return -2;
  }
  else if(dwvnode) {
    std::cerr << "Registration: ERROR: Diffusion weighted volume not supported by this module!" << std::endl;
    return -2;
  }
  else if(vvnode) {
    std::cerr << "Registration: ERROR: Vector volume not supported by this module!" << std::endl;
    return -2;
  }
  else if(svnode)
    std::cerr << "Fixed Image is a Scalar Volume - we can go on." << std::endl;
  else {
    std::cerr << "Fixed Image not recognized!" << std::endl;
    return -1;
  }

  //Check Pixeltype
  int ScalarType;
  imageData = fixedImage->GetImageData();
  ScalarType=imageData->GetScalarType();
  std::cerr << "PixelType Fixed: " << ScalarType << std::endl;
  if(ScalarType!=4) {
    std::cerr << "Fixed image is not of PixelType ,,short''. Please change the PixelType to short by using ";
    std::cerr << "e.g. the module ,,Cast Scalar Volume''." << std::endl;
    return -2;
  }
  imageData = movingImage->GetImageData();
  ScalarType=imageData->GetScalarType();
  std::cerr << "PixelType Moving: " << ScalarType << std::endl;
  if(ScalarType!=4) {
    std::cerr << "Moving image is not of PixelType ,,short''. Please change the PixelType to short by using ";
    std::cerr << "e.g. the module ,,Cast Scalar Volume''." << std::endl;
    return -2;
  }

  //Convert images to ITK
  ImageType::Pointer fixedImageItk = ImageType::New();
  ImageType::Pointer movingImageItk = ImageType::New();
  ConvertVolumeNodeToItkImage(fixedImage, fixedImageItk);
  ConvertVolumeNodeToItkImage(movingImage, movingImageItk);
  if(debugMode) {
    std::cerr << "FixedITKSpacing:  " << fixedImageItk->GetSpacing()[0] << fixedImageItk->GetSpacing()[1];
    std::cerr << fixedImageItk->GetSpacing()[2] << std::endl;
    std::cerr << "MovingITKSpacing: " << movingImageItk->GetSpacing()[0] << movingImageItk->GetSpacing()[1];
    std::cerr << movingImageItk->GetSpacing()[2] << std::endl;
  }

  ImageType::Pointer deformationImageItk = ImageType::New();



  /////////////////////////////////////////////////////////////////////////////////
  ///                                                                           ///
  ///                   Start of actual Featurelets Code                        ///
  ///                                                                           ///
  /////////////////////////////////////////////////////////////////////////////////

  //int x = 0;
  //x = ResampleVolumesToBeIsotropic(fixedImageItk);
  //x = ResampleVolumesToBeIsotropic(movingImageItk);

  //Initializing the Result Storage
  DeformationFieldGenerator resultStorage;

  ImageType::SizeType imageSizeFixed = fixedImageItk->GetLargestPossibleRegion().GetSize();
  ImageType::SizeType imageSizeMoving = movingImageItk->GetLargestPossibleRegion().GetSize();
  double iSize[3];
  iSize[0] = imageSizeMoving[0];
  iSize[1] = imageSizeMoving[1];
  iSize[2] = imageSizeMoving[2];
  double iSizeFix[3];
  iSizeFix[0] = imageSizeFixed[0];
  iSizeFix[1] = imageSizeFixed[1];
  iSizeFix[2] = imageSizeFixed[2];

  double originFix[Dimension];
  originFix[0] = fixedImageItk->GetOrigin()[0];
  originFix[1] = fixedImageItk->GetOrigin()[1];
  originFix[2] = fixedImageItk->GetOrigin()[2];

  ImageType::DirectionType directionFix = fixedImageItk->GetDirection();

  ImageType::SpacingType imageSpaceFixed = fixedImageItk->GetSpacing();
  double spaceFix[Dimension];
  spaceFix[0] = imageSpaceFixed[0];
  spaceFix[1] = imageSpaceFixed[1];
  spaceFix[2] = imageSpaceFixed[2];

  ImageType::SizeType FeatureletSize, SearchRegionSize, TempSize;
  if(rigid) {
    for(unsigned int j=0; j<Dimension; j++) {
      FeatureletSize[j] = imageSizeMoving[j];
      SearchRegionSize[j] = imageSizeFixed[j];
    }
  }
  else {
    for(unsigned int j=0; j<Dimension; j++) {
      FeatureletSize[j] = pnode->GetFeatureletsSize();
      SearchRegionSize[j] = pnode->GetSearchRegionSize();
    }
  }
  if(!rigid && zdifferent) {
    FeatureletSize[2] = pnode->GetFeatureletsSizeZ();
    SearchRegionSize[2] = pnode->GetSearchRegionSizeZ();
  }

  int FeatureletInicialSize[3];
  FeatureletInicialSize[0] = FeatureletSize.GetSize()[0];
  FeatureletInicialSize[1] = FeatureletSize.GetSize()[1];
  FeatureletInicialSize[2] = FeatureletSize.GetSize()[2];

  unsigned int numFeaturelets;
  numFeaturelets = ( imageSizeFixed[0]*imageSizeFixed[1]*imageSizeFixed[2] ) /
    (FeatureletSize[0]*FeatureletSize[1]*FeatureletSize[2]) *2;

  FeatureletRegistrationResultPointer *regResults;
  regResults = new FeatureletRegistrationResultPointer [numFeaturelets];

  int iFeatureletGridPos[3];
  iFeatureletGridPos[0] = iFeatureletGridPos[1] = iFeatureletGridPos[2] = 0;
  int noreg=0, iCurrentFeaturelet=0, regis=0, fixedreg=0;
  int totalnumber = (imageSizeFixed[0]*imageSizeFixed[1]*imageSizeFixed[2]) /
    (FeatureletSize[0]*FeatureletSize[1]*FeatureletSize[2]);

  ImageType::IndexType ImageIndex;
  bool fixed = true;
  Featurelet::Status statusFixed, statusMoving;
  FeatureletRegistrationResultPointer regResult;
  regResult = FeatureletRegistrationResultType::New();
  FeatureletRegistrationResultType *regPointer;

  int progress = 0;

  for (unsigned int i=0; i<imageSizeFixed[0]; i=i+FeatureletSize[0]) {
    iFeatureletGridPos[0]++;
    iFeatureletGridPos[1]=0;
    for (unsigned int j=0; j<imageSizeFixed[1]; j=j+FeatureletSize[1]) {
      iFeatureletGridPos[1]++;
      iFeatureletGridPos[2]=0;
      for (unsigned int k=0; k<imageSizeFixed[2]; k=k+FeatureletSize[2]) {
        this->UpdateFromMRMLScene();
        iFeatureletGridPos[2]++;
        iCurrentFeaturelet++;
        ImageIndex[0] = i;
        ImageIndex[1] = j;
        ImageIndex[2] = k;
        TempSize = FeatureletSize;      //for featurelets on the edges
        if (ImageIndex[0] + FeatureletSize[0] > imageSizeFixed[0])
            TempSize[0] = imageSizeFixed[0] - ImageIndex[0];
        if (ImageIndex[1] + FeatureletSize[1] > imageSizeFixed[1])
            TempSize[1] = imageSizeFixed[1] - ImageIndex[1];
        if (ImageIndex[2] + FeatureletSize[2] > imageSizeFixed[2])
            TempSize[2] = imageSizeFixed[2] - ImageIndex[2];

        // Output for additional Information
        if(debugMode) {
          std::cerr << std::endl << std::endl << std::endl;
          std::cerr << "+++++++   N E X T     F E A T U R E L E T   +++++++" << std::endl;
          std::cerr << "i-j-k: " <<i<< "-" <<j<< "-" <<k<< " -> Imagesize: " << imageSizeFixed << std::endl;
          std::cout << "registered " << regis << " featurelets so far " << std::endl;
          std::cout << "out of (approx.): " << totalnumber << " featurelets" << std::endl;
          std::cout << "Now registering featurelet no. " << iCurrentFeaturelet << " at ";
          std::cout << ImageIndex << " " << TempSize << std::endl;
          std::cout << std::endl;
        }

        SubsampleVolume(fixedImageItk, TempSize, SearchRegionSize, ImageIndex, fixed=true);
        SubsampleVolume(movingImageItk, TempSize, SearchRegionSize, ImageIndex, fixed=false);
        statusMoving = CheckFeaturelet(fixed=false, UseFiducialPoints, fiducialPoints);
        statusFixed  = CheckFeaturelet(fixed=true, UseFiducialPoints, fiducialPoints);

        if( statusFixed==Featurelet::OK && statusMoving==Featurelet::OK ) {
          if((!correl)&&(linear)) {
            regResult = RegisterVolumesI(debugMode, rigid);
            regis++;
          }
          if((correl)&&(!linear)) {
            regResult = RegisterVolumes1(fixedImageItk, movingImageItk, debugMode, rigid);
            regis++;
          }
          if((!correl)&&(!linear)) {
            regResult = RegisterVolumesII(debugMode, rigid);
            regis++;
          }
          if((correl)&&(linear)) {
            regResult = RegisterVolumes2(fixedImageItk, movingImageItk, debugMode, rigid);
            regis++;
          }
        }
        else if( statusFixed==Featurelet::fixedFeaturelet || statusMoving==Featurelet::fixedFeaturelet ) {
          fixedreg++;
          regResult = FeatureletRegistrationResultType::New();
        }
        else {
          noreg++;
          regResult = FeatureletRegistrationResultType::New();
        }

        regResult->SetStatusFixed(statusFixed);
        regResult->SetStatusMoving(statusMoving);
        regResult->SetFeatureletIndex(ImageIndex[0],ImageIndex[1],ImageIndex[2]);
        regResult->SetFeatureletSize(TempSize[0],TempSize[1],TempSize[2]);
        regResult->SetFeatureletGridPosition(iFeatureletGridPos[0],iFeatureletGridPos[1],iFeatureletGridPos[2]);
        regResult->SetImageSizeM(iSize[0], iSize[1], iSize[2]);
        regResult->SetFeatureletNumber(iCurrentFeaturelet);

        regResults[iCurrentFeaturelet] = regResult;
        regPointer = regResult.GetPointer();

        resultStorage.AddFeaturelet(regPointer);
      }
    }
    progress = (i*100)/imageSizeFixed[0];
    pnode->Setprogress(progress);
    //std::cerr << "Slice Number " << i << " of " << imageSizeFixed[0] << std::endl;
    std::cerr << "Progress: " << progress << "%" << std::endl;
    this->UpdateFromMRMLScene();
  }
  progress = 100;
  pnode->Setprogress(progress);
  this->UpdateFromMRMLScene();
  progress = 0;
  pnode->Setprogress(progress);

  resultStorage.SetImageSizeFix(iSizeFix);
  resultStorage.SetImageSpacingFix(spaceFix);
  resultStorage.SetImageOriginFix(originFix);
  resultStorage.SetImageDirectionFix(directionFix);
  resultStorage.SetNumberOfFeaturelets(iFeatureletGridPos);
  resultStorage.SetSizeOfFeaturelets(FeatureletInicialSize);
  resultStorage.GenerateDeformationField();
  resultStorage.RescaleDeformationField();
  resultStorage.TotalDeformationFieldCreator();
  resultStorage.WarpImagebyDeformationField(movingImageItk, fixedImageItk, deformationImageItk);

  vtkAbstractTransform *deformationFieldTransform =
    vtkAbstractTransform::SafeDownCast(resultStorage.VTKDeformationFieldCreator());

  std::cout << "Registration finished! " << std::endl;
  std::cout << "Registered " << regis << " featurelet(s)! " << std::endl;
  std::cout << "  " << fixedreg << " featurelet(s) not moving." << std::endl;
  std::cout << "  " << noreg << " featurelet(s) not registered." << std::endl;

  /// End of actual Featurelets Code
  /////////////////////////////////////////////////////////////////////////////////



  //Checking the ITK images
  if(debugMode) {
    double originFixed[Dimension];
    originFixed[0] = fixedImageItk->GetOrigin()[0];
    originFixed[1] = fixedImageItk->GetOrigin()[1];
    originFixed[2] = fixedImageItk->GetOrigin()[2];
    ImageType::SpacingType imageSpaceFixed = fixedImageItk->GetSpacing();
    double spaceFixed[Dimension];
    spaceFixed[0] = imageSpaceFixed[0];
    spaceFixed[1] = imageSpaceFixed[1];
    spaceFixed[2] = imageSpaceFixed[2];
    double originMoving[Dimension];
    originMoving[0] = movingImageItk->GetOrigin()[0];
    originMoving[1] = movingImageItk->GetOrigin()[1];
    originMoving[2] = movingImageItk->GetOrigin()[2];
    ImageType::SpacingType imageSpaceMoving = movingImageItk->GetSpacing();
    double spaceMoving[Dimension];
    spaceMoving[0] = imageSpaceMoving[0];
    spaceMoving[1] = imageSpaceMoving[1];
    spaceMoving[2] = imageSpaceMoving[2];
    std::cerr << "Spacing Fixed:  "<< spaceFixed[0] << "/" << spaceFixed[1] << "/" << spaceFixed[2] << std::endl;
    std::cerr << "Spacing Moving: "<< spaceMoving[0] << "/" << spaceMoving[1] << "/" << spaceMoving[2] << std::endl;
    std::cerr << "Origin Fixed:   "<< originFixed[0] << "/" << originFixed[1] << "/" << originFixed[2] << std::endl;
    std::cerr << "Origin Moving:  "<< originMoving[0] << "/" << originMoving[1] << "/" << originMoving[2] << std::endl;

    //Checking if the Featurelet- and SearchRegion-Size is the same as in de UI
    std::cerr << "FeatureletsSize: " << FeatureletSize[0] << " x " << FeatureletSize[1];
    std::cerr <<" x " << FeatureletSize[2] << std::endl;
    std::cerr << "SearchRegionSize: " << SearchRegionSize[0] << " x " << SearchRegionSize[1];
    std::cerr << " x " << SearchRegionSize[2] << std::endl;
  }

  vtkSmartPointer<vtkImageData> deformationImage = vtkSmartPointer<vtkImageData>::New();
  ImageType::RegionType region = deformationImageItk->GetBufferedRegion();
  ImageType::SizeType imageSize = region.GetSize();
  int extent[6]={0, (int) imageSize[0]-1, 0, (int) imageSize[1]-1, 0, (int) imageSize[2]-1};
  deformationImage->SetExtent(extent);
#if (VTK_MAJOR_VERSION <= 5)
  deformationImage->SetScalarType(VTK_FLOAT);
  deformationImage->SetNumberOfScalarComponents(1);
  deformationImage->AllocateScalars();
#else
  deformationImage->AllocateScalars(VTK_SHORT, 1);
  //deformationImage->AllocateScalars(VTK_INT, 1);
#endif

  PixelType* deformationPtr = (PixelType*)deformationImage->GetScalarPointer();
  itk::ImageRegionIteratorWithIndex< itk::Image<PixelType, Dimension> > itDeformationItk(
    deformationImageItk, deformationImageItk->GetLargestPossibleRegion() );
  for ( itDeformationItk.GoToBegin(); !itDeformationItk.IsAtEnd(); ++itDeformationItk ) {
    ImageType::IndexType i = itDeformationItk.GetIndex();
    (*deformationPtr) = deformationImageItk->GetPixel(i);
    deformationPtr++;
  }

  deformationField->SetAndObserveTransformFromParent(deformationFieldTransform);

  deformedImage->CopyOrientation(fixedImage);
  deformedImage->SetAndObserveImageData(deformationImage);

  if(debugMode) {
    imageData = deformedImage->GetImageData();
    std::cerr << "Deformed Pixeltype: " << imageData->GetScalarType() << " ";
    std::cerr << imageData->GetScalarTypeAsString() << std::endl;
    std::cerr << std::endl;
  }
  return 0;
}

//-----------------------------------------------------------------------------
int vtkSlicerRegistrationLogic::ShowVolume(vtkMRMLRegistrationNode* pnode, bool fixedImage) {
  vtkMRMLScene *scene = this->GetMRMLScene();

  vtkMRMLVolumeNode *fixedImageNode =
    vtkMRMLVolumeNode::SafeDownCast(scene->GetNodeByID(pnode->GetFixedImageNodeID()));
  vtkMRMLVolumeNode *movingImageNode =
    vtkMRMLVolumeNode::SafeDownCast(scene->GetNodeByID(pnode->GetMovingImageNodeID()));

  if(!fixedImageNode && !movingImageNode) {
    std::cerr << "Failed to look up input volume!" << std::endl;
    return -1;
  }
  if(fixedImage) {
    std::cerr << "Show Fixed Image Node" << std::endl;
    pnode->SetDeformedImageNodeID(fixedImageNode->GetID());
  }
  else {
    std::cerr << "Show Moving Image Node" << std::endl;
    pnode->SetDeformedImageNodeID(movingImageNode->GetID());
  }
  return 0;
}

//-----------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::RegisterNodes() {
  if(!this->GetMRMLScene()) {
    return;
  }
  vtkMRMLRegistrationNode* pNode = vtkMRMLRegistrationNode::New();
  this->GetMRMLScene()->RegisterNodeClass(pNode);
  pNode->Delete();
}

//---------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::UpdateFromMRMLScene() {
  assert(this->GetMRMLScene() != 0);
  this->Modified();
}

//---------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::OnMRMLSceneNodeAdded(vtkMRMLNode* node) {
  if(!node || !this->GetMRMLScene())
    return;

  if(node->IsA("vtkMRMLVolumeNode"))
    this->Modified();
}

//---------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::OnMRMLSceneNodeRemoved(vtkMRMLNode* node) {
  if(!node || !this->GetMRMLScene())
    return;

  if(node->IsA("vtkMRMLVolumeNode"))
    this->Modified();
}

//---------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::OnMRMLSceneEndImport() {
  // If we have a parameter node select it
  vtkMRMLRegistrationNode *paramNode = NULL;
  vtkMRMLNode *node = this->GetMRMLScene()->GetNthNodeByClass(0, "vtkMRMLRegistrationNode");
  if (node) {
    paramNode = vtkMRMLRegistrationNode::SafeDownCast(node);
    vtkSetAndObserveMRMLNodeMacro(this->RegistrationNode, paramNode);
  }
  this->Modified();
}

//---------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::OnMRMLSceneEndClose() {
  this->Modified();
}


// From Slicer-RT-commons
//---------------------------------------------------------------------------
bool vtkSlicerRegistrationLogic::ConvertVolumeNodeToItkImage(vtkMRMLVolumeNode* inVolumeNode,
                                                             ImageType::Pointer outItkVolume) {
  if ( inVolumeNode == NULL ) {
    std::cerr << "Failed to convert volume node to itk image - input MRML volume node is NULL!" << std::endl;
    return false;
  }
  vtkImageData* inVolume = inVolumeNode->GetImageData();
  if ( inVolume == NULL ) {
    std::cerr << "Failed to convert volume node to itk image - image in input MRML volume node is NULL!" << std::endl;
    return false;
  }
  if ( outItkVolume.IsNull() ) {
    std::cerr << "Failed to convert volume node to itk image - output image is NULL!" << std::endl;
    return false;
  }

  // Convert vtkImageData to itkImage
  vtkSmartPointer<vtkImageExport> imageExport = vtkSmartPointer<vtkImageExport>::New();
#if (VTK_MAJOR_VERSION <= 5)
  imageExport->SetInput(inVolume);
#else
  imageExport->SetInputData(inVolume);
#endif
  imageExport->Update();

  // Determine input volume to world transform
  vtkSmartPointer<vtkMatrix4x4> rasToWorldTransformMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
  vtkMRMLTransformNode* inTransformNode = inVolumeNode->GetParentTransformNode();
  if (inTransformNode!=NULL) {
    if (inTransformNode->IsTransformToWorldLinear()==0) {
      std::cerr << "Non-linear transform assigned to an input volume. Only linear transforms supported!" << std::endl;
      return false;
    }
    inTransformNode->GetMatrixTransformToWorld(rasToWorldTransformMatrix);
  }

  vtkSmartPointer<vtkMatrix4x4> inVolumeToRasTransformMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
  inVolumeNode->GetIJKToRASMatrix(inVolumeToRasTransformMatrix);

  vtkSmartPointer<vtkTransform> inVolumeToWorldTransform = vtkSmartPointer<vtkTransform>::New();
  inVolumeToWorldTransform->Identity();
  inVolumeToWorldTransform->PostMultiply();
  inVolumeToWorldTransform->Concatenate(inVolumeToRasTransformMatrix);
  inVolumeToWorldTransform->Concatenate(rasToWorldTransformMatrix);

  // Set ITK image properties
  double outputSpacing[3] = {0.0, 0.0, 0.0};
  inVolumeToWorldTransform->GetScale(outputSpacing);
  outItkVolume->SetSpacing(outputSpacing);
  //std::cerr << "OutITKSpacing: " << outputSpacing[0] << outputSpacing[1] << outputSpacing[2] << std::endl;
  //std::cerr << "OutITKSpacing: " << outItkVolume->GetSpacing()[0] << outItkVolume->GetSpacing()[1];
  //std::cerr << outItkVolume->GetSpacing()[2] << std::endl;

  double outputOrigin[3] = {0.0, 0.0, 0.0};
  inVolumeToWorldTransform->GetPosition(outputOrigin);
  outItkVolume->SetOrigin(outputOrigin);

  double outputOrienationAngles[3] = {0.0, 0.0, 0.0};
  inVolumeToWorldTransform->GetOrientation(outputOrienationAngles);
  vtkSmartPointer<vtkTransform> inVolumeToWorldOrientationTransform = vtkSmartPointer<vtkTransform>::New();
  inVolumeToWorldOrientationTransform->Identity();
  inVolumeToWorldOrientationTransform->RotateX(outputOrienationAngles[0]);
  inVolumeToWorldOrientationTransform->RotateY(outputOrienationAngles[1]);
  inVolumeToWorldOrientationTransform->RotateZ(outputOrienationAngles[2]);
  vtkSmartPointer<vtkMatrix4x4> inVolumeToWorldOrientationTransformMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
  inVolumeToWorldOrientationTransform->GetMatrix(inVolumeToWorldOrientationTransformMatrix);
  itk::Matrix<double,3,3> outputDirectionMatrix;
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      outputDirectionMatrix[i][j] = inVolumeToWorldOrientationTransformMatrix->GetElement(i,j);
      //std::cerr << "Output Direction: " << i << "," << j << " = " << outputDirectionMatrix[i][j] << std::endl;
    }
  }
  outItkVolume->SetDirection(outputDirectionMatrix);

  int inputExtent[6]={0,0,0,0,0,0};
  inVolume->GetExtent(inputExtent);
  ImageType::SizeType inputSize;
  inputSize[0] = inputExtent[1] - inputExtent[0] + 1;
  inputSize[1] = inputExtent[3] - inputExtent[2] + 1;
  inputSize[2] = inputExtent[5] - inputExtent[4] + 1;

  ImageType::IndexType start;
  start[0] = start[1] = start[2] = 0.0;

  ImageType::RegionType region;
  region.SetSize(inputSize);
  region.SetIndex(start);
  outItkVolume->SetRegions(region);

  try {
    outItkVolume->Allocate();
  }
  catch(itk::ExceptionObject & err) {
    std::cerr << "Failed to allocate memory for the image conversion: " << err.GetDescription() << std::endl;
    return false;
  }
  imageExport->Export( outItkVolume->GetBufferPointer() );
  return true;
}


int vtkSlicerRegistrationLogic::ResampleVolumesToBeIsotropic(ImageType* Image) {
  typedef float InternalPixelType;
  typedef itk::Image< InternalPixelType, Dimension > InternalImageType;
  typedef itk::RecursiveGaussianImageFilter< ImageType, InternalImageType > GaussianFilterType;
  typedef itk::RecursiveGaussianImageFilter< InternalImageType, InternalImageType > GaussianFilterType1;
  GaussianFilterType::Pointer smootherX = GaussianFilterType::New();
  GaussianFilterType1::Pointer smootherY = GaussianFilterType1::New();
    smootherX->SetInput( Image );
    smootherY->SetInput( smootherX->GetOutput() );
  ImageType::Pointer inputImage = Image;
  const ImageType::SpacingType& inputSpacing = inputImage->GetSpacing();
  const double isoSpacing = sqrt( inputSpacing[2] * inputSpacing[0] );
    smootherX->SetSigma( isoSpacing );
    smootherY->SetSigma( isoSpacing );
    smootherX->SetDirection( 0 );
    smootherY->SetDirection( 1 );
    smootherX->SetNormalizeAcrossScale( true );
    smootherY->SetNormalizeAcrossScale( true );
  typedef itk::ResampleImageFilter< InternalImageType, ImageType > ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  typedef itk::IdentityTransform< double, Dimension > TransformType;
  TransformType::Pointer transform = TransformType::New();
    transform->SetIdentity();
    resampler->SetTransform( transform );
  typedef itk::MinimumMaximumImageCalculator< ImageType > CalculatorType;
  CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetImage( Image );
  try {
    calculator->Compute();
  }
  catch (itk::ExceptionObject& e) {
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
  }

  typedef itk::LinearInterpolateImageFunction< InternalImageType, double > InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
    resampler->SetInterpolator( interpolator );
    resampler->SetDefaultPixelValue( calculator->GetMinimum() ); // highlight regions without source
  ImageType::SpacingType spacing;
  spacing[0] = isoSpacing;
  spacing[1] = isoSpacing;
  spacing[2] = isoSpacing;
  resampler->SetOutputSpacing( spacing );
  resampler->SetOutputOrigin( inputImage->GetOrigin() );
  ImageType::SizeType inputSize = inputImage->GetLargestPossibleRegion().GetSize();

  typedef ImageType::SizeType::SizeValueType SizeValueType;
  const double dx = inputSize[0] * inputSpacing[0] / isoSpacing;
  const double dy = inputSize[1] * inputSpacing[1] / isoSpacing;
  const double dz = (inputSize[2] - 1 ) * inputSpacing[2] / isoSpacing;

  ImageType::SizeType size;
  size[0] = static_cast < SizeValueType>( dx );
  size[1] = static_cast < SizeValueType>( dy );
  size[2] = static_cast < SizeValueType>( dz );

  resampler->SetSize( size );
  resampler->SetInput( smootherY->GetOutput() );
  try {
    std::cout << "Resample volume to be isotropic" << std::endl;
    resampler->Update();
  }
  catch ( itk::ExceptionObject& e ) {
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
  }
  Image = resampler->GetOutput();
  return EXIT_SUCCESS;
}

int vtkSlicerRegistrationLogic::SubsampleVolume(const ImageType::Pointer Image,
                                                ImageType::SizeType FeatureletSize,
                                                ImageType::SizeType SearchRegionSize,
                                                ImageType::IndexType ImageIndex,
                                                bool fixed) {
  ImageType::RegionType desiredRegion;
    desiredRegion.SetSize(FeatureletSize);
    desiredRegion.SetIndex(ImageIndex);

  // Region of the Featurelet should be in the center of the Search Region
  /*ImageType::IndexType ImageIndexSearchRegion;
    ImageIndexSearchRegion = ImageIndex;
  int indexDifference[3];
    indexDifference[0] = (SearchRegionSize[0] - FeatureletSize[0])/2;
    indexDifference[1] = (SearchRegionSize[1] - FeatureletSize[1])/2;
    indexDifference[2] = (SearchRegionSize[2] - FeatureletSize[2])/2;
  //std::cerr << indexDifference[0] << "|" << indexDifference[1] << "|" << indexDifference[2] << std::endl;
  if( (ImageIndex[0] - indexDifference[0]) > 0 )
    ImageIndexSearchRegion[0] = ImageIndex[0] - indexDifference[0];
  if( (ImageIndex[1] - indexDifference[1]) > 0 )
    ImageIndexSearchRegion[1] = ImageIndex[1] - indexDifference[1];
  if( (ImageIndex[2] - indexDifference[2]) > 0 )
    ImageIndexSearchRegion[2] = ImageIndex[2] - indexDifference[2];*/

  ImageType::RegionType SearchingRegion;
    SearchingRegion.SetSize(SearchRegionSize);
    SearchingRegion.SetIndex(ImageIndex); //ImageIndexSearchRegion

  typedef itk::ExtractImageFilter<ImageType, ImageType> FilterType;
  FilterType::Pointer filterFeaturelet = FilterType::New();
  FilterType::Pointer filterSearchRegion = FilterType::New();
    filterFeaturelet->SetExtractionRegion(desiredRegion);
    filterFeaturelet->SetInput(Image);
    filterSearchRegion->SetExtractionRegion(SearchingRegion);
    filterSearchRegion->SetInput(Image);
  try {
    filterFeaturelet->Update();
    filterSearchRegion->Update();
  }
  catch (itk::ExceptionObject& e ) {
    //std::cerr << e << std::endl;
    return EXIT_FAILURE;
  }
  if (fixed){
    FeatureletPointerFixed = filterFeaturelet->GetOutput();
    SearchRegionPointerFixed = filterSearchRegion->GetOutput();
    //std::cerr << "SubsampleVolume: FeatureletPointer_f:   " << FeatureletPointerFixed << std::endl;
    //std::cerr << "SubsampleVolume: SearchRegionPointer_f: " << SearchRegionPointerFixed << std::endl;
  }
  else{
    FeatureletPointerMoving = filterFeaturelet->GetOutput();
    SearchRegionPointerMoving = filterSearchRegion->GetOutput();
    //std::cerr << "SubsampleVolume: FeatureletPointer_m:   " << FeatureletPointerMoving << std::endl;
    //std::cerr << "SubsampleVolume: SearchRegionPointer_m: " << SearchRegionPointerMoving << std::endl;
  }
  return EXIT_SUCCESS;
}


Featurelet::Status vtkSlicerRegistrationLogic::CheckFeaturelet(bool fixed, bool FiducialPoints,
                                                               vtkMRMLMarkupsFiducialNode *Fiducials) {
  ImageType::SizeType FeatureletSize;
  if (fixed) {
    FeatureletSize = FeatureletPointerFixed->GetLargestPossibleRegion().GetSize();
    //std::cerr << "CheckFeaturelet: FeatureletSize_f: " << FeatureletSize << std::endl;
  }
  else {
    FeatureletSize = FeatureletPointerMoving->GetLargestPossibleRegion().GetSize();
    //std::cerr << "CheckFeaturelet: FeatureletSize_m: " << FeatureletSize << std::endl;
  }
  //typedef itk::NormalizedCorrelationImageToImageMetric<ImageType, ImageType> MetricType;
    //MetricType::Pointer metric = MetricType::New();       //Not used
    if (FeatureletSize[0] < 4)
      return Featurelet::SizeError;
    if (FeatureletSize[1] < 4)
      return Featurelet::SizeError;
    if (FeatureletSize[2] < 2)
      return Featurelet::SizeError;

    typedef itk::MinimumMaximumImageCalculator<ImageType> CalculatorType;
    CalculatorType::Pointer calculator = CalculatorType::New();
    if (fixed) {
      calculator->SetImage(FeatureletPointerFixed);
    }
    else {
      calculator->SetImage(FeatureletPointerMoving);
    }
    try {
      calculator->Compute();
    }
    catch (itk::ExceptionObject& e) {
      std::cerr << e << std::endl;
      return Featurelet::internalError;
    }

    if ( calculator->GetMaximum() == calculator->GetMinimum() )
      return Featurelet::sameColor;

    if (!FiducialPoints)
      return Featurelet::OK;
    else {
      typedef itk::Point<double, 3> PointType;
      ImageType::Pointer FiducialImage = ImageType::New();
      if (fixed)
        FiducialImage = FeatureletPointerFixed;
      else
        FiducialImage = FeatureletPointerMoving;
      int TotalNumber = (int) Fiducials->GetNumberOfFiducials();
      double PositionRAS[3] = {0.0, 0.0, 0.0};
      PointType PositionLPS;
      ImageType::IndexType pixelIndex;
      for (int n=0; n<TotalNumber; n++) {
        Fiducials->GetNthFiducialPosition(n, PositionRAS);
        PositionLPS[0] = PositionRAS[0];
        PositionLPS[1] = PositionRAS[1];
        PositionLPS[2] = PositionRAS[2];
        const bool isInside = FiducialImage->TransformPhysicalPointToIndex(PositionLPS, pixelIndex);
        if(isInside)
          return Featurelet::fixedFeaturelet;
      }
    }
    return Featurelet::OK;
}


// Registration for Correlation and linear interpolation
FeatureletRegistrationResult::Pointer vtkSlicerRegistrationLogic::RegisterVolumes1
          (ImageType* FixedImage, ImageType* MovingImage, bool debugMode, bool rigid) {
  FeatureletRegistrationResultPointer regResult;
  regResult = FeatureletRegistrationResultType::New();

  typedef itk::TranslationTransform< double > TransformType;
  TransformType::Pointer transform = TransformType::New();

  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();

  typedef itk::NormalizedCorrelationImageToImageMetric<ImageType, ImageType> MetricType;
  MetricType::Pointer metric = MetricType::New();

  typedef itk::LinearInterpolateImageFunction< ImageType, double > InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  typedef itk::ImageRegistrationMethod< ImageType, ImageType > RegistrationType;
  RegistrationType::Pointer registration = RegistrationType::New();

  optimizer->MinimizeOn();
  registration->SetTransform(transform);
  registration->SetMetric(metric);
  registration->SetOptimizer(optimizer);
  registration->SetInterpolator(interpolator);

  ImageType::SizeType size = MovingImage->GetLargestPossibleRegion().GetSize();
  ImageType::RegionType region;
  region.SetSize(size);
  registration->SetFixedImage(FixedImage);
  registration->SetMovingImage(FeatureletPointerMoving);
  registration->SetFixedImageRegion(SearchRegionPointerFixed->GetBufferedRegion());

  typedef RegistrationType::ParametersType ParametersType;
  ParametersType initialParameters(transform->GetNumberOfParameters());
  initialParameters[0] = 0.0;  //Initial offset in mm along X
  initialParameters[1] = 0.0;  //Initial offset in mm along Y
  initialParameters[2] = 0.0;  //Initial offset in mm along Z
  registration->SetInitialTransformParameters(initialParameters);

  optimizer->SetMaximumStepLength(MaxStepLength);
  optimizer->SetMinimumStepLength(MinStepLength);
  optimizer->SetNumberOfIterations(NumberIterations);

  if(debugMode || rigid) {
    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    std::cout << "*** Optimization starts ***" << std::endl;
    optimizer->AddObserver( itk::IterationEvent(), observer);
  }

  try {
    registration->Update();
  }
  catch(itk::ExceptionObject& e) {
    std::cerr << e << std::endl;
    regResult->SetStatusRegistration(EXIT_FAILURE);
    return regResult;  //Identificar este error como de regstration error!!!
  }

  OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();
  transform->SetParameters( finalParameters );
  //std::cout << "*** Optimization done ***" << std::endl;
  //TransformType::OutputVectorType offset = transform->GetOffset();
  //const double bestValue = optimizer->GetValue();
  //std::cout << "  Offset1 = " << offset << std::endl;
  //std::cout << "Metric Value= "<< bestValue << std::endl;

  if(true) {
    regResult->SetTransform(transform);
    regResult->SetOptimizer(optimizer);
  }
  return regResult;
}

// Registration for Correlation and neighbour interpolation
FeatureletRegistrationResult::Pointer vtkSlicerRegistrationLogic::RegisterVolumes2
          (ImageType* FixedImage, ImageType* MovingImage, bool debugMode, bool rigid) {
  FeatureletRegistrationResultPointer regResult;
  regResult = FeatureletRegistrationResultType::New();

  typedef itk::TranslationTransform< double > TransformType;
  TransformType::Pointer transform = TransformType::New();

  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();

  typedef itk::NormalizedCorrelationImageToImageMetric< ImageType, ImageType > MetricType;
  MetricType::Pointer metric = MetricType::New();

  typedef itk::NearestNeighborInterpolateImageFunction< ImageType, double > InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  typedef itk::ImageRegistrationMethod< ImageType, ImageType > RegistrationType;
  RegistrationType::Pointer registration = RegistrationType::New();
  optimizer->MinimizeOn();
  registration->SetTransform( transform );
  registration->SetMetric( metric );
  registration->SetOptimizer( optimizer );
  registration->SetInterpolator( interpolator );

  ImageType::SizeType size = MovingImage->GetLargestPossibleRegion().GetSize();
  ImageType::RegionType region;
  region.SetSize(size);
  registration->SetFixedImage(FixedImage);
  registration->SetMovingImage(FeatureletPointerMoving);
  registration->SetFixedImageRegion(SearchRegionPointerFixed->GetBufferedRegion());

  typedef RegistrationType::ParametersType ParametersType;
  ParametersType initialParameters(transform->GetNumberOfParameters());
  initialParameters[0] = 0.0;  //Initial offset in mm along X
  initialParameters[1] = 0.0;  //Initial offset in mm along Y
  initialParameters[2] = 0.0;  //Initial offset in mm along Z

  registration->SetInitialTransformParameters(initialParameters);
  optimizer->SetMaximumStepLength(MaxStepLength);
  optimizer->SetMinimumStepLength(MinStepLength);
  optimizer->SetNumberOfIterations(NumberIterations);

  if(debugMode || rigid) {
    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    std::cout << "*** Optimization starts ***" << std::endl;
    optimizer->AddObserver(itk::IterationEvent(), observer);
  }

  try {
    registration->Update();
  }
  catch(itk::ExceptionObject& e) {
    std::cerr << e << std::endl;
    regResult->SetStatusRegistration(EXIT_FAILURE);
    return regResult;  //Identificar este error como de regstration error!!!
  }

  OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();
  transform->SetParameters(finalParameters);
  //std::cout << "*** Optimization done ***" << std::endl;
  //TransformType::OutputVectorType offset = transform->GetOffset();
  //const double bestValue = optimizer->GetValue();
  //std::cout << "  Offset1 = " << offset << std::endl;
  //std::cout << "Metric Value= "<< bestValue<<std::endl;
  if(true) {
    regResult->SetTransform(transform);
    regResult->SetOptimizer(optimizer);
  }
  return regResult;
}

// Registration for Mutual Information and neighbour interpolation
FeatureletRegistrationResult::Pointer vtkSlicerRegistrationLogic::RegisterVolumesI
          (bool debugMode, bool rigid) {
  FeatureletRegistrationResultPointer regResult;
  regResult = FeatureletRegistrationResultType::New();
  typedef double InternalPixelType;
  typedef itk::Image< InternalPixelType, 3 > InternalImageType;
  typedef itk::TranslationTransform< double, 3 > TransformType;
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef itk::NearestNeighborInterpolateImageFunction <InternalImageType, double> InterpolatorType;
  typedef itk::ImageRegistrationMethod<InternalImageType, InternalImageType> RegistrationType;
  typedef itk::MutualInformationImageToImageMetric<InternalImageType, InternalImageType> MetricType;

  TransformType::Pointer      transform     = TransformType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();
        registration->SetOptimizer(     optimizer     );
        registration->SetTransform(     transform     );
        registration->SetInterpolator(  interpolator  );
  MetricType::Pointer         metric        = MetricType::New();
        registration->SetMetric(metric);
        metric->SetFixedImageStandardDeviation(  0.4 );
        metric->SetMovingImageStandardDeviation( 0.4 );

  typedef itk::NormalizeImageFilter<ImageType, InternalImageType> FixedNormalizeFilterType;
  typedef itk::NormalizeImageFilter<ImageType,InternalImageType> MovingNormalizeFilterType;
  FixedNormalizeFilterType::Pointer fixedNormalizer = FixedNormalizeFilterType::New();
  MovingNormalizeFilterType::Pointer movingNormalizer = MovingNormalizeFilterType::New();
  typedef itk::DiscreteGaussianImageFilter<InternalImageType, InternalImageType> GaussianFilterType;
  GaussianFilterType::Pointer fixedSmoother  = GaussianFilterType::New();
  GaussianFilterType::Pointer movingSmoother = GaussianFilterType::New();
  fixedSmoother->SetVariance(  2.0 );
  movingSmoother->SetVariance( 2.0 );
  fixedNormalizer->SetInput(FeatureletPointerFixed);
  movingNormalizer->SetInput(FeatureletPointerMoving);
  fixedSmoother->SetInput(fixedNormalizer->GetOutput() );
  movingSmoother->SetInput(movingNormalizer->GetOutput() );

       registration->SetFixedImage(    fixedSmoother->GetOutput()    );
       registration->SetMovingImage(   movingSmoother->GetOutput()   );

       fixedNormalizer->Update();
       ImageType::RegionType fixedImageRegion = fixedNormalizer->GetOutput()->GetBufferedRegion();
       registration->SetFixedImageRegion(SearchRegionPointerFixed->GetBufferedRegion());

  typedef RegistrationType::ParametersType ParametersType;
  ParametersType initialParameters( transform->GetNumberOfParameters() );
       initialParameters[0] = 0.0;  //Initial offset in mm along X
       initialParameters[1] = 0.0;  //Initial offset in mm along Y
       initialParameters[2] = 0.0;  //Initial offset in mm along Z

       registration->SetInitialTransformParameters(initialParameters);

  const unsigned int numberOfPixels = fixedImageRegion.GetNumberOfPixels();
  const unsigned int numberOfSamples = static_cast< unsigned int >( numberOfPixels * 0.08 );
       metric->SetNumberOfSpatialSamples( numberOfSamples );
       optimizer->SetNumberOfIterations(NumberIterations);
       optimizer->MaximizeOn();
       optimizer->SetMaximumStepLength(MaxStepLength);
       optimizer->SetMinimumStepLength(MinStepLength);

  if(debugMode || rigid) {
    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    std::cout << "*** Optimization starts ***" << std::endl;
    optimizer->AddObserver( itk::IterationEvent(), observer);
  }

  try {
       registration->Update();
  }
  catch(itk::ExceptionObject& e) {
       std::cerr << e << std::endl;
       regResult->SetStatusRegistration(EXIT_FAILURE);
       return regResult;
  }

  OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();
       transform->SetParameters(finalParameters);
  //std::cout << "*** Optimization done ***" << std::endl;
  //TransformType::OutputVectorType offset = transform->GetOffset();
  //const unsigned int numberOfIterations = optimizer->GetCurrentIteration();     //not used
  //const double bestValue = optimizer->GetValue();
  //std::cout << "  Offset1 = " << offset << std::endl;
  //std::cout << "Metric Value = " << bestValue << std::endl;
  if(true) {
       regResult->SetTransform(transform);
       regResult->SetOptimizer(optimizer);
  }
  return regResult;
}

// Registration for Mutual Information and linear interpolation
FeatureletRegistrationResult::Pointer vtkSlicerRegistrationLogic::RegisterVolumesII
          (bool debugMode, bool rigid) {
  FeatureletRegistrationResultPointer regResult;
  regResult = FeatureletRegistrationResultType::New();
  typedef double InternalPixelType;
  typedef itk::Image<InternalPixelType, 3> InternalImageType;
  typedef itk::TranslationTransform< double, 3 > TransformType;
  typedef itk::RegularStepGradientDescentOptimizer  OptimizerType;
  typedef itk::Function::CosineWindowFunction<4, double, double> cosine;
  typedef itk::ConstantBoundaryCondition< InternalImageType > boundaries;
  typedef itk::WindowedSincInterpolateImageFunction< InternalImageType, 4,
          cosine, boundaries, double > InterpolatorType;
  typedef itk::ImageRegistrationMethod<InternalImageType, InternalImageType > RegistrationType;
  typedef itk::MutualInformationImageToImageMetric<InternalImageType, InternalImageType > MetricType;
  TransformType::Pointer      transform     = TransformType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();
      registration->SetOptimizer(     optimizer     );
      registration->SetTransform(     transform     );
      registration->SetInterpolator(  interpolator  );
  MetricType::Pointer metric = MetricType::New();
      registration->SetMetric( metric  );
  typedef itk::NormalizeImageFilter<ImageType, InternalImageType> FixedNormalizeFilterType;
  typedef itk::NormalizeImageFilter<ImageType, InternalImageType> MovingNormalizeFilterType;
  FixedNormalizeFilterType::Pointer fixedNormalizer = FixedNormalizeFilterType::New();
  MovingNormalizeFilterType::Pointer movingNormalizer = MovingNormalizeFilterType::New();

  typedef itk::DiscreteGaussianImageFilter<InternalImageType, InternalImageType> GaussianFilterType;
  GaussianFilterType::Pointer fixedSmoother  = GaussianFilterType::New();
  GaussianFilterType::Pointer movingSmoother = GaussianFilterType::New();
      fixedSmoother->SetVariance( 2.0 );
      movingSmoother->SetVariance( 2.0 );
      fixedNormalizer->SetInput(FeatureletPointerFixed);
      movingNormalizer->SetInput(FeatureletPointerMoving);
      fixedSmoother->SetInput(fixedNormalizer->GetOutput() );
      movingSmoother->SetInput(movingNormalizer->GetOutput() );
      registration->SetFixedImage(fixedSmoother->GetOutput() );
      registration->SetMovingImage(movingSmoother->GetOutput() );
      fixedNormalizer->Update();

  //ImageType::RegionType fixedImageRegion = fixedNormalizer->GetOutput()->GetBufferedRegion();

    registration->SetFixedImageRegion(SearchRegionPointerFixed->GetBufferedRegion() );

  typedef RegistrationType::ParametersType ParametersType;
  ParametersType initialParameters( transform->GetNumberOfParameters() );
      initialParameters[0] = 0.0;  //Initial offset in mm along X
      initialParameters[1] = 0.0;  //Initial offset in mm along Y
      initialParameters[2] = 0.0;  //Initial offset in mm along Z
      registration->SetInitialTransformParameters( initialParameters );

      optimizer->SetNumberOfIterations(NumberIterations);
      optimizer->MaximizeOn();
      optimizer->SetMaximumStepLength(MaxStepLength);
      optimizer->SetMinimumStepLength(MinStepLength);

  if(debugMode || rigid) {
    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    std::cout << "*** Optimization starts ***" << std::endl;
    optimizer->AddObserver( itk::IterationEvent(), observer );
  }

  try {
      registration->Update();
  }
  catch(itk::ExceptionObject& e) {
      std::cerr << e << std::endl;
      regResult->SetStatusRegistration(EXIT_FAILURE);
      return regResult;
  }

  OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();
      transform->SetParameters( finalParameters );
  //std::cout << "*** Optimization done ***" << std::endl;
  //TransformType::OutputVectorType offset = transform->GetOffset();
  //const unsigned int numberOfIterations = optimizer->GetCurrentIteration();     //not used
  //const double bestValue = optimizer->GetValue();
  //std::cout << "  Offset1 = " << offset << std::endl;
  //std::cout << "Metric Value= "<< bestValue<<std::endl;

  if(true) {
      regResult->SetTransform(transform);
      regResult->SetOptimizer(optimizer);
  }
  return regResult;
}
