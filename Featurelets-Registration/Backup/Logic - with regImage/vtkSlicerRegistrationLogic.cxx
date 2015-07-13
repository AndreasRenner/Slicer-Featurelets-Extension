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
#include "vtkSlicerRegistrationLogic.h"
#include "vtkMRMLRegistrationNode.h"
#include "vtkSlicerVolumesLogic.h"

// MRML includes
#include <vtkMRMLAnnotationROINode.h>
#include <vtkMRMLDiffusionWeightedVolumeDisplayNode.h>
#include <vtkMRMLDiffusionWeightedVolumeNode.h>
#include <vtkMRMLDiffusionTensorVolumeNode.h>
#include <vtkMRMLLinearTransformNode.h>
#include <vtkMRMLRegistrationNode.h>
#include <vtkMRMLScalarVolumeDisplayNode.h>
#include <vtkMRMLScene.h>
#include <vtkMRMLSelectionNode.h>
#include <vtkMRMLVectorVolumeDisplayNode.h>
#include <vtkMRMLVectorVolumeNode.h>
#include <vtkMRMLVolumeNode.h>
#include <vtkMRMLTransformNode.h>

// VTK includes
#include <vtkImageData.h>
#include <vtkImageClip.h>
#include <vtkImageExport.h>
#include <vtkIntArray.h>
#include <vtkMatrix4x4.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>
#include <vtkVersion.h>
#include <vtkTransform.h>

// ITK includes
#include <itkImage.h>
#include <itkImageBase.h>
#include <itkImageDuplicator.h>
#include <itkImageRegionIteratorWithIndex.h>

// Featurelets includes
#include "DeformationFieldGenerator.h"
#include "FeatureletRegistrationResult.h"
#include "regImage.h"
#include "RegisterVolumes.h"
#include "ResampleVolume.h"
//#include "subsampleFeut.h"

// STD includes
#include <cassert>
#include <iostream>

//----------------------------------------------------------------------------
class vtkSlicerRegistrationLogic::vtkInternal
{
public:
  vtkInternal();

  vtkSlicerVolumesLogic* VolumesLogic;
};

//----------------------------------------------------------------------------
vtkSlicerRegistrationLogic::vtkInternal::vtkInternal()
{
  this->VolumesLogic = 0;
}

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerRegistrationLogic);

//----------------------------------------------------------------------------
vtkSlicerRegistrationLogic::vtkSlicerRegistrationLogic()
{
    this->Internal = new vtkInternal;
    this->RegistrationNode = NULL;
}

//----------------------------------------------------------------------------
vtkSlicerRegistrationLogic::~vtkSlicerRegistrationLogic()
{
    delete this->Internal;
    this->SetAndObserveRegistrationNode(NULL); //release the node object to avoid memory leaks
}

//----------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::SetVolumesLogic(vtkSlicerVolumesLogic* logic)
{
  this->Internal->VolumesLogic = logic;
}

//----------------------------------------------------------------------------
vtkSlicerVolumesLogic* vtkSlicerRegistrationLogic::GetVolumesLogic()
{
  return this->Internal->VolumesLogic;
}

//----------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::SetAndObserveRegistrationNode(vtkMRMLRegistrationNode *node)
{
  vtkSetAndObserveMRMLNodeMacro(this->RegistrationNode, node);
}

//---------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndImportEvent);
  events->InsertNextValue(vtkMRMLScene::EndCloseEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEvents(newScene, events.GetPointer());
}

//----------------------------------------------------------------------------
int vtkSlicerRegistrationLogic::RunClicked(vtkMRMLRegistrationNode* pnode)
{
  std::cerr << "*** Entered RunClicked ***" << std::endl;

  vtkMRMLScene *scene = this->GetMRMLScene();

  std::cerr << "Fixed: " << pnode->GetFixedImageNodeID() << std::endl;
  std::cerr << "Moving: " << pnode->GetMovingImageNodeID() << std::endl;
  std::cerr << "Deformed: " << pnode->GetDeformedImageNodeID() << std::endl;

  vtkMRMLVolumeNode *fixedImage =
    vtkMRMLVolumeNode::SafeDownCast(scene->GetNodeByID(pnode->GetFixedImageNodeID()));
  vtkMRMLVolumeNode *movingImage =
    vtkMRMLVolumeNode::SafeDownCast(scene->GetNodeByID(pnode->GetMovingImageNodeID()));
  vtkMRMLVolumeNode *deformedImage = NULL;
    //      vtkMRMLVolumeNode::SafeDownCast(scene->GetNodeByID(pnode->GetDeformedImageNodeID()));

  vtkImageData* imageData = fixedImage->GetImageData();
  std::cerr << "Fixed - Pixeltype: " << imageData->GetScalarType() << " " << imageData->GetScalarTypeAsString() << std::endl;
  imageData = movingImage->GetImageData();
  std::cerr << "Moving - Pixeltype: " << imageData->GetScalarType() << " " << imageData->GetScalarTypeAsString() << std::endl;

  if(!fixedImage || !movingImage)
    {
    std::cerr << "Failed to look up fixed/moving image!" << std::endl;
    return -1;
    }

  // make sure inputs are initialized
  if(!fixedImage || !movingImage)
    {
    std::cerr << "Registration: Inputs are not initialized" << std::endl;
    return -1;
    }

  // check the volume type
  vtkMRMLDiffusionTensorVolumeNode *dtvnode= vtkMRMLDiffusionTensorVolumeNode::SafeDownCast(movingImage);
  vtkMRMLDiffusionWeightedVolumeNode *dwvnode= vtkMRMLDiffusionWeightedVolumeNode::SafeDownCast(movingImage);
  vtkMRMLVectorVolumeNode *vvnode= vtkMRMLVectorVolumeNode::SafeDownCast(movingImage);
  vtkMRMLScalarVolumeNode *svnode = vtkMRMLScalarVolumeNode::SafeDownCast(movingImage);

  if(!this->Internal->VolumesLogic)
    {
      std::cerr << "Registration: ERROR: failed to get hold of Volumes logic" << std::endl;
      return -2;
    }

  std::ostringstream outSS;
  outSS << movingImage->GetName() << "_deformed";

  if(dtvnode)
    {
    std::cerr << "Registration: ERROR: Diffusion tensor volumes are not supported by this module!" << std::endl;
    return -2;
    }
  // need to create clones and display nodes here, since
  // VolumesLogic::CloneVolume() handles only ScalarVolumeNode's
  else if(dwvnode)
    {
    vtkNew<vtkMRMLDiffusionWeightedVolumeNode> outputDWVNode;
    outputDWVNode->CopyWithScene(dwvnode);
    vtkNew<vtkMRMLDiffusionWeightedVolumeDisplayNode> dwiDisplayNode;
    dwiDisplayNode->CopyWithScene(dwvnode->GetDisplayNode());
    scene->AddNode(dwiDisplayNode.GetPointer());

    vtkNew<vtkImageData> outputImageData;
    outputImageData->DeepCopy(dwvnode->GetImageData());
    outputDWVNode->SetAndObserveImageData(outputImageData.GetPointer());

    outputDWVNode->SetAndObserveDisplayNodeID(dwiDisplayNode->GetID());
    outputDWVNode->SetAndObserveStorageNodeID(NULL);
    scene->AddNode(outputDWVNode.GetPointer());

    deformedImage = outputDWVNode.GetPointer();
    }
  else if(vvnode)
    {
    vtkNew<vtkMRMLVectorVolumeNode> outputVVNode;
    outputVVNode->CopyWithScene(dwvnode);
    vtkNew<vtkMRMLVectorVolumeDisplayNode> vvDisplayNode;
    vvDisplayNode->CopyWithScene(vvnode->GetDisplayNode());
    scene->AddNode(vvDisplayNode.GetPointer());

    vtkNew<vtkImageData> outputImageData;
    outputImageData->DeepCopy(vvnode->GetImageData());
    outputVVNode->SetAndObserveImageData(outputImageData.GetPointer());

    outputVVNode->SetAndObserveDisplayNodeID(vvDisplayNode->GetID());
    outputVVNode->SetAndObserveStorageNodeID(NULL);
    scene->AddNode(outputVVNode.GetPointer());

    deformedImage = outputVVNode.GetPointer();
    }
  else if(svnode)
    {
    deformedImage = vtkMRMLVolumeNode::SafeDownCast(scene->GetNodeByID(pnode->GetDeformedImageNodeID()));
    //        this->Internal->VolumesLogic->CloneVolume(this->GetMRMLScene(), movingImage, outSS.str().c_str());
    }
  else
    {
    std::cerr << "Moving Image not recognized!" << std::endl;
    return -1;
    }

  deformedImage->SetName(outSS.str().c_str());


  //Convert images to ITK
  regImage* fixedImageItk = new regImage;
  regImage* movingImageItk = new regImage;
  ConvertVolumeNodeToItkImage(fixedImage, fixedImageItk);
  ConvertVolumeNodeToItkImage(movingImage, movingImageItk);

  ///If CloneVolume is used: not needed!
  //Duplicate the moving image to the deformation image
  regImage* deformationImageItk = new regImage;
  /*typedef itk::ImageDuplicator< ImageType > DuplicatorType;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(movingImageItk);
  duplicator->Update();
  deformationImageItk = duplicator->GetOutput();*/


  /////////////////////////////////////////////////////////////////////////////////
  ///                                                                           ///
  ///                   Start of actual Featurelets Code                        ///
  ///                                                                           ///
  /////////////////////////////////////////////////////////////////////////////////

  //Initializing the Result Storage
  DeformationFieldGenerator resultStorage;
  resultStorage.OpenStream();

  regImage::ImageType::SizeType imageSizeFixed = fixedImageItk->GetImagePointer()->GetLargestPossibleRegion().GetSize();
  regImage::ImageType::SizeType imageSizeMoving = movingImageItk->GetImagePointer()->GetLargestPossibleRegion().GetSize();
  PixelType iSize[Dimension];
  iSize[0] = imageSizeMoving[0];
  iSize[1] = imageSizeMoving[1];
  iSize[2] = imageSizeMoving[2];

  //Not used
  /*PixelType origin[Dimension];
  origin[0] = movingImageItk->GetOrigin()[0];
  origin[1] = movingImageItk->GetOrigin()[1];
  origin[2] = movingImageItk->GetOrigin()[2];

  ImageType::SpacingType imageSpaceMoving = movingImageItk->GetSpacing();
  PixelType space[Dimension];
  space[0] = imageSpaceMoving[0];
  space[1] = imageSpaceMoving[1];
  space[2] = imageSpaceMoving[2];*/

  regImage::ImageType::SizeType FeatureletSize, SearchRegionSize, TempSize;
  for(unsigned int j=0; j<Dimension; j++)
  {
    FeatureletSize[j] = pnode->GetFeatureletsSize();
    SearchRegionSize[j] = pnode->GetSearchRegionSize();
  }

  //Checking if the Size is the same as in de UI
  std::cerr << "FeatureletsSize: " << FeatureletSize[0] << " x " << FeatureletSize[1] << " x " << FeatureletSize[2] << std::endl;
  std::cerr << "SearchRegionSize: " << SearchRegionSize[0] << " x " << SearchRegionSize[1] << " x " << SearchRegionSize[2] << std::endl;

  bool correl = pnode->GetUseCorrelationForSimilarity();
  bool linear = pnode->GetUseLinearCorrelation();

  unsigned int numFeaturelets;
  numFeaturelets = ( imageSizeFixed[0]*imageSizeFixed[1]*imageSizeFixed[2] ) / (FeatureletSize[0]*FeatureletSize[1]*FeatureletSize[2]) *2;

  FeatureletRegistrationResultPointer *regResults;
  regResults = new FeatureletRegistrationResultPointer [numFeaturelets];

  int iFeatureletGridPosition[3];
  iFeatureletGridPosition[0]=0;
  iFeatureletGridPosition[1]=0;
  iFeatureletGridPosition[2]=0;
  int noreg=0, iCurrentFeaturelet=0, regis=0;

  for ( unsigned int i=0; i < imageSizeFixed[0]; i = i + FeatureletSize[0] )
  {
    iFeatureletGridPosition[0]++;
    iFeatureletGridPosition[1]=0;
    for ( unsigned int j=0; j < imageSizeFixed[1]; j = j + FeatureletSize[1] )
    {
      iFeatureletGridPosition[1]++;
      iFeatureletGridPosition[2]=0;
      for ( unsigned int k=0; k < imageSizeFixed[2]; k = k + FeatureletSize[2] )
      {
        iFeatureletGridPosition[2]++;
        iCurrentFeaturelet++;
        regImage::ImageType::IndexType ImageIndex;
        ImageIndex[0] = i;
        ImageIndex[1] = j;
        ImageIndex[2] = k;
        TempSize = FeatureletSize;
        if ( ImageIndex[0] + FeatureletSize[0] > imageSizeFixed[0] )
            TempSize[0] = imageSizeFixed[0] - ImageIndex[0];
        if ( ImageIndex[1] + FeatureletSize[1] > imageSizeFixed[1] )
            TempSize[1] = imageSizeFixed[1] - ImageIndex[1];
        if ( ImageIndex[2] + FeatureletSize[2] > imageSizeFixed[2] )
            TempSize[2] = imageSizeFixed[2] - ImageIndex[2];
        std::cout << "registered " << regis << " featurelets so far " << std::endl;
        std::cout << "out of (approx.): " << (imageSizeFixed[0]*imageSizeFixed[1]*imageSizeFixed[2]) / (FeatureletSize[0]*FeatureletSize[1]*FeatureletSize[2]) << " featurelets" << std::endl;
        std::cout << "Now registering featurelet no. " << iCurrentFeaturelet << " at " << ImageIndex << " " << TempSize << std::endl;

        FeatureletRegistrationResultPointer regResult;
        regResult = FeatureletRegistrationResultType::New();

        SubsampleVolume( fixedImageItk, TempSize, SearchRegionSize, ImageIndex );
        SubsampleVolume( movingImageItk, TempSize, SearchRegionSize, ImageIndex );
        int statusFixed;
        int statusMoving;
        statusMoving = CheckFeaturelet( movingImageItk );
        statusFixed = CheckFeaturelet( fixedImageItk );

        if( statusFixed==0 && statusMoving==0 )
        {
          if((!correl)&&(linear)) //Hier muss noch eine Fallunterscheidung eingefügt werden (mutual information usw.)
          {
            regResult = RegisterVolumesII2( fixedImageItk, movingImageItk );
            regis++;
          }
          if((correl)&&(!linear))
          {
            regResult = RegisterVolumes( fixedImageItk, movingImageItk );
            regis++;
          }
          if((!correl)&&(!linear))
          {
            regResult = RegisterVolumesII( fixedImageItk, movingImageItk );
            regis++;
          }
          if((correl)&&(linear))
          {
            regResult = RegisterVolumes2( fixedImageItk, movingImageItk );
            regis++;
          }
          else
          {
            noreg++;
            regResult = FeatureletRegistrationResultType::New();
          }

          regResult->SetStatusFixed(statusFixed);
          regResult->SetStatusMoving(statusMoving);
          regResult->SetFeatureletIndex(ImageIndex[0],ImageIndex[1],ImageIndex[2]);
          regResult->SetFeatureletSize(TempSize[0],TempSize[1],TempSize[2]);
          regResult->SetFeatureletGridPosition(iFeatureletGridPosition[0],iFeatureletGridPosition[1],iFeatureletGridPosition[2]);
          regResult->SetImageSizeM(iSize[0], iSize[1], iSize[2]);
          regResult->SetFeatureletNumber(iCurrentFeaturelet);

          regResults[iCurrentFeaturelet] = regResult;

          FeatureletRegistrationResultType *regPointer;
          regPointer = regResult.GetPointer();

          resultStorage.AddFeaturelet(regPointer);
          //resultStorage.SetImageSizeMov(iSize);       //Problems with Type-Conversion
          resultStorage.SetNumberOfFeaturelets(iFeatureletGridPosition);

          //resultStorage.SetSizeOfFeaturelets(FeatureletSize);         //Problems with Type-Conversion
          resultStorage.GenerateDeformationField();
          resultStorage.RescaleDeformationField();
          resultStorage.TotaldeformationFieldCreator();

          resultStorage.SaveAll();
        }
      }
    }
  }

  /// End of actual Featurelets Code
  /////////////////////////////////////////////////////////////////////////////////


  vtkSmartPointer<vtkImageData> deformationImage = vtkSmartPointer<vtkImageData>::New();
  regImage::ImageType::RegionType region = deformationImageItk->GetImagePointer()->GetBufferedRegion();
  regImage::ImageType::SizeType imageSize = region.GetSize();
  int extent[6]={0, (int) imageSize[0]-1, 0, (int) imageSize[1]-1, 0, (int) imageSize[2]-1};
  deformationImage->SetExtent(extent);
  //deformationImage->SetScalarType(VTK_FLOAT);             //Old Version
  //deformationImage->SetNumberOfScalarComponents(1);       //Changed with VTK 6
#if (VTK_MAJOR_VERSION <= 5)
  deformationImage->AllocateScalars();
#else
  deformationImage->AllocateScalars(VTK_SHORT, 1); //Allocates memory for the image pixel data
#endif

  /*PixelType* deformationPtr = (PixelType*)deformationImage->GetScalarPointer();
  itk::ImageRegionIteratorWithIndex< itk::Image<PixelType, Dimension> > itDeformationItk(
    deformationImageItk, deformationImageItk->GetImagePointer()->GetLargestPossibleRegion() );
  for ( itDeformationItk.GoToBegin(); !itDeformationItk.IsAtEnd(); ++itDeformationItk )
  {
    itk::Image<PixelType, Dimension>::IndexType i = itDeformationItk.GetIndex();
    (*deformationPtr) = deformationImageItk->GetImagePointer()->GetPixel(i);
    //std::cerr << deformationImageItk->GetPixel(i) << std::endl;
    deformationPtr++;
  }*/

  deformedImage->CopyOrientation(fixedImage);
  deformedImage->SetAndObserveImageData(deformationImage);

  imageData = deformedImage->GetImageData();
  std::cerr << "Deformed - Pixeltype: " << imageData->GetScalarType() << " " << imageData->GetScalarTypeAsString() << std::endl;
  std::cerr << std::endl;

  return 0;
}

//-----------------------------------------------------------------------------
int vtkSlicerRegistrationLogic::ShowVolume(vtkMRMLRegistrationNode* pnode, bool fixedImage)
{
    std::cerr << "*** Entered ShowVolume ***" << std::endl;

  vtkMRMLScene *scene = this->GetMRMLScene();
  //vtkMRMLVolumeNode *deformedImage = NULL;
  vtkMRMLVolumeNode *fixedImageNode = vtkMRMLVolumeNode::SafeDownCast(scene->GetNodeByID(pnode->GetFixedImageNodeID()));
  vtkMRMLVolumeNode *movingImageNode = vtkMRMLVolumeNode::SafeDownCast(scene->GetNodeByID(pnode->GetMovingImageNodeID()));
/*
  if(fixedImage)
      deformedImage->SetAndObserveTransformNodeID(fixedImageNode->GetParentTransformNode()->GetID());
  else
      deformedImage->SetAndObserveTransformNodeID(movingImageNode->GetParentTransformNode()->GetID());
*/
  if(!fixedImageNode && !movingImageNode)
    {
    std::cerr << "Failed to look up input volume!" << std::endl;
    return -1;
    }
  if(fixedImage)
  {
      std::cerr << "Show Fixed Image Node" << std::endl;
      pnode->SetDeformedImageNodeID(fixedImageNode->GetID());
  }
  else
  {
      std::cerr << "Show Moving Image Node" << std::endl;
      pnode->SetDeformedImageNodeID(movingImageNode->GetID());
  }

  //pnode->SetDeformedImageNodeID(deformedImage->GetID());
  return 0;
}

//-----------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::RegisterNodes()
{
    if(!this->GetMRMLScene())
      {
      return;
      }
    vtkMRMLRegistrationNode* pNode = vtkMRMLRegistrationNode::New();
    this->GetMRMLScene()->RegisterNodeClass(pNode);
    pNode->Delete();
}

/*//-----------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::SnapROIToVoxelGrid(vtkMRMLVolumeNode* inputVolume)
{
    if (!inputVolume)
      {
        vtkDebugMacro(
            "vtkSlicerRegistrationLogic: Inversion failed (input Volume is invalid)");
        return;
      }

    vtkSmartPointer<vtkMRMLScene> scene = this->GetMRMLScene();

    if (!scene)
      {
        vtkDebugMacro(
            "vtkSlicerRegistrationLogic: Inversion failed (MRML Scene is not available)");
        return;
      }

    vtkNew<vtkMatrix4x4> rotationMat;

    bool ok = vtkSlicerRegistrationLogic::ComputeIJKToRASRotationOnlyMatrix(
        inputVolume, rotationMat.GetPointer());

    if (!ok)
      {
        vtkDebugMacro(
            "vtkSlicerRegistrationLogic: Inversion failed (IJK to RAS rotation only matrix computation failed)");
        return;
      }

    vtkNew<vtkMRMLLinearTransformNode> InversionTransformNode;
    InversionTransformNode->ApplyTransformMatrix(rotationMat.GetPointer());
    InversionTransformNode->SetScene(this->GetMRMLScene());

    this->GetMRMLScene()->AddNode(InversionTransformNode.GetPointer());

    inputVolume->SetAndObserveTransformNodeID(InversionTransformNode->GetID());
}

//----------------------------------------------------------------------------
bool vtkSlicerRegistrationLogic::ComputeIJKToRASRotationOnlyMatrix(vtkMRMLVolumeNode* inputVolume, vtkMatrix4x4* outputMatrix)
{
  if(inputVolume == NULL || outputMatrix == NULL)
    return false;

  vtkNew<vtkMatrix4x4> rotationMat;
  rotationMat->Identity();

  vtkNew<vtkMatrix4x4> inputIJKToRASMat;
  inputIJKToRASMat->Identity();

  inputVolume->GetIJKToRASMatrix(inputIJKToRASMat.GetPointer());
  const char* scanOrder = vtkMRMLVolumeNode::ComputeScanOrderFromIJKToRAS(inputIJKToRASMat.GetPointer());

  vtkNew<vtkMatrix4x4> orientMat;
  orientMat->Identity();

  bool orientation = vtkSlicerRegistrationLogic::ComputeOrientationMatrixFromScanOrder(scanOrder,orientMat.GetPointer());

  if(!orientation)
    return false;

  orientMat->Invert();

  vtkNew<vtkMatrix4x4> directionsMat;
  directionsMat->Identity();

  inputVolume->GetIJKToRASDirectionMatrix(directionsMat.GetPointer());

  vtkMatrix4x4::Multiply4x4(directionsMat.GetPointer(),orientMat.GetPointer(),rotationMat.GetPointer());

  outputMatrix->DeepCopy(rotationMat.GetPointer());

  return true;
}

//----------------------------------------------------------------------------
bool vtkSlicerRegistrationLogic::IsVolumeTiltedInRAS( vtkMRMLVolumeNode* inputVolume, vtkMatrix4x4* rotationMatrix)
{
  assert(inputVolume);

  vtkNew<vtkMatrix4x4> iJKToRASMat;
  vtkNew<vtkMatrix4x4> directionMat;

  inputVolume->GetIJKToRASMatrix(iJKToRASMat.GetPointer());
  inputVolume->GetIJKToRASDirectionMatrix(directionMat.GetPointer());

  const char* scanOrder = vtkMRMLVolumeNode::ComputeScanOrderFromIJKToRAS(iJKToRASMat.GetPointer());

  vtkNew<vtkMatrix4x4> orientMat;
  vtkSlicerRegistrationLogic::ComputeOrientationMatrixFromScanOrder(scanOrder,orientMat.GetPointer());

  bool same = true;

  for(int i=0; i < 4; ++i)
    {
      for(int j=0; j < 4; ++j)
        {
          if(directionMat->GetElement(i,j) != orientMat->GetElement(i,j))
            {
              same = false;
              break;
            }
        }
      if(same == false)
        break;
    }

  if(!rotationMatrix)
    rotationMatrix = vtkSmartPointer<vtkMatrix4x4>::New();

  if(same)
    {
      rotationMatrix->Identity();
    }
  else
    {
      vtkSlicerRegistrationLogic::ComputeIJKToRASRotationOnlyMatrix(inputVolume,rotationMatrix);
      return true;
    }

  return false;
}

//----------------------------------------------------------------------------
bool
vtkSlicerRegistrationLogic::ComputeOrientationMatrixFromScanOrder(
    const char *order, vtkMatrix4x4 *outputMatrix)
{
  vtkNew<vtkMatrix4x4> orientMat;
  orientMat->Identity();

  if (!strcmp(order, "IS") || !strcmp(order, "Axial IS")
      || !strcmp(order, "Axial"))
    {
      double elems[] =
        { -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };
      orientMat->DeepCopy(elems);
    }
  else if (!strcmp(order, "SI") || !strcmp(order, "Axial SI"))
    {
      double elems[] =
        { -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1 };
      orientMat->DeepCopy(elems);
    }
  else if (!strcmp(order, "RL") || !strcmp(order, "Sagittal RL")
      || !strcmp(order, "Sagittal"))
    {
      double elems[] =
        { 0, 0, -1, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1 };
      orientMat->DeepCopy(elems);
    }
  else if (!strcmp(order, "LR") || !strcmp(order, "Sagittal LR"))
    {
      double elems[] =
        { 0, 0, 1, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1 };
      orientMat->DeepCopy(elems);
    }
  else if (!strcmp(order, "PA") || !strcmp(order, "Coronal PA")
      || !strcmp(order, "Coronal"))
    {
      double elems[] =
        { -1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 1 };
      orientMat->DeepCopy(elems);
    }
  else if (!strcmp(order, "AP") || !strcmp(order, "Coronal AP"))
    {
      double elems[] =
        { -1, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 1 };
      orientMat->DeepCopy(elems);
    }
  else
    {
      return false;
    }


  outputMatrix->DeepCopy(orientMat.GetPointer());
  return true;
}
*/

//---------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);

  this->Modified();
}

//---------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::OnMRMLSceneNodeAdded(vtkMRMLNode* node)
{
    if(!node || !this->GetMRMLScene())
        return;

    if(node->IsA("vtkMRMLVolumeNode"))
        this->Modified();
}

//---------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::OnMRMLSceneNodeRemoved(vtkMRMLNode* node)
{
    if(!node || !this->GetMRMLScene())
        return;

    if(node->IsA("vtkMRMLVolumeNode"))
        this->Modified();
}

//---------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::OnMRMLSceneEndImport()
{
  // If we have a parameter node select it
  vtkMRMLRegistrationNode *paramNode = NULL;
  vtkMRMLNode *node = this->GetMRMLScene()->GetNthNodeByClass(0, "vtkMRMLRegistrationNode");
  if (node)
  {
    paramNode = vtkMRMLRegistrationNode::SafeDownCast(node);
    vtkSetAndObserveMRMLNodeMacro(this->RegistrationNode, paramNode);
  }

  this->Modified();
}

//---------------------------------------------------------------------------
void vtkSlicerRegistrationLogic::OnMRMLSceneEndClose()
{
  this->Modified();
}

// From Slicer-RT-commons
//---------------------------------------------------------------------------
bool vtkSlicerRegistrationLogic::ConvertVolumeNodeToItkImage(vtkMRMLVolumeNode* inVolumeNode, regImage* outItkVolume)
{
  if ( inVolumeNode == NULL )
  {
    std::cerr << "Failed to convert volume node to itk image - input MRML volume node is NULL!" << std::endl;
    return false;
  }

  vtkImageData* inVolume = inVolumeNode->GetImageData();
  if ( inVolume == NULL )
  {
    std::cerr << "Failed to convert volume node to itk image - image in input MRML volume node is NULL!" << std::endl;
    return false;
  }

  if ( outItkVolume->GetImagePointer().IsNull() )
  {
    std::cerr << "Failed to convert volume node to itk image - output image is NULL!" << std::endl;
    return false;
  }

  // Convert vtkImageData to itkImage
  vtkSmartPointer<vtkImageExport> imageExport = vtkSmartPointer<vtkImageExport>::New();
#if (VTK_MAJOR_VERSION <= 5)
  imageExport->SetInput(inVolume);
#else
  imageExport->SetInputData(inVolume); // Was ->SetInput(inVolume) before VTK 6
#endif
  imageExport->Update();

  // Determine input volume to world transform
  vtkSmartPointer<vtkMatrix4x4> rasToWorldTransformMatrix=vtkSmartPointer<vtkMatrix4x4>::New();
  vtkMRMLTransformNode* inTransformNode=inVolumeNode->GetParentTransformNode();
  if (inTransformNode!=NULL)
  {
    if (inTransformNode->IsTransformToWorldLinear() == 0)
    {
      std::cerr << "There is a non-linear transform assigned to an input dose volume. Only linear transforms are supported!" << std::endl;
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
  double outputSpacing[Dimension] = {0.0, 0.0, 0.0};
  inVolumeToWorldTransform->GetScale(outputSpacing);
  outItkVolume->GetImagePointer()->SetSpacing(outputSpacing);

  double outputOrigin[Dimension] = {0.0, 0.0, 0.0};
  inVolumeToWorldTransform->GetPosition(outputOrigin);
  outItkVolume->GetImagePointer()->SetOrigin(outputOrigin);

  double outputOrienationAngles[Dimension] = {0.0, 0.0, 0.0};
  inVolumeToWorldTransform->GetOrientation(outputOrienationAngles);
  vtkSmartPointer<vtkTransform> inVolumeToWorldOrientationTransform = vtkSmartPointer<vtkTransform>::New();
  inVolumeToWorldOrientationTransform->Identity();
  inVolumeToWorldOrientationTransform->RotateX(outputOrienationAngles[0]);
  inVolumeToWorldOrientationTransform->RotateY(outputOrienationAngles[1]);
  inVolumeToWorldOrientationTransform->RotateZ(outputOrienationAngles[2]);
  vtkSmartPointer<vtkMatrix4x4> inVolumeToWorldOrientationTransformMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
  inVolumeToWorldOrientationTransform->GetMatrix(inVolumeToWorldOrientationTransformMatrix);
  itk::Matrix<double,3,3> outputDirectionMatrix;
  for(int i=0; i<3; i++)
  {
    for(int j=0; j<3; j++)
    {
      outputDirectionMatrix[i][j] = inVolumeToWorldOrientationTransformMatrix->GetElement(i,j);
    }
  }
  outItkVolume->GetImagePointer()->SetDirection(outputDirectionMatrix);

  int inputExtent[6]={0,0,0,0,0,0};
  inVolume->GetExtent(inputExtent);
  itk::Image<PixelType, Dimension>::SizeType inputSize;
  inputSize[0] = inputExtent[1] - inputExtent[0] + 1;
  inputSize[1] = inputExtent[3] - inputExtent[2] + 1;
  inputSize[2] = inputExtent[5] - inputExtent[4] + 1;

  itk::Image<PixelType, Dimension>::IndexType start;
  start[0]=start[1]=start[2]=0.0;

  itk::Image<PixelType, Dimension>::RegionType region;
  region.SetSize(inputSize);
  region.SetIndex(start);
  outItkVolume->GetImagePointer()->SetRegions(region);

  try
  {
    outItkVolume->GetImagePointer()->Allocate();
  }
  catch(itk::ExceptionObject & err)
  {
    std::cerr << "Failed to allocate memory for the image conversion: " << err.GetDescription() << std::endl;
    return false;
  }

  imageExport->Export( outItkVolume->GetImagePointer()->GetBufferPointer() );

  return true;
}


///	This function performs subsampling of a volume using ITK classes. In order
///	to avoid aliasing artifacts, the volume must be processed by a low-pass
///	filter before resampling.
int vtkSlicerRegistrationLogic::SubsampleVolume( regImage* Image,
                     regImage::ImageType::SizeType FeatureletSize,
                     regImage::ImageType::SizeType SearchRegionSize,
                     regImage::ImageType::IndexType ImageIndex)
{
  typedef regImage::ImageType ImageType;
  regImage::ImageType::RegionType desiredRegion;
        desiredRegion.SetSize( FeatureletSize );
        desiredRegion.SetIndex( ImageIndex );
  regImage::ImageType::RegionType SearchingRegion;
        SearchingRegion.SetSize( SearchRegionSize );
        SearchingRegion.SetIndex( ImageIndex );

  typedef itk::ExtractImageFilter< ImageType, ImageType > FilterType;
  FilterType::Pointer filter = FilterType::New();
  FilterType::Pointer filter2 = FilterType::New();
        filter->SetExtractionRegion( desiredRegion );
        filter->SetInput(Image->GetImagePointer());
        filter2->SetExtractionRegion(SearchingRegion);
        filter2->SetInput(Image->GetImagePointer());
  try
  {
        filter->Update();
        filter2->Update();
  }
  catch ( itk::ExceptionObject& e )
  {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
  }
  Image->SetFeatureletPointer( filter->GetOutput() );
  Image->SetSearchRegionPointer(filter2->GetOutput());
  return EXIT_SUCCESS;
}


///	This function checks a featurelet if it is necessary to be registered to
///	avoid registration errors or a bad result.
int vtkSlicerRegistrationLogic::CheckFeaturelet( regImage* Image )
{
  typedef regImage::ImageType ImageType;
  ImageType::SizeType FeatureletSize = Image->GetFeatureletPointer()->GetLargestPossibleRegion().GetSize();

  typedef itk::NormalizedCorrelationImageToImageMetric<ImageType, ImageType> MetricType;
    //MetricType::Pointer metric = MetricType::New();
    if ( FeatureletSize[0] <= 4 )
        return 1;
    if ( FeatureletSize[1] <= 4 )
        return  1;
    if ( FeatureletSize[2] <= 4 )
        return 1;

    typedef itk::MinimumMaximumImageCalculator< ImageType > CalculatorType;
    CalculatorType::Pointer calculator = CalculatorType::New();
        calculator->SetImage( Image->GetFeatureletPointer() );
    try {
        calculator->Compute();
    }
    catch ( itk::ExceptionObject& e ) {
        std::cerr << e << std::endl;
        return 3;
        }
    if ( calculator->GetMaximum() == calculator->GetMinimum() ){
        return 2;
    }
    return 0;
}
