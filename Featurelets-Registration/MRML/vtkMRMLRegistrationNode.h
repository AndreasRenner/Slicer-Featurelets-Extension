/*=auto=========================================================================

  Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: vtkMRMLVolumeRenderingParametersNode.h,v $
  Date:      $Date: 2006/03/19 17:12:29 $
  Version:   $Revision: 1.3 $

=========================================================================auto=*/
// .NAME vtkMRMLVolumeRenderingParametersNode - MRML node for storing a slice through RAS space
// .SECTION Description
// This node stores the information about the currently selected volume
//
//

#ifndef __vtkMRMLRegistrationNode_h
#define __vtkMRMLRegistrationNode_h

#include "vtkMRML.h"
#include "vtkMRMLScene.h"
#include "vtkMRMLNode.h"
#include "vtkSlicerRegistrationModuleMRMLExport.h"

class vtkMRMLAnnotationROINode;
class vtkMRMLVolumeNode;

/// \ingroup Slicer_QtModules_Registration
class VTK_SLICER_REGISTRATION_MODULE_MRML_EXPORT vtkMRMLRegistrationNode : public vtkMRMLNode
{
  public:

  static vtkMRMLRegistrationNode *New();
  vtkTypeMacro(vtkMRMLRegistrationNode,vtkMRMLNode);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual vtkMRMLNode* CreateNodeInstance();

  // Set node attributes
  virtual void ReadXMLAttributes( const char** atts);

  // Write this node's information to a MRML file in XML format.
  virtual void WriteXML(ostream& of, int indent);

  // Copy the node's attributes to this object
  virtual void Copy(vtkMRMLNode *node);

  // Get node XML tag name (like Volume, Model)
  virtual const char* GetNodeTagName() {return "Registration";};

  void SetAndObserveFixedImageNodeID(const char* ID);
  void SetAndObserveMovingImageNodeID(const char* ID);
  void SetAndObserveDeformedImageNodeID(const char* ID);
  void SetAndObserveDeformationFieldID(const char* ID);
  void SetAndObserveFiducialPointsID(const char* ID);
  void SetAndObserveROINodeID(const char* ID);          //Not yet implemented

  // Description:
  vtkSetStringMacro(FixedImageNodeID);
  vtkGetStringMacro (FixedImageNodeID);
  vtkSetStringMacro(MovingImageNodeID);
  vtkGetStringMacro (MovingImageNodeID);
  vtkSetStringMacro(DeformedImageNodeID);
  vtkGetStringMacro (DeformedImageNodeID);
  vtkSetStringMacro(DeformationFieldID);
  vtkGetStringMacro (DeformationFieldID);
  vtkSetStringMacro(FiducialPointsID);
  vtkGetStringMacro (FiducialPointsID);

  //Not yet implemented
  vtkSetStringMacro(ROINodeID);
  vtkGetStringMacro (ROINodeID);
  vtkSetMacro(ROIcheckBox,bool);
  vtkGetMacro(ROIcheckBox,bool);
  vtkBooleanMacro(ROIcheckBox,bool);

  vtkGetMacro(FeatureletsSize, int);
  vtkSetMacro(FeatureletsSize, int);
  vtkGetMacro(SearchRegionSize, int);
  vtkSetMacro(SearchRegionSize, int);
  vtkGetMacro(FeatureletsSizeZ, int);
  vtkSetMacro(FeatureletsSizeZ, int);
  vtkGetMacro(SearchRegionSizeZ, int);
  vtkSetMacro(SearchRegionSizeZ, int);


  vtkGetMacro(MaxStepLength, double);
  vtkSetMacro(MaxStepLength, double);
  vtkGetMacro(MinStepLength, double);
  vtkSetMacro(MinStepLength, double);
  vtkGetMacro(NumberIterations, int);
  vtkSetMacro(NumberIterations, int);
  vtkGetMacro(progress, int);
  vtkSetMacro(progress, int);

   // Bool Similarity
  vtkGetMacro(UseCorrelationForSimilarity, bool);
  vtkSetMacro(UseCorrelationForSimilarity, bool);
  vtkBooleanMacro(UseCorrelationForSimilarity, bool);
   // Bool Correlation
  vtkGetMacro(UseLinearCorrelation, bool);
  vtkSetMacro(UseLinearCorrelation, bool);
  vtkBooleanMacro(UseLinearCorrelation, bool);
   // Bool Fiducial
  vtkGetMacro(checkBoxFiducial, bool);
  vtkSetMacro(checkBoxFiducial, bool);
  vtkBooleanMacro(checkBoxFiducial, bool);
   // Bool Debug
  vtkGetMacro(checkBoxDebug, bool);
  vtkSetMacro(checkBoxDebug, bool);
  vtkBooleanMacro(checkBoxDebug, bool);
   // Bool Rigid
  vtkGetMacro(checkBoxRigid, bool);
  vtkSetMacro(checkBoxRigid, bool);
  vtkBooleanMacro(checkBoxRigid, bool);
   // Bool Z Different
  vtkGetMacro(checkBoxZDifferent, bool);
  vtkSetMacro(checkBoxZDifferent, bool);
  vtkBooleanMacro(checkBoxZDifferent, bool);

  // Update the stored reference to another node in the scene
  virtual void UpdateReferenceID(const char *oldID, const char *newID);

protected:
  vtkMRMLRegistrationNode();
  ~vtkMRMLRegistrationNode();

  vtkMRMLRegistrationNode(const vtkMRMLRegistrationNode&);
  void operator=(const vtkMRMLRegistrationNode&);

  char *FixedImageNodeID;
  char *MovingImageNodeID;
  char *DeformedImageNodeID;
  char *DeformationFieldID;
  char *FiducialPointsID;

  //Not yet implemented
  char *ROINodeID;
  bool ROIcheckBox;

  int FeatureletsSize;
  int SearchRegionSize;
  int FeatureletsSizeZ;
  int SearchRegionSizeZ;
  double MaxStepLength;
  double MinStepLength;
  int NumberIterations;
  int progress;
  bool UseCorrelationForSimilarity;
  bool UseLinearCorrelation;
  bool checkBoxFiducial;
  bool checkBoxDebug;
  bool checkBoxRigid;
  bool checkBoxZDifferent;
};

#endif

