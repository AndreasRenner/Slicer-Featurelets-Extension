/*=auto=========================================================================

Portions (c) Copyright 2005 Brigham and Women\"s Hospital (BWH) All Rights Reserved.

See COPYRIGHT.txt
or http://www.slicer.org/copyright/copyright.txt for details.

Program:   3D Slicer
Module:    $RCSfile: vtkMRMLRegistrationNode.cxx,v $
Date:      $Date: 2006/03/17 15:10:10 $
Version:   $Revision: 1.2 $

=========================================================================auto=*/

// VTK includes
#include <vtkCommand.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>

// MRML includes
#include "vtkMRMLVolumeNode.h"

// CropModuleMRML includes
#include "vtkMRMLRegistrationNode.h"

// AnnotationModuleMRML includes
#include "vtkMRMLAnnotationROINode.h"

// STD includes
#include <sstream>

//----------------------------------------------------------------------------
vtkMRMLNodeNewMacro(vtkMRMLRegistrationNode);

//----------------------------------------------------------------------------
vtkMRMLRegistrationNode::vtkMRMLRegistrationNode() {
  this->HideFromEditors = 1;

  this->FixedImageNodeID = NULL;
  this->MovingImageNodeID = NULL;
  this->DeformedImageNodeID = NULL;
  this->DeformationFieldID = NULL;
  this->FiducialPointsID = NULL;

  this->ROINodeID = NULL;       //Not yet implemented
  this->ROIcheckBox = false;    //Not yet implemented

  this->FeatureletsSize = 15;
  this->SearchRegionSize = 30;
  this->FeatureletsSizeZ = 15;
  this->SearchRegionSizeZ = 30;
  this->MaxStepLength = 0.5;
  this->MinStepLength = 0.001;
  this->NumberIterations = 2000;
  this->progress = 0;
  this->UseCorrelationForSimilarity = true;
  this->UseLinearCorrelation = true;
  this->checkBoxFiducial = false;
  this->checkBoxDebug = false;
  this->checkBoxRigid = false;
  this->checkBoxZDifferent = false;
}

//----------------------------------------------------------------------------
vtkMRMLRegistrationNode::~vtkMRMLRegistrationNode() {
  if (this->FixedImageNodeID) {
    this->SetFixedImageNodeID(NULL);
  }
  if (this->MovingImageNodeID) {
    this->SetMovingImageNodeID(NULL);
  }
  if (this->DeformedImageNodeID) {
    this->SetDeformedImageNodeID(NULL);
  }
  if (this->DeformationFieldID) {
    this->SetDeformationFieldID(NULL);
  }
  if (this->FiducialPointsID) {
    this->SetFiducialPointsID(NULL);
  }
  if (this->ROINodeID) {         //Not yet implemented
    this->SetROINodeID(NULL);
  }
}

//----------------------------------------------------------------------------
void vtkMRMLRegistrationNode::ReadXMLAttributes(const char** atts) {
  std::cerr << "Reading Registration param node!" << std::endl;
  Superclass::ReadXMLAttributes(atts);

  const char* attName;
  const char* attValue;

  while (*atts != NULL)   {
    attName = *(atts++);
    attValue = *(atts++);

    if (!strcmp(attName, "fixedImageNodeID")) {
      std::stringstream ss;
      ss << attValue;
      this->SetAndObserveFixedImageNodeID(ss.str().c_str());
      continue;
    }
    if (!strcmp(attName, "movingImageNodeID")) {
      std::stringstream ss;
      ss << attValue;
      this->SetAndObserveMovingImageNodeID(ss.str().c_str());
      continue;
    }
    if (!strcmp(attName, "deformedImageNodeID")) {
      std::stringstream ss;
      ss << attValue;
      this->SetAndObserveDeformedImageNodeID(ss.str().c_str());
      continue;
    }
    if (!strcmp(attName, "deformationFieldID")) {
      std::stringstream ss;
      ss << attValue;
      this->SetAndObserveDeformationFieldID(ss.str().c_str());
      continue;
    }
    if (!strcmp(attName, "fiducialPointsID")) {
      std::stringstream ss;
      ss << attValue;
      this->SetAndObserveFiducialPointsID(ss.str().c_str());
      continue;
    }
    if (!strcmp(attName, "ROINodeID")) {      //Not yet implemented
      std::stringstream ss;
      ss << attValue;
      this->SetAndObserveROINodeID(ss.str().c_str());
      continue;
    }
    if (!strcmp(attName, "FeatureletsSize")) {
      std::stringstream ss;
      ss << attValue;
      int intAttValue;
      ss >> intAttValue;
      this->FeatureletsSize = intAttValue;
      continue;
    }
    if (!strcmp(attName, "SearchRegionSize")) {
      std::stringstream ss;
      ss << attValue;
      int intAttValue;
      ss >> intAttValue;
      this->SearchRegionSize = intAttValue;
      continue;
    }
    if (!strcmp(attName, "FeatureletsSizeZ")) {
      std::stringstream ss;
      ss << attValue;
      int intAttValue;
      ss >> intAttValue;
      this->FeatureletsSizeZ = intAttValue;
      continue;
    }
    if (!strcmp(attName, "SearchRegionSizeZ")) {
      std::stringstream ss;
      ss << attValue;
      int intAttValue;
      ss >> intAttValue;
      this->SearchRegionSizeZ = intAttValue;
      continue;
    }
    if (!strcmp(attName, "MaxStepLength")) {
      std::stringstream ss;
      ss << attValue;
      double intAttValue;
      ss >> intAttValue;
      this->MaxStepLength = intAttValue;
      continue;
    }
    if (!strcmp(attName, "MinStepLength")) {
      std::stringstream ss;
      ss << attValue;
      double intAttValue;
      ss >> intAttValue;
      this->MinStepLength = intAttValue;
      continue;
    }
    if (!strcmp(attName, "NumberIterations")) {
      std::stringstream ss;
      ss << attValue;
      int intAttValue;
      ss >> intAttValue;
      this->NumberIterations = intAttValue;
      continue;
    }
    if (!strcmp(attName, "progress")) {
      std::stringstream ss;
      ss << attValue;
      int intAttValue;
      ss >> intAttValue;
      this->progress = intAttValue;
      continue;
    }
    if (!strcmp(attName, "UseCorrelationForSimilarity")) {
      this->UseCorrelationForSimilarity = (strcmp(attValue,"true") ? false : true);
      continue;
    }
    if (!strcmp(attName, "UseLinearCorrelation")) {
      this->UseLinearCorrelation = (strcmp(attValue,"true") ? false : true);
      continue;
    }
    if (!strcmp(attName, "checkBoxFiducial")) {
      this->checkBoxFiducial = (strcmp(attValue,"false") ? true : false);
      continue;
    }
    if (!strcmp(attName, "checkBoxDebug")) {
      this->checkBoxDebug = (strcmp(attValue,"false") ? true : false);
      continue;
    }
    if (!strcmp(attName, "checkBoxRigid")) {
      this->checkBoxRigid = (strcmp(attValue,"false") ? true : false);
      continue;
    }
    if (!strcmp(attName, "checkBoxZDifferent")) {
      this->checkBoxZDifferent = (strcmp(attValue,"false") ? true : false);
      continue;
    }
  }
  this->WriteXML(std::cout,1);
}

//----------------------------------------------------------------------------
void vtkMRMLRegistrationNode::WriteXML(ostream& of, int nIndent) {
  Superclass::WriteXML(of, nIndent);

  vtkIndent indent(nIndent);

  of << indent << " fixedImageNodeID=\"" << (this->FixedImageNodeID ? this->FixedImageNodeID : "NULL") << "\"";
  of << indent << " movingImageNodeID=\"" << (this->MovingImageNodeID ? this->MovingImageNodeID : "NULL") << "\"";
  of << indent << " deformedImageNodeID=\"" << (this->DeformedImageNodeID ? this->DeformedImageNodeID : "NULL") << "\"";
  of << indent << " deformationFieldID=\"" << (this->DeformationFieldID ? this->DeformationFieldID : "NULL") << "\"";
  of << indent << " fiducialPointsID=\"" << (this->FiducialPointsID ? this->FiducialPointsID : "NULL") << "\"";
  of << indent << " ROINodeID=\"" << (this->ROINodeID ? this->ROINodeID : "NULL") << "\"";        //Not yet implemented
  of << indent << " FeatureletsSize=\"" << this->FeatureletsSize << "\"";
  of << indent << " SearchRegionSize=\"" << this->SearchRegionSize << "\"";
  of << indent << " FeatureletsSizeZ=\"" << this->FeatureletsSizeZ << "\"";
  of << indent << " SearchRegionSizeZ=\"" << this->SearchRegionSizeZ << "\"";
  of << indent << " MaxStepLength=\"" << this->MaxStepLength << "\"";
  of << indent << " MinStepLength=\"" << this->MinStepLength << "\"";
  of << indent << " NumberIterations=\"" << this->NumberIterations << "\"";
  of << indent << " progress=\"" << this->progress << "\"";
  of << indent << " UseCorrelationForSimilarity=\"" << (this->UseCorrelationForSimilarity ? "true" : "false") << "\"";
  of << indent << " UseLinearCorrelation=\"" << (this->UseLinearCorrelation ? "true" : "false") << "\"";
  of << indent << " checkBoxFiducial=\"" << (this->checkBoxFiducial ? "false" : "true") << "\"";
  of << indent << " checkBoxDebug=\"" << (this->checkBoxDebug ? "false" : "true") << "\"";
  of << indent << " checkBoxRigid=\"" << (this->checkBoxRigid ? "false" : "true") << "\"";
  of << indent << " checkBoxZDifferent=\"" << (this->checkBoxZDifferent ? "false" : "true") << "\"";
}

//----------------------------------------------------------------------------
// Copy the node\"s attributes to this object.
// Does NOT copy: ID, FilePrefix, Name, SliceID
void vtkMRMLRegistrationNode::Copy(vtkMRMLNode *anode) {
  Superclass::Copy(anode);
  this->DisableModifiedEventOn();

  vtkMRMLRegistrationNode *node = vtkMRMLRegistrationNode::SafeDownCast(anode);

  this->SetAndObserveFixedImageNodeID(node->FixedImageNodeID);
  this->SetAndObserveMovingImageNodeID(node->MovingImageNodeID);
  this->SetAndObserveDeformedImageNodeID(node->DeformedImageNodeID);
  this->SetAndObserveDeformationFieldID(node->DeformationFieldID);
  this->SetAndObserveFiducialPointsID(node->FiducialPointsID);
  this->SetAndObserveROINodeID(node->ROINodeID);         //Not yet implemented

  this->FeatureletsSize = node->FeatureletsSize;
  this->SearchRegionSize = node->SearchRegionSize;
  this->FeatureletsSizeZ = node->FeatureletsSizeZ;
  this->SearchRegionSizeZ = node->SearchRegionSizeZ;
  this->MaxStepLength = node->MaxStepLength;
  this->MinStepLength = node->MinStepLength;
  this->NumberIterations = node->NumberIterations;
  this->progress = node->progress;
  this->UseCorrelationForSimilarity = node->UseCorrelationForSimilarity;
  this->UseLinearCorrelation = node->UseLinearCorrelation;
  this->checkBoxFiducial = node->checkBoxFiducial;
  this->checkBoxDebug = node->checkBoxDebug;
  this->checkBoxRigid = node->checkBoxRigid;
  this->checkBoxZDifferent = node->checkBoxZDifferent;

  this->DisableModifiedEventOff();
  this->InvokePendingModifiedEvent();
}

//----------------------------------------------------------------------------
void vtkMRMLRegistrationNode::PrintSelf(ostream& os, vtkIndent indent) {
  Superclass::PrintSelf(os,indent);

  os << "FixedImageNodeID: " << ( (this->FixedImageNodeID) ? this->FixedImageNodeID : "None" ) << "\n";
  os << "MovingImageNodeID: " << ( (this->MovingImageNodeID) ? this->MovingImageNodeID : "None" ) << "\n";
  os << "DeformedImageNodeID: " << ( (this->DeformedImageNodeID) ? this->DeformedImageNodeID : "None" ) << "\n";
  os << "DeformationFieldID: " << ( (this->DeformationFieldID) ? this->DeformationFieldID : "None" ) << "\n";
  os << "FiducialPointsID: " << ( (this->FiducialPointsID) ? this->FiducialPointsID : "None" ) << "\n";
  os << "ROINodeID: " << ( (this->ROINodeID) ? this->ROINodeID : "None" ) << "\n";    //Not yet implemented

  os << indent << "FeatureletsSize:   " << this->FeatureletsSize << "\n";
  os << indent << "SearchRegionSize:   " << this->SearchRegionSize << "\n";
  os << indent << "FeatureletsSizeZ:   " << this->FeatureletsSizeZ << "\n";
  os << indent << "SearchRegionSizeZ:   " << this->SearchRegionSizeZ << "\n";
  os << indent << "MaxStepLength:   " << this->MaxStepLength << "\n";
  os << indent << "MinStepLength:   " << this->MinStepLength << "\n";
  os << indent << "NumberIterations:   " << this->NumberIterations << "\n";
  os << indent << "Progress:   " << this->progress << "\n";
  os << indent << "UseCorrelationForSimilarity:   " << (this->UseCorrelationForSimilarity ? "true" : "false") << "\n";
  os << indent << "UseLinearCorrelation:   " << (this->UseLinearCorrelation ? "true" : "false") << "\n";
  os << indent << "checkBoxFiducial:   " << (this->checkBoxFiducial ? "false" : "true") << "\n";
  os << indent << "checkBoxDebug:   " << (this->checkBoxDebug ? "false" : "true") << "\n";
  os << indent << "checkBoxRigid:   " << (this->checkBoxRigid ? "false" : "true") << "\n";
  os << indent << "checkBoxZDifferent:   " << (this->checkBoxZDifferent ? "false" : "true") << "\n";
}

//----------------------------------------------------------------------------
void vtkMRMLRegistrationNode::SetAndObserveFixedImageNodeID(const char* ID) {
  if (this->FixedImageNodeID) {
    this->Scene->RemoveReferencedNodeID(this->FixedImageNodeID, this);
  }
  this->SetFixedImageNodeID(ID);

  if (ID) {
    this->Scene->AddReferencedNodeID(this->FixedImageNodeID, this);
  }
}

//----------------------------------------------------------------------------
void vtkMRMLRegistrationNode::SetAndObserveMovingImageNodeID(const char* ID) {
  if (this->MovingImageNodeID) {
    this->Scene->RemoveReferencedNodeID(this->MovingImageNodeID, this);
  }
  this->SetMovingImageNodeID(ID);

  if (ID) {
    this->Scene->AddReferencedNodeID(this->MovingImageNodeID, this);
  }
}

//----------------------------------------------------------------------------
void vtkMRMLRegistrationNode::SetAndObserveDeformedImageNodeID(const char* ID) {
  if (this->DeformedImageNodeID) {
    this->Scene->RemoveReferencedNodeID(this->DeformedImageNodeID, this);
  }
  this->SetDeformedImageNodeID(ID);

  if (ID) {
    this->Scene->AddReferencedNodeID(this->DeformedImageNodeID, this);
  }
}

//----------------------------------------------------------------------------
void vtkMRMLRegistrationNode::SetAndObserveDeformationFieldID(const char* ID) {
  if (this->DeformationFieldID) {
    this->Scene->RemoveReferencedNodeID(this->DeformationFieldID, this);
  }
  this->SetDeformationFieldID(ID);

  if (ID) {
    this->Scene->AddReferencedNodeID(this->DeformationFieldID, this);
  }
}

//----------------------------------------------------------------------------
void vtkMRMLRegistrationNode::SetAndObserveFiducialPointsID(const char* ID) {
  if (this->FiducialPointsID) {
    this->Scene->RemoveReferencedNodeID(this->FiducialPointsID, this);
  }
  this->SetFiducialPointsID(ID);

  if (ID) {
    this->Scene->AddReferencedNodeID(this->FiducialPointsID, this);
  }
}

//----------------------------------------------------------------------------
void vtkMRMLRegistrationNode::SetAndObserveROINodeID(const char* ID) {
  if (this->ROINodeID) {
    this->Scene->RemoveReferencedNodeID(this->ROINodeID, this);
  }
  this->SetROINodeID(ID);

  if (ID) {
    this->Scene->AddReferencedNodeID(this->ROINodeID, this);
  }
}

//----------------------------------------------------------------------------
void vtkMRMLRegistrationNode::UpdateReferenceID(const char *oldID, const char *newID) {
  if (this->FixedImageNodeID && !strcmp(oldID, this->FixedImageNodeID)) {
    this->SetAndObserveFixedImageNodeID(newID);
  }
  if (this->MovingImageNodeID && !strcmp(oldID, this->MovingImageNodeID)) {
    this->SetAndObserveMovingImageNodeID(newID);
  }
  if (this->DeformedImageNodeID && !strcmp(oldID, this->DeformedImageNodeID)) {
    this->SetAndObserveDeformedImageNodeID(newID);
  }
  if (this->DeformationFieldID && !strcmp(oldID, this->DeformationFieldID)) {
    this->SetAndObserveDeformationFieldID(newID);
  }
  if (this->FiducialPointsID && !strcmp(oldID, this->FiducialPointsID)) {
    this->SetAndObserveFiducialPointsID(newID);
  }
  if (this->ROINodeID && !strcmp(oldID, this->ROINodeID)) {
    this->SetAndObserveROINodeID(newID);
  }
}

// End
