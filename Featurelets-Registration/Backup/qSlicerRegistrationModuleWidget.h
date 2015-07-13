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

#ifndef __qSlicerRegistrationModuleWidget_h
#define __qSlicerRegistrationModuleWidget_h

// SlicerQt includes
#include "qSlicerAbstractModuleWidget.h"

#include "qSlicerRegistrationModuleExport.h"

class qSlicerRegistrationModuleWidgetPrivate;
class vtkMRMLNode;
class vtkMRMLRegistrationParametersNode;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_REGISTRATION_EXPORT qSlicerRegistrationModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:

  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerRegistrationModuleWidget(QWidget *parent=0);
  virtual ~qSlicerRegistrationModuleWidget();

public slots:


protected:
  QScopedPointer<qSlicerRegistrationModuleWidgetPrivate> d_ptr;

  virtual void setup();
  virtual void enter();
  virtual void setMRMLScene(vtkMRMLScene*);

  void initializeParameterNode(vtkMRMLScene*);


protected slots:
  //void initializeNode(vtkMRMLNode*);
  void onFixedImageNodeChanged();
  void onMovingImageNodeChanged();
  void onInputVolumeAdded(vtkMRMLNode*);
  void onApply(); //for the connection to the vtkSlicerRegistrationLogic
  void onShowVolume();
  void onShowVolume2();
  void updateWidget();
  void updateParameters();
  void onEndCloseEvent();

private:
  Q_DECLARE_PRIVATE(qSlicerRegistrationModuleWidget);
  Q_DISABLE_COPY(qSlicerRegistrationModuleWidget);

  vtkMRMLRegistrationParametersNode *parametersNode;
};

#endif
