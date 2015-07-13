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

// Qt includes
#include <QDebug>
#include <QMessageBox>

// SlicerQt includes
#include "qSlicerRegistrationModuleWidget.h"
#include "ui_qSlicerRegistrationModuleWidget.h"

// Registration Logic includes
#include <vtkSlicerRegistrationLogic.h>

// SlicerQt includes
#include <qSlicerAbstractCoreModule.h>

// qMRML includes
#include <qMRMLNodeFactory.h>

// MRML includes
#include <vtkMRMLVolumeNode.h>
#include <vtkMRMLSelectionNode.h>
#include <vtkMRMLApplicationLogic.h>
#include <vtkMRMLRegistrationParametersNode.h>
#include <vtkMRMLAnnotationROINode.h>

#include <vtkNew.h>
#include <vtkMatrix4x4.h>
#include <vtkMRMLLinearTransformNode.h>

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerRegistrationModuleWidgetPrivate: public Ui_qSlicerRegistrationModuleWidget
{
    Q_DECLARE_PUBLIC(qSlicerRegistrationModuleWidget);
  protected:
    qSlicerRegistrationModuleWidget* const q_ptr;

  public:
    qSlicerRegistrationModuleWidgetPrivate(qSlicerRegistrationModuleWidget& object);
    ~qSlicerRegistrationModuleWidgetPrivate();

    vtkSlicerRegistrationLogic* logic() const;

    void performROIVoxelGridAlignment();
    bool checkForVolumeParentTransform() const;
    void showUnsupportedTransVolumeVoxelCroppingDialog() const;

    /// Using this flag prevents overriding the parameter set node contents when the
    ///   QMRMLCombobox selects the first instance of the specified node type when initializing
    bool ModuleWindowInitialized;
};

//-----------------------------------------------------------------------------
// qSlicerRegistrationModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerRegistrationModuleWidgetPrivate::qSlicerRegistrationModuleWidgetPrivate(qSlicerRegistrationModuleWidget& object)
    : q_ptr(&object)
    , ModuleWindowInitialized(false)
{
}

//-----------------------------------------------------------------------------
qSlicerRegistrationModuleWidgetPrivate::~qSlicerRegistrationModuleWidgetPrivate()
{
}

//-----------------------------------------------------------------------------
vtkSlicerRegistrationLogic* qSlicerRegistrationModuleWidgetPrivate::logic() const
{
  Q_Q(const qSlicerRegistrationModuleWidget);
  return vtkSlicerRegistrationLogic::SafeDownCast(q->logic());
}

/*//-----------------------------------------------------------------------------
void qSlicerRegistrationModuleWidgetPrivate::showUnsupportedTransVolumeVoxelCroppingDialog() const
{
  QMessageBox::information(NULL,"Registration","The selected volume is under a transform. Voxel based Registration is only supported for non transformed volumes!");
}

//-----------------------------------------------------------------------------
bool qSlicerRegistrationModuleWidgetPrivate::checkForVolumeParentTransform() const
{
  Q_ASSERT(this->FixedImageComboBox);


  vtkSmartPointer<vtkMRMLVolumeNode> inputVolume = vtkMRMLVolumeNode::SafeDownCast(this->FixedImageComboBox->currentNode());

  if(!inputVolume)
    return false;

   vtkSmartPointer<vtkMRMLLinearTransformNode> volTransform  = vtkMRMLLinearTransformNode::SafeDownCast(inputVolume->GetParentTransformNode());

   if(volTransform)
       return true;


   return false;
}

//-----------------------------------------------------------------------------
void qSlicerRegistrationModuleWidgetPrivate::performROIVoxelGridAlignment()
{
    Q_ASSERT(this->FixedImageComboBox);

    vtkSmartPointer<vtkMRMLVolumeNode> inputVolume = vtkMRMLVolumeNode::SafeDownCast(this->FixedImageComboBox->currentNode());

    if( !inputVolume )
        return;

    vtkNew<vtkMatrix4x4> volRotMat;
    bool volumeTilted = vtkSlicerRegistrationLogic::IsVolumeTiltedInRAS(inputVolume,volRotMat.GetPointer());

    if(volumeTilted)
      {
        vtkSlicerRegistrationLogic* logic = this->logic();
        Q_ASSERT(logic);
        logic->SnapROIToVoxelGrid(inputVolume);
      }
}*/

//-----------------------------------------------------------------------------
// qSlicerRegistrationModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerRegistrationModuleWidget::qSlicerRegistrationModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerRegistrationModuleWidgetPrivate(*this) )
{
    this->parametersNode = NULL;
}

//-----------------------------------------------------------------------------
qSlicerRegistrationModuleWidget::~qSlicerRegistrationModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerRegistrationModuleWidget::setup()
{
  Q_D(qSlicerRegistrationModuleWidget);
  d->setupUi(this);
  this->Superclass::setup();

  /*connect(d->ROIComboBox->nodeFactory(), SIGNAL(nodeInitialized(vtkMRMLNode*)),
          this, SLOT(initializeNode(vtkMRMLNode*)));*/

  connect(d->RunButton, SIGNAL(clicked()),
          this, SLOT(onApply()) );
  connect(d->ShowButton, SIGNAL(clicked()),
          this, SLOT(onShowVolume()) );
  connect(d->ShowButton2, SIGNAL(clicked()),
          this, SLOT(onShowVolume2()) );

  connect(d->FixedImageComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)),
          this, SLOT(onFixedImageNodeChanged()));
  connect(d->FixedImageComboBox, SIGNAL(nodeAdded(vtkMRMLNode*)),
          this, SLOT(onInputVolumeAdded(vtkMRMLNode*)));
  connect(d->MovingImageComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)),
          this, SLOT(onMovingImageNodeChanged()));
  connect(d->MovingImageComboBox, SIGNAL(nodeAdded(vtkMRMLNode*)),
          this, SLOT(onInputVolumeAdded(vtkMRMLNode*)));

}

//-----------------------------------------------------------------------------
void qSlicerRegistrationModuleWidget::enter()
{
  // make sure that there's a parameters node so if there are already some
  // volumes in the scene, they can be set up for use
  this->onFixedImageNodeChanged();
  this->onMovingImageNodeChanged();

  this->Superclass::enter();
}

//-----------------------------------------------------------------------------
void qSlicerRegistrationModuleWidget::setMRMLScene(vtkMRMLScene* scene)
{
  Q_D(qSlicerRegistrationModuleWidget);

  this->Superclass::setMRMLScene(scene);

  if(scene == NULL)
    {
    return;
    }

  this->initializeParameterNode(scene);

  this->updateWidget();
}

//-----------------------------------------------------------------------------
void qSlicerRegistrationModuleWidget::initializeParameterNode(vtkMRMLScene* scene)
{
  vtkCollection* parameterNodes = scene->GetNodesByClass("vtkMRMLRegistrationParametersNode");

  if(parameterNodes->GetNumberOfItems() > 0)
    {
    this->parametersNode = vtkMRMLRegistrationParametersNode::SafeDownCast(parameterNodes->GetItemAsObject(0));
    if(!this->parametersNode)
      {
      qCritical() << "FATAL ERROR: Cannot instantiate RegistrationParameterNode";
      Q_ASSERT(this->parametersNode);
      }
    }
  else
    {
    qDebug() << "No Registration parameter nodes found!";
    this->parametersNode = vtkMRMLRegistrationParametersNode::New();
    scene->AddNode(this->parametersNode);
    this->parametersNode->Delete();
    }

  parameterNodes->Delete();
}

//-----------------------------------------------------------------------------
void qSlicerRegistrationModuleWidget::onApply(){

  Q_D(const qSlicerRegistrationModuleWidget);
  vtkSlicerRegistrationLogic *logic = d->logic();

  if( !this->parametersNode || !d->FixedImageComboBox->currentNode() || !d->MovingImageComboBox->currentNode() )
    return;

  this->parametersNode->SetFixedImageNodeID(d->FixedImageComboBox->currentNode()->GetID());
  this->parametersNode->SetMovingImageNodeID(d->MovingImageComboBox->currentNode()->GetID());

  if(!logic->Apply(this->parametersNode))
    {
    std::cerr << "Propagating to the selection node" << std::endl;
    vtkSlicerApplicationLogic *appLogic = this->module()->appLogic();
    vtkMRMLSelectionNode *selectionNode = appLogic->GetSelectionNode();
    selectionNode->SetReferenceActiveVolumeID(this->parametersNode->GetOutputVolumeNodeID());
    appLogic->PropagateVolumeSelection();
    }
}

//-----------------------------------------------------------------------------
void qSlicerRegistrationModuleWidget::onShowVolume(){
    Q_D(const qSlicerRegistrationModuleWidget);
    vtkSlicerRegistrationLogic *logic = d->logic();

    if( !this->parametersNode || !d->FixedImageComboBox->currentNode() )
      return;

    this->parametersNode->SetFixedImageNodeID(d->FixedImageComboBox->currentNode()->GetID());
    bool fixedImage = true;

    if(!logic->ShowVolume(this->parametersNode, fixedImage))
      {
      std::cerr << "Propagating to the selection node" << std::endl;
      vtkSlicerApplicationLogic *appLogic = this->module()->appLogic();
      vtkMRMLSelectionNode *selectionNode = appLogic->GetSelectionNode();
      selectionNode->SetReferenceActiveVolumeID(this->parametersNode->GetOutputVolumeNodeID());
      appLogic->PropagateVolumeSelection();
      }
}

//-----------------------------------------------------------------------------
void qSlicerRegistrationModuleWidget::onShowVolume2(){
    Q_D(const qSlicerRegistrationModuleWidget);
    vtkSlicerRegistrationLogic *logic = d->logic();

    if(!this->parametersNode || !d->MovingImageComboBox->currentNode())
      return;

    this->parametersNode->SetMovingImageNodeID(d->MovingImageComboBox->currentNode()->GetID());
    bool fixedImage = false;

    if(!logic->ShowVolume(this->parametersNode, fixedImage))
      {
      std::cerr << "Propagating to the selection node" << std::endl;
      vtkSlicerApplicationLogic *appLogic = this->module()->appLogic();
      vtkMRMLSelectionNode *selectionNode = appLogic->GetSelectionNode();
      selectionNode->SetReferenceActiveVolumeID(this->parametersNode->GetOutputVolumeNodeID());
      appLogic->PropagateVolumeSelection();
      }
}

//-----------------------------------------------------------------------------
void qSlicerRegistrationModuleWidget::onFixedImageNodeChanged()
{
  Q_D(qSlicerRegistrationModuleWidget);
  Q_ASSERT(d->FixedImageComboBox);

  vtkMRMLNode* node = d->FixedImageComboBox->currentNode();
  if(node)
    {
    if(d->checkForVolumeParentTransform())
      {
        d->showUnsupportedTransVolumeVoxelCroppingDialog();
        d->FixedImageComboBox->setCurrentNode(NULL);
      }
    else
      {
        d->performROIVoxelGridAlignment();
      }
    }
}

//-----------------------------------------------------------------------------
void qSlicerRegistrationModuleWidget::onMovingImageNodeChanged()
{
  Q_D(qSlicerRegistrationModuleWidget);
  Q_ASSERT(d->MovingImageComboBox);

  vtkMRMLNode* node = d->MovingImageComboBox->currentNode();
  if(node)
    {
    if(d->checkForVolumeParentTransform())
      {
        d->showUnsupportedTransVolumeVoxelCroppingDialog();
        d->MovingImageComboBox->setCurrentNode(NULL);
      }
    else
      {
        d->performROIVoxelGridAlignment();
      }
    }
}

//-----------------------------------------------------------------------------
void qSlicerRegistrationModuleWidget::onInputVolumeAdded(vtkMRMLNode *mrmlNode)
{
  Q_D(qSlicerRegistrationModuleWidget);
  if (!mrmlNode)
    {
    return;
    }

  if (d->FixedImageComboBox->currentNode() != NULL || d->MovingImageComboBox->currentNode() != NULL)
    {
    // there's already a selected node, don't reset it
    return;
    }
  d->FixedImageComboBox->setCurrentNode(mrmlNode);
  d->MovingImageComboBox->setCurrentNode(mrmlNode);
}

//-----------------------------------------------------------------------------
void qSlicerRegistrationModuleWidget::onEndCloseEvent()
{
  this->initializeParameterNode(this->mrmlScene());
}

//-----------------------------------------------------------------------------
void qSlicerRegistrationModuleWidget::updateParameters()
{
  Q_D(qSlicerRegistrationModuleWidget);
  if(!this->parametersNode)
    return;
  vtkMRMLRegistrationParametersNode *pNode = this->parametersNode;

  vtkMRMLNode *fixedImageNode = d->FixedImageComboBox->currentNode();
  vtkMRMLNode *movingImageNode = d->MovingImageComboBox->currentNode();

  if(fixedImageNode)
    pNode->SetFixedImageNodeID(fixedImageNode->GetID());
  else if(movingImageNode)
    pNode->SetMovingImageNodeID(movingImageNode->GetID());
  else
    pNode->SetFixedImageNodeID(NULL);
    pNode->SetMovingImageNodeID(NULL);

}

//-----------------------------------------------------------------------------
void qSlicerRegistrationModuleWidget::updateWidget()
{
  Q_D(qSlicerRegistrationModuleWidget);
  if (!this->parametersNode || !this->mrmlScene())
    {
    return;
    }
  vtkMRMLRegistrationParametersNode *parameterNode = this->parametersNode;

  vtkMRMLNode *fixedImageNode = this->mrmlScene()->GetNodeByID(parameterNode->GetFixedImageNodeID());
  vtkMRMLNode *movingImageNode = this->mrmlScene()->GetNodeByID(parameterNode->GetMovingImageNodeID());

  if(fixedImageNode || movingImageNode)
      d->FixedImageComboBox->setCurrentNode(fixedImageNode);
      d->MovingImageComboBox->setCurrentNode(movingImageNode);
}
