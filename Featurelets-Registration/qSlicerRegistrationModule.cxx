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
#include <QtPlugin>

// Slicer includes
#include <qSlicerCoreApplication.h>
#include <qSlicerModuleManager.h>

// Registration Logic includes
#include <vtkSlicerRegistrationLogic.h>
#include <vtkSlicerVolumesLogic.h>

// Registration includes
#include "qSlicerRegistrationModule.h"
#include "qSlicerRegistrationModuleWidget.h"

//-----------------------------------------------------------------------------
Q_EXPORT_PLUGIN2(qSlicerRegistrationModule, qSlicerRegistrationModule);

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerRegistrationModulePrivate
{
public:
  qSlicerRegistrationModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerRegistrationModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerRegistrationModulePrivate::qSlicerRegistrationModulePrivate() {
}

//-----------------------------------------------------------------------------
// qSlicerRegistrationModule methods

//-----------------------------------------------------------------------------
qSlicerRegistrationModule::qSlicerRegistrationModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerRegistrationModulePrivate) {
}

//-----------------------------------------------------------------------------
qSlicerRegistrationModule::~qSlicerRegistrationModule() {
}

//-----------------------------------------------------------------------------
QString qSlicerRegistrationModule::helpText() const {
  return "<p>The ,,Featurelets-Registration'' module uses pixeltype short. To convert an image from another pixeltype use e.g. the module ,,Cast Scalar Volume''.</p> <p>To threshold an image use e.g. the module ,,Threshold Scalar Volume''.</p> <p>For initial alignment of fixed and moving image e.g. the module ,,Expert Automated Registration'' can be used.</p> <p>If only a subvolume of the whole image should be deformed use the module ,,Crop Volume''.</p> <p>The output deformation field can be visualized using the module ,,Transforms''.</p>";
}

//-----------------------------------------------------------------------------
QString qSlicerRegistrationModule::acknowledgementText() const {
  return "<p>The ,,Featurelets-Registration'' module was developed as a Master Thesis at AKH Vienna under supervision of Univ.-Prof. Dietmar Georg and Hugo Furtado.</p> <p>The idea of the Featurelet algorithm is based on Soehn, M. et al. (2008). The ,,Featurelets-Registration'' module is an implementation of the algorithm developed by Fabri, D. et al. (2014).</p>";
}

//-----------------------------------------------------------------------------
QStringList qSlicerRegistrationModule::contributors() const {
  QStringList moduleContributors;
  moduleContributors << QString("Andreas Renner (Vienna UT)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerRegistrationModule::icon() const {
  return QIcon(":/Icons/Registration.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerRegistrationModule::categories() const {
  return QStringList() << "Registration";
}

//-----------------------------------------------------------------------------
QStringList qSlicerRegistrationModule::dependencies() const {
  return QStringList() << "Volumes";
}

//-----------------------------------------------------------------------------
void qSlicerRegistrationModule::setup() {
  this->Superclass::setup();

  vtkSlicerRegistrationLogic* RegistrationLogic =
    vtkSlicerRegistrationLogic::SafeDownCast(this->logic());

  qSlicerAbstractCoreModule* volumesModule =
    qSlicerCoreApplication::application()->moduleManager()->module("Volumes");
  if (volumesModule) {
    vtkSlicerVolumesLogic* volumesLogic =
      vtkSlicerVolumesLogic::SafeDownCast(volumesModule->logic());
    RegistrationLogic->SetVolumesLogic(volumesLogic);
  }
  else {
    qWarning() << "Volumes module is not found";
  }
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerRegistrationModule::createWidgetRepresentation() {
  return new qSlicerRegistrationModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerRegistrationModule::createLogic() {
  return vtkSlicerRegistrationLogic::New();
}
