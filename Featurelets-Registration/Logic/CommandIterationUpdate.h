#ifndef COMMANDITERATIONUPDATE_H
#define COMMANDITERATIONUPDATE_H


#include "itkCommand.h"
#include "itkObject.h"
#include "itkRegularStepGradientDescentOptimizer.h"

class CommandIterationUpdate : public itk::Command {
public:
  typedef CommandIterationUpdate   Self;
  typedef itk::Command             Superclass; //Used to monitor the evolution of the registratioin process
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() {}

public:
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef const OptimizerType *OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject &event) {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object *object, const itk::EventObject &event) {
    if( !itk::IterationEvent().CheckEvent(&event) ) {
      return;
    }
    OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >( object );
    std::cout << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetValue() << "   ";
    std::cout << optimizer->GetCurrentPosition() << std::endl;
  }
};

#endif // COMMANDITERATIONUPDATE_H
