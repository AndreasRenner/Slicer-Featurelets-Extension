/*
 * FeatureletRegistrationResult.h
 *
 *  Created on: Dec 3, 2009
 *      Author: danifabri
 */

#ifndef FEATURELETREGISTRATIONRESULT_H_
#define FEATURELETREGISTRATIONRESULT_H_

#include "itkObject.h"

namespace Featurelet
{
  enum Status
  {
    OK=0,
    SizeError=1,
    sameColor=2,
    internalError=3
  };
}

class ITK_EXPORT FeatureletRegistrationResult:public
itk::Object
{
public:
  /** Standard class typedefs. */
typedef FeatureletRegistrationResult   Self;
typedef itk::Object                           Superclass;
typedef itk::SmartPointer<Self>               Pointer;
typedef itk::SmartPointer<const Self>         ConstPointer;

typedef itk::TranslationTransform<double> TransformType;
typedef TransformType::Pointer TransformPointer;
typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
typedef	OptimizerType::Pointer OptimizerPointer;

 /** Method for creation through the object factory. */
  itkNewMacro(Self);

 /** Run-time type information (and related methods). */
  itkTypeMacro(FeatureletRegistrationResult, itk::Object);

virtual void PrintSelf(std::ostream& os,itk::Indent indent)
const {
  Superclass::PrintSelf(os, indent);
os << indent
<< "Offset: " << m_Transform->GetOffset() << "\n";
os <<indent
<< "Metric Value II: "<< m_Optimizer->GetValue()<<"\n";
os << indent
<< "status Fixed featurelet: " << m_StatusFixed << "\n";
os << indent
<< "status Moving featurelet: " << m_StatusMoving << "\n";
os << indent
<< "Number of Featurelet: " << m_FeatureletNumber<< "\n";

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
	virtual ~FeatureletRegistrationResult() {};

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

#endif /* FEATURELETREGISTRATIONRESULT_H_ */
