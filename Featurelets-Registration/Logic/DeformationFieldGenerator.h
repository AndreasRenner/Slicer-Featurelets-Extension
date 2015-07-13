/*
 * DeformationFieldGenerator.h
 *
 *  Created on: Jan 20, 2010
 *      Author: danifabri
 */

#ifndef DEFORMATIONFIELDGENERATOR_H_
#define DEFORMATIONFIELDGENERATOR_H_


#endif /* DEFORMATIONFIELDGENERATOR_H_ */

#include "vtkSlicerRegistrationLogic.h"

#include "itkExceptionObject.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkObject.h"
#include "itkResampleImageFilter.h"
#include "itkVectorResampleImageFilter.h"
#include "itkWarpImageFilter.h"

#include <vtkGridTransform.h>
#include <vtkImageData.h>
#include <vtkMatrix4x4.h>
#include <vtkNew.h>
#include <vtkOrientedGridTransform.h>


/** \class NumericTraits< std::complex<double> >
 *  \brief Define traits for type std::complex<double>.
 *  \ingroup DataRepresentation
 */


class DeformationFieldGenerator {
public:
    typedef std::list<FeatureletRegistrationResultPointer> FeatureletRegistrationResultListType;
    static const unsigned int VectorDimension = 1;
	typedef double PixelType;
    typedef itk::Image<PixelType, Dimension> DeformationFieldType;
    typedef itk::Vector<double, Dimension> PixelTypeVector;
    typedef itk::Image<PixelTypeVector, Dimension> DeformationFieldTypeVector;
    typedef itk::TranslationTransform<double> TransformType;

	DeformationFieldGenerator() {}
	virtual ~DeformationFieldGenerator() {}


    void AddFeaturelet(FeatureletRegistrationResult *s) {
      //std::cerr << "DeformationFieldGenerator - Entered Add Featurelet" << std::endl;
      FeatureletRegistrationResult::Pointer regResult;
      regResult = s;
      m_RegistrationResultList.push_back(regResult);
	}

    void GenerateDeformationField() {
      //std::cerr << "DeformationFieldGenerator - Entered GenerateDeformationField" << std::endl;
      FeatureletRegistrationResultListType::iterator it;

      DeformationFieldType::IndexType start;
      start[0] = start[1] = start[2] = 0.0;

      DeformationFieldType::SizeType size;
      size[0] = m_NumberOfFeaturelets[0]+1; //final FeatureletGridPosition
      size[1] = m_NumberOfFeaturelets[1]+1;
      size[2] = m_NumberOfFeaturelets[2]+1;

      DeformationFieldType::RegionType region0; // a cambiar??????????????????????
      region0.SetSize(size);
      region0.SetIndex(start);

      double spacing1[Dimension];
      spacing1[0] = m_SizeOfFeaturelets[0]*m_ImageSpacingFix[0];
      spacing1[1] = m_SizeOfFeaturelets[1]*m_ImageSpacingFix[1];
      spacing1[2] = m_SizeOfFeaturelets[2]*m_ImageSpacingFix[2];

      m_DeformationFieldRawX = DeformationFieldType::New();
      m_DeformationFieldRawX->SetRegions(region0);
      m_DeformationFieldRawX->Allocate();
      m_DeformationFieldRawX->SetSpacing(spacing1);
      m_DeformationFieldRawX->SetOrigin(m_ImageOriginFix);   ///TEST _Fix

      PixelType zeroVector = 0;
      m_DeformationFieldRawX->FillBuffer(zeroVector);

      // the other two components:
      m_DeformationFieldRawY = DeformationFieldType::New();
      m_DeformationFieldRawY->SetRegions(region0);
      m_DeformationFieldRawY->Allocate();
      m_DeformationFieldRawY->SetSpacing(spacing1);
      m_DeformationFieldRawY->SetOrigin(m_ImageOriginFix);   ///TEST _Fix
      m_DeformationFieldRawY->FillBuffer(zeroVector);

      m_DeformationFieldRawZ = DeformationFieldType::New();
      m_DeformationFieldRawZ->SetRegions(region0);
      m_DeformationFieldRawZ->SetSpacing(spacing1);
      m_DeformationFieldRawZ->SetOrigin(m_ImageOriginFix);   ///TEST _Fix
      m_DeformationFieldRawZ->Allocate();
      m_DeformationFieldRawZ->FillBuffer(zeroVector);

      m_MetricValue = DeformationFieldType::New();
      m_MetricValue->SetRegions(region0);
      m_MetricValue->SetSpacing(spacing1);
      m_MetricValue->Allocate();
      m_MetricValue->FillBuffer(zeroVector);

      //fill the image
      double bestValue;
      metricthreshold = 1000;

      PixelType transformVectorX, transformVectorY, transformVectorZ;

      int underthreshold = 0;
      int abovethreshold = 0;

      for (it = m_RegistrationResultList.begin(); it != m_RegistrationResultList.end(); it++) {
        DeformationFieldType::IndexType pixelIndex ;
        TransformType::OutputVectorType offset;
        offset[0]= (*it)->GetTransform()->GetOffset()[0];
        offset[1]= (*it)->GetTransform()->GetOffset()[1];
        offset[2]= (*it)->GetTransform()->GetOffset()[2];
        pixelIndex[0]= (*it)->GetFeatureletGridPosition()[0];
        pixelIndex[1]= (*it)->GetFeatureletGridPosition()[1];
        pixelIndex[2]= (*it)->GetFeatureletGridPosition()[2];

        transformVectorX = offset[0];
        transformVectorY = offset[1];
        transformVectorZ = offset[2];

        bestValue = (*it)->GetOptimizer()->GetValue();
        //std::cout << "best metric value " << bestValue << std::endl;
        if (bestValue < metricthreshold) {
          abovethreshold++;
          m_DeformationFieldRawX->SetPixel(pixelIndex, transformVectorX);
          m_DeformationFieldRawY->SetPixel(pixelIndex, transformVectorY);
          m_DeformationFieldRawZ->SetPixel(pixelIndex, transformVectorZ);
          //std::cout<<"pixel "<<pixelIndex<<" vector: ["<<transformVector<<", "<<transformVector1<<", "<<transformVector2<< "] "<<std::endl;
        }
        else {
          underthreshold++;
          m_DeformationFieldRawX->SetPixel(pixelIndex, 0);
          m_DeformationFieldRawY->SetPixel(pixelIndex, 0);
          m_DeformationFieldRawZ->SetPixel(pixelIndex, 0);
        }
      }
    }

    void RescaleDeformationField() {
      //std::cerr << "DeformationFieldGenerator - Entered RescaleDeformationField" << std::endl;
      typedef itk::ResampleImageFilter<DeformationFieldType,DeformationFieldType> FilterType;
      FilterType::Pointer filterX = FilterType::New();
      FilterType::Pointer filterY = FilterType::New();
      FilterType::Pointer filterZ = FilterType::New();
      TransformType::Pointer transform = TransformType::New();

      typedef itk::LinearInterpolateImageFunction< DeformationFieldType, double > InterpolatorType;
      InterpolatorType::Pointer interpolator = InterpolatorType::New();

      DeformationFieldType::IndexType startII;
      startII[0] = startII[1] = startII[2] = 0.0;

      DeformationFieldType::SizeType sizeII;
      sizeII[0] = m_ImageSizeFix[0];
      sizeII[1] = m_ImageSizeFix[1];
      sizeII[2] = m_ImageSizeFix[2];

      DeformationFieldType::SizeType sizeFeaturelet;
      sizeFeaturelet[0] = m_NumberOfFeaturelets[0]+1;
      sizeFeaturelet[1] = m_NumberOfFeaturelets[1]+1;
      sizeFeaturelet[2] = m_NumberOfFeaturelets[2]+1;

      m_DeformationFieldRawXI = DeformationFieldType::New();
      m_DeformationFieldRawYI = DeformationFieldType::New();
      m_DeformationFieldRawZI = DeformationFieldType::New();

      DeformationFieldType::RegionType regionII;
      regionII.SetSize(sizeII);
      regionII.SetIndex(startII);

      m_DeformationFieldRawXI->SetRegions(regionII);
      m_DeformationFieldRawYI->SetRegions(regionII);
      m_DeformationFieldRawZI->SetRegions(regionII);

      m_DeformationFieldRawXI->Allocate();
      m_DeformationFieldRawYI->Allocate();
      m_DeformationFieldRawZI->Allocate();

      PixelType zeroVector = 0;
      m_DeformationFieldRawXI->FillBuffer(zeroVector);
      m_DeformationFieldRawYI->FillBuffer(zeroVector);
      m_DeformationFieldRawZI->FillBuffer(zeroVector);

      filterX->SetTransform(transform);
      filterX->SetInterpolator(interpolator);
      filterY->SetTransform(transform);
      filterY->SetInterpolator(interpolator);
      filterZ->SetTransform(transform);
      filterZ->SetInterpolator(interpolator);

      double origin[Dimension];
      origin[0] = m_ImageOriginFix[0];
      origin[1] = m_ImageOriginFix[1];
      origin[2] = m_ImageOriginFix[2];

      double spacingFeaturelet[Dimension];
      spacingFeaturelet[0] = m_SizeOfFeaturelets[0]*m_ImageSpacingFix[0];
      spacingFeaturelet[1] = m_SizeOfFeaturelets[1]*m_ImageSpacingFix[1];
      spacingFeaturelet[2] = m_SizeOfFeaturelets[2]*m_ImageSpacingFix[2];

      DeformationFieldType::SpacingType spacing;
      spacing[0] = m_ImageSpacingFix[0];
      spacing[1] = m_ImageSpacingFix[1];
      spacing[2] = m_ImageSpacingFix[2];

      double o2[3];
      for(unsigned int i=0; i<3; i++) {
        double center = origin[i] + (sizeFeaturelet[i] * spacingFeaturelet[i]) / 2.0;
        o2[i] = center - sizeII[i] * spacing[i] / 2.0;
      }

      filterX->SetOutputOrigin(o2);
      filterX->SetSize( sizeII );
      filterX->SetOutputSpacing(spacing);
      filterY->SetOutputOrigin(o2);
      filterY->SetSize( sizeII );
      filterY->SetOutputSpacing(spacing);
      filterZ->SetOutputOrigin(o2);
      filterZ->SetSize( sizeII );
      filterZ->SetOutputSpacing(spacing);

      DeformationFieldType::DirectionType direction;
      direction.SetIdentity();

      filterX->SetOutputDirection( direction );
      filterY->SetOutputDirection( direction );
      filterZ->SetOutputDirection( direction );

      filterX->SetInput(m_DeformationFieldRawX);
      filterY->SetInput(m_DeformationFieldRawY);
      filterZ->SetInput(m_DeformationFieldRawZ);

      filterX->Update();
      filterX->UpdateLargestPossibleRegion();
      filterY->Update();
      filterY->UpdateLargestPossibleRegion();
      filterZ->Update();
      filterZ->UpdateLargestPossibleRegion();

      m_DeformationFieldRawXI = filterX->GetOutput();
      m_DeformationFieldRawYI = filterY->GetOutput();
      m_DeformationFieldRawZI = filterZ->GetOutput();

      m_DeformationFieldRawXI->Update();
      m_DeformationFieldRawYI->Update();
      m_DeformationFieldRawZI->Update();
    }

    void TotalDeformationFieldCreator() {
      //std::cerr << "DeformationFieldGenerator - Entered TotalDeformationField" << std::endl;
      DeformationFieldTypeVector::IndexType startII;
      startII[0] = startII[1] = startII[2] = 0.0;

      DeformationFieldTypeVector::SizeType sizeII;
      sizeII[0] = m_ImageSizeFix[0];
      sizeII[1] = m_ImageSizeFix[1];
      sizeII[2] = m_ImageSizeFix[2];
      DeformationFieldTypeVector::RegionType regionII;
      regionII.SetSize(sizeII);
      regionII.SetIndex(startII);

      m_DeformationFieldRawII = DeformationFieldTypeVector::New();
      m_DeformationFieldRawII->SetRegions(regionII);

      DeformationFieldTypeVector::SpacingType spacingII;
      spacingII[0] = m_ImageSpacingFix[0];
      spacingII[1] = m_ImageSpacingFix[1];
      spacingII[2] = m_ImageSpacingFix[2];
      m_DeformationFieldRawII->SetSpacing(spacingII);

      m_DeformationFieldRawII->SetOrigin(m_ImageOriginFix);

      try {
        m_DeformationFieldRawII->Allocate();
      }
      catch(itk::ExceptionObject &e) {
        std::cerr << "Failed to allocate memory for the ITK Deformation Field: " << e.GetDescription() << std::endl;
        return;
      }

      PixelTypeVector zeroVector;
      zeroVector.Fill(0);
      m_DeformationFieldRawII->FillBuffer(zeroVector);

      PixelTypeVector TotalDeformationVector;
      DeformationFieldTypeVector::IndexType pixelIndexII;

      for (int i=1; i<m_ImageSizeFix[0]-1; i++) {     // Index starts at [1] and ends at [m_ImageSizeMov-1]
        for (int j=1; j<m_ImageSizeFix[1]-1; j++) {   // because on the borders the Deformationfield should
          for (int k=1; k<m_ImageSizeFix[2]-1; k++) { // be equal to the zeroVector.
            pixelIndexII[0] = i;
            pixelIndexII[1] = j;
            pixelIndexII[2] = k;

            TotalDeformationVector[0] = m_DeformationFieldRawXI->GetPixel(pixelIndexII);
            TotalDeformationVector[1] = m_DeformationFieldRawYI->GetPixel(pixelIndexII);
            TotalDeformationVector[2] = m_DeformationFieldRawZI->GetPixel(pixelIndexII);

            m_DeformationFieldRawII->SetPixel(pixelIndexII, TotalDeformationVector);
          }
        }
      }
	}

    //Taken from vtkITKTransformConverter.h
    vtkSmartPointer<vtkOrientedGridTransform> VTKDeformationFieldCreator() {
      //std::cerr << "DeformationFieldGenerator - Entered VTKDeformationFieldCreator" << std::endl;

      grid_Ras = vtkSmartPointer<vtkOrientedGridTransform>::New();
      deformationField = vtkSmartPointer<vtkImageData>::New();

      deformationField->SetOrigin(m_DeformationFieldRawII->GetOrigin()[0], m_DeformationFieldRawII->GetOrigin()[1], m_DeformationFieldRawII->GetOrigin()[2]);
      deformationField->SetSpacing(m_DeformationFieldRawII->GetSpacing()[0], m_DeformationFieldRawII->GetSpacing()[1], m_DeformationFieldRawII->GetSpacing()[2]);
      const vtkSmartPointer<vtkMatrix4x4> gridDirectionMatrix_LPS = vtkSmartPointer<vtkMatrix4x4>::New();
      for (unsigned int row=0; row<3; row++) {
        for (unsigned int column=0; column<3; column++) {
          gridDirectionMatrix_LPS->SetElement(row, column, m_ImageDirectionFix(row,column));
        }
      }
      /*const vtkSmartPointer<vtkMatrix4x4> lpsToRas = vtkSmartPointer<vtkMatrix4x4>::New();
      lpsToRas->SetElement(0,0,-1);
      lpsToRas->SetElement(1,1,-1);
      vtkSmartPointer<vtkMatrix4x4> gridDirectionMatrix_RAS = vtkSmartPointer<vtkMatrix4x4>::New();
      vtkMatrix4x4::Multiply4x4(lpsToRas, gridDirectionMatrix_LPS, gridDirectionMatrix_RAS);*/
      grid_Ras->SetGridDirectionMatrix(gridDirectionMatrix_LPS.GetPointer());
                  //gridDirectionMatrix_RAS.GetPointer());

      DeformationFieldTypeVector::RegionType region = m_DeformationFieldRawII->GetBufferedRegion();
      DeformationFieldTypeVector::SizeType imageSize = region.GetSize();

      int extent[6]={0, (int) imageSize[0]-1, 0, (int) imageSize[1]-1, 0, (int) imageSize[2]-1};
      deformationField->SetExtent(extent);
    #if (VTK_MAJOR_VERSION <= 5)
      deformationField->SetScalarTypeToDouble();
      deformationField->SetNumberOfScalarComponents(3);
      deformationField->AllocateScalars();
    #else
      deformationField->AllocateScalars(VTK_DOUBLE, 3);
    #endif

      double* displacementVectors_Ras = reinterpret_cast<double*>(deformationField->GetScalarPointer());
      itk::ImageRegionConstIterator<DeformationFieldTypeVector> inputIt(m_DeformationFieldRawII, m_DeformationFieldRawII->GetRequestedRegion());
      inputIt.GoToBegin();
      while( !inputIt.IsAtEnd() ) {
        DeformationFieldTypeVector::PixelType displacementVectorLps = inputIt.Get();
        *(displacementVectors_Ras++) = displacementVectorLps[0];
        *(displacementVectors_Ras++) = displacementVectorLps[1];
        *(displacementVectors_Ras++) = displacementVectorLps[2];
        ++inputIt;
      }

    #if (VTK_MAJOR_VERSION <= 5)
      grid_Ras->SetDisplacementGrid( deformationField.GetPointer() );
    #else
      grid_Ras->SetDisplacementGridData( deformationField.GetPointer() );
    #endif
      return grid_Ras;
    }

    int  WarpImagebyDeformationField(ImageType* InputImage, ImageType* FixedImage, ImageType* OutputImage) {
      //std::cerr << "DeformationFieldGenerator - Entered WarpImagebyDeformationField" << std::endl;
      typedef itk::WarpImageFilter<ImageType, ImageType, DeformationFieldTypeVector> WarperType;
      typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
      WarperType::Pointer warper = WarperType::New();
      InterpolatorType::Pointer interpolator = InterpolatorType::New();

      DeformationFieldTypeVector::SizeType sizeII;
      sizeII[0] = m_ImageSizeFix[0];
      sizeII[1] = m_ImageSizeFix[1];
      sizeII[2] = m_ImageSizeFix[2];
      warper->SetOutputSize(sizeII);
      //std::cerr << sizeII[0] << sizeII[1] << sizeII[2] << std::endl;

      warper->SetInput(InputImage);
      warper->SetInterpolator(interpolator);
      warper->SetOutputSpacing(FixedImage->GetSpacing());
      warper->SetOutputOrigin(FixedImage->GetOrigin());
      warper->SetOutputDirection(FixedImage->GetDirection());
      warper->SetDisplacementField(m_DeformationFieldRawII);
      warper->UpdateLargestPossibleRegion();
      warper->Update();
      OutputImage->Graft(warper->GetOutput());
      try {
        OutputImage->Update();
      }
      catch(itk::ExceptionObject &e) {
          std::cerr << e << std::endl;
          return EXIT_FAILURE;
      }

      return EXIT_SUCCESS;
    }

    void SetNumberOfFeaturelets(const int *_NumberOfFeaturelets) {
      //std::cerr << "DeformationFieldGenerator - Entered SetNumberOfFeaturelets" << std::endl;
      m_NumberOfFeaturelets[0] = _NumberOfFeaturelets[0];
      m_NumberOfFeaturelets[1] = _NumberOfFeaturelets[1];
      m_NumberOfFeaturelets[2] = _NumberOfFeaturelets[2];
	}

    void SetSizeOfFeaturelets(const int *_SizeOfFeaturelets) {
      //std::cerr << "DeformationFieldGenerator - Entered SetSizeOfFeaturelets" << std::endl;
      m_SizeOfFeaturelets[0] = _SizeOfFeaturelets[0];
      m_SizeOfFeaturelets[1] = _SizeOfFeaturelets[1];
      m_SizeOfFeaturelets[2] = _SizeOfFeaturelets[2];
	}

    void SetImageSizeMov(const double* _ImageSizeMov) {
      //std::cerr << "DeformationFieldGenerator - Entered SetImageSizeMov" << std::endl;
      m_ImageSizeMov[0] = _ImageSizeMov[0];
      m_ImageSizeMov[1] = _ImageSizeMov[1];
      m_ImageSizeMov[2] = _ImageSizeMov[2];
	}

    void SetImageSpacingMov(const double* _ImageSpacingMov) {
      //std::cerr << "DeformationFieldGenerator - Entered SetImageSpacingMov" << std::endl;
      m_ImageSpacingMov[0] = _ImageSpacingMov[0];
      m_ImageSpacingMov[1] = _ImageSpacingMov[1];
      m_ImageSpacingMov[2] = _ImageSpacingMov[2];
    }

    void SetImageOriginMov(const double* _ImageOriginMov) {
      //std::cerr << "DeformationFieldGenerator - Entered SetImageOriginMov" << std::endl;
      m_ImageOriginMov[0] = _ImageOriginMov[0];
      m_ImageOriginMov[1] = _ImageOriginMov[1];
      m_ImageOriginMov[2] = _ImageOriginMov[2];
    }

    void SetImageSizeFix(const double* _ImageSizeFix) {
      //std::cerr << "DeformationFieldGenerator - Entered SetImageSizeFix" << std::endl;
      m_ImageSizeFix[0] = _ImageSizeFix[0];
      m_ImageSizeFix[1] = _ImageSizeFix[1];
      m_ImageSizeFix[2] = _ImageSizeFix[2];
    }

    void SetImageSpacingFix(const double* _ImageSpacingFix) {
      //std::cerr << "DeformationFieldGenerator - Entered SetImageSpacingFix" << std::endl;
      m_ImageSpacingFix[0] = _ImageSpacingFix[0];
      m_ImageSpacingFix[1] = _ImageSpacingFix[1];
      m_ImageSpacingFix[2] = _ImageSpacingFix[2];
    }

    void SetImageOriginFix(const double* _ImageOriginFix) {
      //std::cerr << "DeformationFieldGenerator - Entered SetImageOriginFix" << std::endl;
      m_ImageOriginFix[0] = _ImageOriginFix[0];
      m_ImageOriginFix[1] = _ImageOriginFix[1];
      m_ImageOriginFix[2] = _ImageOriginFix[2];
    }

    void SetImageDirectionFix(DeformationFieldType::DirectionType _ImageDirectionFix){
      //std::cerr << "DeformationFieldGenerator - Entered SetImageDirectionFix" << std::endl;
      m_ImageDirectionFix = _ImageDirectionFix;
    }

public:
    FeatureletRegistrationResultListType m_RegistrationResultList;
	int m_NumberOfFeaturelets[3];
    int m_SizeOfFeaturelets[3];
	int m_ImageSizeMov[3];
    int m_ImageSizeFix[3];
    double m_ImageSpacingMov[3];
    double m_ImageOriginMov[3];
    double m_ImageSpacingFix[3];
    double m_ImageOriginFix[3];
	double metricthreshold;
    DeformationFieldType::DirectionType m_ImageDirectionFix;
	DeformationFieldType::Pointer m_MetricValue;
	DeformationFieldType::Pointer m_DeformationFieldRawZ;
	DeformationFieldType::Pointer m_DeformationFieldRawY;
	DeformationFieldType::Pointer m_DeformationFieldRawX;
	DeformationFieldType::Pointer m_DeformationFieldRawZI;
	DeformationFieldType::Pointer m_DeformationFieldRawYI;
	DeformationFieldType::Pointer m_DeformationFieldRawXI;
    DeformationFieldTypeVector::Pointer m_DeformationFieldRawII;
    vtkSmartPointer<vtkImageData> deformationField;
    vtkSmartPointer<vtkOrientedGridTransform> grid_Ras;
};
