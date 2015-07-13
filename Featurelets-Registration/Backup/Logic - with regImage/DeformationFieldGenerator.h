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
#include "FeatureletRegistrationResult.h"
#include "regImage.h"

#include "itkObject.h"
#include "itkWarpImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkVectorResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkIdentityTransform.h"

#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include <itkNearestNeighborInterpolateImageFunction.h>


/** \class NumericTraits< std::complex<double> >
 *  \brief Define traits for type std::complex<double>.
 *  \ingroup DataRepresentation
 */


class DeformationFieldGenerator {
public:
	typedef std::list<FeatureletRegistrationResult::Pointer> FeatureletRegistrationResultListType;
    static const unsigned int VectorDimension = 1;

    typedef regImage::ImageType ImageType;

	typedef double PixelType;
    typedef itk::Image < PixelType, Dimension > DeformationFieldType;
    typedef itk::Vector< double, Dimension > PixelTypeVector;
    typedef itk::Image< PixelTypeVector, Dimension > DeformationFieldTypeVector;

	DeformationFieldGenerator() {}
	virtual ~DeformationFieldGenerator() {}

    void AddFeaturelet(FeatureletRegistrationResult *s)
    {
      FeatureletRegistrationResult::Pointer regResult;
      regResult = s;
      m_RegistrationResultList.push_back(regResult);
	}

    void OpenStream(const char *_fileName = "regResults.txt")
    {
      ofs.open(_fileName);
	}

    void SaveAll()
    {
	}

    DeformationFieldType::Pointer GenerateDeformationField()
    {
    FeatureletRegistrationResultListType::iterator it;
    DeformationFieldType::IndexType start;

    start[0]=0;
    start[1]=0;
    start[2]=0;

    DeformationFieldType::SizeType size;

    size[0]=m_NumberOfFeaturelets[0]+1;//final FeatureletGridPosition
    size[1]=m_NumberOfFeaturelets[1]+1;
    size[2]=m_NumberOfFeaturelets[2]+1;

    DeformationFieldType::RegionType region0; // a cambiar??????????????????????????????

    region0.SetSize(size);
    region0.SetIndex(start);

    double spacing1[Dimension];
    spacing1[0]=m_SizeOfFeaturelets[0];//15;
    spacing1[1]=m_SizeOfFeaturelets[1];//15;
    spacing1[2]=m_SizeOfFeaturelets[2];//15;

    m_DeformationFieldRawX = DeformationFieldType::New();
    m_DeformationFieldRawX->SetRegions(region0);
    m_DeformationFieldRawX->Allocate();
    m_DeformationFieldRawX->SetSpacing(spacing1);

    PixelType zeroVector;
    zeroVector=0;
    m_DeformationFieldRawX->FillBuffer(zeroVector);
    m_DeformationFieldRawX->Print(std::cout);

    // the other two components:
    m_DeformationFieldRawY = DeformationFieldType::New();
    m_DeformationFieldRawY->SetRegions(region0);
    m_DeformationFieldRawY->Allocate();
    m_DeformationFieldRawY->SetSpacing(spacing1);
    m_DeformationFieldRawY->FillBuffer(zeroVector);
    m_DeformationFieldRawZ = DeformationFieldType::New();
    m_DeformationFieldRawZ->SetRegions(region0);
    m_DeformationFieldRawZ->SetSpacing(spacing1);
    m_DeformationFieldRawZ->Allocate();
    m_DeformationFieldRawZ->FillBuffer(zeroVector);

    m_MetricValue = DeformationFieldType::New();
    m_MetricValue->SetRegions(region0);
    m_MetricValue->SetSpacing(spacing1);
    m_MetricValue->Allocate();
    m_MetricValue->FillBuffer(zeroVector);

    //fill the image
	double bestValue;
    metricthreshold = 1000; ////////////////////////////////////////////////////////////////////////////

	PixelType transformVector;

    PixelType transformVector1;
    PixelType transformVector2;

    int underthreshold=0;
    int	abovethreshold=0;

    for (it = m_RegistrationResultList.begin(); it != m_RegistrationResultList.end(); it++)
    {
      DeformationFieldType::IndexType pixelIndex ;
      FeatureletRegistrationResult::TransformType::OutputVectorType offset;
      offset[0]= (*it)->GetTransform()->GetOffset()[0];
      offset[1]= (*it)->GetTransform()->GetOffset()[1];
      offset[2]= (*it)->GetTransform()->GetOffset()[2];
      pixelIndex[0]= (*it)->GetFeatureletGridPosition()[0];
      pixelIndex[1]= (*it)->GetFeatureletGridPosition()[1];
      pixelIndex[2]= (*it)->GetFeatureletGridPosition()[2];

      transformVector = offset[0];
      transformVector1 = offset[1];
      transformVector2 = offset[2];

      bestValue=(*it)->GetOptimizer()->GetValue();
      std::cout <<"best metric value "<< bestValue << std::endl;
      if (bestValue < metricthreshold)
      {
        abovethreshold++;
        m_DeformationFieldRawX->SetPixel(pixelIndex, transformVector);
        m_DeformationFieldRawY->SetPixel(pixelIndex, transformVector1);
        m_DeformationFieldRawZ->SetPixel(pixelIndex, transformVector2);
        std::cout<<"pixel "<<pixelIndex<<" vector: ["<<transformVector<<", "<<transformVector1<<", "<<transformVector2<< "] "<<std::endl;
      }
      else
      {
        underthreshold++;
        m_DeformationFieldRawX->SetPixel(pixelIndex, 0);
        m_DeformationFieldRawY->SetPixel(pixelIndex, 0);
        m_DeformationFieldRawZ->SetPixel(pixelIndex, 0);
      }
    }
    return m_DeformationFieldRawX, m_DeformationFieldRawY, m_DeformationFieldRawZ;
    }

    void RescaleDeformationField()
    {
	typedef itk::ResampleImageFilter<DeformationFieldType,DeformationFieldType> FilterType;
	FilterType::Pointer filter= FilterType::New();
	FilterType::Pointer filter3= FilterType::New();
	FilterType::Pointer filter2= FilterType::New();

	typedef itk::IdentityTransform< double, Dimension > TransformType;
	TransformType::Pointer transform = TransformType::New();

    typedef itk::LinearInterpolateImageFunction< DeformationFieldType, double > InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();

    DeformationFieldType::IndexType startII;
    startII[0]= 0;
    startII[1]= 0;
    startII[2]= 0;

    DeformationFieldType::SizeType sizeII;
    sizeII[0]=m_ImageSizeMov[0];
    sizeII[1]=m_ImageSizeMov[1];
    sizeII[2]=m_ImageSizeMov[2];

    DeformationFieldType::SizeType size;
    size[0]=m_NumberOfFeaturelets[0]+1; //+1
    size[1]=m_NumberOfFeaturelets[1]+1; //+1
    size[2]=m_NumberOfFeaturelets[2]+1; //+1

    m_DeformationFieldRawXI= DeformationFieldType::New();
    m_DeformationFieldRawYI= DeformationFieldType::New();
    m_DeformationFieldRawZI= DeformationFieldType::New();

    DeformationFieldType::RegionType regionII;
    regionII.SetSize(sizeII);
    regionII.SetIndex(startII);

    m_DeformationFieldRawXI->SetRegions(regionII);
    m_DeformationFieldRawYI->SetRegions(regionII);
    m_DeformationFieldRawZI->SetRegions(regionII);

    m_DeformationFieldRawXI->Allocate();
    m_DeformationFieldRawYI->Allocate();
    m_DeformationFieldRawZI->Allocate();

    PixelType zeroVector;
    zeroVector=0;
    m_DeformationFieldRawXI->FillBuffer(zeroVector);
    m_DeformationFieldRawYI->FillBuffer(zeroVector);
    m_DeformationFieldRawZI->FillBuffer(zeroVector);

    filter->SetTransform(transform);
    filter->SetInterpolator(interpolator);
    filter2->SetTransform(transform);
    filter2->SetInterpolator(interpolator);
    filter3->SetTransform(transform);
    filter3->SetInterpolator(interpolator);

    double origin[ Dimension ];
    origin[0] = 0.0;  // X space coordinate of origin
    origin[1] = 0.0;  // Y space coordinate of origin
    origin[2] = 0.0;  // Z space coordinate of origin

    double spacing1[Dimension];
    spacing1[0]=m_SizeOfFeaturelets[0];
    spacing1[1]=m_SizeOfFeaturelets[1];
    spacing1[2]=m_SizeOfFeaturelets[2];

    DeformationFieldType::SpacingType spacing;
	spacing[0] = 1;
	spacing[1] = 1;
	spacing[2] = 1;

    double o2[Dimension];
    for( unsigned int i=0; i<3; i++)
    {
      double center= origin[i] + (size[i] * spacing1[i]) / 2.0;
      o2[i] = center - sizeII[i] * spacing[i] / 2.0;
    }

    filter->SetOutputOrigin(o2);
    filter->SetSize( sizeII );
    filter->SetOutputSpacing(spacing);
    filter2->SetOutputOrigin(o2);
    filter2->SetSize( sizeII );
    filter2->SetOutputSpacing(spacing);
    filter3->SetOutputOrigin(o2);
    filter3->SetSize( sizeII );
    filter3->SetOutputSpacing(spacing);

    DeformationFieldType::DirectionType direction;
    direction.SetIdentity();

    filter->SetOutputDirection( direction );
    filter2->SetOutputDirection( direction );
    filter3->SetOutputDirection( direction );

    filter->SetInput(m_DeformationFieldRawX);
    filter2->SetInput(m_DeformationFieldRawY);
    filter3->SetInput(m_DeformationFieldRawZ);

    filter->Update();
    filter->UpdateLargestPossibleRegion();
    filter2->Update();
    filter2->UpdateLargestPossibleRegion();
    filter3->Update();
    filter3->UpdateLargestPossibleRegion();

    m_DeformationFieldRawXI=filter->GetOutput();
    m_DeformationFieldRawYI=filter2->GetOutput();
    m_DeformationFieldRawZI=filter3->GetOutput();

    m_DeformationFieldRawXI->Update();
    m_DeformationFieldRawYI->Update();
    m_DeformationFieldRawZI->Update();
    }

    int ResampleVolumesToBeIsotropic2(DeformationFieldType::Pointer Image)
    {
      typedef itk::ResampleImageFilter<DeformationFieldType, DeformationFieldType> ResampleFilterType;
      ResampleFilterType::Pointer resampler2 = ResampleFilterType::New();

      typedef itk::IdentityTransform< double, Dimension > TransformType;
      TransformType::Pointer transform = TransformType::New();
      transform->SetIdentity();
      resampler2->SetTransform( transform );

      typedef itk::LinearInterpolateImageFunction< DeformationFieldType, double > InterpolatorType;
      InterpolatorType::Pointer interpolator = InterpolatorType::New();
      resampler2->SetInterpolator( interpolator );

      std::cout<<"error1"<<std::endl;
      DeformationFieldType::SpacingType spacing;
      spacing[0] = (m_NumberOfFeaturelets[0]+1)*(m_SizeOfFeaturelets[0])/(m_ImageSizeMov[0]);
      spacing[0] = (m_NumberOfFeaturelets[1]+1)*(m_SizeOfFeaturelets[1])/(m_ImageSizeMov[1]);
      spacing[0] = (m_NumberOfFeaturelets[2]+1)*(m_SizeOfFeaturelets[2])/(m_ImageSizeMov[2]);
      std::cout<<"error2"<<std::endl;
      resampler2->SetOutputSpacing( spacing );

      typedef DeformationFieldType::SizeType::SizeValueType SizeValueType;

      const double dx2 = m_ImageSizeMov[0];
      const double dy2 = m_ImageSizeMov[1];
      const double dz2 = m_ImageSizeMov[2];

      DeformationFieldType::SizeType size2;
      size2[0] = static_cast < SizeValueType>( dx2 );
      size2[1] = static_cast < SizeValueType>( dy2 );
      size2[2] = static_cast < SizeValueType>( dz2 );
      std::cout<<"error3"<<std::endl;
      resampler2->SetSize( size2 );
      resampler2->SetInput(Image);
      std::cout<<"error4"<<std::endl;
      resampler2->UpdateLargestPossibleRegion();
      std::cout<<"error4"<<std::endl;
      Image=resampler2->GetOutput() ;
      Image->Update();
      return EXIT_SUCCESS;
	}

    void RescaleDeformationField2()
    {
      ResampleVolumesToBeIsotropic2(m_DeformationFieldRawXI);
      ResampleVolumesToBeIsotropic2(m_DeformationFieldRawYI);
      ResampleVolumesToBeIsotropic2(m_DeformationFieldRawZI);
    }

    void TotaldeformationFieldCreator()
    {
      DeformationFieldTypeVector::IndexType startII;
      startII[0]= 0;
      startII[1]= 0;
      startII[2]= 0;

      DeformationFieldTypeVector::SizeType sizeII;
      sizeII[0]=m_ImageSizeMov[0];
      sizeII[1]=m_ImageSizeMov[1];
      sizeII[2]=m_ImageSizeMov[2];

      DeformationFieldTypeVector::RegionType regionII;
      regionII.SetSize(sizeII);
      regionII.SetIndex(startII);

      m_DeformationFieldRawII= DeformationFieldTypeVector::New();
      m_DeformationFieldRawII->SetRegions(regionII);
      m_DeformationFieldRawII->Allocate();

      PixelTypeVector zeroVector;
      zeroVector.Fill(0);
      m_DeformationFieldRawII->FillBuffer(zeroVector);

      for (int i=0; i< m_ImageSizeMov[0]; i++)
      {
        for (int j=0; j< m_ImageSizeMov[1]; j++)
        {
          for (int k=0; k< m_ImageSizeMov[2]; k++)
          {
            PixelTypeVector TotalDeformationVector;

            DeformationFieldTypeVector::IndexType pixelIndexII;
            pixelIndexII[0]=i;
            pixelIndexII[1]=j;
            pixelIndexII[2]=k;

            TotalDeformationVector[0]=m_DeformationFieldRawXI->GetPixel(pixelIndexII);
            TotalDeformationVector[1]=m_DeformationFieldRawYI->GetPixel(pixelIndexII);
            TotalDeformationVector[2]=m_DeformationFieldRawZI->GetPixel(pixelIndexII);

            m_DeformationFieldRawII->SetPixel(pixelIndexII,TotalDeformationVector);
          }
        }
      }
      m_DeformationFieldRawII->Print(std::cout);
	}

    int WriteImageToFile(std::string outputFilename)
    {
      typedef itk::ImageFileWriter< DeformationFieldTypeVector> FieldWriterType;
      FieldWriterType::Pointer fwriter = FieldWriterType::New();
      fwriter->SetFileName( outputFilename );////por definir
      fwriter->SetInput(m_DeformationFieldRawII);
      fwriter->UpdateLargestPossibleRegion();

      try
      {
        fwriter->Update();
      }
      catch ( itk::ExceptionObject& err )
      {
        std::cerr << "ExceptionObject caught !" <<std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
      }
      return EXIT_SUCCESS;
    }

    int WriteImageToFile2(ImageType*InputImage,std::string outputFilename2)
    {
      typedef itk::ImageFileWriter<ImageType> FieldWriterType;
      FieldWriterType::Pointer fwriter = FieldWriterType::New();
      fwriter->SetFileName( outputFilename2 );////por definir
      fwriter->SetInput(InputImage);
      fwriter->UpdateLargestPossibleRegion();

      try
      {
        fwriter->Update();
      }
      catch ( itk::ExceptionObject& err )
      {
        std::cerr << "ExceptionObject caught !" <<std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
      }
      return EXIT_SUCCESS;
    }

    int WarpImagebyDeformationField(ImageType* InputImage, std::string outputFilename)
    {
      typedef itk::WarpImageFilter<ImageType,ImageType, DeformationFieldTypeVector > WarperType;
      typedef itk::LinearInterpolateImageFunction<ImageType,double> InterpolatorType;
      typedef  signed short  OutputPixelType;
      typedef itk::Image< OutputPixelType, 3 > OutputImageType;
      typedef itk::ImageFileWriter< OutputImageType >  WriterType;
      WriterType::Pointer      writer =  WriterType::New();
      WarperType::Pointer warper = WarperType::New();
      InterpolatorType::Pointer interpolator = InterpolatorType::New();

      DeformationFieldTypeVector::SizeType sizeII;

      sizeII[0]=m_ImageSizeMov[0];
      sizeII[1]=m_ImageSizeMov[1];
      sizeII[2]=m_ImageSizeMov[2];

      warper->SetOutputSize(sizeII);
      warper->SetInput(InputImage);
      warper->SetInterpolator( interpolator );
      warper->SetOutputSpacing( InputImage->GetSpacing() );
      warper->SetOutputOrigin( InputImage->GetOrigin() );
      warper->SetOutputDirection( InputImage->GetDirection() ); /////////////////////
      warper->SetDisplacementField(m_DeformationFieldRawII);///por generar a partir del anterior y reescalando
      warper->UpdateLargestPossibleRegion();
      warper->Update();
      writer->SetFileName(outputFilename); //
      writer->SetInput( warper->GetOutput()); //
      writer->UpdateLargestPossibleRegion(); //
      writer->Update(); //
      return EXIT_SUCCESS; //
      std::cout<<"Metric Threshold"<<metricthreshold<<std::endl;
    }

    void SetNumberOfFeaturelets(const int *_NumberOfFeaturelets)
    {
      m_NumberOfFeaturelets[0] = _NumberOfFeaturelets[0];
      m_NumberOfFeaturelets[1] = _NumberOfFeaturelets[1];
      m_NumberOfFeaturelets[2] = _NumberOfFeaturelets[2];
	}

    void SetSizeOfFeaturelets(const int *_SizeOfFeaturelets)
    {
      m_SizeOfFeaturelets[0] = _SizeOfFeaturelets[0];
      m_SizeOfFeaturelets[1] = _SizeOfFeaturelets[1];
      m_SizeOfFeaturelets[2] = _SizeOfFeaturelets[2];
	}

    void SetImageSizeMov(const double* _ImageSizeMov)
    {
      m_ImageSizeMov[0] = _ImageSizeMov[0];
      m_ImageSizeMov[1] = _ImageSizeMov[1];
      m_ImageSizeMov[2] = _ImageSizeMov[2];
	}

public:
	FeatureletRegistrationResultListType m_RegistrationResultList;
	std::ofstream ofs;
	int m_NumberOfFeaturelets[3];
	int m_SizeOfFeaturelets[3];
	int m_ImageSizeMov[3];
	double metricthreshold;
	DeformationFieldType::Pointer m_MetricValue;
	DeformationFieldType::Pointer m_DeformationFieldRawZ;
	DeformationFieldType::Pointer m_DeformationFieldRawY;
	DeformationFieldType::Pointer m_DeformationFieldRawX;
	DeformationFieldType::Pointer m_DeformationFieldRawZI;
	DeformationFieldType::Pointer m_DeformationFieldRawYI;
	DeformationFieldType::Pointer m_DeformationFieldRawXI;
    DeformationFieldTypeVector::Pointer m_DeformationFieldRawII;
};
