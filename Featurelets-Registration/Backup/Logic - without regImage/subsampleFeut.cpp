#include "subsampleFeut.h"

///	\fn			SubsampleVolume( regImage* Image, ImageType::SizeType FeatureletSize , ImageType::IndexType ImageIndex )
///	\brief		This function subsamples a volume
///	\param		Image The name of the #regImage object
///	\param		FeatureletSize The size of one featurelet
///	\param		ImageIndex The current index of the image
///	\return		This function returns an integer
///	This function performs subsampling of a volume using ITK classes. In order
///	to avoid aliasing artifacts, the volume must be processed by a low-pass
///	filter before resampling.

int SubsampleVolume( ImageType::Pointer Image,
                     ImageType::SizeType FeatureletSize,
                     ImageType::SizeType SearchRegionSize,
                     ImageType::IndexType ImageIndex)
{
  ImageType::RegionType desiredRegion;
        desiredRegion.SetSize( FeatureletSize );
        desiredRegion.SetIndex( ImageIndex );
  ImageType::RegionType SearchingRegion;
        SearchingRegion.SetSize( SearchRegionSize );
        SearchingRegion.SetIndex( ImageIndex );

  typedef itk::ExtractImageFilter< ImageType, ImageType > FilterType;
  FilterType::Pointer filter = FilterType::New();
  FilterType::Pointer filter2 = FilterType::New();
        filter->SetExtractionRegion( desiredRegion );
        filter->SetInput( Image->GetImagePointer() );
        filter2->SetExtractionRegion(SearchingRegion);
        filter2->SetInput(Image->GetImagePointer());
  try
  {
        filter->Update();
        filter2->Update();
  }
  catch ( itk::ExceptionObject& e )
  {
        if (details) std::cerr << e << std::endl;
        return EXIT_FAILURE;
  }
  Image->SetFeatureletPointer( filter->GetOutput() );
  Image->SetSearchRegionPointer(filter2->GetOutput());
  return EXIT_SUCCESS;
}

///	\fn			CheckFeaturelet( regImage* Image )
///	\brief		This function checks a featurelet if it is necessary to be registered
///	\param		Image The name of the #regImage object
///	\return		This function returns an integer
///	This function checks a featurelet if it is necessary to be registered to
///	avoid registration errors or a bad result.

Featurelet::Status CheckFeaturelet( ImageType::Pointer Image )
{
  ImageType::SizeType FeatureletSize=
            Image->GetLargestPossibleRegion().GetSize();

  typedef itk::NormalizedCorrelationImageToImageMetric<ImageType, ImageType> MetricType;
    MetricType::Pointer metric = MetricType::New();
    if ( FeatureletSize[0] <= 4 )
        return Featurelet:: SizeError;
    if ( FeatureletSize[1] <= 4 )
        return  Featurelet:: SizeError;
    if ( FeatureletSize[2] <= 4 )
        return Featurelet:: SizeError;

    typedef itk::MinimumMaximumImageCalculator< ImageType > CalculatorType;
    CalculatorType::Pointer calculator = CalculatorType::New();
        calculator->SetImage( Image->GetFeatureletPointer() );
    try {
        calculator->Compute();
    }
    catch ( itk::ExceptionObject& e ) {
        if (details) std::cerr << e << std::endl;
        return Featurelet::internalError;
        }
    if ( calculator->GetMaximum() == calculator->GetMinimum() ){
        return Featurelet::sameColor;
    }
    return Featurelet:: OK;
}
