#include "ResampleVolume.h"

///	\fn			ResampleVolumesToBeIsotropic( regImage* Image )
///	\brief		This function resamples a voulme to be isotropic
///	\param		Image The name of the #regImage object
///	\return		This function returns an integer
///	It is unfortunate that it is still very common to find medical image
///	datasets that have been acquired with large inter-sclice spacings that
///	result in voxels with anisotropic shapes. In many cases these voxels have
///	ratios of [1:5] or even [1:10] between the resolution in the plane (x,y) and
///	the resolution along the z axis. Such datasets are close to useless for the
///	purpose of computer assisted image analysis.

int ResampleVolumesToBeIsotropic(regImage* Image)
{
    typedef regImage::ImageType ImageType;
    typedef float InternalPixelType;
    typedef itk::Image< InternalPixelType, Dimension > InternalImageType;
    typedef itk::RecursiveGaussianImageFilter< ImageType, InternalImageType >
            GaussianFilterType;
    GaussianFilterType::Pointer smootherX = GaussianFilterType::New();
    typedef itk::RecursiveGaussianImageFilter< InternalImageType, InternalImageType > GaussianFilterType1;
    GaussianFilterType1::Pointer smootherY = GaussianFilterType1::New();
        smootherX->SetInput( Image->GetImagePointer() );
        smootherY->SetInput( smootherX->GetOutput() );
    ImageType::Pointer inputImage = Image->GetImagePointer();
    const ImageType::SpacingType& inputSpacing = inputImage->GetSpacing();
    const double isoSpacing = sqrt( inputSpacing[2] * inputSpacing[0] );
        smootherX->SetSigma( isoSpacing );
        smootherY->SetSigma( isoSpacing );
        smootherX->SetDirection( 0 );
        smootherY->SetDirection( 1 );
        smootherX->SetNormalizeAcrossScale( true );
        smootherY->SetNormalizeAcrossScale( true );
    typedef itk::ResampleImageFilter< InternalImageType, ImageType > ResampleFilterType;
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    typedef itk::IdentityTransform< double, Dimension > TransformType;
    TransformType::Pointer transform = TransformType::New();
        transform->SetIdentity();
        resampler->SetTransform( transform );
    typedef itk::MinimumMaximumImageCalculator< ImageType > CalculatorType;
    CalculatorType::Pointer calculator = CalculatorType::New();
        calculator->SetImage( Image->GetImagePointer() );
    try
    {
        calculator->Compute();
    }
    catch ( itk::ExceptionObject& e )
    {
        if (details) std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }
    typedef itk::LinearInterpolateImageFunction< InternalImageType, double > InterpolatorType;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
        resampler->SetInterpolator( interpolator );
        resampler->SetDefaultPixelValue( calculator->GetMinimum() ); // highlight regions without source
    ImageType::SpacingType spacing;
        spacing[0] = isoSpacing;
        spacing[1] = isoSpacing;
        spacing[2] = isoSpacing;
    resampler->SetOutputSpacing( spacing );
    resampler->SetOutputOrigin( inputImage->GetOrigin() );
    ImageType::SizeType inputSize = inputImage->GetLargestPossibleRegion().GetSize();
    typedef ImageType::SizeType::SizeValueType SizeValueType;
    const double dx = inputSize[0] * inputSpacing[0] / isoSpacing;
    const double dy = inputSize[1] * inputSpacing[1] / isoSpacing;
    const double dz = (inputSize[2] - 1 ) * inputSpacing[2] / isoSpacing;
    ImageType::SizeType size;
        size[0] = static_cast < SizeValueType>( dx );
        size[1] = static_cast < SizeValueType>( dy );
        size[2] = static_cast < SizeValueType>( dz );
    resampler->SetSize( size );
    resampler->SetInput( smootherY->GetOutput() );
    try
    {
        std::cout << "Resample volume to be isotropic" << std::endl;
        resampler->Update();
    }
    catch ( itk::ExceptionObject& e )
    {
        if ( details ) std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }
    Image->SetImagePointer(resampler->GetOutput());
    return EXIT_SUCCESS;
}


int ResampleVolumes(regImage* Image)
{
    typedef regImage::ImageType ImageType;
    typedef float InternalPixelType;
    typedef itk::Image< InternalPixelType, Dimension > InternalImageType;

    typedef itk::ResampleImageFilter< InternalImageType, ImageType > ResampleFilterType;
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();

    typedef itk::IdentityTransform< double, Dimension > TransformType;
    TransformType::Pointer transform = TransformType::New();
        transform->SetIdentity();
        resampler->SetTransform( transform );

    typedef itk::MinimumMaximumImageCalculator< ImageType > CalculatorType;
    CalculatorType::Pointer calculator = CalculatorType::New();
        calculator->SetImage( Image->GetImagePointer() );

    try
    {
        calculator->Compute();
    }
    catch (itk::ExceptionObject& e)
    {
        if (details) std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::LinearInterpolateImageFunction< InternalImageType, double > InterpolatorType;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
        resampler->SetInterpolator( interpolator );
        resampler->SetDefaultPixelValue( calculator->GetMinimum() ); // highlight regions without source

    ImageType::SpacingType spacing;
        spacing[0] = 1;
        spacing[1] = 1;
        spacing[2] = 1;

    double oriGyn[3];
        oriGyn[0] = 0.0;
        oriGyn[1] = 0.0;
        oriGyn[2] = 0.0;
    resampler->SetOutputSpacing( spacing );
    resampler->SetOutputOrigin( oriGyn/*inputImage->GetOrigin()*/ );
    try
    {
        std::cout << "Resample volume" << std::endl;
        resampler->Update();
    }
    catch ( itk::ExceptionObject& e )
    {
        if ( details ) std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }
    Image->SetImagePointer(resampler->GetOutput());
    return EXIT_SUCCESS;
}
