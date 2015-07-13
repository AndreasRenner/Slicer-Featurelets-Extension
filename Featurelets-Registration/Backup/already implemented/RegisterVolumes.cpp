#include "RegisterVolumes.h"

///	\fn			RegisterVolumes( ImageType* FixedImage, ImageType* MovingImage )
///	\brief		This function registers two volumes
///	\param		FixedImage The name of the fixed #ImageType object
///	\param		MovingImage The name of the moving #ImageType object
///	\return		This function returns an integer
///	This function uses of the VersorRigid3DTransform class for performing
///	registration of two 3D images. The class CenteredTransformInitializer is
///	used to initialize the center and translation of the transform. The case of
///	rigid registration of 3D images is probably one of the most commonly found
///	cases of image registration.

FeatureletRegistrationResult::Pointer RegisterVolumes(ImageType* FixedImage, ImageType* MovingImage )
{
  FeatureletRegistrationResultPointer regResult;
  regResult = FeatureletRegistrationResultType::New();

  typedef itk::TranslationTransform< double > TransformType;
  TransformType::Pointer transform = TransformType::New();

  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();

  typedef itk::NormalizedCorrelationImageToImageMetric<ImageType, ImageType> MetricType;
  MetricType::Pointer metric = MetricType::New();

  typedef itk::LinearInterpolateImageFunction< ImageType, double > InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  typedef itk::ImageRegistrationMethod< ImageType, ImageType > RegistrationType;
  RegistrationType::Pointer registration = RegistrationType::New();

  optimizer->MinimizeOn();
  registration->SetTransform( transform );
  registration->SetMetric( metric );
  registration->SetOptimizer( optimizer );
  registration->SetInterpolator( interpolator );

  ImageType::SizeType size =
                    MovingImage->GetImagePointer()->GetLargestPossibleRegion().GetSize();
  ImageType::RegionType region;
  region.SetSize(size);
  registration->SetFixedImage( FixedImage->GetImagePointer() );
  registration->SetMovingImage( MovingImage->GetFeatureletPointer() );
  registration->SetFixedImageRegion( FixedImage->GetSearchRegionPointer()->GetBufferedRegion()); //******GetBufferedRegion()

  typedef RegistrationType::ParametersType ParametersType;
  ParametersType initialParameters( transform->GetNumberOfParameters() );
  initialParameters[0] = 0.0;  // Initial offset in mm along X
  initialParameters[1] = 0.0;  // Initial offset in mm along Y
  initialParameters[2] = 0.0; // Initial offset in mm along Z

  registration->SetInitialTransformParameters( initialParameters );

  optimizer->SetMaximumStepLength(0.05);
  optimizer->SetMinimumStepLength( 0.0001 );
  optimizer->SetNumberOfIterations( 2000 );

  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();

  optimizer->AddObserver( itk::IterationEvent(), observer );

  try
  {
    registration->Update();
  }
  catch( itk::ExceptionObject& e )
  {
    if ( details ) std::cerr << e << std::endl;
    regResult->SetStatusRegistration(EXIT_FAILURE);
    return regResult;  //Identificar este error como de regstration error!!!
  }

  OptimizerType::ParametersType finalParameters =
                registration->GetLastTransformParameters();
  transform->SetParameters( finalParameters );
  TransformType::OutputVectorType offset = transform->GetOffset();
  const double bestValue = optimizer->GetValue();
  std::cout << "  Offset1 = " << offset << std::endl;
  std::cout << "Metric Value= "<< bestValue << std::endl;

  if(true)
  {
    regResult->SetTransform(transform);
    regResult->SetOptimizer(optimizer);
  }
  return regResult;
}

FeatureletRegistrationResult::Pointer RegisterVolumes2(ImageType* FixedImage, ImageType* MovingImage)
{
  FeatureletRegistrationResultPointer regResult;
  regResult = FeatureletRegistrationResultType::New();

  typedef itk::TranslationTransform< double > TransformType;
  TransformType::Pointer transform = TransformType::New();

  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();

  typedef itk::NormalizedCorrelationImageToImageMetric< ImageType, ImageType > MetricType;
  MetricType::Pointer metric = MetricType::New();

  typedef itk::NearestNeighborInterpolateImageFunction< ImageType, double > InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  typedef itk::ImageRegistrationMethod< ImageType, ImageType > RegistrationType;
  RegistrationType::Pointer registration = RegistrationType::New();
  optimizer->MinimizeOn();
  registration->SetTransform( transform );
  registration->SetMetric( metric );
  registration->SetOptimizer( optimizer );
  registration->SetInterpolator( interpolator );

  ImageType::SizeType size = MovingImage->GetImagePointer()->GetLargestPossibleRegion().GetSize();
  ImageType::RegionType region;
  region.SetSize(size);
  registration->SetFixedImage( FixedImage->GetImagePointer() );
  registration->SetMovingImage( MovingImage->GetFeatureletPointer() );
  registration->SetFixedImageRegion( FixedImage->GetSearchRegionPointer()->GetBufferedRegion()); //******GetBufferedRegion()

  typedef RegistrationType::ParametersType ParametersType;
  ParametersType initialParameters( transform->GetNumberOfParameters() );
  initialParameters[0] = 0.0;  // Initial offset in mm along X
  initialParameters[1] = 0.0;  // Initial offset in mm along Y
  initialParameters[2] = 0.0; // Initial offset in mm along Z

  registration->SetInitialTransformParameters( initialParameters );
  optimizer->SetMaximumStepLength(0.05);
  optimizer->SetMinimumStepLength( 0.0001 );
  optimizer->SetNumberOfIterations( 2000 );

  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();

  optimizer->AddObserver( itk::IterationEvent(), observer );

  try
  {
    registration->Update();
  }
  catch( itk::ExceptionObject& e )
  {
    if ( details ) std::cerr << e << std::endl;
    regResult->SetStatusRegistration(EXIT_FAILURE);
    return regResult;  //Identificar este error como de regstration error!!!
  }
  OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();

  transform->SetParameters( finalParameters );
  TransformType::OutputVectorType offset = transform->GetOffset();
  const double bestValue = optimizer->GetValue();
  std::cout << "  Offset1 = " << offset << std::endl;
  std::cout << "Metric Value= "<< bestValue<<std::endl;
  if(true)
  {
    regResult->SetTransform(transform);
    regResult->SetOptimizer(optimizer);
  }
  return regResult;
}


FeatureletRegistrationResult::Pointer RegisterVolumesII2( ImageType* FixedImage,
                                                          ImageType* MovingImage ) {
    FeatureletRegistrationResultPointer regResult;
    regResult = FeatureletRegistrationResultType::New();
    typedef   float  InternalPixelType;
    typedef itk::Image< InternalPixelType, 3 > InternalImageType;

  ////////// Software Guide : BeginCodeSnippet /////////////

    typedef itk::TranslationTransform< double, 3 > TransformType;
    typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
    typedef itk::NearestNeighborInterpolateImageFunction <InternalImageType, double> InterpolatorType;
    typedef itk::ImageRegistrationMethod<InternalImageType, InternalImageType> RegistrationType;
    typedef itk::MutualInformationImageToImageMetric<InternalImageType, InternalImageType> MetricType;

    TransformType::Pointer      transform     = TransformType::New();
    OptimizerType::Pointer      optimizer     = OptimizerType::New();
    InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
    RegistrationType::Pointer   registration  = RegistrationType::New();
        registration->SetOptimizer(     optimizer     );
        registration->SetTransform(     transform     );
        registration->SetInterpolator(  interpolator  );
    MetricType::Pointer         metric        = MetricType::New();

        registration->SetMetric( metric  );

        metric->SetFixedImageStandardDeviation(  0.4 );
         metric->SetMovingImageStandardDeviation( 0.4 );

         // Software Guide : EndCodeSnippet

    typedef itk::NormalizeImageFilter<ImageType, InternalImageType> FixedNormalizeFilterType;
    typedef itk::NormalizeImageFilter<ImageType,InternalImageType> MovingNormalizeFilterType;
    FixedNormalizeFilterType::Pointer fixedNormalizer = FixedNormalizeFilterType::New();
    MovingNormalizeFilterType::Pointer movingNormalizer = MovingNormalizeFilterType::New();
    typedef itk::DiscreteGaussianImageFilter<InternalImageType, InternalImageType> GaussianFilterType;
    GaussianFilterType::Pointer fixedSmoother  = GaussianFilterType::New();
    GaussianFilterType::Pointer movingSmoother = GaussianFilterType::New();
    fixedSmoother->SetVariance( 2.0 );
    movingSmoother->SetVariance( 2.0 );
    fixedNormalizer->SetInput( FixedImage->GetFeatureletPointer() );
    movingNormalizer->SetInput( MovingImage->GetFeatureletPointer() );
    fixedSmoother->SetInput( fixedNormalizer->GetOutput() );
    movingSmoother->SetInput( movingNormalizer->GetOutput() );

        registration->SetFixedImage(    fixedSmoother->GetOutput()    );
        registration->SetMovingImage(   movingSmoother->GetOutput()   );



        fixedNormalizer->Update();
        ImageType::RegionType fixedImageRegion =
        fixedNormalizer->GetOutput()->GetBufferedRegion();
        registration->SetFixedImageRegion(
                    FixedImage->GetSearchRegionPointer()->GetBufferedRegion() );

    typedef RegistrationType::ParametersType ParametersType;
    ParametersType initialParameters( transform->GetNumberOfParameters() );
        initialParameters[0] = 0.0;  // Initial offset in mm along X
        initialParameters[1] = 0.0;  // Initial offset in mm along Y
        initialParameters[2] = 0.0;

        registration->SetInitialTransformParameters( initialParameters );

    const unsigned int numberOfPixels = fixedImageRegion.GetNumberOfPixels();
    const unsigned int numberOfSamples =
                        static_cast< unsigned int >( numberOfPixels * 0.08 );

        metric->SetNumberOfSpatialSamples( numberOfSamples );

        optimizer->SetNumberOfIterations( 2000 );
        optimizer->MaximizeOn();
        optimizer->SetMaximumStepLength(0.5);
        optimizer->SetMinimumStepLength( 0.001);

    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
        optimizer->AddObserver( itk::IterationEvent(), observer );

        try {
            registration->Update();
        }
        catch( itk::ExceptionObject& e ) {
            if ( details ) std::cerr << e << std::endl;
            regResult->SetStatusRegistration(EXIT_FAILURE);
            return regResult;
        }

    OptimizerType::ParametersType finalParameters =
              registration->GetLastTransformParameters();
        transform->SetParameters( finalParameters );
    TransformType::OutputVectorType offset = transform->GetOffset();
    const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
    const double bestValue = optimizer->GetValue();
    std::cout << "  Offset1 = " << offset << std::endl;
    std::cout << "Metric Value= "<< bestValue<<std::endl;
    if(true) {
        regResult->SetTransform(transform);
        regResult->SetOptimizer(optimizer);
    }
    return regResult;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////MODIFYING &playing around!!! //////////////////////////////////////
FeatureletRegistrationResult::Pointer RegisterVolumesII(
        ImageType* FixedImage, ImageType* MovingImage ) {

    FeatureletRegistrationResultPointer regResult;
    regResult = FeatureletRegistrationResultType::New();
    typedef   float                                    InternalPixelType;
    typedef itk::Image< InternalPixelType, 3 > InternalImageType;

  ////////// Software Guide : BeginCodeSnippet  ////////////////////
    typedef itk::TranslationTransform< double, 3 > TransformType;
    typedef itk::RegularStepGradientDescentOptimizer  OptimizerType;
    typedef itk::Function::CosineWindowFunction<4, double, double> cosine;
    typedef itk::ConstantBoundaryCondition< InternalImageType > boundaries;
    typedef itk::WindowedSincInterpolateImageFunction< InternalImageType, 4,
            cosine, boundaries, double > InterpolatorType;
    typedef itk::ImageRegistrationMethod<
                                    InternalImageType,
                                    InternalImageType >  RegistrationType;
    typedef itk::MutualInformationImageToImageMetric<
                                          InternalImageType,
                                          InternalImageType >    MetricType;
    TransformType::Pointer      transform     = TransformType::New();
    OptimizerType::Pointer      optimizer     = OptimizerType::New();
    InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
    RegistrationType::Pointer   registration  = RegistrationType::New();

        registration->SetOptimizer(     optimizer     );
        registration->SetTransform(     transform     );
        registration->SetInterpolator(  interpolator  );
    MetricType::Pointer         metric        = MetricType::New();

        registration->SetMetric( metric  );
    typedef itk::NormalizeImageFilter<
                                ImageType,
                                InternalImageType
                                        > FixedNormalizeFilterType;
    typedef itk::NormalizeImageFilter<
                                ImageType,
                                InternalImageType
                                              > MovingNormalizeFilterType;
    FixedNormalizeFilterType::Pointer fixedNormalizer =
                                            FixedNormalizeFilterType::New();
    MovingNormalizeFilterType::Pointer movingNormalizer =
                                            MovingNormalizeFilterType::New();
    typedef itk::DiscreteGaussianImageFilter<
                                      InternalImageType,
                                      InternalImageType
                                                    > GaussianFilterType;
    GaussianFilterType::Pointer fixedSmoother  = GaussianFilterType::New();
    GaussianFilterType::Pointer movingSmoother = GaussianFilterType::New();
        fixedSmoother->SetVariance( 2.0 );
        movingSmoother->SetVariance( 2.0 );
        fixedNormalizer->SetInput( FixedImage->GetFeatureletPointer() );
        movingNormalizer->SetInput( MovingImage->GetFeatureletPointer() );
        fixedSmoother->SetInput( fixedNormalizer->GetOutput() );
        movingSmoother->SetInput( movingNormalizer->GetOutput() );
        registration->SetFixedImage( fixedSmoother->GetOutput() );
        registration->SetMovingImage( movingSmoother->GetOutput() );
        fixedNormalizer->Update();

    ImageType::RegionType fixedImageRegion =
       fixedNormalizer->GetOutput()->GetBufferedRegion();

        registration->SetFixedImageRegion(
                    FixedImage->GetSearchRegionPointer()->GetBufferedRegion() );

    typedef RegistrationType::ParametersType ParametersType;
    ParametersType initialParameters( transform->GetNumberOfParameters() );

        initialParameters[0] = 0.0;  // Initial offset in mm along X
        initialParameters[1] = 0.0;  // Initial offset in mm along Y
        initialParameters[2] = 0.0;
        registration->SetInitialTransformParameters( initialParameters );

    const unsigned int numberOfPixels = fixedImageRegion.GetNumberOfPixels();
    const unsigned int numberOfSamples =
                        static_cast< unsigned int >( numberOfPixels * 0.08 );

        optimizer->SetNumberOfIterations( 2000 );
        optimizer->MaximizeOn();
        optimizer->SetMaximumStepLength(0.5);
        optimizer->SetMinimumStepLength( 0.001);
    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();

        optimizer->AddObserver( itk::IterationEvent(), observer );
  try {
            registration->Update();
        }
        catch( itk::ExceptionObject& e ) {
            if ( details ) std::cerr << e << std::endl;
            regResult->SetStatusRegistration(EXIT_FAILURE);
            return regResult;
        }

    OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();

        transform->SetParameters( finalParameters );

    TransformType::OutputVectorType offset = transform->GetOffset();
    const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
    const double bestValue = optimizer->GetValue();

    std::cout << "  Offset1 = " << offset << std::endl;
    std::cout << "Metric Value= "<< bestValue<<std::endl;

    if(true) {
        regResult->SetTransform(transform);
        regResult->SetOptimizer(optimizer);
    }
    return regResult;
  }
//////////////////////////////////////////////////////////////////////////////////
//////////////////////END of MODIFYING &playing around!!! //////////////////////////////////////
