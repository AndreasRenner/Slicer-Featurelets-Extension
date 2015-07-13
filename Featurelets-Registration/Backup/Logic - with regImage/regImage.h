/*=========================================================================
Program:	registration
Module:		regImage.h
Language:	C++
Created on:	April 23rd, 2009

Author:		Amon Bhatia

Written at the University of Applied Sciences Technikum-Vienna
Bachelor degree course Biomedical Engineering
=========================================================================*/
#ifndef REGIMAGE_H_
#define REGIMAGE_H_

// include standard header files
#include "stdio.h"
#include <string>
#include <vector>
#include <iostream>

// include stl header files
#include <map>

// include standard itk header files
#include "itkImage.h"
#include "itkImageFileReader.h"

///	\typedef	PixelType
///	\brief		Definition of the dicom data type for #ImageType
typedef short PixelType;

///	\var		Dimension
///	\brief		Definition of three dimensions for #ImageType
const unsigned int Dimension = 3;

///	\typedef	map
///	\brief		Definition of a #map as associative array with strings
typedef std::map<std::string, std::string> map;

///	\class		regImage
///	\brief		This class is the core image class
///
///	This class is the core image class, which contains of a map and two
///	#ImageType pointer, one for the origin and a temporary one for caching
///	featurelet volumes

class regImage
{
	public:
	///	\typedef	ImageType
	///	\brief		Definition of the #ImageType with #PixelType and #Dimension
	typedef itk::Image< PixelType, Dimension > ImageType;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	typedef	ReaderType::Pointer ReaderPointer;

    regImage();
    ~regImage();
    ImageType::Pointer GetImagePointer();
    void SetImagePointer( ImageType::Pointer m_Image );
    ImageType::Pointer GetFeatureletPointer();
    void SetFeatureletPointer( ImageType::Pointer m_Image );
    ImageType::Pointer GetSearchRegionPointer();
    void SetSearchRegionPointer( ImageType::Pointer m_Image );
    std::string GetDicomHeaderTag( std::string m_DicomHeaderTagKey );
    void SetDicomHeaderTag( std::string m_DicomHeaderTagKey, std::string m_DicomHeaderTag );
    itkGetObjectMacro(Reader, ReaderType);

	private:
		ReaderPointer m_Reader;
		///	\var		DicomHeaderMap
		///	\brief		This #map contains all the DICOM header information
		map DicomHeaderMap;
		///	\var		Image
		///	\brief		This pointer contains the original volumetric image
		ImageType::Pointer Image;

		///	\var		Featurelet
		///	\brief		This pointer contains the current featurelet
		ImageType::Pointer Featurelet;
		ImageType::Pointer SearchRegion;
		///	\var		FeatureletSize
		///	\brief		This array contains the calculated size of one featurelet
		ImageType::SizeType FeatureletSize;
		ImageType::SizeType SearchRegionSize;

};


#endif /* REGIMAGE_H_ */
