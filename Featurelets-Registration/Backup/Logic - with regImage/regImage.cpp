/*=========================================================================
Program:	registration
Module:		regImage.h
Language:	C++
Created on:	April 23rd, 2009

Author:		Amon Bhatia

Written at the University of Applied Sciences Technikum-Vienna
Bachelor degree course Biomedical Engineering
=========================================================================*/

#include "stdio.h"
#include <string>
#include <vector>
#include <iostream>
// include stl header files
#include <map>

// include standard itk header files
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "regImage.h"


typedef short PixelType;

///	\typedef	map
///	\brief		Definition of a #map as associative array with strings
typedef std::map<std::string, std::string> map;

///	\class		regImage
///	\brief		This class is the core image class
///
///	This class is the core image class, which contains of a map and two
///	#ImageType pointer, one for the origin and a temporary one for caching
///	featurelet volumes
regImage::ImageType::Pointer regImage::GetImagePointer() {
    return Image;
}

///	\fn			regImage::regImage()
///	\brief		This function is the default constructor
regImage::regImage() {
}

///	\fn			regImage::~regImage()
///	\brief		This function is the default destructor
regImage::~regImage() {
}

///	\fn			regImage::GetImagePointer()
///	\brief		This function returns the current pointer of an image
///	\return		This function returns the #ImageType pointer of the #Image


///	\fn			regImage::SetImagePointer( ImageType::Pointer m_Image )
///	\brief		This function sets the current pointer of an image
///	\param		m_Image The #ImageType pointer of the #Image
void regImage::SetImagePointer( regImage::ImageType::Pointer m_Image ) {
    Image = m_Image;
}



///	\fn			regImage::GetFeatureletPointer()
///	\brief		This function returns the current pointer of a featurelet
///	\return		This function returns the #ImageType pointer of the #Image
regImage::ImageType::Pointer regImage::GetFeatureletPointer() {
    return Featurelet;
}

///	\fn			regImage::SetFeatureletPointer( ImageType::Pointer m_Featurelet )
///	\brief		This function sets the current pointer of a featurelet
///	\param		m_Featurelet The #ImageType pointer of the #Image
void regImage::SetFeatureletPointer( regImage::ImageType::Pointer m_Featurelet ) {
    Featurelet = m_Featurelet;
}


regImage::ImageType::Pointer regImage::GetSearchRegionPointer() {
    return SearchRegion;
}

///	\fn			regImage::SetFeatureletPointer( ImageType::Pointer m_Featurelet )
///	\brief		This function sets the current pointer of a featurelet
///	\param		m_Featurelet The #ImageType pointer of the #Image
void regImage::SetSearchRegionPointer( regImage::ImageType::Pointer m_SearchRegion ) {
    SearchRegion = m_SearchRegion;
}




///	\fn			regImage::GetDicomHeaderTag( std::string m_DicomHeaderTagKey )
///	\brief		This function returns the value of the DicomHeaderMap
///	\param		m_DicomHeaderTagKey The key of the DICOM header tag
///	\return		This function returns a string
std::string regImage::GetDicomHeaderTag( std::string m_DicomHeaderTagKey ) {
    return DicomHeaderMap[m_DicomHeaderTagKey];
}

///	\fn			regImage::SetDicomHeaderTag( std::string m_DicomHeaderTagKey, std::string m_DicomHeaderTag )
///	\brief		This function sets the value of the DicomHeaderMap
///	\param		m_DicomHeaderTagKey The key of the DICOM header tag
///	\param		m_DicomHeaderTag The tag of the DICOM header
void regImage::SetDicomHeaderTag( std::string m_DicomHeaderTagKey, std::string m_DicomHeaderTag ) {
    DicomHeaderMap[m_DicomHeaderTagKey] = m_DicomHeaderTag;
}



