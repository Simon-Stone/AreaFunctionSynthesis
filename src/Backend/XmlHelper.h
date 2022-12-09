// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2008, Peter Birkholz, Hamburg, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************

#ifndef __XML_HELPER_H__
#define __XML_HELPER_H__

#include <string>

#include "XmlNode.h"

// Disable compiler warnings in Visual Studio, because the compiler
// does not yet support the throw() statement behind function declarations.
#ifdef WIN32
#pragma warning( disable : 4290 )
#endif

// ****************************************************************************
/// Static xml support class for the xml parser (XmlNode).
// ****************************************************************************

class XmlHelper
{
public:
  static XmlNode *getChildNode(XmlNode *node, const char *childName, int index = 0) throw (std::string);

  static void readAttribute(XmlNode *node, const char *attrName, double &attrValue) throw (std::string);
  static void readAttribute(XmlNode *node, const char *attrName, int &attrValue) throw (std::string);
  static void readAttribute(XmlNode *node, const char *attrName, std::string &attrValue) throw (std::string);
};

// ****************************************************************************

#endif
