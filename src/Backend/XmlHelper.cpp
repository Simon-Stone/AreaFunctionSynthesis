// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2008, Peter Birkholz, Hamburg, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************

#include "XmlHelper.h"


// ****************************************************************************
// ****************************************************************************

XmlNode *XmlHelper::getChildNode(XmlNode *node, const char *childName, int index) throw (std::string)
{
  if ((node == NULL) || (childName == NULL))
  { 
    throw std::string("Invalid parameters for getChildNode(...).");
  }

  XmlNode *childNode = node->getChildElement(childName, index);

  if (childNode == NULL)
  {
    char st[512];
    sprintf(st, "The child element <%s> of the node <%s> at position %d"
      " does not exist!", childName, node->name.c_str(), index);
    throw std::string(st);
  }

  return childNode;
}


// ****************************************************************************
// ****************************************************************************

void XmlHelper::readAttribute(XmlNode *node, const char *attrName, double &attrValue) throw (std::string)
{
  if ((node == NULL) || (attrName == NULL))
  {
    throw std::string("Invalid parameters for readAttribute(...).");
  }

  if (node->hasAttribute(attrName) == false)
  {
    throw std::string("The attribute '") + attrName + "' for the element <" + 
      node->name + "> does not exist!"; 
  }

  attrValue = node->getAttributeDouble(attrName);
}

// ****************************************************************************
// ****************************************************************************

void XmlHelper::readAttribute(XmlNode *node, const char *attrName, int &attrValue) throw (std::string)
{
  if ((node == NULL) || (attrName == NULL))
  {
    throw std::string("Invalid parameters for readAttribute(...).");
  }

  if (node->hasAttribute(attrName) == false)
  {
    throw std::string("The attribute '") + attrName + "' for the element <" + 
      node->name + "> does not exist!"; 
  }

  attrValue = node->getAttributeInt(attrName);
}

// ****************************************************************************
// ****************************************************************************

void XmlHelper::readAttribute(XmlNode *node, const char *attrName, std::string &attrValue) throw (std::string)
{
  if ((node == NULL) || (attrName == NULL))
  {
    throw std::string("Invalid parameters for readAttribute(...).");
  }

  if (node->hasAttribute(attrName) == false)
  {
    throw std::string("The attribute '") + attrName + "' for the element <" + 
      node->name + "> does not exist!"; 
  }

  attrValue = node->getAttributeString(attrName);
}

// ****************************************************************************

