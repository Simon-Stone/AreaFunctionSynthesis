// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2010, Peter Birkholz, Hamburg, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************

#ifndef __XML_NODE_H__
#define __XML_NODE_H__

#include <string>
#include <vector>

using namespace std;

// ****************************************************************************
// This file implements a simple XML parser using the STL.
// ****************************************************************************

struct XmlError;
struct XmlAttribute;
class XmlNode;

// ****************************************************************************
// Some public functions.
// Use xmlParseString(...) or xmlParseFile(...) to get the root node of the
// XML-document. After processing, the caller must delete the returned node
// with "delete" to free the memory of the XML-tree again.
// If the return value of the parsing functions is NULL, use xmlPrintErrors(...)
// to print the parsing errors.
// ****************************************************************************

XmlNode *xmlParseString(string &input, string tag, vector<XmlError> *errors = NULL);
XmlNode *xmlParseFile(string fileName, string tag, vector<XmlError> *errors = NULL);
void xmlPrintErrors(vector<XmlError> &errors);
void xmlTest();

// ****************************************************************************
// Structure for an error detected during parsing.
// ****************************************************************************

struct XmlError
{
  int line;
  int column;
  string text;
};

// ****************************************************************************
// Structure for an attribute of an XML element.
// ****************************************************************************

struct XmlAttribute
{
  string name;
  string value;
};

// ****************************************************************************
/// An XML node.
// ****************************************************************************

class XmlNode
{
  // ****************************************************************
  // Public data.
  // ****************************************************************

public:
  enum NodeType
  {
    ELEMENT,
    TEXT,
    OTHER       ///< declaration-nodes, comments, CDATA-nodes, etc.
  };

  // Common members in all types of nodes.
  XmlNode *parent;
  NodeType type;

  // Exclusively used in ELEMENT nodes.
  string name;
  vector<XmlNode*> child;           ///< All child nodes.
  vector<XmlNode*> childElement;    ///< Only element child nodes.
  vector<XmlAttribute> attribute;

  // Exclusively used in TEXT, COMMENT and OTHER nodes.
  // These are leaf nodes without children or attributes.
  string text;

  // ****************************************************************
  // Public functions.
  // ****************************************************************

public:
  XmlNode(NodeType t, XmlNode *parent);
  ~XmlNode();

  int numChildElements(string name);
  XmlNode *getChildElement(string name, int index = 0);
  bool hasAttribute(string name);
  int getAttributeInt(string name);
  double getAttributeDouble(string name);
  string getAttributeString(string name);

  string toXmlString();

  // ****************************************************************
  // Private functions.
  // ****************************************************************

private:
  void toXmlString(ostream &os, int indent);
};


#endif
