#ifndef __TIME_FUNCTION_H__
#define __TIME_FUNCTION_H__

#include <vector>

using namespace std;

// ****************************************************************************
/// This class represents a time function that is defined by a set of time-
/// value nodes. Between the nodes, the function is linearly interpolated.
// ****************************************************************************

class TimeFunction
{
  // **************************************************************************
  // Public data.
  // **************************************************************************

public:
  struct Node
  {
    double time;
    double value;
  };

  // **************************************************************************
  /// Public functions.
  // **************************************************************************

public:
  TimeFunction();
  bool setNodes(const Node n[], const int numNodes);
  bool setNodes(const vector<Node> n);
  void getNodes(vector<Node> &nodes);
  double getValue(double t);
  
  static void test();

  // **************************************************************************
  // Private data.
  // **************************************************************************

private:
  vector<Node> nodes;
};


#endif
