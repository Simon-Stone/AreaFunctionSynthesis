#ifndef _AREA_FUNCTION_PICTURE_H__
#define _AREA_FUNCTION_PICTURE_H__

#include "BasicPicture.h"
#include "Graph.h"


/// ****************************************************************************
/// This class represents the area function picture.
/// ****************************************************************************
enum ctrlPoint { CTRL_NONE, CTRL_LARYNX, CTRL_MAX_POST, CTRL_CONSTRICT, CTRL_MAX_ANT, CTRL_INCISOR, CTRL_LIP};

static const double VT_LEN_MAX = 25.0;

class AreaFunctionPicture : public BasicPicture
{
	// **************************************************************************
	// Public data.
	// **************************************************************************

public:
	Graph *graph;
	bool showTubeFunction;
	bool showAreaFunction;

	// **************************************************************************
	// Public functions.
	// **************************************************************************

public:
	AreaFunctionPicture(wxWindow *parent);
	virtual void draw(wxDC &dc);
	

	// **************************************************************************
	// Private data.
	// **************************************************************************

private:
	/* Area function control points */
	wxPoint ctrlPointLar;
	wxPoint ctrlPointP;
	wxPoint ctrlPointC;
	wxPoint ctrlPointA;
	wxPoint ctrlPointIn;
	wxPoint ctrlPointLip;
	ctrlPoint moveCtrlPointIdx;
	bool isDragging;

	

	/* Last mouse event coordinates */
	int lastMx;
	int lastMy;

	// **************************************************************************
	// Private functions.
	// **************************************************************************
	void OnMouseEvent(wxMouseEvent &event);
private:


	// ****************************************************************************
	// Declare the event table right at the end
	// ****************************************************************************

	DECLARE_EVENT_TABLE()
};


#endif

// ****************************************************************************