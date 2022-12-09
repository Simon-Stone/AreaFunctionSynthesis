#include "AreaFunctionPicture.h"
#include "../Backend/OneDimAreaFunction.h"
#include "Data.h"

// ****************************************************************************
// The event table.
// ****************************************************************************

BEGIN_EVENT_TABLE(AreaFunctionPicture, BasicPicture)
	EVT_MOUSE_EVENTS(AreaFunctionPicture::OnMouseEvent)
END_EVENT_TABLE()

// ****************************************************************************
/// Constructor. Passes the parent parameter.
// ****************************************************************************

AreaFunctionPicture::AreaFunctionPicture(wxWindow *parent) : BasicPicture(parent)
{
	// ****************************************************************
	// Init variables.
	// ****************************************************************
	isDragging = false;
	showTubeFunction = false;
	showAreaFunction = true;

	// ****************************************************************
	// Init. the scales.
	// ****************************************************************

	const int LEFT_MARGIN = 80;
	const int RIGHT_MARGIN = 0;
	const int TOP_MARGIN = 0;
	const int BOTTOM_MARGIN = 25;


	graph = new Graph();
	graph->init(this, LEFT_MARGIN, RIGHT_MARGIN, TOP_MARGIN, BOTTOM_MARGIN);
	graph->initAbscissa(PQ_LENGTH, 0.0, 1.0,
		0.0, 0.0, 0.0, 18.0, 18.0, 18.0,
		1, 1, true, true, true);
	graph->initLinearOrdinate(PQ_AREA, 0.0, 1.0,
		0.0, -2.0, -2.0,
		1.0, 10.0, 8.0,
		1, 2, true, true, true);

	
	//Data *data = Data::getInstance();
	OneDimAreaFunction::Param currentParams[OneDimAreaFunction::NUM_AF_PARAMS];
	//data->synthesizer->vtAreaFunction->getParameters(currentParams);
	//data->refreshAreaFunction();
	this->Refresh();
}


// ****************************************************************************
/// Draw the picture.
// ****************************************************************************

void AreaFunctionPicture::draw(wxDC &dc)
{
	Data *data = Data::getInstance();

	int graphX, graphY, graphW, graphH;
	graph->getDimensions(graphX, graphY, graphW, graphH);

	int width, height;
	this->GetSize(&width, &height);

	dc.SetBackground(wxBrush(wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOW) ));
	dc.Clear();

	// Set abscissa and ordinate dimensions to fit current area function(s)
	double maxA_cm2 = 0.0;

	for (int i = 0; i < data->synthesizer->vtAreaFunction->NUM_AREA_SAMPLES; i++)
	{
		if (data->synthesizer->vtAreaFunction->areaFunction_cm2[i] > maxA_cm2)
		{
			maxA_cm2 = data->synthesizer->vtAreaFunction->areaFunction_cm2[i];
		}
		graph->linearOrdinate.positiveLimit = maxA_cm2 + 5.0;
	}
	OneDimAreaFunction::Param currentParams[OneDimAreaFunction::NUM_AF_PARAMS];
	data->synthesizer->vtAreaFunction->getParameters(currentParams);
	graph->abscissa.positiveLimit = (currentParams[OneDimAreaFunction::LVT_CM].x + 5.0);
	
	
	// ****************************************************************
	// Paint the scales.
	// ****************************************************************
	graph->paintAbscissa(dc);
	graph->paintOrdinate(dc);
	wxPoint p0, p1;

	if (showAreaFunction)
	{
		// ****************************************************************
		// Draw the vocal tract area function.
		// ****************************************************************
		dc.SetPen(wxPen(wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOWTEXT), 2));

		p0.x = graph->getXPos(0.0);
		p0.y = graph->getYPos(data->synthesizer->vtAreaFunction->areaFunction_cm2[0]);
		for (int i = 1; i < data->synthesizer->vtAreaFunction->NUM_AREA_SAMPLES; i++)
		{
			p1.x = graph->getXPos(i*data->synthesizer->vtAreaFunction->deltaX_cm);
			p1.y = graph->getYPos(data->synthesizer->vtAreaFunction->areaFunction_cm2[i]);
			dc.DrawLine(p0, p1);
			p0 = p1;
		}
	}

	if (showTubeFunction)
	{
		dc.SetPen(wxPen(wxSystemSettings::GetColour(wxSYS_COLOUR_WINDOWTEXT) , 2));

		for (int i = 0; i < data->synthesizer->tubeFunction->NUM_PHARYNX_MOUTH_SECTIONS; i++)
		{
			p0.x = graph->getXPos(data->synthesizer->tubeFunction->pharynxMouthSection[i].pos_cm);
			p0.y = graph->getYPos(data->synthesizer->tubeFunction->pharynxMouthSection[i].area_cm2);
			dc.DrawLine(p0.x, p0.y, p0.x, graph->getYPos(0.0));
			p1.x = graph->getXPos(data->synthesizer->tubeFunction->pharynxMouthSection[i].pos_cm
				+ data->synthesizer->tubeFunction->pharynxMouthSection[i].length_cm);
			p1.y = p0.y;
			dc.DrawLine(p0, p1);

			dc.DrawLine(p1.x, p1.y, p1.x, graph->getYPos(0.0));
		}
	}

	if (showAreaFunction || showTubeFunction)
	{
		// ****************************************************************
		// Draw control points.
		// ****************************************************************
        const auto radius = FromDIP(10);

		OneDimAreaFunction::Param currentParams[OneDimAreaFunction::NUM_AF_PARAMS];
		data->synthesizer->vtAreaFunction->getParameters(currentParams);

		dc.SetBrush(wxBrush(wxSystemSettings::GetColour(wxSYS_COLOUR_HIGHLIGHT)));
		dc.SetPen(wxPen(wxSystemSettings::GetColour(wxSYS_COLOUR_HIGHLIGHT), 2));
		// Draw larynx tube control point
		ctrlPointLar.x = graph->getXPos(currentParams[OneDimAreaFunction::LLAR_CM].x);
		ctrlPointLar.y = graph->getYPos(currentParams[OneDimAreaFunction::ALAR_CM2].x);
		dc.DrawCircle(ctrlPointLar, radius);

		dc.SetBrush(wxBrush(wxSystemSettings::GetColour(wxSYS_COLOUR_HIGHLIGHT)));
		dc.SetPen(wxPen(wxSystemSettings::GetColour(wxSYS_COLOUR_HIGHLIGHT), 2));

		// Draw posterior maximum control point
		ctrlPointP.x = graph->getXPos(currentParams[OneDimAreaFunction::XP_CM].x);
		ctrlPointP.y = graph->getYPos(currentParams[OneDimAreaFunction::AP_CM2].x);
		dc.DrawCircle(ctrlPointP, radius);

		// Draw constriction control point
		ctrlPointC.x = graph->getXPos(currentParams[OneDimAreaFunction::XC_CM].x);
		ctrlPointC.y = graph->getYPos(currentParams[OneDimAreaFunction::AC_CM2].x);
		dc.DrawCircle(ctrlPointC, radius);

		// Draw anterior maximum control point
		ctrlPointA.x = graph->getXPos(currentParams[OneDimAreaFunction::XA_CM].x);
		ctrlPointA.y = graph->getYPos(currentParams[OneDimAreaFunction::AA_CM2].x);
		dc.DrawCircle(ctrlPointA, radius);

		// Draw incisor control point
		ctrlPointIn.x = graph->getXPos(currentParams[OneDimAreaFunction::XIN_CM].x);
		ctrlPointIn.y = graph->getYPos(currentParams[OneDimAreaFunction::AIN_CM2].x);
		dc.DrawCircle(ctrlPointIn, radius);

		// Draw vocal tract length and lip opening control poin
		ctrlPointLip.x = graph->getXPos(currentParams[OneDimAreaFunction::LVT_CM].x);
		ctrlPointLip.y = graph->getYPos(currentParams[OneDimAreaFunction::ALIP_CM2].x);
		dc.DrawCircle(ctrlPointLip, radius);
	}
}

void AreaFunctionPicture::OnMouseEvent(wxMouseEvent &event)
{
	Data *data = Data::getInstance();
	if (data->synthesizer->modelControlMode == Synthesizer::MANUAL)
	{
		int mx = event.GetX();
		int my = event.GetY();

		OneDimAreaFunction::Param currentParams[OneDimAreaFunction::NUM_AF_PARAMS];

		data->synthesizer->vtAreaFunction->getParameters(currentParams);

		int graphX, graphY, graphW, graphH;
		graph->getDimensions(graphX, graphY, graphW, graphH);
		if (graphW < 1)
		{
			graphW = 1;
		}

		// ****************************************************************
		// The mouse is entering the window.
		// ****************************************************************

		if (event.Entering())
		{
			lastMx = mx;
			lastMy = my;
			return;
		}

		// ****************************************************************
		// Check if any control point is under the cursor
		// ****************************************************************
		if ((mx >= ctrlPointLar.x - 4) && (mx <= ctrlPointLar.x + 4) && (my >= ctrlPointLar.y - 4) && (my <= ctrlPointLar.y + 4))
		{
			this->SetCursor(wxCursor(wxCURSOR_HAND));
			moveCtrlPointIdx = CTRL_LARYNX;
		}
		else if ((mx >= ctrlPointP.x - 4) && (mx <= ctrlPointP.x + 4) && (my >= ctrlPointP.y - 4) && (my <= ctrlPointP.y + 4))
		{
			this->SetCursor(wxCursor(wxCURSOR_HAND));
			moveCtrlPointIdx = CTRL_MAX_POST;
		}
		else if ((mx >= ctrlPointC.x - 4) && (mx <= ctrlPointC.x + 4) && (my >= ctrlPointC.y - 4) && (my <= ctrlPointC.y + 4))
		{
			this->SetCursor(wxCursor(wxCURSOR_HAND));
			moveCtrlPointIdx = CTRL_CONSTRICT;
		}
		else if ((mx >= ctrlPointA.x - 4) && (mx <= ctrlPointA.x + 4) && (my >= ctrlPointA.y - 4) && (my <= ctrlPointA.y + 4))
		{
			this->SetCursor(wxCursor(wxCURSOR_HAND));
			moveCtrlPointIdx = CTRL_MAX_ANT;
		}
		else if ((mx >= ctrlPointIn.x - 4) && (mx <= ctrlPointIn.x + 4) && (my >= ctrlPointIn.y - 4) && (my <= ctrlPointIn.y + 4))
		{
			this->SetCursor(wxCursor(wxCURSOR_HAND));
			moveCtrlPointIdx = CTRL_INCISOR;
		}
		else if ((mx >= ctrlPointLip.x - 4) && (mx <= ctrlPointLip.x + 4) && (my >= ctrlPointLip.y - 4) && (my <= ctrlPointLip.y + 4))
		{
			this->SetCursor(wxCursor(wxCURSOR_HAND));
			moveCtrlPointIdx = CTRL_LIP;
		}
		else
		{
			this->SetCursor(wxCursor(wxCURSOR_ARROW));
		}


		// ****************************************************************
		// Stop dragging a control point when it is released or the area
		// is left.
		// ****************************************************************

		if ((event.ButtonUp()) || (event.Leaving()))
		{
			if (lastMx == mx && lastMy == my && !isDragging)
			{
				/* Mouse event was a single click (not a drag) */
				if (event.Button(wxMOUSE_BTN_LEFT))
				{
					switch (moveCtrlPointIdx)
					{
					case CTRL_NONE:
						break;
					case CTRL_LARYNX:
						currentParams[OneDimAreaFunction::POWLAR].x += 0.1;
						break;
					case CTRL_MAX_POST:
						currentParams[OneDimAreaFunction::POWP].x += 0.1;
						break;
					case CTRL_CONSTRICT:
						currentParams[OneDimAreaFunction::POWC].x += 0.1;
						break;
					case CTRL_MAX_ANT:
						currentParams[OneDimAreaFunction::POWA].x += 0.1;
						break;
					case CTRL_INCISOR:
						break;
					case CTRL_LIP:
						break;
					default:
						break;
					}
				}
				else if (event.Button(wxMOUSE_BTN_RIGHT))
				{
					switch (moveCtrlPointIdx)
					{
					case CTRL_NONE:
						break;
					case CTRL_LARYNX:
						currentParams[OneDimAreaFunction::POWLAR].x -= 0.1;
						if (currentParams[OneDimAreaFunction::POWLAR].x < 0)
						{
							currentParams[OneDimAreaFunction::POWLAR].x = 0;
						}
						break;
					case CTRL_MAX_POST:
						currentParams[OneDimAreaFunction::POWP].x -= 0.1;
						if (currentParams[OneDimAreaFunction::POWP].x < 0)
						{
							currentParams[OneDimAreaFunction::POWP].x = 0;
						}
						break;
					case CTRL_CONSTRICT:
						currentParams[OneDimAreaFunction::POWC].x -= 0.1;
						if (currentParams[OneDimAreaFunction::POWC].x < 0)
						{
							currentParams[OneDimAreaFunction::POWC].x = 0;
						}
						break;
					case CTRL_MAX_ANT:
						currentParams[OneDimAreaFunction::POWA].x -= 0.1;
						if (currentParams[OneDimAreaFunction::POWA].x < 0)
						{
							currentParams[OneDimAreaFunction::POWA].x = 0;
						}
						break;
					case CTRL_INCISOR:
						break;
					case CTRL_LIP:
						break;
					default:
						break;
					}
				}
				data->synthesizer->vtAreaFunction->setParameters(currentParams);
				data->refreshAreaFunction();
				this->Refresh();
			}
			isDragging = false;
			lastMx = mx;
			lastMy = my;
			moveCtrlPointIdx = CTRL_NONE;
			this->SetCursor(wxCursor(wxCURSOR_ARROW));

			return;
		}

		// ****************************************************************
		// The left mouse button just changed to down.
		// ****************************************************************

		if (event.ButtonDown(wxMOUSE_BTN_LEFT))
		{
			lastMx = mx;
			lastMy = my;
			return;
		}

		// ****************************************************************
		// The user is dragging the mouse.
		// ****************************************************************

		if (event.Dragging())
		{
			if ((event.LeftIsDown()) && (mx >= graphX) && (mx < graphX + graphW))
			{
				isDragging = true;
				double deltaPosX_cm, deltaPosY_cm2;
				// Start moving the cut position.
				switch (moveCtrlPointIdx)
				{
				case CTRL_NONE:
					this->SetCursor(wxCursor(wxCURSOR_ARROW));
					break;
				case CTRL_LARYNX:
					/* Larynx tube length and area is fixed for now */
					deltaPosX_cm = graph->getAbsXValue(mx) - graph->getAbsXValue(lastMx);
					currentParams[OneDimAreaFunction::LLAR_CM].x += deltaPosX_cm;
					if (currentParams[OneDimAreaFunction::LLAR_CM].x > currentParams[OneDimAreaFunction::XP_CM].x)
					{
						currentParams[OneDimAreaFunction::LLAR_CM].x = currentParams[OneDimAreaFunction::XP_CM].x;
					}
					if (currentParams[OneDimAreaFunction::LLAR_CM].x < 0.0)
					{
						currentParams[OneDimAreaFunction::LLAR_CM].x = 0.0;
					}
					deltaPosY_cm2 = graph->getAbsYValue(my) - graph->getAbsYValue(lastMy);
					currentParams[OneDimAreaFunction::ALAR_CM2].x += deltaPosY_cm2;
					break;
				case CTRL_MAX_POST:
					deltaPosX_cm = graph->getAbsXValue(mx) - graph->getAbsXValue(lastMx);
					currentParams[OneDimAreaFunction::XP_CM].x += deltaPosX_cm;
					if (currentParams[OneDimAreaFunction::XP_CM].x > currentParams[OneDimAreaFunction::XC_CM].x)
					{
						currentParams[OneDimAreaFunction::XP_CM].x = currentParams[OneDimAreaFunction::XC_CM].x;
					}
					if (currentParams[OneDimAreaFunction::XP_CM].x < currentParams[OneDimAreaFunction::LLAR_CM].x)
					{
						currentParams[OneDimAreaFunction::XP_CM].x = currentParams[OneDimAreaFunction::LLAR_CM].x;
					}
					deltaPosY_cm2 = graph->getAbsYValue(my) - graph->getAbsYValue(lastMy);
					currentParams[OneDimAreaFunction::AP_CM2].x += deltaPosY_cm2;
					if (currentParams[OneDimAreaFunction::AP_CM2].x < -2.0)
					{
						currentParams[OneDimAreaFunction::AP_CM2].x = -2.0;
					}
					break;
				case CTRL_CONSTRICT:
					deltaPosX_cm = graph->getAbsXValue(mx) - graph->getAbsXValue(lastMx);
					currentParams[OneDimAreaFunction::XC_CM].x += deltaPosX_cm;
					if (currentParams[OneDimAreaFunction::XC_CM].x > currentParams[OneDimAreaFunction::XA_CM].x)
					{
						currentParams[OneDimAreaFunction::XC_CM].x = currentParams[OneDimAreaFunction::XA_CM].x;
					}
					if (currentParams[OneDimAreaFunction::XC_CM].x < currentParams[OneDimAreaFunction::XP_CM].x)
					{
						currentParams[OneDimAreaFunction::XC_CM].x = currentParams[OneDimAreaFunction::XP_CM].x;
					}
					deltaPosY_cm2 = graph->getAbsYValue(my) - graph->getAbsYValue(lastMy);
					currentParams[OneDimAreaFunction::AC_CM2].x += deltaPosY_cm2;
					if (currentParams[OneDimAreaFunction::AC_CM2].x < -2.0)
					{
						currentParams[OneDimAreaFunction::AC_CM2].x = -2.0;
					}
					break;
				case CTRL_MAX_ANT:
					deltaPosX_cm = graph->getAbsXValue(mx) - graph->getAbsXValue(lastMx);
					currentParams[OneDimAreaFunction::XA_CM].x += deltaPosX_cm;
					if (currentParams[OneDimAreaFunction::XA_CM].x > currentParams[OneDimAreaFunction::XIN_CM].x)
					{
						currentParams[OneDimAreaFunction::XA_CM].x = currentParams[OneDimAreaFunction::XIN_CM].x;
					}
					if (currentParams[OneDimAreaFunction::XA_CM].x < currentParams[OneDimAreaFunction::XC_CM].x)
					{
						currentParams[OneDimAreaFunction::XA_CM].x = currentParams[OneDimAreaFunction::XC_CM].x;
					}
					deltaPosY_cm2 = graph->getAbsYValue(my) - graph->getAbsYValue(lastMy);
					currentParams[OneDimAreaFunction::AA_CM2].x += deltaPosY_cm2;
					if (currentParams[OneDimAreaFunction::AA_CM2].x < -2.0)
					{
						currentParams[OneDimAreaFunction::AA_CM2].x = -2.0;
					}
					break;
				case CTRL_INCISOR:
					deltaPosX_cm = graph->getAbsXValue(mx) - graph->getAbsXValue(lastMx);
					currentParams[OneDimAreaFunction::XIN_CM].x += deltaPosX_cm;
					if (currentParams[OneDimAreaFunction::XIN_CM].x > currentParams[OneDimAreaFunction::LVT_CM].x)
					{
						currentParams[OneDimAreaFunction::XIN_CM].x = currentParams[OneDimAreaFunction::LVT_CM].x;
					}
					if (currentParams[OneDimAreaFunction::XIN_CM].x < currentParams[OneDimAreaFunction::XA_CM].x)
					{
						currentParams[OneDimAreaFunction::XIN_CM].x = currentParams[OneDimAreaFunction::XA_CM].x;
					}
					deltaPosY_cm2 = graph->getAbsYValue(my) - graph->getAbsYValue(lastMy);
					currentParams[OneDimAreaFunction::AIN_CM2].x += deltaPosY_cm2;
					if (currentParams[OneDimAreaFunction::AIN_CM2].x < -2.0)
					{
						currentParams[OneDimAreaFunction::AIN_CM2].x = -2.0;
					}
					break;
				case CTRL_LIP:
					deltaPosX_cm = graph->getAbsXValue(mx) - graph->getAbsXValue(lastMx);
					currentParams[OneDimAreaFunction::LVT_CM].x += deltaPosX_cm;
					if (currentParams[OneDimAreaFunction::LVT_CM].x > VT_LEN_MAX)
					{
						currentParams[OneDimAreaFunction::LVT_CM].x = VT_LEN_MAX;
					}
					if (currentParams[OneDimAreaFunction::LVT_CM].x < currentParams[OneDimAreaFunction::XIN_CM].x)
					{
						currentParams[OneDimAreaFunction::LVT_CM].x = currentParams[OneDimAreaFunction::XIN_CM].x;
					}
					deltaPosY_cm2 = graph->getAbsYValue(my) - graph->getAbsYValue(lastMy);
					currentParams[OneDimAreaFunction::ALIP_CM2].x += deltaPosY_cm2;
					if (currentParams[OneDimAreaFunction::ALIP_CM2].x < -2.0)
					{
						currentParams[OneDimAreaFunction::ALIP_CM2].x = -2.0;
					}
					break;
				default:
					break;
				}

				data->synthesizer->vtAreaFunction->setParameters(currentParams);
				data->refreshAreaFunction();
				// Update the picture.
				this->Refresh();
			}
			lastMx = mx;
			lastMy = my;
			return;
		}
	}
}
