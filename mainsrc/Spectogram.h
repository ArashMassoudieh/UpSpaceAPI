#ifdef QT_version
#include <qwt_plot.h>
#include <qwt_plot_spectrogram.h>
#endif // QT_version
#include "Grid.h"

#ifdef QT_version
class Plot : public QwtPlot
{
	//Q_OBJECT

public:
	enum ColorMap
	{
		RGBMap,
		IndexMap,
		HueMap,
		AlphaMap
	};

	Plot(CGrid *g, QWidget * = NULL);

	public Q_SLOTS:
	void showContour(bool on);
	void showSpectrogram(bool on);
	void setColorMap(int);
	void setAlpha(int);
	CGrid *grid;

#ifndef QT_NO_PRINTER
	void printPlot();
#endif

private:
	QwtPlotSpectrogram *d_spectrogram;

	int d_mapType;
	int d_alpha;
};
#endif
