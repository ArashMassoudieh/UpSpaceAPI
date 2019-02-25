#ifndef HETEVAL_H
#define HETEVAL_H

#include <QtWidgets/QMainWindow>
#include "ui_heteval.h"
#include <qpushbutton.h>
#include "Grid.h"

namespace Ui {

	class HETEVAL;
}

class HETEVAL : public QMainWindow
{
	Q_OBJECT

public:
	explicit HETEVAL(QWidget *parent = 0);
	~HETEVAL();
	CGrid *grid;
	Ui::HETEVALClass* get_ui()
	{
		return &ui;
	}

private:
	Ui::HETEVALClass ui;

private slots:
	void handleButton();
	void handleRunButton();
};

#endif // HETEVAL_H
