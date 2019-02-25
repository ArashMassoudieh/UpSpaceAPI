#include "heteval.h"
#include "qdebug.h"
#include <QCoreApplication>
#include <qfiledialog.h>
#include "StringOP.h"
#include <iostream>
#include <fstream>
#include <Grid.h>


using namespace std;

HETEVAL::HETEVAL(QWidget *parent)
	: QMainWindow(parent)
{
    
        ui.setupUi(this);
	//m_button = new QPushButton("Click here", this);
	// set size and location of the button
	//m_button->setGeometry(QRect(QPoint(200, 100), QSize(200, 50)));
	connect(ui.OpenFile, SIGNAL(released()), this, SLOT(handleButton()));
	connect(ui.RunButton, SIGNAL(released()), this, SLOT(handleRunButton()));
	
}

HETEVAL::~HETEVAL()
{

}

void HETEVAL::handleRunButton()
{
	if (grid->commands.size() > 0)
	{
		grid->main_window = this;
		grid->runcommands_qt();
	}
	else
		ui.ShowOutput->append("No input file has been opened");

}

void HETEVAL::handleButton()
{
	QString fileName = QFileDialog::getOpenFileName(this,tr("Open Input File"), "", tr("Input Files (*.txt *.inp *.*)"));
	qDebug() << fileName;
	ifstream file(fileName.toStdString());
	vector<string> s;
        qDebug() << 1;
	grid = new CGrid(fileName.toStdString());
        qDebug() << 2;
	while (!file.eof())
	{
		string line;
		getline(file, line);
		QString qline = QString::fromStdString(line);
		ui.ShowInput->append(qline);
		qDebug() << qline;
		//qApp->processEvents();
	}
	ui.ShowInput->update();
}
