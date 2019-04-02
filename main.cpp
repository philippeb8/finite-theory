/**
    Finite Theory Simulator
    Copyright (C) 2011 Phil Bouchard <phil@fornux.com>

    Based on:
        "Galaxy Rotation with Dark Matter" by Dr. Armando Pisani
        https://www.compadre.org/OSP/document/ServeFile.cfm?ID=11512&DocID=2444

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#define EDITION "5.6"

#include "main.h"

#include <cmath>
#include <limits>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <typeinfo>

#include <QDebug>
#include <QtWidgets/QApplication>
#include <qevent.h>
#include <qpainter.h>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QToolButton>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QToolTip>
#include <qrect.h>
#include <qpoint.h>
#include <QtWidgets/QColorDialog>
#include <QtWidgets/QFileDialog>
#include <qcursor.h>
#include <qimage.h>
#include <QStringList>
#include <QtWidgets/QMenu>
#include <QHash>
#include <QtWidgets/QLayout>
#include <qobject.h>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QStyleFactory>
//Added by qt3to4:
#include <QtWidgets/QHBoxLayout>
#include <QTimerEvent>
#include <QtWidgets/QLabel>
#include <QPixmap>
#include <QResizeEvent>
#include <QMouseEvent>
#include <QtWidgets/QVBoxLayout>
#include <QPaintEvent>
#include <QVector>
#include <QtGui/QPainter>
#include <QtWidgets/QDesktopWidget>
#include <QtCore/QProcess>

using namespace std;


void bitBlt( QPaintDevice * dst, int x, int y, const QPixmap* src, int sx, int sy, int sw, int sh )
{
    QPainter p( dst );
    p.drawPixmap( x, y, *src, sx, sy, sw, sh );
}


/** 
	@brief			Calculates the next position of the planet or photon
	@param planet	Planets that will affect the movement of the planet that is moving (this)
    @param dt       Time interval
*/

inline void Planet::operator () (const real & dt)
{
    alpha = alpha + omega * dt;

    p = vector3(r * cos(alpha), r * sin(alpha), 0);

/*
    real massf = 1.0;
    real totalmass = pow(planet.back().vel, 2) * planet.back().r;
    real starmass = totalmass / (8.0 * planet.size()); // here Mtot/Mvisible=8 as example

    vvis = sqrt(i * starmass * massf / r);
    omegavis = vvis / r ;
    massvis = vvis * vvis * r ;

    real rdm0 = 1.0;
    real dmf = 1.0;
    real md = totalmass - starmass * massf * planet.size();
    real mdk = md * dmf / (planet.back().r / rdm0 - atan(planet.back().r / rdm0));

    vdark = sqrt(mdk * (r / rdm0 - atan(r / rdm0)) / r );
    massdk = vdark * vdark * r;
*/
}

Dual::Dual(Canvas * pParent) : p(pParent)
{
    start();
}

void Dual::run()
{
	Scribble * q = static_cast<Scribble *>(p->topLevelWidget());
	
    while (true)
	{
		while (! p->isVisible())
			QThread::msleep(100);
		
        for (size_t j = 0; j < q->planet.size(); ++ j)
            for (size_t i = 1; i < q->planet[j].size(); ++ i)
                q->planet[j][i](q->pTime->value());
    }
}

Canvas::Canvas( int type, QWidget *parent )
    : QWidget( parent ),
      type(type), pen( Qt::red, 3 ), polyline(3), first(true), mousePressed( false ), buffer( width(), height() )
{
	Scribble * q = static_cast<Scribble *>(topLevelWidget());

#ifndef QT_NO_CURSOR
    setCursor( Qt::CrossCursor );
#endif

    QPalette p = palette();
    p.setColor(QPalette::Base, Qt::black);
    setPalette(p);

	startTimer(100);

    // launch a thread for each set of planets or photons
    new Dual(this);
}

Canvas::~Canvas()
{
	exit(-1);
}

void Canvas::timerEvent(QTimerEvent *)
{
	Scribble * q = static_cast<Scribble *>(topLevelWidget());

    if (! isVisible())
        return;

    switch (type)
    {
    case 0:
        {
            QRect r(q->planet[0][0].p[0] / scale[0] - 5 + width()/2, q->planet[0][0].p[1] / scale[0] - 5 + height()/2, 10, 10);
            QPainter painter;

            painter.begin( &buffer );
            painter.setBrush(Qt::yellow);
            painter.drawEllipse(r);
            painter.end();

            for (size_t j = 0; j < q->planet.size(); ++ j)
            {
                for (size_t i = 1; i < q->planet[j].size(); ++ i)
                {
                    QRect e(q->planet[j][i].o[0] / scale[0] - 2 + width()/2, q->planet[j][i].o[1] / scale[0] - 2 + height()/2, 4+1, 4+1);

                    q->planet[j][i].o = q->planet[j][i].p;

                    QRect r(q->planet[j][i].o[0] / scale[0] - 2 + width()/2, q->planet[j][i].o[1] / scale[0] - 2 + height()/2, 4, 4);
                    QPainter painter;

                    painter.begin( &buffer );
                    painter.setPen(q->planet[j][i].c);
                    painter.setBrush(q->planet[j][i].c);
                    painter.fillRect(e, palette().base());

                    if (q->pTheory[j]->checkState() == Qt::Checked)
                        painter.drawEllipse(r);

                    painter.end();
                    r |= e;
                }

                update();
            }
        }
        break;

    case 1:
        if (first)
        {
            first = false;

            real vmax = 0.0;

            for (size_t j = 0; j < q->planet.size(); ++ j)
                for (size_t i = 1; i < q->planet[j].size(); ++ i)
                    if (vmax < q->planet[j][i].vel)
                        vmax = q->planet[j][i].vel;

            for (size_t j = 0; j < q->planet.size(); ++ j)
            {
                for (size_t i = 1; i < q->planet[j].size(); ++ i)
                {
                    QRect r(q->planet[j][i].r / scale[0] - 2, height() - q->planet[j][i].vel / vmax * height() - 2, 4, 4);
                    QPainter painter;

                    painter.begin( &buffer );
                    painter.setPen(q->planet[j][i].c);
                    painter.setBrush(q->planet[j][i].c);

                    if (q->pTheory[j]->checkState() == Qt::Checked)
                        painter.drawEllipse(r);

                    painter.end();
                }
            }

            for (real y = 0; y < vmax; y += vmax / 5.0)
            {
                QPainter painter;

                painter.begin( &buffer );
                painter.setPen(Qt::white);
                painter.setBrush(Qt::white);
                painter.drawText(0, height() - y / vmax * height(), QString::number(y / Scribble::zoom) + " km/s");
                painter.end();
            }

            update();
        }

        break;
    };
}

void Canvas::clearScreen()
{
    buffer.fill( palette().base().color() );
    repaint();
}

void Canvas::mousePressEvent( QMouseEvent *e )
{
    mousePressed = true;
    polyline[2] = polyline[1] = polyline[0] = e->pos();
}

void Canvas::mouseReleaseEvent( QMouseEvent * )
{
    mousePressed = false;
}

void Canvas::mouseMoveEvent( QMouseEvent *e )
{
    if ( mousePressed )
    {
        QPainter painter;
        painter.begin( &buffer );
        painter.setPen( Qt::gray );
        polyline[2] = polyline[1];
        polyline[1] = polyline[0];
        polyline[0] = e->pos();
        painter.drawPolyline( polyline );
        painter.end();

        QRect r = polyline.boundingRect();
        r = r.normalized();
        r.setLeft( r.left() - pen.width() );
        r.setTop( r.top() - pen.width() );
        r.setRight( r.right() + pen.width() );
        r.setBottom( r.bottom() + pen.width() );

        bitBlt( this, r.x(), r.y(), &buffer, r.x(), r.y(), r.width(), r.height() );
    }
}

void Canvas::resizeEvent( QResizeEvent *e )
{
    QWidget::resizeEvent( e );

    int w = width() > buffer.width() ?
            width() : buffer.width();
    int h = height() > buffer.height() ?
            height() : buffer.height();

    buffer = QPixmap(w, h);

    buffer.fill( palette().base().color() );
}

void Canvas::paintEvent( QPaintEvent *e )
{
    QWidget::paintEvent( e );

    QVector<QRect> rects = e->region().rects();
    for ( uint i = 0; i < rects.count(); ++ i ) 
	{
        QRect r = rects[(int)i];
        bitBlt( this, r.x(), r.y(), &buffer, r.x(), r.y(), r.width(), r.height() );
    }
}

//------------------------------------------------------

char * const Scribble::theory[Scribble::nt] = {"Initial", "Visible", "Dark Matter", "Total", "Finite Theory"};

Scribble::Scribble( QWidget *parent, const char *name )
    : QMainWindow( parent )
{
    ntime = 0.0000005;

    QMenu *file = new QMenu( "&File", this );
    file->addAction( "&Restart", this, SLOT(slotRestart()), Qt::CTRL+Qt::Key_R );
    file->addSeparator();
    file->addAction( "E&xit", qApp, SLOT(quit()), Qt::CTRL+Qt::Key_Q );

    QMenu *help = new QMenu( "&Help", this );
    help->addAction( "&About", this, SLOT(slotAbout()));

	QMenuBar * menu = new QMenuBar( this );
    menu->addMenu( file );
    menu->addSeparator();
    menu->addMenu( help );
    setMenuBar(menu);

    QToolBar *tools = new QToolBar( this );

    pTime = new QDoubleSpinBox( tools );
    pTime->setRange(0.0000005, 50000000000);
    pTime->setDecimals(7);
    pTime->setSingleStep(10);
    pTime->setToolTip("Time Interval (s)");
    pTime->setValue( ntime );

    tools->addWidget(pTime);
    tools->addSeparator();

    for (size_t i = 0; i < nt; ++ i)
    {
        pTheory[i] = new QCheckBox(theory[i], tools);
        tools->addWidget(pTheory[i]);
    }

    pTheory[1]->setCheckState(Qt::Checked);
    pTheory[3]->setCheckState(Qt::Checked);
    pTheory[4]->setCheckState(Qt::Checked);

    addToolBar(tools);
	
	pTabWidget = new QTabWidget(this);

    QPalette* palette = new QPalette();
    palette->setColor(QPalette::WindowText, Qt::darkRed);

    for (size_t i = 0; i < ntabs; ++ i)
	{
		pTab[i] = new QWidget(pTabWidget);

        canvas[i] = new Canvas(i, pTab[i]);
		canvas[i]->setSizePolicy(QSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding));

        QBoxLayout * l = new QVBoxLayout(pTab[i]);
	    l->addWidget( canvas[i] );
    }

    pTabWidget->addTab(pTab[0], "Galactic Rotation");
    pTabWidget->addTab(pTab[1], "Galactic Rotation Curve");

    setCentralWidget( pTabWidget );

    // galactic rotation
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<real> dis(0.0, 1.0);

    planet.resize(nt);

    for (size_t j = 0; j < planet.size(); ++ j)
    {
        planet[j].reserve(np + 1);
        planet[j].push_back(Planet("Nucleus", Qt::black, 2E+12L));

        for (int i = 0; i < np; ++ i)
            planet[j].push_back(Planet("Star1", Qt::red, 50000L));
    }

    for (size_t i = 1; i < np + 1; ++ i)
    {
        planet[0][i].alpha = PI / 2 * ceil(dis(gen) * 4.0);
        planet[0][i].r = h * sqrt(tan(PI * 0.5 * emax * i / (np + 1)));
        planet[0][i].vel = v0 * 2.0 / PI * atan(planet[0][i].r / r0);
        planet[0][i].omega = planet[0][i].vel / planet[0][i].r;
        planet[0][i].mass = planet[0][i].vel * planet[0][i].vel * planet[0][i].r;
    }

    totalmass = pow(planet[0][np].vel, 2) * planet[0][np].r;
    starmass = totalmass / (8.0 * (np + 1));

    for (size_t i = 1; i < np + 1; ++ i)
    {
        planet[1][i].alpha = planet[0][i].alpha;
        planet[1][i].r = planet[0][i].r;
        planet[1][i].vel = sqrt(i * starmass * massf / planet[1][i].r);
        planet[1][i].omega = planet[1][i].vel / planet[1][i].r;
        planet[1][i].mass = planet[1][i].vel * planet[1][i].vel * planet[1][i].r;
    }

    md = totalmass - starmass * massf * (np + 1);
    mdk = md * dmf / (planet[0][np].r / rdm0 - atan(planet[0][np].r / rdm0));

    for (size_t i = 1; i < np + 1; i++)
    {
        planet[2][i].alpha = planet[0][i].alpha;
        planet[2][i].r = planet[0][i].r;
        planet[2][i].vel = sqrt(mdk * (planet[2][i].r / rdm0 - atan(planet[2][i].r / rdm0)) / planet[2][i].r);
        planet[2][i].omega = planet[2][i].vel / planet[2][i].r;
        planet[2][i].mass = planet[2][i].vel * planet[2][i].vel * planet[2][i].r;
    }

    for (size_t i = 1; i < np + 1; i++)
    {
        planet[3][i].alpha = planet[0][i].alpha;
        planet[3][i].r = planet[0][i].r;
        planet[3][i].vel = sqrt(planet[1][i].vel * planet[1][i].vel + planet[2][i].vel * planet[2][i].vel);
        planet[3][i].omega = planet[3][i].vel / planet[3][i].r;
        planet[3][i].mass = planet[1][i].mass + planet[2][i].mass;
    }

    for (size_t i = 1; i < np + 1; i++)
    {
        planet[4][i].alpha = planet[0][i].alpha;
        planet[4][i].r = planet[0][i].r;
        planet[4][i].vel = sqrt(mdk * (planet[4][i].r / rdm0 - atan(planet[4][i].r / rdm0)) / planet[4][i].r);
        planet[4][i].omega = planet[4][i].vel / planet[4][i].r + 2.83273668e-16;
        planet[4][i].mass = planet[4][i].vel * planet[4][i].vel * planet[4][i].r;
    }

    for (size_t j = planet.size(); j > 1; -- j)
        for (size_t i = 1; i < np + 1; ++ i)
            planet[j - 2][i].c = planet[j - 1][i].c.dark();
}

void Scribble::slotRestart()
{
    qApp->quit();
    QProcess::startDetached(qApp->arguments()[0], qApp->arguments());
}

void Scribble::slotAbout()
{
    QMessageBox::about( this, "Finite Theory of the Universe " EDITION, "\nCopyright (c) 2011-2019\n\nPhil Bouchard <phil@fornux.com>\n");
}
	

int main( int argc, char **argv )
{
    QApplication a( argc, argv );

    a.setStyle("windows");

    Scribble scribble;

    scribble.setWindowTitle("Finite Theory of the Universe " EDITION);
    scribble.showMaximized();
	
    return a.exec();
}
