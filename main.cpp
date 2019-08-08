/**
    Finite Theory Simulator
    Copyright (C) 2011 Phil Bouchard <phil@fornux.com>

    Based on:
        "Galaxy Rotation with Dark Matter" simulator by Dr. Armando Pisani
        https://www.compadre.org/OSP/document/ServeFile.cfm?ID=11512&DocID=2444

        "Finite Theory of the Universe [...]" debate on Cosmoquest.org by Phil Bouchard
        http://forum.cosmoquest.org/showthread.php?137103-Finite-Theory-of-the-Universe-Dark-Matter-Disproof-and-Faster-Than-Light-Speed&p=2464825#post2464825

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

#define EDITION "6.1"

#include "main.h"

#include <cmath>
#include <random>
#include <limits>
#include <fstream>
#include <iostream>

#include <QtWidgets/QApplication>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QToolButton>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QToolTip>
#include <QtWidgets/QMenu>
#include <QtWidgets/QLayout>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QStyleFactory>
#include <QtWidgets/QHBoxLayout>
#include <QTimerEvent>
#include <QtWidgets/QLabel>
#include <QPixmap>
#include <QResizeEvent>
#include <QMouseEvent>
#include <QtWidgets/QVBoxLayout>
#include <QPaintEvent>
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

Dual::Dual(Scribble * pParent) : QThread(pParent)
{
    start();
}

void Dual::run()
{
    Scribble * q = static_cast<Scribble *>(parent());

    while (true)
        for (size_t j = 0; j < q->planet.size(); ++ j)
            for (size_t i = 0; i < q->planet[j].size(); ++ i)
                q->planet[j][i](q->pTime->value());
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
    p.setColor(QPalette::Base, Qt::white);
    setPalette(p);

    for (size_t i = 0; i < Scribble::nt; ++ i)
        connect(q->pTheory[i], SIGNAL(stateChanged(int)), SLOT(slotForceUpdate()));

	startTimer(100);
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
            QRect r(q->planet[0][0].p[0] / Scribble::rmax * width()/2 - 5 + width()/2, q->planet[0][0].p[1] / Scribble::rmax * height()/2 - 5 + height()/2, 10, 10);
            QPainter painter;

            painter.begin( &buffer );
            painter.setBrush(Qt::yellow);
            painter.drawEllipse(r);
            painter.end();

            for (size_t j = 0; j < q->planet.size(); ++ j)
            {
                for (size_t i = 0; i < q->planet[j].size(); ++ i)
                {
                    QRect e(q->planet[j][i].o[0] / Scribble::rmax * width()/2 - 2 + width()/2, q->planet[j][i].o[1] / Scribble::rmax * height()/2 - 2 + height()/2, 4+1, 4+1);

                    q->planet[j][i].o = q->planet[j][i].p;

                    QRect r(q->planet[j][i].o[0] / Scribble::rmax * width()/2 - 2 + width()/2, q->planet[j][i].o[1] / Scribble::rmax * height()/2 - 2 + height()/2, 4, 4);
                    QPainter painter;

                    painter.begin( &buffer );
                    painter.setPen(q->planet[j][i].c);
                    painter.setBrush(q->planet[j][i].c);
                    painter.fillRect(e, palette().base());

                    if (q->pTheory[j]->checkState() == Qt::Checked)
                        painter.drawEllipse(r);

                    painter.end();
                }
            }

            update();
        }
        break;

    case 1:
        if (first)
        {
            first = false;

            QRect e(0, 0, width(), height());
            QPainter painter;

            painter.begin( &buffer );
            painter.fillRect(e, palette().base());
            painter.end();

            for (size_t j = 0; j < q->planet.size(); ++ j)
            {
                for (size_t i = 2; i < q->planet[j].size(); ++ i)
                {
                    QPoint r(q->planet[j][i].r / Scribble::rmax * width(), height() - q->planet[j][i].vel / q->vmax * height());
                    QPoint s(q->planet[j][i - 1].r / Scribble::rmax * width(), height() - q->planet[j][i - 1].vel / q->vmax * height());
                    QPainter painter;

                    painter.begin( &buffer );

                    if (q->pTheory[j]->checkState() == Qt::Checked)
                    {
                        if (q->planet[j][i].err == 0)
                        {
                            painter.setPen(QPen(q->planet[j][i].c, 2));
                            painter.drawLine(s.x(), s.y(), r.x(), r.y());
                        }
                        else
                        {
                            painter.setBrush(QBrush(q->planet[j][i].c, Qt::SolidPattern));
                            painter.drawEllipse(r, 4, 4);
                            painter.setPen(QPen(q->planet[j][i].c, 1));
                            painter.drawLine(r.x(), r.y() - q->planet[j][i].err / q->vmax * height(), r.x(), r.y() + q->planet[j][i].err / q->vmax * height());
                        }
                    }

                    painter.end();
                }
            }

            for (real x = 0; x < Scribble::rmax; x += Scribble::rmax / 5.0)
            {
                QPainter painter;

                painter.begin( &buffer );
                painter.setPen(Qt::black);
                painter.setBrush(Qt::black);
                painter.drawText(x / Scribble::rmax * width(), height() - painter.boundingRect(QRectF(), 0, QString::number(x)).height(), QString::number(x) + " kpc");
                painter.end();
            }

            for (real y = 0; y < q->vmax; y += q->vmax / 5.0)
            {
                QPainter painter;

                painter.begin( &buffer );
                painter.setPen(Qt::black);
                painter.setBrush(Qt::black);
                painter.drawText(0, height() - y / q->vmax * height(), QString::number(y) + " km/s");
                painter.end();
            }

            update();
        }

        break;
    };
}

void Canvas::slotForceUpdate()
{
    first = true;
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
    first = true;

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

char * const Scribble::theory[Scribble::nt] = {"Initial", "Visible", "Dark Matter", "Total", "Finite Theory", "Observed"};

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
    pTime->setDecimals(abs(log10(ntime)) + 1);
    pTime->setSingleStep(ntime / (abs(log10(ntime)) + 1));
    pTime->setToolTip("Time Interval (s)");
    pTime->setValue(ntime);

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
    pTheory[5]->setCheckState(Qt::Checked);

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
        planet[j].reserve(np);

        for (int i = 0; i < np; ++ i)
            planet[j].push_back(Planet("Star", QColor(Qt::GlobalColor(Qt::darkRed + j))));
    }

    for (size_t i = 0; i < np; ++ i)
    {
        planet[0][i].alpha = PI / 2 * ceil(dis(gen) * 4.0);
        planet[0][i].r = h * sqrt(tan(PI * 0.5 * emax * i / np));
        planet[0][i].vel = v0 * 2.0 / PI * atan(planet[0][i].r / r0);
        planet[0][i].omega = planet[0][i].vel / planet[0][i].r;
        planet[0][i].mass = planet[0][i].vel * planet[0][i].vel * planet[0][i].r;
    }

    real const totalmass = planet[0][np - 1].vel * planet[0][np - 1].vel * planet[0][np - 1].r;
    real const starmass = totalmass / (8.0 * np);

    for (size_t i = 0; i < np; ++ i)
    {
        planet[1][i].alpha = planet[0][i].alpha;
        planet[1][i].r = planet[0][i].r;
        planet[1][i].vel = sqrt(i * starmass * massf / planet[1][i].r);
        planet[1][i].omega = planet[1][i].vel / planet[1][i].r;
        planet[1][i].mass = planet[1][i].vel * planet[1][i].vel * planet[1][i].r;
    }

    real const md = totalmass - starmass * massf * np;
    real const mdk = md * dmf / (planet[0][np - 1].r / rdm0 - atan(planet[0][np - 1].r / rdm0));

    for (size_t i = 0; i < np; ++ i)
    {
        planet[2][i].alpha = planet[0][i].alpha;
        planet[2][i].r = planet[0][i].r;
        planet[2][i].vel = sqrt(mdk * (planet[2][i].r / rdm0 - atan(planet[2][i].r / rdm0)) / planet[2][i].r);
        planet[2][i].omega = planet[2][i].vel / planet[2][i].r;
        planet[2][i].mass = planet[2][i].vel * planet[2][i].vel * planet[2][i].r;
    }

    for (size_t i = 0; i < np; ++ i)
    {
        planet[3][i].alpha = planet[0][i].alpha;
        planet[3][i].r = planet[0][i].r;
        planet[3][i].vel = sqrt(planet[1][i].vel * planet[1][i].vel + planet[2][i].vel * planet[2][i].vel);
        planet[3][i].omega = planet[3][i].vel / planet[3][i].r;
        planet[3][i].mass = planet[1][i].mass + planet[2][i].mass;
    }

    for (size_t i = 0; i < np; ++ i)
    {
        planet[4][i].alpha = planet[0][i].alpha;
        planet[4][i].r = planet[0][i].r;
        planet[4][i].vel = planet[1][i].vel + sping * planet[4][i].r;
        planet[4][i].omega = planet[4][i].vel / planet[4][i].r;
        planet[4][i].mass = planet[4][i].vel * planet[4][i].vel * planet[4][i].r;
    }

    ifstream f("../finite-theory/data/rc-ngc_2403.dat");

    planet[5].clear();
    planet[5].reserve(1000);

    for (size_t i = 0; f; ++ i)
    {
        planet[5].push_back(Planet("Star", Qt::red));

        f >> planet[5][i].r >> planet[5][i].vel >> planet[5][i].err;

        planet[5][i].alpha = planet[0][i].alpha;
        planet[5][i].omega = planet[5][i].vel / planet[5][i].r;
        planet[5][i].mass = planet[5][i].vel * planet[5][i].vel * planet[5][i].r;
    }

    planet[5].pop_back();

    vmax = 0.0;

    for (size_t j = 0; j < planet.size(); ++ j)
        for (size_t i = 0; i < planet[j].size(); ++ i)
            if (vmax < planet[j][i].vel)
                vmax = planet[j][i].vel;

    new Dual(this);
}

Scribble::~Scribble()
{
    exit(-1);
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
    scribble.resize( 500, 360 );
    scribble.showMaximized();
	
    return a.exec();
}
