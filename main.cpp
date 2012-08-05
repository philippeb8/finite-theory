/**
    Finite Theory Simulator
    Copyright (C) 2011 Phil Bouchard <philippeb8@gmail.com>

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

#define EDITION "3.3"

#include "main.h"

//#include <unistd.h>
#include <stdlib.h>
//#include <sys/sysinfo.h>

#include <cmath>
#include <limits>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <typeinfo>

#include <qapplication.h>
#include <qevent.h>
#include <qpainter.h>
#include <q3toolbar.h>
#include <qtoolbutton.h>
#include <qspinbox.h>
#include <qtooltip.h>
#include <qrect.h>
#include <qpoint.h>
#include <qcolordialog.h>
#include <q3filedialog.h>
#include <qcursor.h>
#include <qimage.h>
#include <q3strlist.h>
#include <q3popupmenu.h>
#include <q3intdict.h>
#include <qlayout.h>
#include <qobject.h>
#include <qmenubar.h>
#include <qmessagebox.h>
#include <qstylefactory.h>
#include <qapplication.h>
//Added by qt3to4:
#include <Q3HBoxLayout>
#include <QTimerEvent>
#include <QLabel>
#include <QPixmap>
#include <QResizeEvent>
#include <QMouseEvent>
#include <Q3VBoxLayout>
#include <QPaintEvent>
#include <Q3MemArray>
#include <Q3Painter>
#include <QDesktopWidget>

using namespace std;

const real scale = 1e9L;
const real upper = 30.L;

inline real Planet::FR(real m, real d)
{
	return 1.L / ((m / d + h) / h);
}

inline real Planet::GR(real m, real d)
{
	return 1.L / ((m / d + h) / h);
}

inline real Planet::NW(real, real)
{
	return 1.L;
}

inline void Planet::operator () (const vector<Planet> &planet, const real & upper)
{
	vector3 va(0.L, 0.L, 0.L);
	
	for (size_t i = 0; i < planet.size(); i ++)
	{
		if (planet[i].m == m || planet[i].m == 0.0L)
			continue;
		
		const vector3 normal(p[0] - planet[i].p[0], p[1] - planet[i].p[1], p[2] - planet[i].p[2]);
		const real norm2 = pow(normal[0], 2) + pow(normal[1], 2) + pow(normal[2], 2);
		const real norm = sqrt(norm2);
		
		// a = Gm/r^2 decomposed in scalar
		va[0] += - G * planet[i].m / norm2 * normal[0] / norm;
		va[1] += - G * planet[i].m / norm2 * normal[1] / norm;
		va[2] += - G * planet[i].m / norm2 * normal[2] / norm;

		if (i == 0 && pd > norm)
		{
			pd = norm;

			pp[0] = p;
		}
	}
	
	// save
	const vector3 q(p[0], p[1], p[2]);
	
	// sun - planet
	const real nnorm_p = p.norm();
	vector3 vnn = p / nnorm_p;

	if (t[1] != 0.L)
	{
		t[1] = t[0];
		t[0] = upper * f(planet[0].m, nnorm_p);
	}
	else
	{
		t[0] = upper * f(planet[0].m, nnorm_p);
		t[1] = t[0];
	}

	vector3 vi(v * t[0] + va * (t[0] * t[0] / 2.));

	p += vi;
	v += va * t[0];

	switch (eType)
	{
	case PP:
		if (q[1] < 0 && p[1] >= 0)
		{
			//const real mu = - q[1] / (p[1] - q[1]);
		
			//i = q + (p - q) * mu;
			pp[1] = pp[0];
			
			pd = std::numeric_limits<real>::max();

			updated = true;
		}
		break;

	case LB:
		if (q[0] >= -200000000000.L && p[0] < -200000000000.L)
		{
			//const real mu = - q[0] / (p[0] - q[0]);
		
			//i = q + (p - q) * mu;
			pp[1] = pp[0];
			
			pd = std::numeric_limits<real>::max();

			updated = true;
		}
		break;
	}
}

Dual::Dual(Canvas * pParent, int id) : p(pParent), i(id)
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
			
		for (size_t j = 0; j < p->planet.size(); j ++)
			p->planet[j][i](p->planet[j], q->pTime->value());
	}
}

const bool no_writing = FALSE;

Canvas::Canvas( Type eType, QWidget *parent, const char *name )
    : QWidget( parent, name, Qt::WStaticContents ), 
      eType(eType), pen( Qt::red, 3 ), polyline(3), mousePressed( FALSE ), buffer( width(), height() )
{
	setAttribute(Qt::WA_PaintOutsidePaintEvent, true);
	
	Scribble * q = static_cast<Scribble *>(topLevelWidget());
	
	static const real pos[12][3] =
	{
		{0.L, 0.L, 0.L},
		{-57025548112.2453L, 3197006916.08582L, 5283916036.50742L},
		{26317771130.7392L, 105373484164.43L, 481049442.321637L}, 
		{-40584904469.4072L, -146162841483.741L, 582517208.913105L}, 
		{192608888576.284L, -72078449728.0548L, -5537406864.12226L}, 
		{-230068941192.889L, -766153804794.071L, 9039825087.87588L}, 
		{1359034179077.08L, -555461097003.149L, -48376702948.4567L}, 
		{1905563957085.85L, 2247912953966.77L, -16532490448.7952L}, 
		{1788083649521.39L, 4079380837677.57L, -125881827325.591L}, 
		{-4043923627184.17L, 3575690969311.01L, 795204553555.504L}, 
		
		{250000000000.L, -40000000000.L, 0.L}, 
		{250000000000.L, -40000000000.L, 0.L}
	};
	
	static const real vel[12][3] =
	{
		{0.L, 0.L, 0.L},
		{-13058.0445420602L, -46493.5791091285L, -2772.42900405547L},
		{-33720.199494784L, 8727.97495192353L, 2044.70922687897L},
		{28173.5639447033L, -8286.58463896112L, 13.3258392757908L},
		{9453.24519302534L, 24875.9047777036L, 333.149595901334L},
		{12310.4853583322L, -3126.10777330552L, -250.842129088533L},
		{3203.18660260855L, 8810.22721786771L, -260.876357307397L},
		{-5198.23543233994L, 4090.32678482699L, 78.6156634354517L},
		{-5005.88142339012L, 2215.0599004751L, 70.452880649377L},
		{-2122.7269723267L, -4538.25658137665L, 1101.51599904528L},
		
		{-300000.L, 0.L, 0.L}, 
		{-300000.L, 0.L, 0.L}
	};

	static const Planet Sun 	  ("Sun", 		Qt::yellow, 1.98911e30L, pos[0], vel[0]);
	static const Planet Mercury   ("Mercury", 	Qt::red, 3.302E+23L, pos[1], vel[1]);
	static const Planet Venus 	  ("Venus", 	Qt::cyan, 4.8685E+24L, pos[2], vel[2]);
	static const Planet Earth 	  ("Earth", 	Qt::blue, 5.9736E+24L, pos[3], vel[3]);
	static const Planet Mars 	  ("Mars", 		Qt::yellow, 6.41850000000001E+23L, pos[4], vel[4]);
	static const Planet Jupiter   ("Jupiter", 	Qt::magenta, 1.8986E+27L, pos[5], vel[5]);
	static const Planet Saturn 	  ("Saturn", 	Qt::darkRed, 5.6842928E+26L, pos[6], vel[6]);
	static const Planet Uranus 	  ("Uranus", 	Qt::green, 8.68320000000002E+25L, pos[7], vel[7]);
	static const Planet Neptune   ("Neptune", 	Qt::darkBlue, 1.0243E+26L, pos[8], vel[8]);
	static const Planet Pluto 	  ("Pluto", 	Qt::darkGray, 1.27E+22L, pos[9], vel[9]);

	static const Planet Photon1   ("Photon1", 	Qt::darkGreen, 0.0L, pos[10], vel[10], Planet::FR, Planet::LB);
	static const Planet Photon2   ("Photon2", 	Qt::darkRed, 0.0L, pos[11], vel[11], Planet::NW, Planet::LB);

	switch (eType)
	{
	case PP:
		planet.resize(2);
		planet[0].reserve(10);
		planet[0].push_back(Sun);
		planet[0].push_back(Mercury);
		planet[0].push_back(Venus);
		planet[0].push_back(Earth);
		planet[0].push_back(Mars);
		planet[0].push_back(Jupiter);
		planet[0].push_back(Saturn);
		planet[0].push_back(Uranus);
		planet[0].push_back(Neptune);
		planet[0].push_back(Pluto);
		
		//for (int i = 0; i < 9 - 2 /*get_nprocs()*/ * 4; i ++)
		//	planet[0].pop_back();

		planet[1] = planet[0];
		for (size_t i = 0; i < planet[1].size(); i ++)
		{
			planet[1][i].f = Planet::FR;
			planet[1][i].c = planet[1][i].c.dark();
		}

		stats.resize(planet[0].size());

		for (size_t i = 1; i < planet[0].size(); i ++)
		{
			QPixmap p(12, 8);
			p.fill(planet[0][i].c);
			q->pPlanet->insertItem(p, planet[0][i].n);
		}

		connect(q->pPlanet, SIGNAL(activated(int)), SLOT(slotPlanet(int)));
		break;
	
	case LB:
		planet.resize(2);
		planet[0].reserve(10);
		planet[0].push_back(Sun);
		planet[0].push_back(Photon1);
		planet[1].push_back(Sun);
		planet[1].push_back(Photon2);

		stats.resize(planet[0].size());
		break;
	}
	
    if ((qApp->argc() > 0) && !buffer.load(qApp->argv()[1]))
        buffer.fill( colorGroup().base() );
    setBackgroundMode( Qt::PaletteBase );
#ifndef QT_NO_CURSOR
    setCursor( Qt::crossCursor );
#endif

	startTimer(100);
	//QApplication::postEvent(this, new QTimerEvent(0));

	for (size_t i = 1; i < planet[0].size(); i ++)
		new Dual(this, i);
}

Canvas::~Canvas()
{
	exit(-1);
}

void Canvas::slotPlanet(int i)
{
	Scribble * p = static_cast<Scribble *>(topLevelWidget());
	
	++ i;
	
	// precession
	for (size_t j = 0; j < planet.size(); ++ j)
		for (int x = 0; x < 3; ++ x)
		{
			ostringstream s;
			
			if (stats[i].mean[0].size() > 0)
			{
				s.setf(ios::scientific, ios::floatfield);
				s << std::setprecision(numeric_limits<real>::digits10);
				s << stats[i].precession[j][x];
			}

			p->pLabel[eType][j][x]->setText(s.str().c_str());
		}

	// anomaly
	for (int x = 0; x < 3; ++ x)
	{
		ostringstream s;
		
		if (stats[i].mean[0].size() > 0)
		{
			s.setf(ios::scientific, ios::floatfield);
			s << std::setprecision(numeric_limits<real>::digits10);
			s << stats[i].precession[1][x] - stats[i].precession[0][x];
		}

		p->pLabel[eType][2][x]->setText(s.str().c_str());
	}
	
	// mean
	for (int x = 0; x < 3; ++ x)
	{
		ostringstream s[4];
		
		// median
		if (stats[i].mean[0].size() > 0)
		{
			s[0].setf(ios::scientific, ios::floatfield);
			s[0] << std::setprecision(numeric_limits<real>::digits10);
			
			set<real>::iterator k = stats[i].mean[x].begin();
			set<real>::reverse_iterator l = stats[i].mean[x].rbegin();
			advance(k, stats[i].mean[x].size() / 2);
			advance(l, stats[i].mean[x].size() / 2);
			
			const real median = (* k + * l) / 2;
			
			s[0] << median;

			// median absolute deviation
			if (stats[i].mean[0].size() > 1)
			{
				s[1].setf(ios::scientific, ios::floatfield);
				s[1] << std::setprecision(numeric_limits<real>::digits10);
				
				set<real> dev;
				for (set<real>::iterator m = stats[i].mean[x].begin(); m != stats[i].mean[x].end(); ++ m)
					dev.insert(abs(* m - median));
					
				set<real>::iterator m = dev.begin();
				set<real>::reverse_iterator n = dev.rbegin();
				advance(m, dev.size() / 2);
				advance(n, dev.size() / 2);
				
				const real mad = (* m + * n) / 2;
				
				s[1] << mad;
				
				if (stats[i].best[1][x] > mad)
				{
					stats[i].best[0][x] = median;
					stats[i].best[1][x] = mad;
				}
				
				s[2].setf(ios::scientific, ios::floatfield);
				s[2] << std::setprecision(numeric_limits<real>::digits10);
				
				s[3].setf(ios::scientific, ios::floatfield);
				s[3] << std::setprecision(numeric_limits<real>::digits10);
					
				s[2] << stats[i].best[0][x];
				s[3] << stats[i].best[1][x];
			}
		}

		p->pLabel[eType][3][x]->setText(s[0].str().c_str());
		p->pLabel[eType][4][x]->setText(s[1].str().c_str());
		p->pLabel[eType][5][x]->setText(s[2].str().c_str());
		p->pLabel[eType][6][x]->setText(s[3].str().c_str());
	}
}

void Canvas::timerEvent(QTimerEvent *)
{
	Scribble * q = static_cast<Scribble *>(topLevelWidget());

	if (! isVisible())
		return;

	{
		QRect r(planet[0][0].p[0] / scale - 5 + width()/2, planet[0][0].p[1] / scale - 5 + height()/2, 10, 10);
		Q3Painter painter;
		painter.begin( &buffer );
		painter.setBrush(Qt::yellow);
		painter.drawEllipse(r);
		painter.end();
		bitBlt( this, r.x(), r.y(), &buffer, r.x(), r.y(), r.width(), r.height() );

		update(r);
	}

	for (size_t i = 1; i < planet[0].size(); ++ i)
	{
		for (size_t j = 0; j < planet.size(); ++ j)
			if (planet[j][i].updated)
			{
				planet[j][i].ps[1][0] = planet[j][i].ps[0][0];
				planet[j][i].ps[1][1] = planet[j][i].ps[0][1];
				planet[j][i].ps[1][2] = planet[j][i].ps[0][2];

				planet[j][i].ps[0][0] = sqrt(pow(planet[j][i].pp[1][0], 2) + pow(planet[j][i].pp[1][1], 2) + pow(planet[j][i].pp[1][2], 2));
				planet[j][i].ps[0][1] = atan2(planet[j][i].pp[1][1], planet[j][i].pp[1][0]);
				planet[j][i].ps[0][2] = acos(planet[j][i].pp[1][2] / planet[j][i].ps[0][0]);
				
				for (int x = 0; x < 3; ++ x)
					stats[i].precession[j][x] = planet[j][i].ps[1][x] - planet[j][i].ps[0][x];
			}

		if (planet[0][i].updated && planet[1][i].updated)
		{
			for (int x = 0; x < 3; ++ x)
				stats[i].mean[x].insert(stats[i].precession[1][x] - stats[i].precession[0][x]);
			
			if (size_t(q->pPlanet->currentItem() + 1) == i)
				slotPlanet(i - 1);

			planet[0][i].updated = false;
			planet[1][i].updated = false;
		}

		for (size_t j = 0; j < planet.size(); ++ j)
		{
			QRect e(planet[j][i].o[0] / scale - 2 + width()/2, planet[j][i].o[1] / scale - 2 + height()/2, 4, 4);

			planet[j][i].o[0] = planet[j][i].p[0];
			planet[j][i].o[1] = planet[j][i].p[1];
			planet[j][i].o[2] = planet[j][i].p[2];
			
			QRect r(planet[j][i].o[0] / scale - 2 + width()/2, planet[j][i].o[1] / scale - 2 + height()/2, 4, 4);
			Q3Painter painter;
			painter.begin( &buffer );
			painter.setPen(planet[j][i].c);
			painter.setBrush(planet[j][i].c);
			painter.eraseRect(e);
			painter.drawEllipse(r);
			painter.end();
			r |= e;
			bitBlt( this, r.x(), r.y(), &buffer, r.x(), r.y(), r.width(), r.height() );

			update(r);
		}
	}
}

void Canvas::clearScreen()
{
    buffer.fill( colorGroup().base() );
    repaint( FALSE );
}

void Canvas::mousePressEvent( QMouseEvent *e )
{
    mousePressed = TRUE;
    polyline[2] = polyline[1] = polyline[0] = e->pos();
}

void Canvas::mouseReleaseEvent( QMouseEvent * )
{
    mousePressed = FALSE;
}

void Canvas::mouseMoveEvent( QMouseEvent *e )
{
    if ( mousePressed ) {
        Q3Painter painter;
        painter.begin( &buffer );
        painter.setPen( Qt::gray );
        polyline[2] = polyline[1];
        polyline[1] = polyline[0];
        polyline[0] = e->pos();
        painter.drawPolyline( polyline );
        painter.end();

        QRect r = polyline.boundingRect();
        r = r.normalize();
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

    QPixmap tmp( buffer );
    buffer.resize( w, h );
    buffer.fill( colorGroup().base() );
    bitBlt( &buffer, 0, 0, &tmp, 0, 0, tmp.width(), tmp.height() );
}

void Canvas::paintEvent( QPaintEvent *e )
{
    QWidget::paintEvent( e );

    Q3MemArray<QRect> rects = e->region().rects();
    for ( uint i = 0; i < rects.count(); ++ i ) 
	{
        QRect r = rects[(int)i];
        bitBlt( this, r.x(), r.y(), &buffer, r.x(), r.y(), r.width(), r.height() );
    }
}

//------------------------------------------------------

Scribble::Scribble( QWidget *parent, const char *name )
    : Q3MainWindow( parent, name ), nc(0)
{
	ntime[0] = upper;
	ntime[1] = 1;
	
	Q3PopupMenu *file = new Q3PopupMenu( this );
	file->insertItem( "E&xit", qApp, SLOT(quit()), Qt::CTRL+Qt::Key_Q );

	Q3PopupMenu *view = new Q3PopupMenu( this );
	view->insertItem( "&Perihelion Precession", this, SLOT(slotPP()));
	view->insertItem( "&Bending of Light", this, SLOT(slotLB()));

	Q3PopupMenu *help = new Q3PopupMenu( this );
	help->insertItem( "&About", this, SLOT(slotAbout()));

	QMenuBar * menu = new QMenuBar( this );
	menu->insertItem( "&File", file );
	menu->insertItem( "&View", view );
	menu->insertSeparator();
	menu->insertItem( "&Help", help );
	menu->setSeparator( QMenuBar::InWindowsStyle );

    Q3ToolBar *tools = new Q3ToolBar( this );

    bClear = new QToolButton( QIcon(), "Clear Screen", "Clear Screen", this, SLOT( slotClear() ), tools );
    bClear->setText( "Clear Screen" );

    tools->addSeparator();

    pTime = new QSpinBox( 1, 36000, 10, tools );
    QToolTip::add( pTime, "Time Interval (s)" );
    pTime->setValue( ntime[nc] );

    tools->addSeparator();

	pPlanet = new QComboBox(tools);
    QToolTip::add( pPlanet, "Planet" );
	connect(pPlanet, SIGNAL(activated(int)), SLOT(slotPlanet(int)));

	
	pTabWidget = new QTabWidget(this);
	connect(pTabWidget, SIGNAL(currentChanged(QWidget *)), SLOT(slotChanged(QWidget *)));

	for (unsigned i = 0; i < 2; ++ i)
	{
		pTab[i] = new QWidget(pTabWidget);

		pLabel[i][0][0] = new QLabel(pTab[i]);
		pLabel[i][0][1] = new QLabel(pTab[i]);
		pLabel[i][0][2] = new QLabel(pTab[i]);
		pLabel[i][1][0] = new QLabel(pTab[i]);
		pLabel[i][1][1] = new QLabel(pTab[i]);
		pLabel[i][1][2] = new QLabel(pTab[i]);
		pLabel[i][2][0] = new QLabel(pTab[i]);
		pLabel[i][2][1] = new QLabel(pTab[i]);
		pLabel[i][2][2] = new QLabel(pTab[i]);
		pLabel[i][3][0] = new QLabel(pTab[i]);
		pLabel[i][3][1] = new QLabel(pTab[i]);
		pLabel[i][3][2] = new QLabel(pTab[i]);
		pLabel[i][4][0] = new QLabel(pTab[i]);
		pLabel[i][4][1] = new QLabel(pTab[i]);
		pLabel[i][4][2] = new QLabel(pTab[i]);
		pLabel[i][5][0] = new QLabel(pTab[i]);
		pLabel[i][5][1] = new QLabel(pTab[i]);
		pLabel[i][5][2] = new QLabel(pTab[i]);
		pLabel[i][6][0] = new QLabel(pTab[i]);
		pLabel[i][6][1] = new QLabel(pTab[i]);
		pLabel[i][6][2] = new QLabel(pTab[i]);
		
	    QToolTip::add(pLabel[i][0][0], QString("Newton ") + QChar(0x0394) + QChar(0x03C1));
	    QToolTip::add(pLabel[i][0][1], QString("Newton ") + QChar(0x0394) + QChar(0x03C6));
	    QToolTip::add(pLabel[i][0][2], QString("Newton ") + QChar(0x0394) + QChar(0x03B8));
	    QToolTip::add(pLabel[i][1][0], QString("Finite Theory ") + QChar(0x0394) + QChar(0x03C1));
	    QToolTip::add(pLabel[i][1][1], QString("Finite Theory ") + QChar(0x0394) + QChar(0x03C6));
	    QToolTip::add(pLabel[i][1][2], QString("Finite Theory ") + QChar(0x0394) + QChar(0x03B8));
	    QToolTip::add(pLabel[i][2][0], QString("Finite Theory ") + QChar(0x0394) + QChar(0x03C1) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C1));
	    QToolTip::add(pLabel[i][2][1], QString("Finite Theory ") + QChar(0x0394) + QChar(0x03C6) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C6));
	    QToolTip::add(pLabel[i][2][2], QString("Finite Theory ") + QChar(0x0394) + QChar(0x03B8) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03B8));
	    QToolTip::add(pLabel[i][3][0], QChar(0x03BC) + QString("(Finite Theory ") + QChar(0x0394) + QChar(0x03C1) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C1) + QString(")"));
	    QToolTip::add(pLabel[i][3][1], QChar(0x03BC) + QString("(Finite Theory ") + QChar(0x0394) + QChar(0x03C6) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C6) + QString(")"));
	    QToolTip::add(pLabel[i][3][2], QChar(0x03BC) + QString("(Finite Theory ") + QChar(0x0394) + QChar(0x03B8) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03B8) + QString(")"));
	    QToolTip::add(pLabel[i][4][0], QString("MAD(Finite Theory ") + QChar(0x0394) + QChar(0x03C1) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C1) + QString(")"));
	    QToolTip::add(pLabel[i][4][1], QString("MAD(Finite Theory ") + QChar(0x0394) + QChar(0x03C6) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C6) + QString(")"));
	    QToolTip::add(pLabel[i][4][2], QString("MAD(Finite Theory ") + QChar(0x0394) + QChar(0x03B8) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03B8) + QString(")"));
	    QToolTip::add(pLabel[i][5][0], QString("Best of ") + QChar(0x03BC) + QString("(Finite Theory ") + QChar(0x0394) + QChar(0x03C1) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C1) + QString(")"));
	    QToolTip::add(pLabel[i][5][1], QString("Best of ") + QChar(0x03BC) + QString("(Finite Theory ") + QChar(0x0394) + QChar(0x03C6) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C6) + QString(")"));
	    QToolTip::add(pLabel[i][5][2], QString("Best of ") + QChar(0x03BC) + QString("(Finite Theory ") + QChar(0x0394) + QChar(0x03B8) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03B8) + QString(")"));
	    QToolTip::add(pLabel[i][6][0], QString("Best of ") + QString("MAD(Finite Theory ") + QChar(0x0394) + QChar(0x03C1) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C1) + QString(")"));
	    QToolTip::add(pLabel[i][6][1], QString("Best of ") + QString("MAD(Finite Theory ") + QChar(0x0394) + QChar(0x03C6) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C6) + QString(")"));
	    QToolTip::add(pLabel[i][6][2], QString("Best of ") + QString("MAD(Finite Theory ") + QChar(0x0394) + QChar(0x03B8) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03B8) + QString(")"));

		canvas[i] = new Canvas((Canvas::Type) (i), pTab[i]);
		canvas[i]->setSizePolicy(QSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding));

		Q3BoxLayout * l = new Q3VBoxLayout(pTab[i]);
	    l->addWidget( canvas[i] );

		Q3BoxLayout * h1 = new Q3HBoxLayout(l);
		h1->addWidget(pLabel[i][0][0], 1);
		h1->addWidget(pLabel[i][0][1], 1);
		h1->addWidget(pLabel[i][0][2], 1);
		Q3BoxLayout * h2 = new Q3HBoxLayout(l);
		h2->addWidget(pLabel[i][1][0], 1);
		h2->addWidget(pLabel[i][1][1], 1);
		h2->addWidget(pLabel[i][1][2], 1);
		Q3BoxLayout * h3 = new Q3HBoxLayout(l);
		h3->addWidget(pLabel[i][2][0], 1);
		h3->addWidget(pLabel[i][2][1], 1);
		h3->addWidget(pLabel[i][2][2], 1);
		Q3BoxLayout * h4 = new Q3HBoxLayout(l);
		h4->addWidget(pLabel[i][3][0], 1);
		h4->addWidget(pLabel[i][3][1], 1);
		h4->addWidget(pLabel[i][3][2], 1);
		Q3BoxLayout * h5 = new Q3HBoxLayout(l);
		h5->addWidget(pLabel[i][4][0], 1);
		h5->addWidget(pLabel[i][4][1], 1);
		h5->addWidget(pLabel[i][4][2], 1);
		Q3BoxLayout * h6 = new Q3HBoxLayout(l);
		h6->addWidget(pLabel[i][5][0], 1);
		h6->addWidget(pLabel[i][5][1], 1);
		h6->addWidget(pLabel[i][5][2], 1);
		Q3BoxLayout * h7 = new Q3HBoxLayout(l);
		h7->addWidget(pLabel[i][6][0], 1);
		h7->addWidget(pLabel[i][6][1], 1);
		h7->addWidget(pLabel[i][6][2], 1);

		pLabel[i][5][1]->setPaletteForegroundColor(Qt::darkRed);
	}
	
	pTabWidget->addTab(pTab[0], "Perihelion Precession");
	pTabWidget->addTab(pTab[1], "Light Bending");
	//pTab[1]->hide();

    setCentralWidget( pTabWidget );
}

void Scribble::slotClear()
{
    canvas[nc]->clearScreen();
}

void Scribble::slotPlanet(int)
{
	pLabel[0][0][0]->clear();
	pLabel[0][0][1]->clear();
	pLabel[0][0][2]->clear();
	pLabel[0][1][0]->clear();
	pLabel[0][1][1]->clear();
	pLabel[0][1][2]->clear();
}

void Scribble::slotPP()
{
    ntime[nc] = pTime->value();

	nc = 0;
	pPlanet->setEnabled(true);
	pTabWidget->showPage(pTab[nc]);
    pTime->setValue( ntime[nc] );
}

void Scribble::slotLB()
{
    ntime[nc] = pTime->value();

	nc = 1;
	pPlanet->setEnabled(false);
	pTabWidget->showPage(pTab[nc]);
    pTime->setValue( ntime[nc] );
}

void Scribble::slotChanged(QWidget * p)
{
	if (p == pTab[0])
		slotPP();
	else
		slotLB();
}

void Scribble::slotAbout()
{
    QMessageBox::about( this, "Finite Theory of the Universe", "\nCopyright (c) 2011\n\nPhil Bouchard <philippeb8@gmail.com>\n");
}
	

int main( int argc, char **argv )
{
    QApplication a( argc, argv );

    Scribble scribble;

    scribble.resize( 500, 360 );
    scribble.setCaption("Finite Theory of the Universe " EDITION);
	a.setStyle("windows");
    a.setMainWidget( &scribble );
	
    if ( QApplication::desktop()->width() > 550 && QApplication::desktop()->height() > 366 )
        scribble.show();
    else
        scribble.showMaximized();
	
    return a.exec();
}
