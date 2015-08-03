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

#define EDITION "4.20"

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

const real scale[] = {8e8, 8e9L, 8e9L, 8e9L, 8e10L};
const real upper = 10.L;

// FT time formula
// observer is infinitly far away
real Planet::FR1(real m, real d, real h)
{
    return (h) / (m / d + h);
}

// observer is in a null environment
real Planet::FR2(real m, real d, real h)
{
    return (h + 1.L) / (m / d + h);
}

// Newton time formula
real Planet::NW(real, real, real)
{
	return 1.L;
}


void bitBlt( QPaintDevice * dst, int x, int y, const QPixmap* src, int sx, int sy, int sw, int sh )
{
    QPainter p( dst );
    p.drawPixmap( x, y, *src, sx, sy, sw, sh );
}


/** 
	@brief			Calculates the next position of the planet or photon
	@param planet	Planets that will affect the movement of the planet that is moving (this)
	@param upper	Time interval
*/

inline void Planet::operator () (const vector<Planet> &planet, const real & upper)
{
	// net acceleration vector (with all planets)
    vector3 va(0.L, 0.L, 0.L);
	
	// iterate through all planets
	for (size_t i = 0; i < planet.size(); i ++)
	{
		// if same planet or photon then skip
		if (planet[i].m == m || planet[i].m == 0.0L)
			continue;
		
		// vector and norm between the moving planet and the other one
		const vector3 normal(p[0] - planet[i].p[0], p[1] - planet[i].p[1], p[2] - planet[i].p[2]);
		const real norm2 = pow(normal[0], 2) + pow(normal[1], 2) + pow(normal[2], 2);
		const real norm = sqrt(norm2);
		
		// a = Gm/r^2 decomposed in scalar
        va[0] += - G * planet[i].m / norm2 * normal[0] / norm;
        va[1] += - G * planet[i].m / norm2 * normal[1] / norm;
        va[2] += - G * planet[i].m / norm2 * normal[2] / norm;

		// save position of the planet when perihelion is found
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

	// if this not the first time
	if (t[1] != 0.L)
	{
		// save old time value
		t[1] = t[0];

		// Newton: t = upper
		// FT: t = upper / ((m / d + h) / h)
        t[0] = upper * f(planet[0].m, nnorm_p, planet[0].h);
	}
	else
	{
		// Newton: t = upper
		// FT: t = upper / ((m / d + h) / h)
        t[0] = upper * f(planet[0].m, nnorm_p, planet[0].h);

		// copy time value
		t[1] = t[0];
	}

    vector3 vi(v[0] * t[0] + va * (t[0] * t[0] / 2.));

	// p = p + v*t + (a*t^2)/2
	p += vi;

    // v = v + a*t
    v[0] += va * t[0];
	
	switch (eType)
	{
	// perihelion precession disparity
	case PP:
		// if the planet crossed the x axis (y = 0)
		if (q[1] < 0 && p[1] >= 0)
		{
			pp[1] = pp[0];
			
 			pd = std::numeric_limits<real>::max();

            updated = true;
		}
		break;

	// gravitational light bending
	case LB:
		// if the photon crossed x = -200000000000
		if (q[0] >= -200000000000.L && p[0] < -200000000000.L)
		{
			pp[1] = pp[0];
			
 			pd = std::numeric_limits<real>::max();

            updated = true;
		}
		break;

    // big bang & pioneer 10
    case BB:
    case V1:
        if (floor(q[0] / 1e10) != floor(p[0] / 1e10))
        {
            v[1] = v[0];

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
		// stop the processing until the tab becomes visible
		while (! p->isVisible())
			QThread::msleep(100);
		
		// move the same planet or photon according to Newton & FT
		for (size_t j = 0; j < p->planet.size(); j ++)
			p->planet[j][i](p->planet[j], q->pTime->value());
	}
}

const bool no_writing = false;

Canvas::Canvas( Type eType, QWidget *parent, real scale )
    : QWidget( parent/*, name, Qt::WStaticContents*/ ),
      eType(eType), pen( Qt::red, 3 ), polyline(3), mousePressed( false ), buffer( width(), height() ), scale(scale)
{
//	setAttribute(Qt::WA_PaintOutsidePaintEvent, true);
	
	Scribble * q = static_cast<Scribble *>(topLevelWidget());
	
	// initial position of each planet and photon
    static const real pos[][3] =
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
        {250000000000.L, -40000000000.L, 0.L},

        {-50000000000.L, 50000000000.L, 0.L},
        {0.L, 50000000000.L, 0.L},
        {50000000000.L, 50000000000.L, 0.L},
        {-50000000000.L, 0.L, 0.L},
        {50000000000.L, 0.L, 0.L},
        {-50000000000.L, -50000000000.L, 0.L},
        {0.L, -50000000000.L, 0.L},
        {50000000000.L, -50000000000.L, 0.L},

        {-4000000000000.L, 0.L, 0.L},
        {-3800000000000.L, 0.L, 0.L},
        {-300000000000.L, 0.L, 0.L},
        {-100000000000.L, 0.L, 0.L},
        {100000000000.L, 0.L, 0.L},
        {300000000000.L, 0.L, 0.L},
        {3800000000000.L, 0.L, 0.L},
        {4000000000000.L, 0.L, 0.L},

        {17048116800000.L, 0.L, 0.L},
        {17048116800000.L, 0.L, 0.L},
    };
	
	// initial velocity of each planet and photon
    static const real vel[][3] =
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
        {-300000.L, 0.L, 0.L},

        {-50000.L, 50000.L, 0.L},
        {0.L, 50000.L, 0.L},
        {50000.L, 50000.L, 0.L},
        {-80000L, 0.L, 0.L},
        {80000L, 0.L, 0.L},
        {-50000.L, -50000.L, 0.L},
        {0.L, -50000.L, 0.L},
        {50000.L, -50000.L, 0.L},

        {0.L, 5.7749e-6L, 0.L},
        {0.L, 5.9249e-6L, 0.L},
        {0.L, 2.1087e-5L, 0.L},
        {0.L, 3.6523e-5L, 0.L},
        {0.L, -3.6523e-5L, 0.L},
        {0.L, -2.1087e-5L, 0.L},
        {0.L, -5.9249e-6L, 0.L},
        {0.L, -5.7749e-6L, 0.L},

        {11992.L, 0.L, 0.L},
        {11992.L, 0.L, 0.L},
    };

	// name, color, mass, position and velocity of each moving object
	static const Planet Sun 	  ("Sun", 		Qt::yellow, 1.98911E+30L, pos[0], vel[0]);
	static const Planet Mercury   ("Mercury", 	Qt::red, 3.302E+23L, pos[1], vel[1]);
	static const Planet Venus 	  ("Venus", 	Qt::cyan, 4.8685E+24L, pos[2], vel[2]);
	static const Planet Earth 	  ("Earth", 	Qt::blue, 5.9736E+24L, pos[3], vel[3]);
	static const Planet Mars 	  ("Mars", 		Qt::yellow, 6.41850000000001E+23L, pos[4], vel[4]);
	static const Planet Jupiter   ("Jupiter", 	Qt::magenta, 1.8986E+27L, pos[5], vel[5]);
	static const Planet Saturn 	  ("Saturn", 	Qt::darkRed, 5.6842928E+26L, pos[6], vel[6]);
	static const Planet Uranus 	  ("Uranus", 	Qt::green, 8.68320000000002E+25L, pos[7], vel[7]);
	static const Planet Neptune   ("Neptune", 	Qt::darkBlue, 1.0243E+26L, pos[8], vel[8]);
	static const Planet Pluto 	  ("Pluto", 	Qt::darkGray, 1.27E+22L, pos[9], vel[9]);

    static const Planet Photon1   ("Photon1", 	Qt::darkGreen, 0.0L, pos[10], vel[10], Planet::FR1, Planet::LB);
    static const Planet Photon2   ("Photon2", 	Qt::darkRed, 0.0L, pos[11], vel[11], Planet::NW, Planet::LB);

    static const Planet Pioneer1   ("Pioneer1", 	Qt::darkGreen, 258.8L, pos[28], vel[28], Planet::FR1, Planet::V1);
    static const Planet Pioneer2   ("Pioneer2", 	Qt::darkRed, 258.8L, pos[29], vel[29], Planet::NW, Planet::V1);

    static const Planet Core	  ("Core", 		Qt::black, 2E+11L, pos[0], vel[0], Planet::NW, Planet::BB);
    static const Planet Galaxy1   ("Galaxy1", 	Qt::red, 50000L, pos[12], vel[12], Planet::NW, Planet::BB);
    static const Planet Galaxy2   ("Galaxy2", 	Qt::cyan, 50000L, pos[13], vel[13], Planet::NW, Planet::BB);
    static const Planet Galaxy3   ("Galaxy3", 	Qt::blue, 50000L, pos[14], vel[14], Planet::NW, Planet::BB);
    static const Planet Galaxy4   ("Galaxy4", 	Qt::yellow, 50000L, pos[15], vel[15], Planet::NW, Planet::BB);
    static const Planet Galaxy5   ("Galaxy5", 	Qt::magenta, 50000L, pos[16], vel[16], Planet::NW, Planet::BB);
    static const Planet Galaxy6   ("Galaxy6", 	Qt::darkRed, 50000L, pos[17], vel[17], Planet::NW, Planet::BB);
    static const Planet Galaxy7   ("Galaxy7", 	Qt::green, 50000L, pos[18], vel[18], Planet::NW, Planet::BB);
    static const Planet Galaxy8   ("Galaxy8", 	Qt::darkBlue, 50000L, pos[19], vel[19], Planet::NW, Planet::BB);

    static const Planet Nucleus   ("Nucleus", 	Qt::black, 2E+12L, pos[0], vel[0], Planet::NW, Planet::GR);
    static const Planet Star1     ("Star1", 	Qt::red, 50000L, pos[20], vel[20], Planet::NW, Planet::GR);
    static const Planet Star2     ("Star2", 	Qt::red, 50000L, pos[21], vel[21], Planet::NW, Planet::GR);
    static const Planet Star3     ("Star3", 	Qt::red, 50000L, pos[22], vel[22], Planet::NW, Planet::GR);
    static const Planet Star4     ("Star4", 	Qt::red, 50000L, pos[23], vel[23], Planet::NW, Planet::GR);
    static const Planet Star5     ("Star5", 	Qt::red, 50000L, pos[24], vel[24], Planet::NW, Planet::GR);
    static const Planet Star6     ("Star6", 	Qt::red, 50000L, pos[25], vel[25], Planet::NW, Planet::GR);
    static const Planet Star7     ("Star7", 	Qt::red, 50000L, pos[26], vel[26], Planet::NW, Planet::GR);
    static const Planet Star8     ("Star8", 	Qt::red, 50000L, pos[27], vel[27], Planet::NW, Planet::GR);

    switch (eType)
	{
	// perihelion precession disparity
	case PP:
		planet.resize(2);

		// store each planet using the Newton time formula
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
		
		// copy & change each planet for the FT time formula
		planet[1] = planet[0];
		for (size_t i = 0; i < planet[1].size(); i ++)
		{
            planet[1][i].f = Planet::FR1;
			planet[1][i].c = planet[1][i].c.dark();
		}

		stats.resize(planet[0].size());

		// prepare the UI
		for (size_t i = 1; i < planet[0].size(); i ++)
		{
			QPixmap p(12, 8);
			p.fill(planet[0][i].c);
            q->pPlanet[0]->addItem(p, planet[0][i].n);
		}

        connect(q->pPlanet[0], SIGNAL(activated(int)), SLOT(slotPlanet(int)));
		break;
	
	// gravitational light bending
	case LB:
		planet.resize(2);

		// store the Sun & the photon using the Newton time formula
		planet[0].reserve(2);
		planet[0].push_back(Sun);
		planet[0].push_back(Photon1);

		// store the Sun & the photon using the FT time formula
		planet[1].reserve(2);
		planet[1].push_back(Sun);
		planet[1].push_back(Photon2);

		stats.resize(planet[0].size());
		break;

    // big bang
    case BB:
        planet.resize(2);

        // store the Sun & the planets using FT time formula
        planet[0].reserve(9);
        planet[0].push_back(Core);
        planet[0].push_back(Galaxy1);
        planet[0].push_back(Galaxy2);
        planet[0].push_back(Galaxy3);
        planet[0].push_back(Galaxy4);
        planet[0].push_back(Galaxy5);
        planet[0].push_back(Galaxy6);
        planet[0].push_back(Galaxy7);
        planet[0].push_back(Galaxy8);

        // copy & change each planet for the FT time formula
        planet[1] = planet[0];
        for (size_t i = 0; i < planet[1].size(); i ++)
        {
            planet[1][i].f = Planet::Planet::FR2;
            planet[1][i].h = H[1];
            planet[1][i].c = planet[1][i].c.dark();
        }

        stats.resize(planet[0].size());

        // prepare the UI
        for (size_t i = 1; i < planet[0].size(); i ++)
        {
            QPixmap p(12, 8);
            p.fill(planet[0][i].c);
            q->pPlanet[1]->addItem(p, planet[0][i].n);
        }

        connect(q->pPlanet[1], SIGNAL(activated(int)), SLOT(slotGalaxy(int)));
        break;

        // galactic rotation
    case GR:
        planet.resize(2);

        // store the Sun & the planets using FT time formula
        planet[0].reserve(9);
        planet[0].push_back(Nucleus);
        planet[0].push_back(Star1);
        planet[0].push_back(Star2);
        planet[0].push_back(Star3);
        planet[0].push_back(Star4);
        planet[0].push_back(Star5);
        planet[0].push_back(Star6);
        planet[0].push_back(Star7);
        planet[0].push_back(Star8);

        // copy & change each planet for the FT time formula
        planet[1] = planet[0];
        for (size_t i = 0; i < planet[1].size(); i ++)
        {
            planet[1][i].f = Planet::Planet::FR2;
            planet[1][i].h = H[1];
            planet[1][i].c = planet[1][i].c.dark();
        }

        stats.resize(planet[0].size());
        break;

    // pioneer 10
    case V1:
        planet.resize(2);

        // store the Sun & the photon using the Newton time formula
        planet[0].reserve(2);
        planet[0].push_back(Sun);
        planet[0].push_back(Pioneer1);

        // store the Sun & the photon using the FT time formula
        planet[1].reserve(2);
        planet[1].push_back(Sun);
        planet[1].push_back(Pioneer2);

        stats.resize(planet[0].size());
        break;
    }
	
//    if ((qApp->argc() > 0) && !buffer.load(qApp->argv()[1]))
//        buffer.fill( palette().base().color() );
//    setBackgroundMode( Qt::PaletteBase );
#ifndef QT_NO_CURSOR
    setCursor( Qt::CrossCursor );
#endif

	startTimer(100);

	// launch a thread for each planet or photon
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
	
    switch (eType)
    {
    case PP:
    case LB:
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
        break;
    }
}

void Canvas::slotGalaxy(int i)
{
    Scribble * p = static_cast<Scribble *>(topLevelWidget());

    ++ i;

    switch (eType)
    {
    case BB:
    case V1:
        for (size_t j = 0; j < planet.size(); ++ j)
            for (int x = 0; x < 3; ++ x)
            {
                ostringstream s;

                s.setf(ios::scientific, ios::floatfield);
                s << std::setprecision(numeric_limits<real>::digits10);
                s << planet[j][i].p[x];

                p->pLabel[eType][j][x]->setText(s.str().c_str());
            }

        for (size_t j = 0; j < planet.size(); ++ j)
            for (int x = 0; x < 3; ++ x)
            {
                ostringstream s;

                s.setf(ios::scientific, ios::floatfield);
                s << std::setprecision(numeric_limits<real>::digits10);
                s << planet[j][i].v[1][x];

                p->pLabel[eType][j + 2][x]->setText(s.str().c_str());
            }

        for (int x = 0; x < 3; ++ x)
        {
            ostringstream s;

            s.setf(ios::scientific, ios::floatfield);
            s << std::setprecision(numeric_limits<real>::digits10);
            s << planet[1][i].v[1][x] - planet[0][i].v[1][x];

            p->pLabel[eType][4][x]->setText(s.str().c_str());
        }
        break;
    }
}

void Canvas::timerEvent(QTimerEvent *)
{
	Scribble * q = static_cast<Scribble *>(topLevelWidget());

    if (! isVisible())
        return;

	{
		QRect r(planet[0][0].p[0] / scale - 5 + width()/2, planet[0][0].p[1] / scale - 5 + height()/2, 10, 10);
        QPainter painter;
		painter.begin( &buffer );
		painter.setBrush(Qt::yellow);
		painter.drawEllipse(r);
		painter.end();
		bitBlt( this, r.x(), r.y(), &buffer, r.x(), r.y(), r.width(), r.height() );

		update(r);
	}

    switch (eType)
    {
    case PP:
    case LB:
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
            }
        }
        break;
    }

    for (size_t i = 1; i < planet[0].size(); ++ i)
    {
        for (size_t j = 0; j < planet.size(); ++ j)
		{
            QRect e(planet[j][i].o[0] / scale - 2 + width()/2, planet[j][i].o[1] / scale - 2 + height()/2, 4+1, 4+1);

			planet[j][i].o[0] = planet[j][i].p[0];
			planet[j][i].o[1] = planet[j][i].p[1];
			planet[j][i].o[2] = planet[j][i].p[2];
			
			QRect r(planet[j][i].o[0] / scale - 2 + width()/2, planet[j][i].o[1] / scale - 2 + height()/2, 4, 4);
            QPainter painter;
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

        if (planet[0][i].updated && planet[1][i].updated)
        {
            if (size_t(q->pPlanet[0]->currentIndex() + 1) == i)
                slotPlanet(i - 1);

            if (size_t(q->pPlanet[1]->currentIndex() + 1) == i)
                slotGalaxy(i - 1);

            planet[0][i].updated = false;
            planet[1][i].updated = false;
        }
    }
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
    if ( mousePressed ) {
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
//    bitBlt( &buffer, 0, 0, &tmp, 0, 0, tmp.width(), tmp.height() );
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

Scribble::Scribble( QWidget *parent, const char *name )
    : QMainWindow( parent ), nc(0)
{
	ntime[0] = upper;
	ntime[1] = 1;
    ntime[2] = 1;
    ntime[3] = 50000000000;
    ntime[4] = 1;

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
//	menu->setSeparator( QMenuBar::IPlanet::NWindowsStyle );
    setMenuBar(menu);

    QToolBar *tools = new QToolBar( this );

//    bClear = new QToolButton( QIcon(), "Clear Screen", "Clear Screen", this, SLOT( slotClear() ), tools );
//    bClear->setText( "Clear Screen" );

    pTime = new QDoubleSpinBox( tools );
    pTime->setRange(1, 50000000000);
    pTime->setDecimals(0);
    pTime->setSingleStep(10);
    pTime->setToolTip("Time Interval (s)");
    pTime->setValue( ntime[nc] );

    tools->addWidget(pTime);
    tools->addSeparator();

    pPlanet[0] = new QComboBox(tools);
    pPlanet[0]->setToolTip("Planet" );
    connect(pPlanet[0], SIGNAL(activated(int)), SLOT(slotPlanet(int)));

    pPlanet[1] = new QComboBox(tools);
    pPlanet[1]->setToolTip("Galaxy" );
    connect(pPlanet[1], SIGNAL(activated(int)), SLOT(slotGalaxy(int)));

    tools->addWidget(pPlanet[0]);
    tools->addWidget(pPlanet[1]);
    addToolBar(tools);
	
	pTabWidget = new QTabWidget(this);
    connect(pTabWidget, SIGNAL(currentChanged(int)), SLOT(slotChanged(int)));
    setCentralWidget(pTabWidget);

    QPalette* palette = new QPalette();
    palette->setColor(QPalette::WindowText,Qt::darkRed);

    for (unsigned i = 0; i < ntabs; ++ i)
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
		
        switch ((Canvas::Type) (i))
        {
        case Canvas::PP:
        case Canvas::LB:
            pLabel[i][0][0]->setToolTip(QString("Newton ") + QChar(0x0394) + QChar(0x03C1));
            pLabel[i][0][1]->setToolTip(QString("Newton ") + QChar(0x0394) + QChar(0x03C6));
            pLabel[i][0][2]->setToolTip(QString("Newton ") + QChar(0x0394) + QChar(0x03B8));
            pLabel[i][1][0]->setToolTip(QString("Finite Theory ") + QChar(0x0394) + QChar(0x03C1));
            pLabel[i][1][1]->setToolTip(QString("Finite Theory ") + QChar(0x0394) + QChar(0x03C6));
            pLabel[i][1][2]->setToolTip(QString("Finite Theory ") + QChar(0x0394) + QChar(0x03B8));
            pLabel[i][2][0]->setToolTip(QString("Finite Theory ") + QChar(0x0394) + QChar(0x03C1) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C1));
            pLabel[i][2][1]->setToolTip(QString("Finite Theory ") + QChar(0x0394) + QChar(0x03C6) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C6));
            pLabel[i][2][2]->setToolTip(QString("Finite Theory ") + QChar(0x0394) + QChar(0x03B8) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03B8));
            pLabel[i][3][0]->setToolTip(QChar(0x03BC) + QString("(Finite Theory ") + QChar(0x0394) + QChar(0x03C1) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C1) + QString(")"));
            pLabel[i][3][1]->setToolTip(QChar(0x03BC) + QString("(Finite Theory ") + QChar(0x0394) + QChar(0x03C6) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C6) + QString(")"));
            pLabel[i][3][2]->setToolTip(QChar(0x03BC) + QString("(Finite Theory ") + QChar(0x0394) + QChar(0x03B8) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03B8) + QString(")"));
            pLabel[i][4][0]->setToolTip(QString("MAD(Finite Theory ") + QChar(0x0394) + QChar(0x03C1) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C1) + QString(")"));
            pLabel[i][4][1]->setToolTip(QString("MAD(Finite Theory ") + QChar(0x0394) + QChar(0x03C6) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C6) + QString(")"));
            pLabel[i][4][2]->setToolTip(QString("MAD(Finite Theory ") + QChar(0x0394) + QChar(0x03B8) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03B8) + QString(")"));
            pLabel[i][5][0]->setToolTip(QString("Best of ") + QChar(0x03BC) + QString("(Finite Theory ") + QChar(0x0394) + QChar(0x03C1) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C1) + QString(")"));
            pLabel[i][5][1]->setToolTip(QString("Best of ") + QChar(0x03BC) + QString("(Finite Theory ") + QChar(0x0394) + QChar(0x03C6) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C6) + QString(")"));
            pLabel[i][5][2]->setToolTip(QString("Best of ") + QChar(0x03BC) + QString("(Finite Theory ") + QChar(0x0394) + QChar(0x03B8) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03B8) + QString(")"));
            pLabel[i][6][0]->setToolTip(QString("Best of ") + QString("MAD(Finite Theory ") + QChar(0x0394) + QChar(0x03C1) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C1) + QString(")"));
            pLabel[i][6][1]->setToolTip(QString("Best of ") + QString("MAD(Finite Theory ") + QChar(0x0394) + QChar(0x03C6) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03C6) + QString(")"));
            pLabel[i][6][2]->setToolTip(QString("Best of ") + QString("MAD(Finite Theory ") + QChar(0x0394) + QChar(0x03B8) + QString(" - Newton ") + QChar(0x0394) + QChar(0x03B8) + QString(")"));
            pLabel[i][5][1]->setPalette(*palette);
            break;

        case Canvas::BB:
        case Canvas::V1:
            pLabel[i][0][0]->setToolTip(QString("Newton x"));
            pLabel[i][0][1]->setToolTip(QString("Newton y"));
            pLabel[i][0][2]->setToolTip(QString("Newton z"));
            pLabel[i][1][0]->setToolTip(QString("Finite Theory x"));
            pLabel[i][1][1]->setToolTip(QString("Finite Theory y"));
            pLabel[i][1][2]->setToolTip(QString("Finite Theory z"));
            pLabel[i][2][0]->setToolTip(QString("Newton vx"));
            pLabel[i][2][1]->setToolTip(QString("Newton vy"));
            pLabel[i][2][2]->setToolTip(QString("Newton vz"));
            pLabel[i][3][0]->setToolTip(QString("Finite Theory vx"));
            pLabel[i][3][1]->setToolTip(QString("Finite Theory vy"));
            pLabel[i][3][2]->setToolTip(QString("Finite Theory vz"));
            pLabel[i][4][0]->setToolTip(QString("Finite Theory vx - Newton vx"));
            pLabel[i][4][1]->setToolTip(QString("Finite Theory vy - Newton vy"));
            pLabel[i][4][2]->setToolTip(QString("Finite Theory vz - Newton vz"));
            pLabel[i][4][0]->setPalette(*palette);
            break;
        }

        canvas[i] = new Canvas((Canvas::Type) (i), pTab[i], scale[i]);
		canvas[i]->setSizePolicy(QSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding));

        QBoxLayout * l = new QVBoxLayout(pTab[i]);
	    l->addWidget( canvas[i] );

        QBoxLayout * h1 = new QHBoxLayout();
		h1->addWidget(pLabel[i][0][0], 1);
		h1->addWidget(pLabel[i][0][1], 1);
		h1->addWidget(pLabel[i][0][2], 1);
        QBoxLayout * h2 = new QHBoxLayout();
		h2->addWidget(pLabel[i][1][0], 1);
		h2->addWidget(pLabel[i][1][1], 1);
		h2->addWidget(pLabel[i][1][2], 1);
        QBoxLayout * h3 = new QHBoxLayout();
		h3->addWidget(pLabel[i][2][0], 1);
		h3->addWidget(pLabel[i][2][1], 1);
		h3->addWidget(pLabel[i][2][2], 1);
        QBoxLayout * h4 = new QHBoxLayout();
		h4->addWidget(pLabel[i][3][0], 1);
		h4->addWidget(pLabel[i][3][1], 1);
		h4->addWidget(pLabel[i][3][2], 1);
        QBoxLayout * h5 = new QHBoxLayout();
		h5->addWidget(pLabel[i][4][0], 1);
		h5->addWidget(pLabel[i][4][1], 1);
		h5->addWidget(pLabel[i][4][2], 1);
        QBoxLayout * h6 = new QHBoxLayout();
		h6->addWidget(pLabel[i][5][0], 1);
		h6->addWidget(pLabel[i][5][1], 1);
		h6->addWidget(pLabel[i][5][2], 1);
        QBoxLayout * h7 = new QHBoxLayout();
		h7->addWidget(pLabel[i][6][0], 1);
		h7->addWidget(pLabel[i][6][1], 1);
		h7->addWidget(pLabel[i][6][2], 1);

        l->addLayout(h1);
        l->addLayout(h2);
        l->addLayout(h3);
        l->addLayout(h4);
        l->addLayout(h5);
        l->addLayout(h6);
        l->addLayout(h7);
        pTab[i]->setLayout(l);
        pTab[i]->show();
	}
	
	pTabWidget->addTab(pTab[0], "Perihelion Precession");
	pTabWidget->addTab(pTab[1], "Light Bending");
    pTabWidget->addTab(pTab[2], "Big Bang");
    pTabWidget->addTab(pTab[3], "Galactic Rotation");
    pTabWidget->addTab(pTab[4], "Pioneer 10");
    //pTab[1]->hide();

    setCentralWidget( pTabWidget );
}

void Scribble::slotRestart()
{
    qApp->quit();
    QProcess::startDetached(qApp->arguments()[0], qApp->arguments());
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
    pPlanet[0]->setEnabled(true);
    pPlanet[1]->setEnabled(false);
    pTabWidget->setCurrentWidget(pTab[nc]);
    pTime->setValue( ntime[nc] );
}

void Scribble::slotLB()
{
    ntime[nc] = pTime->value();

	nc = 1;
    pPlanet[0]->setEnabled(false);
    pPlanet[1]->setEnabled(false);
    pTabWidget->setCurrentWidget(pTab[nc]);
    pTime->setValue( ntime[nc] );
}

void Scribble::slotBB()
{
    ntime[nc] = pTime->value();

    nc = 2;
    pPlanet[0]->setEnabled(false);
    pPlanet[1]->setEnabled(true);
    pTabWidget->setCurrentWidget(pTab[nc]);
    pTime->setValue( ntime[nc] );
}

void Scribble::slotGR()
{
    ntime[nc] = pTime->value();

    nc = 3;
    pPlanet[0]->setEnabled(false);
    pPlanet[1]->setEnabled(false);
    pTabWidget->setCurrentWidget(pTab[nc]);
    pTime->setValue( ntime[nc] );
}

void Scribble::slotV1()
{
    ntime[nc] = pTime->value();

    nc = 4;
    pPlanet[0]->setEnabled(false);
    pPlanet[1]->setEnabled(false);
    pTabWidget->setCurrentWidget(pTab[nc]);
    pTime->setValue( ntime[nc] );
}

void Scribble::slotChanged(int i)
{
    switch (i)
    {
    case 0: slotPP(); break;
    case 1: slotLB(); break;
    case 2: slotBB(); break;
    case 3: slotGR(); break;
    case 4: slotV1(); break;
    }
}

void Scribble::slotAbout()
{
    QMessageBox::about( this, "Finite Theory of the Universe", "\nCopyright (c) 2011-2015\n\nPhil Bouchard <pbouchard8@gmail.com>\n");
}
	

int main( int argc, char **argv )
{
    QApplication a( argc, argv );

    Scribble scribble;

    scribble.resize( 500, 360 );
    scribble.setWindowTitle("Finite Theory of the Universe " EDITION);
	a.setStyle("windows");
	
    scribble.showMaximized();
	
    return a.exec();
}
