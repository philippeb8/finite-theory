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

#define EDITION "5.1.3"

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
#include <mutex>

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
//#include <QtWidgets/QDesktopWidget>
#include <QtCore/QProcess>

using namespace std;

const ::real upper = 0.1L;

// FT time formula
// observer is infinitly far away
::real Planet::FT_Time(::real m, ::real d, ::real h)
{
    return abs(m) / abs(d);
}

::real Planet::FT_Force(::real G, ::real m1, ::real m2, ::real d, ::real h)
{
    return (2 * G * m2 * m1 * d * pow(h, 3)) / pow(d * h + abs(m2) + abs(m1), 3) - (G * m2 * m1 * pow(h, 2)) / pow(d * h + abs(m2) + abs(m1), 2);
}

// Newton time formula
::real Planet::NW_Time(::real m, ::real d, ::real h)
{
    return 0;
}

::real Planet::NW_Force(::real G, ::real m1, ::real m2, ::real d, ::real h)
{
    return G * m1 * m2 / (d * d);
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

inline void Planet::operator () (const vector<Planet> &planet, const ::real & upper)
{
    // net acceleration vector (with all planets)
    netforce = vector3(0.L, 0.L, 0.L);

    switch (eType)
    {
    default:
        // same effect:
#if 1
        tg[0] = 0;
        te[0] = 0;
#endif

        // iterate through all planets
        for (size_t i = 0; i < planet.size(); i ++)
        {
            // if same particle then skip
            if (planet[i].n == n)
                continue;

            // vector and norm between the moving particle and the other one
            vector3 normal((p[0] - planet[i].p[0]), (p[1] - planet[i].p[1]), (p[2] - planet[i].p[2]));

            const ::real norm2 = pow(normal[0], 2) + pow(normal[1], 2) + pow(normal[2], 2);
            const ::real norm = sqrt(norm2);

#if 1
            // calculate gravitational and electric forces
            const ::real fg = force(G, planet[i].m, m, norm, hg);
            const ::real fe = force(K, planet[i].q, q, norm, he);

            netforce[0] -= fg * normal[0] / norm;
            netforce[1] -= fg * normal[1] / norm;
            netforce[2] -= fg * normal[2] / norm;

            netforce[0] += fe * normal[0] / norm;
            netforce[1] += fe * normal[1] / norm;
            netforce[2] += fe * normal[2] / norm;

            // calculate gravitational and electric time dilation increments
            tg[0] += time(planet[i].m, norm, hg);
            te[0] += time(planet[i].q, norm, he);
#endif
        }
        break;
    }

#if 0
    // spherical coordinates
    if (! first)
    {
        ps[2] = ps[1];
        ps[1] = ps[0];
    }

    ps[0][0] = sqrt(pow(p[0], 2) + pow(p[1], 2) + pow(p[2], 2));
    ps[0][1] = atan2(p[1], p[0]);
    ps[0][2] = acos(p[2] / ps[0][0]);

    if (first)
    {
        ps[2] = ps[0];
        ps[1] = ps[0];
    }
#endif

    // save
    const vector3 s(p[0], p[1], p[2]);

#if 1
    // calculate net gravitational and electric time dilation factors
    tg[0] = hg / (tg[0] + hg);
    te[0] = he / (te[0] + he);

    // save old time value
    if (! first)
    {
        tg[1] = tg[0];
        te[1] = te[0];
    }

    if (first)
    {
        tg[1] = tg[0];
        te[1] = te[0];
    }

    const ::real dt = upper * tg[0] * te[0];

    // v = v + a*t
    v[0] += netforce / m * dt;

    // p = p + v*t + (a*t^2)/2
    p += v[0] * dt;
#endif

    switch (eType)
	{
    // perihelion precession
    case PP:
        if (first || ps[1][0] < ps[2][0] && ps[1][0] < ps[0][0])
        {
            if (! first)
            {
                ps[4] = ps[3];
            }

            ps[3] = ps[1];

            if (first)
            {
                ps[4] = ps[3];
            }
            else
            {
                updated = true;
            }
        }
        break;

	// gravitational light bending
	case LB:
        if (s[0] >= -200000000000.L && p[0] < -200000000000.L)
		{
			pp[1] = pp[0];
			
            updated = true;
		}
		break;

    // big bang & pioneer 10
    case BB:
    case V1:
        if (floor(s[0] / 1e10) != floor(p[0] / 1e10))
        {
            v[1] = v[0];

            updated = true;
        }
        break;

    case NU:
    case QU:
        updated = true;

        break;
    }

    first = false;
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
		// stop the processing until the tab becomes visible
		while (! p->isVisible())
			QThread::msleep(100);
		
		// move the same planet or photon according to Newton & FT
        {
          vector<Planet> temporary = p->planet;

          for (size_t i = 0; i < p->planet.size(); ++ i)
            temporary[i](p->planet, q->pTime->value());

          p->planet = temporary;
        }
	}
}

const bool no_writing = false;

Canvas::Canvas( Type eType, size_t t, QWidget *parent)
    : QWidget( parent/*, name, Qt::WStaticContents*/ ),
      eType(eType), t(t), pen( Qt::red, 3 ), polyline(3), mousePressed( false ), buffer( width(), height() )
{
//	setAttribute(Qt::WA_PaintOutsidePaintEvent, true);

#if 0
    switch (eType)
    {
    case PP: scale = 8e8L; break;
    case LB: scale = 8e10L; break;
    case BB: scale = 8e9L; break;
    case GR: scale = 8e10L; break;
    case V1: scale = 8e11L; break;
    case NU: scale = 1e-11L; break;
    case QU: scale = 1e-17L; break;
    }
#endif

	Scribble * q = static_cast<Scribble *>(topLevelWidget());
	
	// initial position of each planet and photon
    static const ::real pos[2][75][3] =
	{
        // newton
        {
            {0.L, 0.L, 0.L},
            //{-57025548112.2453L, 3197006916.08582L, 5283916036.50742L},
            //{-57358990223.0187831L, 0.L, 0.L},
            {-45922743041.70308L, 0.L, 0.L},
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

            // galaxy:
            {-2.46e20L * 4/4, 0.L, 0.L},
            {-2.46e20L * 3/4, 0.L, 0.L},
            {-2.46e20L * 2/4, 0.L, 0.L},
            {-2.46e20L * 1/4, 0.L, 0.L},
            {2.46e20L * 1/4, 0.L, 0.L},
            {2.46e20L * 2/4, 0.L, 0.L},
            {2.46e20L * 3/4, 0.L, 0.L},
            {2.46e20L * 4/4, 0.L, 0.L},

            {17048116800000.L, 0.L, 0.L},
            {17048116800000.L, 0.L, 0.L},

            // quarks positions:
            {0e-16L, 0.L, 0.L},
            {2e-16L, 0.L, 0.L},
            {3e-16L, 0.L, 0.L},

            {-1e-13L, 1e-13L, 0.L},
            {0e-13L, 1e-13L, 0.L},
            {1e-13L, 1e-13L, 0.L},

            {1e-13L, -1e-13L, 0.L},
            {0e-13L, -1e-13L, 0.L},
            {-1e-13L, -1e-13L, 0.L},

            // protons, neutrons and electrons:
            {1e-15L, 0.L, 0.L},
            {-1e-15L, 0.L, 0.L},
            {0.L, 1e-15L, 0.L},
            {0.L, -1e-15L, 0.L},
            {0.L, 0.L, 1e-15L},
            {0.L, 0.L, -1e-15L},
            {1e-15L, 1e-15L, 0.L},
            {-1e-15L, 1e-15L, 0.L},
            {0.L, 1e-15L, 1e-15L},
            {0.L, -1e-15L, 1e-15L},
            {1e-15L, 0.L, 1e-15L},
            {1e-15L, 0.L, -1e-15L},
            {5.29177e-11 * 1 * 1, 0.L, 0.L},
            {-5.29177e-11 * 1 * 1, 0.L, 0.L},
            {0.L, 5.29177e-11 * 2 * 2, 0.L},
            {0.L, -5.29177e-11 * 2 * 2, 0.L},
            {0.L, 0.L, 5.29177e-11 * 3 * 3},
            {0.L, 0.L, -5.29177e-11 * 3 * 3},

            {5e-10L + 1e-15L, 0.L, 0.L},
            {5e-10L + -1e-15L, 0.L, 0.L},
            {5e-10L + 0.L, 1e-15L, 0.L},
            {5e-10L + 0.L, -1e-15L, 0.L},
            {5e-10L + 0.L, 0.L, 1e-15L},
            {5e-10L + 0.L, 0.L, -1e-15L},
            {5e-10L + 1e-15L, 1e-15L, 0.L},
            {5e-10L + -1e-15L, 1e-15L, 0.L},
            {5e-10L + 0.L, 1e-15L, 1e-15L},
            {5e-10L + 0.L, -1e-15L, 1e-15L},
            {5e-10L + 1e-15L, 0.L, 1e-15L},
            {5e-10L + 1e-15L, 0.L, -1e-15L},
            {5e-10L + 5.29177e-11 * 1 * 1, 0.L, 0.L},
            {5e-10L + -5.29177e-11 * 1 * 1, 0.L, 0.L},
            {5e-10L + 0.L, 5.29177e-11 * 2 * 2, 0.L},
            {5e-10L + 0.L, -5.29177e-11 * 2 * 2, 0.L},
            {5e-10L + 0.L, 0.L, 5.29177e-11 * 3 * 3},
            {5e-10L + 0.L, 0.L, -5.29177e-11 * 3 * 3},
        },
        // finite theory
        {
            {0.L, 0.L, 0.L},
            //{-57025548112.2453L, 3197006916.08582L, 5283916036.50742L},
            //{-57358990223.0187831L, 0.L, 0.L},
            {-45922743041.70308L, 0.L, 0.L},
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

            // galaxy:
            {-4e20L * 8/4, 0.L, 0.L},
            {-4e20L * 7/4, 0.L, 0.L},
            {-4e20L * 6/4, 0.L, 0.L},
            {-4e20L * 5/4, 0.L, 0.L},
            {4e20L * 5/4, 0.L, 0.L},
            {4e20L * 6/4, 0.L, 0.L},
            {4e20L * 7/4, 0.L, 0.L},
            {4e20L * 8/4, 0.L, 0.L},

            {17048116800000.L, 0.L, 0.L},
            {17048116800000.L, 0.L, 0.L},

            // quarks positions:
            {0e-15L, 0.L, 0.L},
            {0e-15L, 0.L, 0.L},
            {0e-15L, 0.L, 0.L},

            {-1e-13L, 1e-13L, 0.L},
            {0e-13L, 1e-13L, 0.L},
            {1e-13L, 1e-13L, 0.L},

            {1e-13L, -1e-13L, 0.L},
            {0e-13L, -1e-13L, 0.L},
            {-1e-13L, -1e-13L, 0.L},

            // protons, neutrons and electrons:
            {1e-15L, 0.L, 0.L},
            {-1e-15L, 0.L, 0.L},
            {0.L, 1e-15L, 0.L},
            {0.L, -1e-15L, 0.L},
            {0.L, 0.L, 1e-15L},
            {0.L, 0.L, -1e-15L},
            {1e-15L, 1e-15L, 0.L},
            {-1e-15L, 1e-15L, 0.L},
            {0.L, 1e-15L, 1e-15L},
            {0.L, -1e-15L, 1e-15L},
            {1e-15L, 0.L, 1e-15L},
            {1e-15L, 0.L, -1e-15L},
            {5.29177e-11 * 1 * 1, 0.L, 0.L},
            {-5.29177e-11 * 1 * 1, 0.L, 0.L},
            {0.L, 5.29177e-11 * 2 * 2, 0.L},
            {0.L, -5.29177e-11 * 2 * 2, 0.L},
            {0.L, 0.L, 5.29177e-11 * 3 * 3},
            {0.L, 0.L, -5.29177e-11 * 3 * 3},

            {5e-10L + 1e-15L, 0.L, 0.L},
            {5e-10L + -1e-15L, 0.L, 0.L},
            {5e-10L + 0.L, 1e-15L, 0.L},
            {5e-10L + 0.L, -1e-15L, 0.L},
            {5e-10L + 0.L, 0.L, 1e-15L},
            {5e-10L + 0.L, 0.L, -1e-15L},
            {5e-10L + 1e-15L, 1e-15L, 0.L},
            {5e-10L + -1e-15L, 1e-15L, 0.L},
            {5e-10L + 0.L, 1e-15L, 1e-15L},
            {5e-10L + 0.L, -1e-15L, 1e-15L},
            {5e-10L + 1e-15L, 0.L, 1e-15L},
            {5e-10L + 1e-15L, 0.L, -1e-15L},
            {5e-10L + 5.29177e-11 * 1 * 1, 0.L, 0.L},
            {5e-10L + -5.29177e-11 * 1 * 1, 0.L, 0.L},
            {5e-10L + 0.L, 5.29177e-11 * 2 * 2, 0.L},
            {5e-10L + 0.L, -5.29177e-11 * 2 * 2, 0.L},
            {5e-10L + 0.L, 0.L, 5.29177e-11 * 3 * 3},
            {5e-10L + 0.L, 0.L, -5.29177e-11 * 3 * 3},
        }
    };
	
	// initial velocity of each planet and photon
    static const ::real vel[2][75][3] =
	{
        // newton:
        {
            {0.L, 0.L, 0.L},
            //{-13058.0445420602L, -46493.5791091285L, -2772.42900405547L},
            //{0.L, -48372.0145148178242L, 0.L},
            {0.L, -59148.05641434967L, 0.L},
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

            // galaxy:
            {0.L, 7e4L * 4/4, 0.L},
            {0.L, 7e4L * 4/3, 0.L},
            {0.L, 7e4L * 4/2, 0.L},
            {0.L, 7e4L * 4/1, 0.L},
            {0.L, -7e4L * 4/1, 0.L},
            {0.L, -7e4L * 4/2, 0.L},
            {0.L, -7e4L * 4/3, 0.L},
            {0.L, -7e4L * 4/4, 0.L},

            {11992.L, 0.L, 0.L},
            {11992.L, 0.L, 0.L},

            // quarks velocities:
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},

            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},

            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},

            // protons, neutrons and electrons:
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 5e5L, 0.L},
            {0.L, -5e5L, 0.L},
            {5e5L, 0.L, 0.L},
            {-5e5L, 0.L, 0.L},
            {0.L, 0.L, 5e5L},
            {0.L, 0.L, -5e5L},

            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 5e5L, 0.L},
            {-5e5 + 0.L, -5e5L, 0.L},
            {-5e5 + 5e5L, 0.L, 0.L},
            {-5e5 + -5e5L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 5e5L},
            {-5e5 + 0.L, 0.L, -5e5L}
        },
        // finite theory:
        {
            {0.L, 0.L, 0.L},
            //{-13058.0445420602L, -46493.5791091285L, -2772.42900405547L},
            //{0.L, -48372.0145148178242L, 0.L},
            {0.L, -59148.05641434967L, 0.L},
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

            // galaxy:
            {0.L, 2.566161566333395e4L, 0.L},
            {0.L, 2.368653289642146e4L, 0.L},
            {0.L, 2.0657202132228e4L, 0.L},
            {0.L, 1.561716195518291e4L, 0.L},
            {0.L, -1.561716195518291e4L, 0.L},
            {0.L, -2.0657202132228e4L, 0.L},
            {0.L, -2.368653289642146e4L, 0.L},
            {0.L, -2.566161566333395e4L, 0.L},

            {11992.L, 0.L, 0.L},
            {11992.L, 0.L, 0.L},

            // quarks velocities:
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},

            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},

            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},

            // protons, neutrons and electrons:
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 0.L, 0.L},
            {0.L, 5e5L, 0.L},
            {0.L, -5e5L, 0.L},
            {5e5L, 0.L, 0.L},
            {-5e5L, 0.L, 0.L},
            {0.L, 0.L, 5e5L},
            {0.L, 0.L, -5e5L},

            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 0.L},
            {-5e5 + 0.L, 5e5L, 0.L},
            {-5e5 + 0.L, -5e5L, 0.L},
            {-5e5 + 5e5L, 0.L, 0.L},
            {-5e5 + -5e5L, 0.L, 0.L},
            {-5e5 + 0.L, 0.L, 5e5L},
            {-5e5 + 0.L, 0.L, -5e5L}
        }
    };

	// name, color, mass, position and velocity of each moving object
    static const Planet Sun 	  ("Sun", 		Qt::yellow, 1.98911E+30L, 0, pos[0][0], vel[0][0], Planet::NW_Time, Planet::NW_Force, Planet::V1, H[0], Eta);
    static const Planet Mercury   ("Mercury", 	Qt::red, 3.302E+23L, 0, pos[0][1], vel[0][1], Planet::NW_Time, Planet::NW_Force, Planet::V1, H[0], Eta);
    static const Planet Venus 	  ("Venus", 	Qt::cyan, 4.8685E+24L, 0, pos[0][2], vel[0][2], Planet::NW_Time, Planet::NW_Force, Planet::V1, H[0], Eta);
    static const Planet Earth 	  ("Earth", 	Qt::blue, 5.9736E+24L, 0, pos[0][3], vel[0][3], Planet::NW_Time, Planet::NW_Force, Planet::V1, H[0], Eta);
    static const Planet Mars 	  ("Mars", 		Qt::yellow, 6.41850000000001E+23L, 0, pos[0][4], vel[0][4], Planet::NW_Time, Planet::NW_Force, Planet::V1, H[0], Eta);
    static const Planet Jupiter   ("Jupiter", 	Qt::magenta, 1.8986E+27L, 0, pos[0][5], vel[0][5], Planet::NW_Time, Planet::NW_Force, Planet::V1, H[0], Eta);
    static const Planet Saturn 	  ("Saturn", 	Qt::darkRed, 5.6842928E+26L, 0, pos[0][6], vel[0][6], Planet::NW_Time, Planet::NW_Force, Planet::V1, H[0], Eta);
    static const Planet Uranus 	  ("Uranus", 	Qt::green, 8.68320000000002E+25L, 0, pos[0][7], vel[0][7], Planet::NW_Time, Planet::NW_Force, Planet::V1, H[0], Eta);
    static const Planet Neptune   ("Neptune", 	Qt::darkBlue, 1.0243E+26L, 0, pos[0][8], vel[0][8], Planet::NW_Time, Planet::NW_Force, Planet::V1, H[0], Eta);
    static const Planet Pluto 	  ("Pluto", 	Qt::darkGray, 1.27E+22L, 0, pos[0][9], vel[0][9], Planet::NW_Time, Planet::NW_Force, Planet::V1, H[0], Eta);

    static const Planet Photon1   ("Photon1", 	Qt::darkGreen, 3.7e-36L, 0, pos[0][10], vel[0][10], Planet::FT_Time, Planet::FT_Force, Planet::LB, H[0], Eta);
    static const Planet Photon2   ("Photon2", 	Qt::darkRed, 3.7e-36L, 0, pos[0][11], vel[0][11], Planet::NW_Time, Planet::NW_Force, Planet::LB, H[0], Eta);

    static const Planet Pioneer1   ("Pioneer1", 	Qt::darkGreen, 258.8L, 0, pos[0][28], vel[0][28], Planet::FT_Time, Planet::FT_Force, Planet::V1, H[0], Eta);
    static const Planet Pioneer2   ("Pioneer2", 	Qt::darkRed, 258.8L, 0, pos[0][29], vel[0][29], Planet::NW_Time, Planet::NW_Force, Planet::V1, H[0], Eta);

    static const Planet Quark1   ("Quark1", 	Qt::red, 8.38e-30, -Q*1/3, pos[0][30], vel[0][30], Planet::NW_Time, Planet::NW_Force, Planet::QU, H[0], Eta);
    static const Planet Quark2   ("Quark2", 	Qt::blue, 3.92e-30, Q*2/3, pos[0][31], vel[0][31], Planet::NW_Time, Planet::NW_Force, Planet::QU, H[0], Eta);
    static const Planet Quark3   ("Quark3", 	Qt::blue, 3.92e-30, Q*2/3, pos[0][32], vel[0][32], Planet::NW_Time, Planet::NW_Force, Planet::QU, H[0], Eta);

    static const Planet Quark4   ("Quark4", 	Qt::red, 8.38e-30, -Q*1/3, pos[0][33], vel[0][33], Planet::NW_Time, Planet::NW_Force, Planet::QU, H[0], Eta);
    static const Planet Quark5   ("Quark5", 	Qt::blue, 3.92e-30, Q*2/3, pos[0][34], vel[0][34], Planet::NW_Time, Planet::NW_Force, Planet::QU, H[0], Eta);
    static const Planet Quark6   ("Quark6", 	Qt::blue, 3.92e-30, Q*2/3, pos[0][35], vel[0][35], Planet::NW_Time, Planet::NW_Force, Planet::QU, H[0], Eta);

    static const Planet Quark7   ("Quark7", 	Qt::blue, 3.92e-30, Q*2/3, pos[0][36], vel[0][36], Planet::NW_Time, Planet::NW_Force, Planet::QU, H[0], Eta);
    static const Planet Quark8   ("Quark8", 	Qt::red, 8.38e-30, -Q*1/3, pos[0][37], vel[0][37], Planet::NW_Time, Planet::NW_Force, Planet::QU, H[0], Eta);
    static const Planet Quark9   ("Quark9", 	Qt::red, 8.38e-30, -Q*1/3, pos[0][38], vel[0][38], Planet::NW_Time, Planet::NW_Force, Planet::QU, H[0], Eta);

    static const Planet Proton1     ("Proton1",   Qt::red, 1.6726e-27, Q, pos[0][39], vel[0][39], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Proton2     ("Proton2",   Qt::red, 1.6726e-27, Q, pos[0][40], vel[0][40], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Proton3     ("Proton3",   Qt::red, 1.6726e-27, Q, pos[0][41], vel[0][41], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Proton4     ("Proton4",   Qt::red, 1.6726e-27, Q, pos[0][42], vel[0][42], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Proton5     ("Proton5",   Qt::red, 1.6726e-27, Q, pos[0][43], vel[0][43], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Proton6     ("Proton6",   Qt::red, 1.6726e-27, Q, pos[0][44], vel[0][44], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Neutron1     ("Neutron1",   Qt::yellow, 1.6726e-27, 0, pos[0][45], vel[0][45], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Neutron2     ("Neutron2",   Qt::yellow, 1.6726e-27, 0, pos[0][46], vel[0][46], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Neutron3     ("Neutron3",   Qt::yellow, 1.6726e-27, 0, pos[0][47], vel[0][47], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Neutron4     ("Neutron4",   Qt::yellow, 1.6726e-27, 0, pos[0][48], vel[0][48], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Neutron5     ("Neutron5",   Qt::yellow, 1.6726e-27, 0, pos[0][49], vel[0][49], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Neutron6     ("Neutron6",   Qt::yellow, 1.6726e-27, 0, pos[0][50], vel[0][50], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Electron1   ("Electron1", Qt::blue, 9.109e-31, -Q, pos[0][51], vel[0][51], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Electron2   ("Electron2", Qt::blue, 9.109e-31, -Q, pos[0][52], vel[0][52], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Electron3   ("Electron3", Qt::blue, 9.109e-31, -Q, pos[0][53], vel[0][53], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Electron4   ("Electron4", Qt::blue, 9.109e-31, -Q, pos[0][54], vel[0][54], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Electron5   ("Electron5", Qt::blue, 9.109e-31, -Q, pos[0][55], vel[0][55], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Electron6   ("Electron6", Qt::blue, 9.109e-31, -Q, pos[0][56], vel[0][56], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);

    static const Planet Proton7     ("Proton7",   Qt::red, 1.6726e-27, Q, pos[0][57], vel[0][57], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Proton8     ("Proton8",   Qt::red, 1.6726e-27, Q, pos[0][58], vel[0][58], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Proton9     ("Proton9",   Qt::red, 1.6726e-27, Q, pos[0][59], vel[0][59], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Proton10     ("Proton10",   Qt::red, 1.6726e-27, Q, pos[0][60], vel[0][60], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Proton11     ("Proton11",   Qt::red, 1.6726e-27, Q, pos[0][61], vel[0][61], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Proton12     ("Proton12",   Qt::red, 1.6726e-27, Q, pos[0][62], vel[0][62], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Neutron7     ("Neutron7",   Qt::yellow, 1.6726e-27, 0, pos[0][63], vel[0][63], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Neutron8     ("Neutron8",   Qt::yellow, 1.6726e-27, 0, pos[0][64], vel[0][64], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Neutron9     ("Neutron9",   Qt::yellow, 1.6726e-27, 0, pos[0][65], vel[0][65], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Neutron10     ("Neutron10",   Qt::yellow, 1.6726e-27, 0, pos[0][66], vel[0][66], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Neutron11     ("Neutron11",   Qt::yellow, 1.6726e-27, 0, pos[0][67], vel[0][67], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Neutron12     ("Neutron12",   Qt::yellow, 1.6726e-27, 0, pos[0][68], vel[0][68], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Electron7   ("Electron7", Qt::blue, 9.109e-31, -Q, pos[0][69], vel[0][69], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Electron8   ("Electron8", Qt::blue, 9.109e-31, -Q, pos[0][70], vel[0][70], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Electron9   ("Electron9", Qt::blue, 9.109e-31, -Q, pos[0][71], vel[0][71], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Electron10   ("Electron10", Qt::blue, 9.109e-31, -Q, pos[0][72], vel[0][72], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Electron11   ("Electron11", Qt::blue, 9.109e-31, -Q, pos[0][73], vel[0][73], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);
    static const Planet Electron12   ("Electron12", Qt::blue, 9.109e-31, -Q, pos[0][74], vel[0][74], Planet::NW_Time, Planet::NW_Force, Planet::NU, H[0], Eta);


    static const Planet Core	  ("Core", 		Qt::black, 2E+11L, 0, pos[0][0], vel[0][0], Planet::NW_Time, Planet::NW_Force, Planet::BB, H[1], Eta);
    static const Planet Galaxy1   ("Galaxy1", 	Qt::red, 50000L, 0, pos[0][12], vel[0][12], Planet::NW_Time, Planet::NW_Force, Planet::BB, H[1], Eta);
    static const Planet Galaxy2   ("Galaxy2", 	Qt::cyan, 50000L, 0, pos[0][13], vel[0][13], Planet::NW_Time, Planet::NW_Force, Planet::BB, H[1], Eta);
    static const Planet Galaxy3   ("Galaxy3", 	Qt::blue, 50000L, 0, pos[0][14], vel[0][14], Planet::NW_Time, Planet::NW_Force, Planet::BB, H[1], Eta);
    static const Planet Galaxy4   ("Galaxy4", 	Qt::yellow, 50000L, 0, pos[0][15], vel[0][15], Planet::NW_Time, Planet::NW_Force, Planet::BB, H[1], Eta);
    static const Planet Galaxy5   ("Galaxy5", 	Qt::magenta, 50000L, 0, pos[0][16], vel[0][16], Planet::NW_Time, Planet::NW_Force, Planet::BB, H[1], Eta);
    static const Planet Galaxy6   ("Galaxy6", 	Qt::darkRed, 50000L, 0, pos[0][17], vel[0][17], Planet::NW_Time, Planet::NW_Force, Planet::BB, H[1], Eta);
    static const Planet Galaxy7   ("Galaxy7", 	Qt::green, 50000L, 0, pos[0][18], vel[0][18], Planet::NW_Time, Planet::NW_Force, Planet::BB, H[1], Eta);
    static const Planet Galaxy8   ("Galaxy8", 	Qt::darkBlue, 50000L, 0, pos[0][19], vel[0][19], Planet::NW_Time, Planet::NW_Force, Planet::BB, H[1], Eta);

    static const Planet Buldge1    ("Buldge1", 	Qt::black, 2e10L * 2e30L, 0, pos[0][0], vel[0][0], Planet::NW_Time, Planet::NW_Force, Planet::GR, H[2], Eta);
    static const Planet Star11     ("Star11", 	Qt::red, 2e30L, 0, pos[0][20], vel[0][20], Planet::NW_Time, Planet::NW_Force, Planet::GR, H[2], Eta);
    static const Planet Star12     ("Star12", 	Qt::red, 2e30L, 0, pos[0][21], vel[0][21], Planet::NW_Time, Planet::NW_Force, Planet::GR, H[2], Eta);
    static const Planet Star13     ("Star13", 	Qt::red, 2e30L, 0, pos[0][22], vel[0][22], Planet::NW_Time, Planet::NW_Force, Planet::GR, H[2], Eta);
    static const Planet Star14     ("Star14", 	Qt::red, 2e30L, 0, pos[0][23], vel[0][23], Planet::NW_Time, Planet::NW_Force, Planet::GR, H[2], Eta);
    static const Planet Star15     ("Star15", 	Qt::red, 2e30L, 0, pos[0][24], vel[0][24], Planet::NW_Time, Planet::NW_Force, Planet::GR, H[2], Eta);
    static const Planet Star16     ("Star16", 	Qt::red, 2e30L, 0, pos[0][25], vel[0][25], Planet::NW_Time, Planet::NW_Force, Planet::GR, H[2], Eta);
    static const Planet Star17     ("Star17", 	Qt::red, 2e30L, 0, pos[0][26], vel[0][26], Planet::NW_Time, Planet::NW_Force, Planet::GR, H[2], Eta);
    static const Planet Star18     ("Star18", 	Qt::red, 2e30L, 0, pos[0][27], vel[0][27], Planet::NW_Time, Planet::NW_Force, Planet::GR, H[2], Eta);

    static const Planet Buldge2    ("Buldge2", 	Qt::black, 2e10L * 2e30L, 0, pos[1][0], vel[1][0], Planet::FT_Time, Planet::FT_Force, Planet::GR, H[2], Eta);
    static const Planet Star21     ("Star21", 	Qt::red, 2e30L, 0, pos[1][20], vel[1][20], Planet::FT_Time, Planet::FT_Force, Planet::GR, H[2], Eta);
    static const Planet Star22     ("Star22", 	Qt::red, 2e30L, 0, pos[1][21], vel[1][21], Planet::FT_Time, Planet::FT_Force, Planet::GR, H[2], Eta);
    static const Planet Star23     ("Star23", 	Qt::red, 2e30L, 0, pos[1][22], vel[1][22], Planet::FT_Time, Planet::FT_Force, Planet::GR, H[2], Eta);
    static const Planet Star24     ("Star24", 	Qt::red, 2e30L, 0, pos[1][23], vel[1][23], Planet::FT_Time, Planet::FT_Force, Planet::GR, H[2], Eta);
    static const Planet Star25     ("Star25", 	Qt::red, 2e30L, 0, pos[1][24], vel[1][24], Planet::FT_Time, Planet::FT_Force, Planet::GR, H[2], Eta);
    static const Planet Star26     ("Star26", 	Qt::red, 2e30L, 0, pos[1][25], vel[1][25], Planet::FT_Time, Planet::FT_Force, Planet::GR, H[2], Eta);
    static const Planet Star27     ("Star27", 	Qt::red, 2e30L, 0, pos[1][26], vel[1][26], Planet::FT_Time, Planet::FT_Force, Planet::GR, H[2], Eta);
    static const Planet Star28     ("Star28", 	Qt::red, 2e30L, 0, pos[1][27], vel[1][27], Planet::FT_Time, Planet::FT_Force, Planet::GR, H[2], Eta);

    switch (eType)
	{
	// perihelion precession disparity
	case PP:
		// store each planet using the Newton time formula
        planet.reserve(10);
        planet.push_back(Sun);
        planet.push_back(Mercury);
        planet.push_back(Venus);
        planet.push_back(Earth);
        planet.push_back(Mars);
        planet.push_back(Jupiter);
        planet.push_back(Saturn);
        planet.push_back(Uranus);
        planet.push_back(Neptune);
        planet.push_back(Pluto);
		
		// copy & change each planet for the FT time formula
        if (t == 1)
            for (size_t i = 0; i < planet.size(); i ++)
            {
                planet[i].time = Planet::FT_Time;
                planet[i].force = Planet::FT_Force;
            }

		break;
	
	// gravitational light bending
	case LB:
        // store the Sun & the photon using the Newton time formula
        planet.reserve(2);
        planet.push_back(Sun);

        switch (t)
        {
        case 0: planet.push_back(Photon1); break;
        case 1: planet.push_back(Photon2); break;
        }

        break;

    // big bang
    case BB:
        // store the Sun & the planets using FT time formula
        planet.reserve(9);
        planet.push_back(Core);
        planet.push_back(Galaxy1);
        planet.push_back(Galaxy2);
        planet.push_back(Galaxy3);
        planet.push_back(Galaxy4);
        planet.push_back(Galaxy5);
        planet.push_back(Galaxy6);
        planet.push_back(Galaxy7);
        planet.push_back(Galaxy8);

        // copy & change each planet for the FT time formula
        if (t == 1)
            for (size_t i = 0; i < planet.size(); i ++)
            {
                planet[i].time = Planet::FT_Time;
                planet[i].force = Planet::FT_Force;
            }

        break;

        // galactic rotation
    case GR:
        planet.reserve(9);

        switch (t)
        {
        case 0:
            // store the buldge & the stars using NW time formula
            planet.push_back(Buldge1);
            planet.push_back(Star11);
            planet.push_back(Star12);
            planet.push_back(Star13);
            planet.push_back(Star14);
            planet.push_back(Star15);
            planet.push_back(Star16);
            planet.push_back(Star17);
            planet.push_back(Star18);
            break;

        case 1:
            // store the buldge & the stars using FT time formula
            planet.push_back(Buldge2);
            planet.push_back(Star21);
            planet.push_back(Star22);
            planet.push_back(Star23);
            planet.push_back(Star24);
            planet.push_back(Star25);
            planet.push_back(Star26);
            planet.push_back(Star27);
            planet.push_back(Star28);
            break;
        }

        break;

    // pioneer 10
    case V1:
        // store the Sun & the photon using the Newton time formula
        planet.reserve(2);
        planet.push_back(Sun);

        switch (t)
        {
        case 0: planet.push_back(Pioneer1); break;
        case 1: planet.push_back(Pioneer2); break;
        }

        break;

    // atomic
    case NU:
        planet.reserve(24);

        // store the Sun & the photon using the Newton time formula
        planet.push_back(Proton1);
        planet.push_back(Proton2);
        planet.push_back(Proton3);
        planet.push_back(Proton4);
        planet.push_back(Proton5);
        planet.push_back(Proton6);
        planet.push_back(Neutron1);
        planet.push_back(Neutron2);
        planet.push_back(Neutron3);
        planet.push_back(Neutron4);
        planet.push_back(Neutron5);
        planet.push_back(Neutron6);
        planet.push_back(Electron1);
        planet.push_back(Electron2);
        planet.push_back(Electron3);
        planet.push_back(Electron4);
        planet.push_back(Electron5);
        planet.push_back(Electron6);

        // copy & change each planet for the FT time formula
        if (t == 1)
            for (size_t i = 0; i < planet.size(); i ++)
            {
                planet[i].time = Planet::FT_Time;
                planet[i].force = Planet::FT_Force;
            }

        break;

        // quantum
        case QU:
            // store the Sun & the photon using the Newton time formula
            planet.reserve(12);
#if 0
            planet.push_back(Quark1);
            planet.push_back(Quark2);
            planet.push_back(Quark3);
#endif
            planet.push_back(Quark4);
            planet.push_back(Quark5);
            planet.push_back(Quark6);
            planet.push_back(Quark7);
            planet.push_back(Quark8);
            planet.push_back(Quark9);

            // copy & change each planet for the FT time formula
            if (t == 1)
                for (size_t i = 0; i < planet.size(); i ++)
                {
                    planet[i].time = Planet::FT_Time;
                    planet[i].force = Planet::FT_Force;
                }

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
    new Dual(this);
}

Canvas::~Canvas()
{
}

void Canvas::slotPlanet(int i)
{
	Scribble * p = static_cast<Scribble *>(topLevelWidget());
	
	++ i;

#if 0
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
                    s << std::setprecision(numeric_limits<::real>::digits10);
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
                s << std::setprecision(numeric_limits<::real>::digits10);
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
                s[0] << std::setprecision(numeric_limits<::real>::digits10);

                set<::real>::iterator k = stats[i].mean[x].begin();
                set<::real>::reverse_iterator l = stats[i].mean[x].rbegin();
                advance(k, stats[i].mean[x].size() / 2);
                advance(l, stats[i].mean[x].size() / 2);

                const ::real median = (* k + * l) / 2;

                s[0] << median;

                // median absolute deviation
                if (stats[i].mean[0].size() > 1)
                {
                    s[1].setf(ios::scientific, ios::floatfield);
                    s[1] << std::setprecision(numeric_limits<::real>::digits10);

                    set<::real> dev;
                    for (set<::real>::iterator m = stats[i].mean[x].begin(); m != stats[i].mean[x].end(); ++ m)
                        dev.insert(abs(* m - median));

                    set<::real>::iterator m = dev.begin();
                    set<::real>::reverse_iterator n = dev.rbegin();
                    advance(m, dev.size() / 2);
                    advance(n, dev.size() / 2);

                    const ::real mad = (* m + * n) / 2;

                    s[1] << mad;

                    if (stats[i].best[1][x] > mad)
                    {
                        stats[i].best[0][x] = median;
                        stats[i].best[1][x] = mad;
                    }

                    s[2].setf(ios::scientific, ios::floatfield);
                    s[2] << std::setprecision(numeric_limits<::real>::digits10);

                    s[3].setf(ios::scientific, ios::floatfield);
                    s[3] << std::setprecision(numeric_limits<::real>::digits10);

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

    case NU:
    case QU:
#if 0
        // position
        for (size_t j = 0; j < planet[i].size(); ++ j)
            for (int x = 0; x < 3; ++ x)
            {
                ostringstream s;

                s.setf(ios::scientific, ios::floatfield);
                s << std::setprecision(numeric_limits<::real>::digits10);
                s << planet[i][j].p[x];

                p->pLabel[eType][j][x]->setText(s.str().c_str());
            }
#endif
        break;
    }
#endif
}

void Canvas::slotGalaxy(int i)
{
    Scribble * p = static_cast<Scribble *>(topLevelWidget());

    ++ i;

    switch (eType)
    {
    case BB:
    case V1:
        for (int x = 0; x < 3; ++ x)
        {
            ostringstream s;

            s.setf(ios::scientific, ios::floatfield);
            s << std::setprecision(numeric_limits<::real>::digits10);
            s << planet[i].p[x];

            p->pLabel[eType][t][x]->setText(s.str().c_str());
        }

        for (int x = 0; x < 3; ++ x)
        {
            ostringstream s;

            s.setf(ios::scientific, ios::floatfield);
            s << std::setprecision(numeric_limits<::real>::digits10);
            s << planet[i].v[1][x];

            p->pLabel[eType][t + 2][x]->setText(s.str().c_str());
        }

        for (int x = 0; x < 3; ++ x)
        {
            ostringstream s;

            s.setf(ios::scientific, ios::floatfield);
            s << std::setprecision(numeric_limits<::real>::digits10);
            s << planet[i].v[1][x] - planet[t].v[1][x];

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
        QRect r(0, 0, width(), height());
        QPainter painter;
        painter.begin( &buffer );
        painter.setBrush(Qt::white);
        painter.drawRect(0, 0, width(), height());
        painter.end();
        bitBlt( this, r.x(), r.y(), &buffer, r.x(), r.y(), r.width(), r.height() );

        update(r);
    }

#if 1
    if (q->pCheck->isChecked())
    {
        vector3 max = {numeric_limits<::real>::min(), numeric_limits<::real>::min(), numeric_limits<::real>::min()};

        for (size_t j = 0; j < planet.size(); ++ j)
        {
            if (abs(planet[j].p[0]) > max[0])
                max[0] = abs(planet[j].p[0]);

            if (abs(planet[j].p[1]) > max[1])
                max[1] = abs(planet[j].p[1]);
        }

        ::real new_scale = numeric_limits<::real>::min();

        for (size_t x = 0; x < 2; ++ x)
            if (pow(10, ceil(log10(abs(max[x] / 200)))) > new_scale)
                new_scale = pow(10, ceil(log10(abs(max[x] / 200))));

        scale = new_scale;
    }

    if (initial == 0.L)
    {
        initial = scale;
    }
#endif

    //if (new_scale != scale && ! isinf(new_scale) && ! isnan(new_scale) && new_scale != numeric_limits<::real>::min())
    ostringstream o[2];
    o[0] << "Scale: " << scale;
    o[1] << "Zoom: " << zoom;

    q->pScale->setText(o[0].str().c_str());
    q->pZoom->setText(o[1].str().c_str());

#if 0
    {
        QRect r((planet[0].p[0] / scale - 5 + width()/2), (planet[0].p[1] / scale - 5 + height()/2), 10, 10);
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
        for (size_t i = 0; i < planet.size(); ++ i)
        {
            if (planet[i].updated)
            {
                for (int x = 0; x < 3; ++ x)
                    stats[i].precession[t][x] = planet[i].ps[3][x] - planet[i].ps[4][x];
            }

            if (planet[i].updated && planet[i].updated)
            {
                for (int x = 0; x < 3; ++ x)
                    stats[i].mean[x].insert(stats[i].precession[1][x] - stats[i].precession[0][x]);
            }
        }
        break;
    }
#endif

    for (size_t i = 0; i < planet.size(); ++ i)
    {
        ::real const radius = (planet[i].m / planet[0].m) / (scale * zoom / initial) + 4;

        QRect e((planet[i].o[0] / (scale * zoom) - radius + width()/2), (planet[i].o[1] / (scale * zoom) - radius + height()/2), (2 * radius), (2 * radius));

        planet[i].o[0] = planet[i].p[0];
        planet[i].o[1] = planet[i].p[1];
        planet[i].o[2] = planet[i].p[2];

        vector3 normal(planet[i].netforce[0], planet[i].netforce[1], planet[i].netforce[2]);

        const ::real norm2 = pow(normal[0], 2) + pow(normal[1], 2) + pow(normal[2], 2);
        const ::real norm = sqrt(norm2);

        QRect r((planet[i].o[0] / (scale * zoom) - radius + width()/2), (planet[i].o[1] / (scale * zoom) - radius + height()/2), (2 * radius), (2 * radius));
        QPainter painter;
        painter.begin( &buffer );
        painter.setPen(planet[i].c);
        painter.setBrush(planet[i].c);
        painter.eraseRect(e);
        painter.drawEllipse(r);

#if 0
        {
            QPointF start(planet[i].o[0] / (scale * zoom) + width()/2, planet[i].o[1] / (scale * zoom) + height()/2);
            QPointF end((planet[i].o[0]) / (scale * zoom) + width()/2 + 20 / zoom * normal[0] / norm, (planet[i].o[1]) / (scale * zoom) + height()/2 + 20 / zoom * normal[1] / norm);

            painter.drawLine(start, end);

            double arrowHeadLength = 4 / zoom;
            double arrowHeadAngle = M_PI / 4;

            double angle = std::atan2(end.y() - start.y(), end.x() - start.x());

            QPointF arrowP1 = end - QPointF(arrowHeadLength * std::cos(angle - arrowHeadAngle),
                                            arrowHeadLength * std::sin(angle - arrowHeadAngle));
            QPointF arrowP2 = end - QPointF(arrowHeadLength * std::cos(angle + arrowHeadAngle),
                                            arrowHeadLength * std::sin(angle + arrowHeadAngle));

            QPolygonF arrowHead;
            arrowHead << end << arrowP1 << arrowP2;
            painter.drawPolygon(arrowHead);
        }
#endif

        painter.end();
        r |= e;
        bitBlt( this, r.x(), r.y(), &buffer, r.x(), r.y(), r.width(), r.height() );

        update(r);

        if (planet[i].updated)
        {
            //if (size_t(q->pPlanet[j]->currentIndex() + 1) == i)
                slotPlanet(i - 1);

            planet[i].updated = false;

#if 0
            if (size_t(q->pPlanet[1]->currentIndex() + 1) == i)
                slotGalaxy(i - 1);

            planet[1][i].updated = false;
#endif
        }
    }
}

void Canvas::wheelEvent(QWheelEvent *event)
{
    QPoint d = event->pixelDelta() / 8;

    if (!d.isNull())
    {
        QPoint n = d / 15;

        zoom *= pow(10, signbit(n.y()) ? 0.1 : -0.1);
    }

    event->accept();
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

    QVector<QRect> rects = QVector<QRect>(e->region().begin(), e->region().end());
    for ( uint i = 0; i < rects.count(); ++ i ) 
	{
        QRect r = rects[(int)i];
        bitBlt( this, r.x(), r.y(), &buffer, r.x(), r.y(), r.width(), r.height() );
    }
}

//------------------------------------------------------

DualCanvas::DualCanvas(Canvas::Type eType, QWidget *parent)
    : QWidget(parent)
{
    left = new Canvas(eType, 0, this);
    right = new Canvas(eType, 1, this);

    QLabel * title[] = {new QLabel("Newton", this), new QLabel("Finite Theory", this)};

    QVBoxLayout *vlayout[] = {new QVBoxLayout(), new QVBoxLayout()};

    QHBoxLayout *hlayout = new QHBoxLayout(this);

    vlayout[0]->addWidget(title[0]);
    vlayout[0]->addWidget(left, 1);
    vlayout[1]->addWidget(title[1]);
    vlayout[1]->addWidget(right, 1);

    hlayout->addLayout(vlayout[0]);
    hlayout->addLayout(vlayout[1]);
}

//------------------------------------------------------

Scribble::Scribble( QWidget *parent, const char *name )
    : QMainWindow( parent ), nc(0)
{
    ntime[0] = 100;
	ntime[1] = 1;
    ntime[2] = 1;
    ntime[3] = 50000000000;
    ntime[4] = 1;
    ntime[5] = 1e-19;
    ntime[6] = 1e-19;

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
    pTime->setRange(1e-35, 1e35);
    pTime->setDecimals(35);
    pTime->setSingleStep(1);
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

    pCheck = new QCheckBox( "Auto Zoom", tools );
    pCheck->setChecked(true);

    tools->addWidget(pCheck);
    tools->addSeparator();

    pScale = new QLabel( tools );

    tools->addWidget(pScale);
    tools->addSeparator();

    pZoom = new QLabel( tools );

    tools->addWidget(pZoom);

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

        canvas[i] = new DualCanvas((Canvas::Type) (i), pTab[i]);
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
    pTabWidget->addTab(pTab[5], "Nuclear");
    pTabWidget->addTab(pTab[6], "Quantum");
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

void Scribble::slotNU()
{
    ntime[nc] = pTime->value();

    nc = 5;
    pPlanet[0]->setEnabled(false);
    pPlanet[1]->setEnabled(false);
    pTabWidget->setCurrentWidget(pTab[nc]);
    pTime->setValue( ntime[nc] );
}

void Scribble::slotQU()
{
    ntime[nc] = pTime->value();

    nc = 6;
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
    case 5: slotNU(); break;
    case 6: slotQU(); break;
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

    scribble.showMaximized();
    scribble.setWindowTitle("Finite Theory of the Universe " EDITION);
	a.setStyle("windows");
	
    //scribble.showMaximized();
    scribble.show();
	
    return a.exec();
}
