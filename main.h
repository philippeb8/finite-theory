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

#ifndef SCRIBBLE_H
#define SCRIBBLE_H

#include <set>
#include <cmath>
#include <limits>
#include <vector>
#include <random>

#include <qcolor.h>
#include <QtWidgets/QMainWindow>
#include <qpen.h>
#include <qpoint.h>
#include <QtWidgets/QLabel>
#include <qpixmap.h>
#include <QtWidgets/QWidget>
#include <qstring.h>
#include <QPolygon>
#include <qthread.h>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QTabWidget>
//Added by qt3to4:
#include <QResizeEvent>
#include <QMouseEvent>
#include <QTimerEvent>
#include <QPaintEvent>

class QMouseEvent;
class QResizeEvent;
class QPaintEvent;
class QToolButton;
class QDoubleSpinBox;

typedef long double real;

const real PI = std::acos(real(-1));
const real C = 299792458.L;
const real G = 6.67428e-11L;
const real H[] = {C*C/(2*G), 0., 1e20};

struct vector3
{
	typedef long double T;
	static const size_t N = 3;

	T elem_[N];

	vector3()
	{
	}
	
	vector3(const T & b1, const T & b2, const T & b3)
	{
		elem_[0] = b1;
		elem_[1] = b2;
		elem_[2] = b3;
	}

	void operator = (const vector3 & b)
	{
		for (size_t i = 0; i < N; ++ i)
			elem_[i] = b.elem_[i];
	}

	T & operator [] (const size_t n) 
	{ 
		return elem_[n]; 
	}

	const T & operator [] (const size_t n) const
	{ 
		return elem_[n]; 
	}

	real norm () const
	{
		real nnorm2 = 0.L;
		
		for (size_t i = 0; i < N; ++ i)
			nnorm2 += elem_[i] * elem_[i];

		return sqrt(nnorm2);
	}

	vector3 cross(const vector3 & b) const
	{
		vector3 result;

		result[0] = elem_[1]*b.elem_[2] - elem_[2]*b.elem_[1];
		result[1] = elem_[2]*b.elem_[0] - elem_[0]*b.elem_[2];
		result[2] = elem_[0]*b.elem_[1] - elem_[1]*b.elem_[0];
		
		return result;
	}

	vector3 operator - () const
	{
		vector3 result;

		for (size_t i = 0; i < N; ++ i)
			result[i] = - elem_[i];
		
		return result;
	}

	real operator * (const vector3 & b) const
	{
		real ndot = 0.L;
		
		for (size_t i = 0; i < N; ++ i)
			ndot += elem_[i] * b.elem_[i];

		return ndot;
	}

	vector3 operator * (const real & b) const
	{
		vector3 result;

		for (size_t i = 0; i < N; ++ i)
			result[i] = elem_[i] * b;
		
		return result;
	}

	vector3 operator / (const real & b) const
	{
		vector3 result;

		for (size_t i = 0; i < N; ++ i)
			result[i] = elem_[i] / b;
		
		return result;
	}

	vector3 operator + (const vector3 & b) const
	{
		vector3 result;

		for (size_t i = 0; i < N; ++ i)
			result[i] = elem_[i] + b.elem_[i];
		
		return result;
	}

	vector3 operator - (const vector3 & b) const
	{
		vector3 result;

		for (size_t i = 0; i < N; ++ i)
			result[i] = elem_[i] - b.elem_[i];
		
		return result;
	}

	void operator += (const vector3 & b)
	{
		for (size_t i = 0; i < N; ++ i)
			elem_[i] += b.elem_[i];
	}

	void operator -= (const vector3 & b)
	{
		for (size_t i = 0; i < N; ++ i)
			elem_[i] -= b.elem_[i];
	}

	void operator *= (const vector3 & b)
	{
		for (size_t i = 0; i < N; ++ i)
			elem_[i] *= b.elem_[i];
	}

	void operator /= (const vector3 & b)
	{
		for (size_t i = 0; i < N; ++ i)
			elem_[i] /= b.elem_[i];
	}

	void operator += (const real & b)
	{
		for (size_t i = 0; i < N; ++ i)
			elem_[i] += b;
	}

	void operator -= (const real & b)
	{
		for (size_t i = 0; i < N; ++ i)
			elem_[i] -= b;
	}

	void operator *= (const real & b)
	{
		for (size_t i = 0; i < N; ++ i)
			elem_[i] *= b;
	}

	void operator /= (const real & b)
	{
		for (size_t i = 0; i < N; ++ i)
			elem_[i] /= b;
	}
};

struct Planet
{
    static real FR1(real m, real d, real h);
    static real FR2(real m, real d, real h);
    static real NW(real m, real d, real h);

    char const * n;						// name
	QColor c;							// color
	real m;								// mass
	vector3 p;							// position
    vector3 v[2];						// current & saved velocity
    vector3 o;							// old position
    real t[2];							// current & old time intervals according to Newton or FT
    bool first;                         // first cycle
	bool updated;						// the cycle of the planet or the photon arrival line has been completed
    vector3 pp[2];						// current & old saved positions on the perihelion
    vector3 ps[5];						// current & old polar coordinates of pp
    real (* f)(real, real, real);   	// function pointer to Newton time formula or FT time formula
    real h;                             // fudge factor

    real alpha;
    real alphavis;
    real alphatot;
    real r;
    real x;
    real y;
    real xvis;
    real yvis;
    real xtot;
    real ytot;
    real vel;
    real omega;
    real omegavis;
    real omegatot;
    real vvis;
    real vdark;
    real vtot;
    real mass;
    real massvis;
    real massdk;
    real masst;

    enum Type {PP, LB, BB, GR, V1} eType;		// is for the perihelion precession disparity or the gravitational light bending

    Planet(char const * n, const QColor & c, real m, const real pp[3], const real pv[3], real (* f)(real, real, real) = NW, Type eType = PP, real h = H[0])
        : n(n), c(c), m(m), p(pp[0], pp[1], pp[2]), v({vector3(pv[0], pv[1], pv[2])}), first(true), updated(false), f(f), eType(eType), h(h)
	{
	}
	
    void operator () (const real & dt);
};

class Canvas;
	
class Dual : public QThread
{
public:
    Dual(Canvas *);
	virtual void run();
	
protected:
	Canvas * p;
};


class Canvas : public QWidget
{
    Q_OBJECT
	friend class Dual;

public:
    enum Type {PP, LB, BB, GR, V1} eType;

    Canvas( Type eType, QWidget *parent = 0, real scale = 8e9L );
    ~Canvas();
    void clearScreen();
	
protected slots:
	void slotPlanet(int);
    void slotGalaxy(int);

protected:
    void mousePressEvent( QMouseEvent *e );
    void mouseReleaseEvent( QMouseEvent *e );
    void mouseMoveEvent( QMouseEvent *e );
    void resizeEvent( QResizeEvent *e );
    void paintEvent( QPaintEvent *e );
	void timerEvent( QTimerEvent *e );

    QPen pen;
    QPolygon polyline;

    bool mousePressed;

    QPixmap buffer;

	std::vector< std::vector<Planet> > planet;
	
    real scale;

    static int const nt = 3;
    static int const np = 200;
    real const zoom = 1e12L;
    real const t = 0.0;
    real const dt = 0.0005;
    real const rmax = 20.0 * zoom;
    real const rmin = 0.5 * zoom;
    real const h = 4.5 * zoom;
    real const r0 = 1.0 * zoom;
    real const v0 = 140 * zoom;
    real const rdm0 = 1.0 * zoom;
    real const massf = 1.0;
    real const dmf = 1.0;
    real const fit = 0.0;
    real const emax = 2.0 / PI * atan(pow((rmax / h), 2));
    real totalmass;
    real starmass;
    real md;
    real mdk;

    struct Stats
	{
		vector3 precession[2];
		std::set<real> mean[3];
		vector3 best[2];
		
		Stats()
		{
			best[1][0] = std::numeric_limits<real>::max();
			best[1][1] = std::numeric_limits<real>::max();
			best[1][2] = std::numeric_limits<real>::max();
		}
	};
	
	std::vector<Stats> stats;
};

class Scribble : public QMainWindow
{
    Q_OBJECT
	friend class Canvas;

public:
    Scribble( QWidget *parent = 0, const char *name = 0 );

protected slots:
    void slotRestart();
    void slotClear();
    void slotPlanet(int);
	void slotPP();
	void slotLB();
    void slotBB();
    void slotGR();
    void slotV1();
    void slotChanged(int);
	void slotAbout();
	
public:
    static const int ntabs = 5;

    unsigned nc;
    real ntime[ntabs];

	QTabWidget *pTabWidget;
    Canvas* canvas[ntabs];
    QWidget * pTab[ntabs];
    QLabel *pLabel[ntabs][8][3];
    QDoubleSpinBox *pTime;
    QComboBox *pPlanet[3];
    QToolButton *bPColor, *bSave, *bClear;
};

#endif
