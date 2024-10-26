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
#define BOOST_DISABLE_ASSERTS
#include <boost/multiprecision/cpp_bin_float.hpp>

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

typedef double real;

const real C = 299792458.L;
const real G = 6.67428e-11L;
const real K = 8.987551e9L;
const real Eta = 1e-8; //0.0013342L;
const real Q = 1.602176634e-19L;
const real H[] = {C*C/(G), 0., 1e20};

struct vector3
{
    typedef real T;
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

    vector3 & operator = (const vector3 & b)
	{
		for (size_t i = 0; i < N; ++ i)
			elem_[i] = b.elem_[i];

        return * this;
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
    static real FR(real m, real d, real h);
    static real NW(real m, real d, real h);

    char const * n;						// name
	QColor c;							// color
	real m;								// mass
    real q;                             // charge
	vector3 p;							// position
    vector3 v[2];						// current & saved velocity
    vector3 o;							// old position
    vector3 a;							// acceleration or force
    real t[2];							// current & old time intervals according to Newton or FT
    bool first;                         // first cycle
	bool updated;						// the cycle of the planet or the photon arrival line has been completed
    vector3 pp[2];						// current & old saved positions on the perihelion
    vector3 ps[5];						// current & old polar coordinates of pp
    real (* f)(real, real, real);   	// function pointer to Newton time formula or FT time formula
    real hg, he;                             // fudge factor

    enum Type {PP, LB, BB, GR, V1, NU, QU} eType;		// is for the perihelion precession disparity or the gravitational light bending

    Planet(char const * n, const QColor & c, real m, real q, const real pp[3], const real pv[3], real (* f)(real, real, real), Type eType, real hg, real he)
        : n(n), c(c), m(m), q(q), p(pp[0], pp[1], pp[2]), first(true), updated(false), f(f), eType(eType), hg(hg), he(he)
	{
        v[0] = vector3(pv[0], pv[1], pv[2]);
	}
	
	void operator () (const std::vector<Planet> &p, const real & upper);
};

class Canvas;
	
class Dual : public QThread
{
public:
	Dual(Canvas *, int);
	virtual void run();
	
protected:
	Canvas * p;
	int i;
};


class Canvas : public QWidget
{
    Q_OBJECT
	friend class Dual;

public:
    enum Type {PP, LB, BB, GR, V1, NU, QU} eType;

    Canvas( Type eType, QWidget *parent = 0);
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
    void wheelEvent( QWheelEvent *e );

    QPen pen;
    QPolygon polyline;

    bool mousePressed;

    QPixmap buffer;

    std::vector< std::vector<Planet> > planet;

    real scale = 0.L, zoom = 0.2L;

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
    void slotNU();
    void slotQU();
    void slotChanged(int);
	void slotAbout();
	
public:
    static const unsigned ntabs = 7;

    unsigned nc;
    real ntime[ntabs];

	QTabWidget *pTabWidget;
    Canvas* canvas[ntabs];
    QWidget * pTab[ntabs];
    QLabel *pLabel[ntabs][8][3];
    QDoubleSpinBox *pTime;
    QLabel *pScale;
    QLabel *pZoom;
    QComboBox *pPlanet[2];
    QToolButton *bPColor, *bSave, *bClear;
};

#endif
