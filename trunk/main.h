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

#include <qcolor.h>
#include <q3mainwindow.h>
#include <qpen.h>
#include <qpoint.h>
#include <qlabel.h>
#include <qpixmap.h>
#include <qwidget.h>
#include <qstring.h>
#include <q3pointarray.h>
#include <qthread.h>
#include <qcombobox.h>
#include <qtabwidget.h>
//Added by qt3to4:
#include <QResizeEvent>
#include <QMouseEvent>
#include <QTimerEvent>
#include <QPaintEvent>

class QMouseEvent;
class QResizeEvent;
class QPaintEvent;
class QToolButton;
class QSpinBox;

typedef long double real;

const real G = 6.67428e-11;
const real h = 1.3450632e27/2;

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
	static real FR(real m, real d);
	static real GR(real m, real d);
	static real NW(real m, real d);

	char const * n;
	QColor c;
	real m;
	vector3 p;
	vector3 v;
	vector3 a;
	vector3 o;
	real t[2];
	bool updated;
	//vector3 i;
	vector3 s;
	real pd;
	vector3 pp[2];
	vector3 ps[2];
	real (* f)(real, real);	
	
	enum Type {PP, LB} eType;

	Planet(char const * n, const QColor & c, real m, const real pp[3], const real pv[3], real (* f)(real, real) = NW, Type eType = PP)
	: n(n), c(c), m(m), p(pp[0], pp[1], pp[2]), v(pv[0], pv[1], pv[2]), updated(false), pd(std::numeric_limits<real>::max()), f(f), eType(eType)
	{
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
	enum Type {PP, LB} eType;

    Canvas( Type eType, QWidget *parent = 0, const char *name = 0 );
    ~Canvas();
    void clearScreen();
	
protected slots:
	void slotPlanet(int);
	
protected:
    void mousePressEvent( QMouseEvent *e );
    void mouseReleaseEvent( QMouseEvent *e );
    void mouseMoveEvent( QMouseEvent *e );
    void resizeEvent( QResizeEvent *e );
    void paintEvent( QPaintEvent *e );
	void timerEvent( QTimerEvent *e );

    QPen pen;
    Q3PointArray polyline;

    bool mousePressed;

    QPixmap buffer;

	std::vector< std::vector<Planet> > planet;
	
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

class Scribble : public Q3MainWindow
{
    Q_OBJECT
	friend class Canvas;

public:
    Scribble( QWidget *parent = 0, const char *name = 0 );

protected slots:
    void slotClear();
    void slotPlanet(int);
	void slotPP();
	void slotLB();
	void slotChanged(QWidget * p);
	void slotAbout();
	
public:
	unsigned nc, ntime[2];

	QTabWidget *pTabWidget;
    Canvas* canvas[2];
	QWidget * pTab[2];
	QLabel *pLabel[2][7][3];
    QSpinBox *pTime;
    QComboBox *pPlanet;
    QToolButton *bPColor, *bSave, *bClear;
};

#endif
