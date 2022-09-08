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

#ifndef SCRIBBLE_H
#define SCRIBBLE_H

#include <cmath>
#include <limits>

#include <QPen>
#include <QThread>
#include <QResizeEvent>
#include <QMouseEvent>
#include <QTimerEvent>
#include <QPaintEvent>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QLabel>
#include <QtWidgets/QWidget>
#include <QPolygon>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QCheckBox>


class QMouseEvent;
class QResizeEvent;
class QPaintEvent;
class QToolButton;
class QDoubleSpinBox;

typedef double real;

constexpr const real PI = 3.14159265358979323846; //std::acos(real(-1));
constexpr const real C = 299792458.L;
constexpr const real G = 6.67428e-11L;
constexpr const real H[] = {C*C/(2*G), 0., 1e20};


size_t const N = 3;

struct vector3
{
    real elem_[N];

    vector3() noexcept
	{
	}
	
    vector3(const real & b1, const real & b2, const real & b3) noexcept
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

    real & operator [] (const size_t n)
	{ 
		return elem_[n]; 
	}

    const real & operator [] (const size_t n) const
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
    char const * n;						// name
	QColor c;							// color
    vector3 o, p;                       // old & new positions

    real alpha = real();
    real r = real();
    real x = real();
    real y = real();
    real vel = real();
    real err = real();
    real omega = real();
    real mass = real();

    Planet(char const * n, const QColor c)
        : n(n), c(c), o(), p()
	{
	}

    void execute(const real & dt);
};


class Scribble;

class Dual : public QThread
{
public:
    Dual(Scribble *);

    virtual void run();
};


class Canvas : public QWidget
{
    Q_OBJECT

public:
    Canvas( int type, QWidget *parent );

protected slots:
    void slotForceUpdate();

protected:
    void mousePressEvent( QMouseEvent *e );
    void mouseReleaseEvent( QMouseEvent *e );
    void mouseMoveEvent( QMouseEvent *e );
    void resizeEvent( QResizeEvent *e );
    void paintEvent( QPaintEvent *e );
    void timerEvent( QTimerEvent *e );

    int type;

    QPen pen;
    QPolygon polyline;

    bool first, mousePressed;

    QPixmap buffer;
};


int const ntabs = 3;
int const nt = 6;
int const np = 200;

class Scribble : public QMainWindow
{
    Q_OBJECT

public:
    Scribble( QWidget *parent, const char *name );
    ~Scribble();

protected slots:
    void slotRestart();
    void slotAbout();
	
public:
    static char * const theory[nt];

    real ntime;

	QTabWidget *pTabWidget;
    Canvas* canvas[ntabs];
    QWidget * pTab[ntabs];
    QDoubleSpinBox *pTime;
    QCheckBox *pTheory[nt];
    QToolButton *bPColor, *bSave, *bClear;

    static real const t;
    static real const rmax;
    static real const rmin;
    static real const h;
    static real const r0;
    static real const v0;
    static real const rdm0;
    static real const dmf;
    static real const fit;
    static real const emax;
    static real const massf;
    static real const sping;
    real vmax = 0;
    real mmax = 0;

    std::vector< std::vector<Planet> > planet;
};

#endif
