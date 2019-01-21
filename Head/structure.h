#ifndef MINI_STRUC
#define MINI_STRUC

#include "windowDecl.h"

namespace mini
{
	//ÈýÎ¬Ê¸Á¿
	template<typename T>
	struct vector
	{
		T _x, _y, _z;
		vector(){ _x = _y = _z = 0; }
		vector(const vector& vec) :_x(vec._x), _y(vec._y), _z(vec._z){ }
		vector(T x, T y, T z) :_x(x), _y(y), _z(z){ }
		void operator = (const vector &vec)
		{
			_x = vec._x; _y = vec._y;
			_z = vec._z;
		}
		vector operator + (const vector &vec) const
		{
			vector vecr(_x + vec._x, _y + vec._y, _z + vec._z);
			return vecr;
		}
		vector operator - (const vector &vec) const
		{
			vector vecr(_x - vec._x, _y - vec._y, _z - vec._z);
			return vecr;
		}
		vector operator-() const
		{
			return vector(-_x, -_y, -_z);
		}
		vector operator * (const vector &vec) const
		{
			vector vecr(_x * vec._x, _y * vec._y, _z * vec._z);
			return vecr;
		}
		vector operator * (float a)const
		{
			return vector(_x * a, _y * a, _z * a);
		}
		vector operator / (float a)const
		{
			return vector(_x / a, _y / a, _z / a);
		}
	};
	typedef vector<float> color;
	typedef vector<float> normal;
	typedef vector<float> point;
};

#endif
