#ifndef MINI_FUN
#define MINI_FUN

#include "structure.h"
#include <random>
#include<cmath>

namespace mini
{
#pragma region FloatOrVector
	//点乘
	template<typename T>
	T Dot(const vector<T> &vec1, const vector<T> &vec2)
	{
		return vec1._x * vec2._x + vec1._y * vec2._y + vec1._z * vec2._z;
	}

	//范围限制
	template<typename T>
	T Clamp(T t, T x, T y)
	{
		return t < x ? x : (t > y ? y : t);
	}

	//范围限制
	template<typename T>
	vector<T> Clamp(vector<T> t, T x, T y)
	{
		vector<T> res;
		res._x = Clamp(t._x, x, y);
		res._y = Clamp(t._y, x, y);
		res._z = Clamp(t._z, x, y);
		return res;
	}
	
	//clamp to (0-1)
	template<typename T>
	T Saturate(T t)
	{
		return t < 0 ? 0: (t > 1 ? 1 : t);
	}	
	template<typename T>
	vector<T> Saturate(vector<T> t)
	{
			vector<T> res;
			res._x = Saturate(t._x);
			res._y = Saturate(t._y);
			res._z = Saturate(t._z);
			return res;
		}

	//长度计算
	float Length(const vector<float> &vec)
	{
		return sqrt(Dot(vec, vec));
	}

	float Distance(const vector<float> &vec0, const vector<float> &vec1)
	{
		vector<float> vec = vec1 - vec0;
		return Length(vec);
	}

	//长度平方
	float LenSquare(const vector<float> &vec)
	{
		return vec._x* vec._x + vec._y*vec._y + vec._z*vec._z;
	}

	//叉乘
	vector<float> Cross(const vector<float> &vec0, const vector<float> &vec1)
	{
		vector<float> vec;
		vec._x = vec0._y * vec1._z - vec0._z * vec1._y;
		vec._y = vec0._z * vec1._x - vec0._x * vec1._z;
		vec._z = vec0._x * vec1._y - vec0._y * vec1._x;
		return vec;
	}

	//归一化
	vector<float> Normalize(const vector<float>& vec0)
	{
		vector<float> vec;
		float len = Length(vec0);
		if (len != 0.0f)
		{
			float rhw = 1.0f / len;
			vec._x = vec0._x * rhw;
			vec._y = vec0._y * rhw;
			vec._z = vec0._z * rhw;
		}
		return vec;
	}

	float Max(const vector<float>& vec0)
	{
		float e = vec0._x > vec0._y ? vec0._x : vec0._y;
		return e > vec0._z ? e : vec0._z;
	}

	//计算行列式 纵向或横向结果一样
	float Determinant(const vector<float> &po1, const vector<float> &po2, const vector<float> &po3)
	{
		return Dot(Cross(po1, po2), po3);
	}

	//幂运算，用于Gamma校正
	vector<float> Power(const vector<float> & vec0, float po)
	{
		vector<float> vec;
		vec._x = std::pow(vec0._x , po);
		vec._y = std::pow(vec0._y, po);
		vec._z = std::pow(vec0._z, po);

		return vec;
	}
#pragma endregion

#pragma region Other
	std::default_random_engine generator;
	std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
	//return 0-1
	float RandNumber() {
		return distribution(generator);
	}
#pragma endregion
}

#endif // !MINI_FUN
