#ifndef  MINI_CORE
#define  MINI_CORE
#define PI 3.14159265358979323846f
#define E 2.718281828459f
#define DELTA 0.0001f

#include <cfloat>
#include <windows.h>

#include "structure.h"
#include "function.h"

namespace mini
{
	//Material type
	enum RefType
	{
		DIFFUSE,
		MIRROR, TRANSP,		//Alpha distribution
		BLIN, COOK, EMPTY		
	};

	//Ray struct
	struct Ray
	{
		point _origin;
		vector<float> _direction;

		Ray(point origin, vector<float> dir) :_origin(origin), _direction(Normalize(dir)){}
	};

	//Material base class
	class Material
	{
	public:
		RefType _reftype;
		Material(RefType reftype) :_reftype(reftype){}
		virtual ~Material(){}

		virtual color getColor(){ return color(); }
		virtual color getEmmision(){ return color(); }
		virtual float getFresnel(){ return 0; }
		virtual float getReftrectRate(){ return 0; }
		virtual float getRoughness(){ return 0; }
		virtual float getReflectionPower() { return 0; }
	};

	//漫反射材质
	class DiffuseMaterial : public Material
	{
	public:
		DiffuseMaterial(RefType reftype, color emmision, color materialCol) :Material(reftype), _emmision(emmision), _materialCol(materialCol){}
		~DiffuseMaterial(){}

		virtual color getColor(){ return _materialCol; }
		virtual color getEmmision(){ return _emmision; }

	private:
		color _emmision;
		color _materialCol;
	};

	//镜面反射材质
	class MirrorMaterial : public Material
	{
	public:
		MirrorMaterial(RefType reftype, color matCol) :Material(reftype), _materialCol(matCol){}
		~MirrorMaterial(){}

		virtual color getColor() { return _materialCol; }

	private:
		color _materialCol;
	};

	//透明材质
	class TramsparentMaterial :public Material
	{
	public:
		TramsparentMaterial(RefType reftype, float fresnel, float reftrectRate, color matCol) 
			:Material(reftype), _fresnel(fresnel), _reftrectRate(reftrectRate), _materialCol(matCol) {};
		~TramsparentMaterial(){};

		virtual float getFresnel(){ return _fresnel; }
		virtual float getReftrectRate(){ return _reftrectRate; }
		virtual color getColor() { return _materialCol; }

	private:
		float _fresnel;
		float _reftrectRate;
		color _materialCol;
	};

	//使用Beckmann分布的微表面反射材质		//参考https://en.wikipedia.org/wiki/Specular_highlight#Cook.E2.80.93Torrance_model
	class CookTorranceMaterial :public Material
	{
	public:
		CookTorranceMaterial(RefType type, color emission, color matCol, float fresnel, float roughness)
			:Material(type), _emission(emission), _materialCol(matCol), _fresnel(fresnel), _roughness(roughness){};
		~CookTorranceMaterial(){};

		virtual color getColor(){ return _materialCol; }
		virtual float getFresnel(){ return _fresnel; }
		virtual color getEmmision(){ return _emission; }
		virtual float getRoughness(){ return _roughness; }

	private:
		color _emission;
		color _materialCol;
		float _fresnel;
		float _roughness;
	};

	//使用Blin选项的微表面材质	//早期发现的材质，比较常用，因为易计算
	class BlinMaterial : public Material		//diffuse+mirror
	{
	public:
		BlinMaterial(RefType reftype, color emission, color materialCol, float fresnel, float reflectionPower)
			:Material(reftype), _emission(emission), _materialCol(materialCol), _fresnel(fresnel), _reflectionPower(reflectionPower){}
		~BlinMaterial(){}

		virtual color getColor(){ return _materialCol; }
		virtual color getEmmision() { return _emission; }
		virtual float getFresnel() { return _fresnel; }
		virtual float getReflectionPower(){ return _reflectionPower; }

	private:
		color _emission;
		color _materialCol;
		float _fresnel;
		float _reflectionPower;
	};

	//Intersection class
	class Model;	//提前声明
	struct Intersection
	{
		point _pos;
		normal _nor;
		float _t;
		RefType _type;
		Material* _material;
		Model* _model;
		static const  Intersection _empty;

		Intersection(point po, normal nor, float t, RefType type, Material* material, Model* model)
			:_pos(po), _nor(nor), _t(t), _type(type), _material(material), _model(model){}
	};
	const Intersection Intersection::_empty = Intersection(point(), normal(), 0, EMPTY, NULL, NULL);

	//Model base class
	class Model
	{
	public:
		Model(){};
		virtual ~Model(){};

		//virtual bool intersectAndSetT(Ray ray) = 0;
		virtual bool inModel(point po) = 0;
		virtual Intersection getIntersection(Ray ray) = 0;

	};

	//Sphere
	class Sphere : public Model
	{
	public:
		Sphere(point center, float redius, Material*  material) :_center(center), _redius(redius), _material(material){}
		~Sphere(){};

		bool inModel(point po){
			return LenSquare(po - _center) <= (_redius + DELTA)*(_redius + DELTA);
		}

		Intersection getIntersection(Ray ray)
		{
			float _t;
			normal _normal;

			vector<float> orTocen = _center - ray._origin;
			float ocLenSquare = LenSquare(orTocen);
			float pro = Dot(ray._direction, orTocen);
			float inLenSquare = ocLenSquare - pro*pro;

			if (inLenSquare > _redius*_redius) return Intersection::_empty;

			point po;
			if (_material->_reftype == TRANSP && this->inModel(ray._origin))
			{
				_t = pro + pow(_redius*_redius - inLenSquare, 0.5);
				po = ray._origin + ray._direction*_t;
				_normal = Normalize(_center - po);
			}
			else
			{
				_t = pro - pow(_redius*_redius - inLenSquare, 0.5);
				po = ray._origin + ray._direction*_t;
				_normal = Normalize(po - _center);
			}

			if (_t < DELTA){
				return Intersection::_empty;
			}
			else
			{
				return Intersection(ray._origin + ray._direction*_t, _normal, _t, _material->_reftype, _material, this);
			}
		}

	private:
		point _center;
		float _redius;
		Material* _material;
	};

	//axis align Box
	class Box : public Model
	{
	public:
		Box(point minPo, point maxPo, Material* mate) :_min(minPo), _max(maxPo), _material(mate){}

		bool inModel(point po)
		{
			return abs(po._x - _min._x + _max._x - po._x) + abs(po._y - _min._y + _max._y - po._y) + abs(po._z - _min._z + _max._z - po._z)
				< _max._x - _min._x + _max._y - _min._y + _max._z - _min._z + DELTA;
		}

		Intersection getIntersection(Ray ray)
		{
			float t;
			normal nor;

			if (abs(ray._direction._x) < DELTA)
			{
				if (ray._origin._x < _min._x || ray._origin._x > _max._x)
				{
					return  Intersection::_empty;
				}
				else
				{
					float tz1 = (_min._z - ray._origin._z) / ray._direction._z;
					float tz2 = (_max._z - ray._origin._z) / ray._direction._z;

					float ty1 = (_min._y - ray._origin._y) / ray._direction._y;
					float ty2 = (_max._y - ray._origin._y) / ray._direction._y;

					InterAssist(-FLT_MAX / 3, FLT_MAX / 3, ty1, ty2, tz1, tz2, &t, &nor, ray._origin);
				}
			}
			else if (abs(ray._direction._y) < DELTA)
			{
				if (ray._origin._y < _min._y || ray._origin._y > _max._y)
				{
					return  Intersection::_empty;
				}
				else
				{
					float tz1 = (_min._z - ray._origin._z) / ray._direction._z;
					float tz2 = (_max._z - ray._origin._z) / ray._direction._z;

					float tx1 = (_min._x - ray._origin._x) / ray._direction._x;
					float tx2 = (_max._x - ray._origin._x) / ray._direction._x;

					InterAssist(tx1, tx2, -FLT_MAX / 3, FLT_MAX / 3, tz1, tz2, &t, &nor, ray._origin);
				}
			}
			else if (abs(ray._direction._z) < DELTA)
			{
				if (ray._origin._z < _min._z || ray._origin._z > _max._z)
				{
					return  Intersection::_empty;
				}
				else
				{
					float ty1 = (_min._y - ray._origin._y) / ray._direction._y;
					float ty2 = (_max._y - ray._origin._y) / ray._direction._y;

					float tx1 = (_min._x - ray._origin._x) / ray._direction._x;
					float tx2 = (_max._x - ray._origin._x) / ray._direction._x;

					InterAssist(tx1, tx2, ty1, ty2, -FLT_MAX / 3, FLT_MAX / 3, &t, &nor, ray._origin);

				}
			}
			else
			{
				float tx1 = (_min._x - ray._origin._x) / ray._direction._x;
				float tx2 = (_max._x - ray._origin._x) / ray._direction._x;

				float ty1 = (_min._y - ray._origin._y) / ray._direction._y;
				float ty2 = (_max._y - ray._origin._y) / ray._direction._y;

				float tz1 = (_min._z - ray._origin._z) / ray._direction._z;
				float tz2 = (_max._z - ray._origin._z) / ray._direction._z;

				InterAssist(tx1, tx2, ty1, ty2, tz1, tz2, &t, &nor, ray._origin);
			}

			if (t < DELTA){
				return Intersection::_empty;
			}
			else
			{
				return Intersection(ray._origin + ray._direction*t, nor, t, _material->_reftype, _material, this);
			}
		}

	private:
		point _min;
		point _max;
		Material* _material;


		void InterAssist(float nx, float fx, float ny, float fy, float nz, float fz, float*const t, normal*const nor, point& po)
		{
			float tmp = 0;
			normal n1 = normal(-1, 0, 0);
			normal n2 = normal(0, -1, 0);
			normal n3 = normal(0, 0, -1);
			if (nx - fx > 0)
			{
				tmp = nx, nx = fx, fx = tmp;
				n1 = normal(1, 0, 0);
			}
			if (ny - fy > 0)
			{
				tmp = ny, ny = fy, fy = tmp;
				n2 = normal(0, 1, 0);
			}
			if (nz - fz > 0)
			{
				tmp = nz, nz = fz, fz = tmp;
				n3 = normal(0, 0, 1);
			}

			//两者之间最远的
			float n = 0;
			if (nx < ny)
			{
				*nor = n2;
				n = ny;
			}
			else
			{
				*nor = n1;
				n = nx;
			}

			if (n < nz)
			{
				*nor = n3;
				n = nz;
			}

			//两者之间最近的
			float f = 0;
			bool isInTransp = (_material->_reftype == TRANSP && inModel(po));
			if (fx > fy)
			{
				f = fy;
				if (isInTransp) *nor = n2;
			}
			else
			{
				f = fx;
				if (isInTransp) *nor = n1;
			}
			if (f > fz)
			{
				f = fz;
				if (isInTransp) *nor = n3;
			}

			if (n > f)
			{
				*t = -1;
			}
			else if (f < 0)
			{
				*t = -1;
			}
			else
			{
				*t = n;
				if (isInTransp) *t = f;
			}
		}


	};

	//equation: Dot(normal, point) + d = 0	;具有方向性的面
	class Plane :public Model
	{
	public:
		//center: point in plane; nor: plane normal
		Plane(point center, normal nor, Material* mate) : _material(mate)
		{
			_normal = Normalize(nor);
			_d = -Dot(center, _normal);
		}
		Plane(const Plane& plane) : _d(plane._d), _normal(plane._normal){}
		~Plane(){}

		//在平面的负法线方向一侧
		bool inModel(point po)
		{
			return Dot(po, _normal) + _d < 0;
		}

		Intersection getIntersection(Ray ray)
		{
			normal nor = Dot(ray._direction, _normal) > 0 ? -_normal : _normal;
			float t = -(_d + Dot(_normal, ray._origin)) / Dot(_normal, ray._direction);
			return t < DELTA ? Intersection::_empty : Intersection(ray._origin + ray._direction*t, nor, t, _material->_reftype, _material, this);
		}

	private:
		float _d;
		normal _normal;
		Material* _material;
	};

	//以逆时针顺序初始化 法线由右手准则确定 仅支持凸多边形
	class Polygon :public Model
	{
	public:
		Polygon(point po1, point po2, point po3, point po4, Material* mate) :_po1(po1), _po2(po2), _po3(po3), _po4(po4), _material(mate)
		{
			_normal = Normalize(Cross(po2 - po1, po3 - po1));
			_d = -Dot(_po1, _normal);
		}
		~Polygon(){}

		//在面的负法线方向上一侧
		bool inModel(point po)
		{
			return Dot(_normal, _po1 - po) > 0;
		}

		Intersection getIntersection(Ray ray)
		{
			normal nor = Dot(ray._direction, _normal) > 0 ? -_normal : _normal;
			float t = -(_d + Dot(_normal, ray._origin)) / Dot(_normal, ray._direction);
			if (t < DELTA) 
				return Intersection::_empty;
			else
			{
				point po = ray._origin + ray._direction*t;
				vector<float> v1 = Normalize(_po1 - po);
				vector<float> v2 = Normalize(_po2 - po);
				vector<float> v3 = Normalize(_po3 - po);
				vector<float> v4 =Normalize( _po4 - po);
				if (acos(Dot(v1,v2)) + acos(Dot(v2,v3))+acos(Dot(v3,v4))+acos(Dot(v4,v1)) > 2*PI- DELTA)			//根据夹角相加是否等于2PI来判断是否在多边形内部
					return Intersection(ray._origin + ray._direction*t, nor, t, _material->_reftype, _material, this);
				else
					return Intersection::_empty;
			}
		}

	private:
		point _po1;
		point _po2;
		point _po3;
		point _po4;
		float _d;
		normal _normal;
		Material* _material;
	};

	//以逆时针顺序初始化 法线由右手准则确定
	class Triangle : public Model
	{
	public:
		Triangle(point po1, point po2, point po3, Material* mate) :_po1(po1), _po2(po2), _po3(po3), _material(mate)
		{ 
			_normal = Normalize(Cross(po2 - po1, po3 - po1));
		}
		~Triangle(){}

		//在三角面的负法线方向上一侧
		bool inModel(point po)
		{
			return Dot(_normal, _po1 - po) > 0;
		}

		Intersection getIntersection(Ray ray)
		{
			float t;

			vector<float> coeffi1 = _po1 - _po2;
			vector<float> coeffi2 = _po1 - _po3;
			vector<float> coeffi3 = ray._direction;
			vector<float> val = _po1 - ray._origin;
			float det = Determinant(coeffi1, coeffi2, coeffi3);
			float beta = Determinant(val, coeffi2, coeffi3) / det;
			float gamma = Determinant(coeffi1, val, coeffi3) / det;
			t = Determinant(coeffi1, coeffi2, val) / det;

			if (beta + gamma <= 1 && beta >= 0 && gamma >= 0)
			{
				normal nor = Dot(ray._direction, _normal) > 0 ? -_normal : _normal;
				return t < DELTA ? Intersection::_empty : Intersection(ray._origin + ray._direction*t, nor, t, _material->_reftype, _material, this);
			}
			else
			{
				return Intersection::_empty;
			}

		}

	private:
		point _po1;
		point _po2;
		point _po3;
		normal _normal;
		Material* _material;
	};

	//Camera
	class Camera
	{
	public:
		Camera(){};
		virtual~Camera(){};

		//x,y~(-1,1)
		virtual Ray getRay(float x, float y) = 0;
	};

	//Perspective Camera
	class PerspectiveCamera :public Camera
	{
	public:
		PerspectiveCamera(float fov,float aspect, point origin, vector<float> front, vector<float> up) : _aspect(aspect), _fov(fov*PI / 180), _origin(origin)
		{
			_right = Normalize(Cross(front, up));
			_up = Normalize(Cross(_right, front));
			_front = Normalize(Cross(_up, _right));
		}
		~PerspectiveCamera(){}

		//归一化屏幕坐标，范围：0~1
		Ray getRay(float x, float y)
		{
			vector<float> dir = _up * std::tan(_fov / 2) *(y*2.0f-1.0f) + _right * std::tan(_fov / 2)*_aspect*(x*2.0f-1.0f) + _front;
			dir = Normalize(dir);
			return Ray(_origin, dir);
		}

	private:
		float _fov;	//竖直方向上夹角
		float _aspect;		//width/height
		point _origin;
		vector<float> _front;
		vector<float> _up;
		vector<float> _right;
	};

	//Orthogonal Camera
	class OrthoCamera :public Camera
	{
	public:
		OrthoCamera(float fov, float aspect, point origin, vector<float> front, vector<float> up) : _aspect(aspect), _fov(fov*PI / 180), _origin(origin)
		{
			_right = Normalize(Cross(front, up));
			_up = Normalize(Cross(_right, front));
			_front = Normalize(Cross(_up, _right));
		}
		~OrthoCamera() {}

		//归一化屏幕坐标，范围：0~1
		Ray getRay(float x, float y)
		{
			vector<float> dir = _up * std::tan(_fov / 2) *(y*2.0f - 1.0f) + _right * std::tan(_fov / 2)*_aspect*(x*2.0f - 1.0f) + _front;
			dir = Normalize(dir);
			return Ray(_origin, dir);
		}

	private:
		float _fov;	//竖直方向上夹角
		float _aspect;		//width/height
		point _origin;
		vector<float> _front;
		vector<float> _up;
		vector<float> _right;
	};

	//Pinhole Camera
	class PinholeCamera :public Camera
	{
	public:
		PinholeCamera(float fov, float aspect, point origin, vector<float> front, vector<float> up) : _aspect(aspect), _fov(fov*PI / 180), _origin(origin)
		{
			_right = Normalize(Cross(front, up));
			_up = Normalize(Cross(_right, front));
			_front = Normalize(Cross(_up, _right));
		}
		~PinholeCamera() {}

		//归一化屏幕坐标，范围：0~1
		Ray getRay(float x, float y)
		{
			vector<float> dir = _up * std::tan(_fov / 2) *(y*2.0f - 1.0f) + _right * std::tan(_fov / 2)*_aspect*(x*2.0f - 1.0f) + _front;
			dir = Normalize(dir);
			return Ray(_origin, dir);
		}

	private:
		float _fov;	//竖直方向上夹角
		float _aspect;		//width/height
		point _origin;
		vector<float> _front;
		vector<float> _up;
		vector<float> _right;
	};

}

#endif