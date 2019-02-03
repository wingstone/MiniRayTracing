#ifndef  MINI_CORE
#define  MINI_CORE
#define PI 3.14159265358979323846f
#define E 2.718281828459f
#define DELTA 0.0001f

#include <cfloat>
#include <windows.h>
#include <fstream>

#include "structure.h"
#include "function.h"

namespace mini
{
	//Material type
	enum RefType
	{
		DIFFUSE,
		MIRROR, TRANSP,		//Alpha distribution
		GGX,
		REFLECTION, EMPTY		
	};

	//Ray struct
	struct Ray
	{
		point _origin;
		vector<float> _direction;

		Ray(point origin, vector<float> dir) :_origin(origin), _direction(Normalize(dir)){}
	};

	//=========================================
	//material

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

		//indir is view dir, outdir is light(direct or indirect) dir
		virtual vector<float> Simple_BRDF(const vector<float>& inDir, const float rand1, const float rand2, const vector<float>& nor, float *pdf, float *brdf)
		{ return vector<float>(0.f, 0.f, 0.f); }
		virtual float CalculateBRDF(const vector<float>& nor, const vector<float>& inDir, const vector<float>& outDir, const vector<float>& H)
		{
			return 0.f;
		}
		virtual float CalculatePDF(const vector<float>& nor, const vector<float>& inDir, const vector<float>& outDir, const vector<float>& H)
		{
			return 0.f;
		}

		//main for ggx
		virtual void Simple_BSDF(const vector<float>& inDir, const float rand1, const float rand2, const vector<float>& nor,
			float *pdf_r, float *brdf, float *pdf_t, float *btdf, float *fresnel,
			vector<float> *rdir, vector<float> *tdir)
		{
		}
	};

	//漫反射材质
	class DiffuseMaterial : public Material
	{
	public:
		DiffuseMaterial(RefType reftype, color emmision, color materialCol) :Material(reftype), _emmision(emmision), _materialCol(materialCol)
		{
			useCosSample = false;
		}
		~DiffuseMaterial(){}

		virtual color getColor(){ return _materialCol; }
		virtual color getEmmision(){ return _emmision; }
		virtual vector<float> Simple_BRDF(const vector<float>& inDir, const float rand1, const float rand2, const vector<float>& nor, float *pdf, float *brdf)
		{
			vector<float> w = nor;
			w = Normalize(w);
			vector<float> u = w._z <= 0.5f ? Cross(vector<float>(0.f, 0.f, 1.f), w) : Cross(vector<float>(0.f, 1.f, 0.f), w);
			u = Normalize(u);
			vector<float> v = Cross(w, u);
			v = Normalize(v);


			vector<float> dir;
			if (!useCosSample)
			{
				//半球采样 
				float cosTheta = rand1;
				float sinTheta = sqrt(1 - cosTheta * cosTheta);
				float phi = 2.f * PI*rand2;
				dir = u * sinTheta*cos(phi) + v * sinTheta* sin(phi) + w * cosTheta;
				*pdf = 1.f / (2.f*PI);
			}
			else
			{
				//余弦采样
				float redius = sqrt(rand1);
				float sinTheta = redius;
				float cosTheta = sqrt(1 - rand1);
				float phi = 2 * PI*rand2;
				dir = u * sinTheta*cos(phi) + v * sinTheta* sin(phi) + w * cosTheta;
				*pdf = cosTheta / PI;
			}

			*brdf = 1.f / PI;
			return dir;
		}
		virtual float CalculateBRDF(const vector<float>& nor, const vector<float>& inDir, const vector<float>& outDir, const vector<float>& H)
		{
			return 1.f / PI;
		}
		virtual float CalculatePDF(const vector<float>& nor, const vector<float>& inDir, const vector<float>& outDir, const vector<float>& H)
		{
			if (!useCosSample)
			{
				return 1.f / (2.f*PI);
			}
			else
			{
				float cosTheta = Dot(outDir, nor);
				return cosTheta / PI;
			}
		}

	private:
		color _emmision;
		color _materialCol;
		bool useCosSample;
	};

	class OrenNayerMaterial : public Material
	{
	public:
		OrenNayerMaterial(RefType reftype, color emmision, color materialCol, float roughness = 0.f, bool useCosSample = true)
			:Material(reftype), _emmision(emmision), _materialCol(materialCol), _roughness(roughness), _useCosSample(useCosSample)
		{
			float roughSqure = _roughness * _roughness;
			A = 1.f - 0.5f*roughSqure / (roughSqure + 0.33);
			B = 0.45f*roughSqure / (roughSqure + 0.09f);
		}
		~OrenNayerMaterial() {}

		virtual color getColor() { return _materialCol; }
		virtual color getEmmision() { return _emmision; }
		virtual vector<float> Simple_BRDF(const vector<float>& inDir, const float rand1, const float rand2, const vector<float>& nor, float *pdf, float *brdf)
		{
			vector<float> w = nor;
			w = Normalize(w);
			vector<float> u = w._z <= 0.5f ? Cross(vector<float>(0.f, 0.f, 1.f), w) : Cross(vector<float>(0.f, 1.f, 0.f), w);
			u = Normalize(u);
			vector<float> v = Cross(w, u);
			v = Normalize(v);


			vector<float> dir;
			if (!_useCosSample)
			{
				//半球采样 
				float cosTheta = rand1;
				float sinTheta = sqrt(1 - cosTheta * cosTheta);
				float phi = 2.f * PI*rand2;
				dir = u * sinTheta*cos(phi) + v * sinTheta* sin(phi) + w * cosTheta;
				*pdf = 1.f / (2.f*PI);
			}
			else
			{
				//余弦采样
				float redius = sqrt(rand1);
				float sinTheta = redius;
				float cosTheta = sqrt(1 - rand1);
				float phi = 2 * PI*rand2;
				dir = u * sinTheta*cos(phi) + v * sinTheta* sin(phi) + w * cosTheta;
				*pdf = cosTheta / PI;
			}

			*brdf = 1.f / PI;
			return dir;
		}
		virtual float CalculateBRDF(const vector<float>& nor, const vector<float>& inDir, const vector<float>& outDir, const vector<float>& H)
		{
			//reference https://en.wikipedia.org/wiki/Oren%E2%80%93Nayar_reflectance_model
			float theta_i = std::acos(Dot(nor, outDir));
			float theta_r = std::acos(Dot(nor, -inDir));
			float alpha = max(theta_i, theta_r);
			float beta = min(theta_i, theta_r);

			vector<float> v1 = outDir - nor * cos(theta_i);
			vector<float> v2 = -inDir - nor * cos(theta_r);
			float cos_delt_phi = Dot(Normalize(v1), Normalize(v2));

			return 1.f / PI*(A+B*max(0, cos_delt_phi)*sin(alpha)*tan(beta));
		}
		virtual float CalculatePDF(const vector<float>& nor, const vector<float>& inDir, const vector<float>& outDir, const vector<float>& H)
		{
			if (!_useCosSample)
			{
				return 1.f / (2.f*PI);
			}
			else
			{
				float cosTheta = Dot(outDir, nor);
				return cosTheta / PI;
			}
		}

	private:
		color _emmision;
		color _materialCol;
		float _roughness;
		bool _useCosSample;
		float A, B;
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
		virtual vector<float> Simple_BRDF(const vector<float>& inDir, const float rand1, const float rand2, const vector<float>& nor, float *pdf, float *brdf)
		{
			vector<float> w = nor;
			w = Normalize(w);
			vector<float> u = w._z <= 0.5f ? Cross(vector<float>(0.f, 0.f, 1.f), w) : Cross(vector<float>(0.f, 1.f, 0.f), w);
			u = Normalize(u);
			vector<float> v = Cross(w, u);
			v = Normalize(v);

			//use isotropic Beckmann distribution
			float logSimple = rand1 == 0 ? 0 : log(rand1);
			float tan2ThetaH = -_roughness * _roughness*logSimple;
			float cosThetaH = 1.f / sqrt(1 + tan2ThetaH);
			float sinThetaH = sqrt(1 - cosThetaH * cosThetaH);
			float phi = 2.f * PI*rand2;
			vector<float> H = u * sinThetaH*cos(phi) + v * sinThetaH* sin(phi) + w * cosThetaH;
			
			vector<float> dir = inDir - H * Dot(inDir, H) * 2.f;

			float roughness2 = _roughness * _roughness;
			float D = exp(-tan2ThetaH / roughness2) / (PI*roughness2* pow(cosThetaH, 4.f));		//there are more different distribution models in pbrt

			float g1 = 2 * Dot(H, nor)*Dot(nor, -inDir) / Dot(-inDir, H);
			float g2 = 2 * Dot(H, nor)*Dot(nor, dir) / Dot(-inDir, H);
			float G = g1 < g2 ? g1 : g2;																									//there are more different masking models in pbrt
			G = 1 < g1 ? 1 : g1;
			G = G < 0 ? 0 : G;

			float fresnel = _fresnel + (1 - _fresnel)*(pow(1 + Dot(inDir, nor), 5));		//反射项所占比例
			*brdf = D * G*fresnel / (4.f*Dot(-inDir, nor)*Dot(nor, dir));

			*pdf = D * cosThetaH / (4.f* Dot(dir, H));

			if (Dot(dir, nor) < 0.f) *pdf = -1.f;

			return dir;
		}
		virtual float CalculateBRDF(const vector<float>& nor, const vector<float>& inDir, const vector<float>& outDir, const vector<float>& H)
		{
			float roughness2 = _roughness * _roughness;
			float cosThetaH = Dot(H, nor);
			float cosThetaH2 = cosThetaH * cosThetaH;
			float tan2ThetaH = (1.f - cosThetaH2) / cosThetaH2;
			float D = exp(-tan2ThetaH / roughness2) / (PI*roughness2* pow(cosThetaH, 4.f));		//there are more different distribution models in pbrt

			float g1 = 2 * Dot(H, nor)*Dot(nor, -inDir) / Dot(-inDir, H);
			float g2 = 2 * Dot(H, nor)*Dot(nor, outDir) / Dot(-inDir, H);
			float G = g1 < g2 ? g1 : g2;																									//there are more different masking models in pbrt
			G = 1 < g1 ? 1 : g1;
			G = G < 0 ? 0 : G;

			float fresnel = _fresnel + (1 - _fresnel)*(pow(1 + Dot(inDir, nor), 5));		//反射项所占比例
			return D * G*fresnel / (4.f*Dot(-inDir, nor)*Dot(nor, outDir));
		}
		virtual float CalculatePDF(const vector<float>& nor, const vector<float>& inDir, const vector<float>& outDir, const vector<float>& H)
		{
			float roughness2 = _roughness * _roughness;
			float cosThetaH = Dot(H, nor);
			float cosThetaH2 = cosThetaH * cosThetaH;
			float tan2ThetaH = (1.f - cosThetaH2) / cosThetaH2;
			float D = exp(-tan2ThetaH / roughness2) / (PI*roughness2* pow(cosThetaH, 4.f));

			return  D * cosThetaH / (4.f* Dot(outDir, H));
		}

	private:
		color _emission;
		color _materialCol;
		float _fresnel;
		float _roughness;
	};

	//使用Blin选项的微表面材质	//早期发现的材质，比较常用，因为易计算
	class BlinMaterial : public Material
	{
	public:
		BlinMaterial(RefType reftype, color emission, color materialCol, float fresnel, float reflectionPower)
			:Material(reftype), _emission(emission), _materialCol(materialCol), _fresnel(fresnel), _reflectionPower(reflectionPower){}
		~BlinMaterial(){}

		virtual color getColor(){ return _materialCol; }
		virtual color getEmmision() { return _emission; }
		virtual float getFresnel() { return _fresnel; }
		virtual float getReflectionPower(){ return _reflectionPower; }
		virtual vector<float> Simple_BRDF(const vector<float>& inDir, const float rand1, const float rand2, const vector<float>& nor, float *pdf, float *brdf)
		{
			
			float phi = 2.f*PI*rand1;
			float cosThetaH = pow(rand2, 1.f/(_reflectionPower + 1.f));
			float sinThetaH = sqrt(1.f - cosThetaH * cosThetaH);

			vector<float> w = nor;
			w = Normalize(w);
			vector<float> u = w._z <= 0.5f ? Cross(vector<float>(0.f, 0.f, 1.f), w) : Cross(vector<float>(0.f, 1.f, 0.f), w);
			u = Normalize(u);
			vector<float> v = Cross(w, u);
			v = Normalize(v);

			vector<float> H = u * sinThetaH*cos(phi) + v * sinThetaH* sin(phi) + w * cosThetaH;

			*brdf = (_reflectionPower + 2.f)*pow(Dot(H, nor), _reflectionPower) / (2.f*PI);

			vector<float> dir = inDir - H * Dot(inDir, H) * 2.f;

			*pdf = (_reflectionPower + 1.f) * pow(cosThetaH, _reflectionPower) / (2.f * PI * 4.f * Dot(dir, H));

			if (Dot(dir, nor) < 0.f) *pdf = -1.f;

			return dir;
		}
		virtual float CalculateBRDF(const vector<float>& nor, const vector<float>& inDir, const vector<float>& outDir, const vector<float>& H)
		{
			return (_reflectionPower + 2.f)*pow(Dot(H, nor), _reflectionPower) / (2.f*PI);
		}
		virtual float CalculatePDF(const vector<float>& nor, const vector<float>& inDir, const vector<float>& outDir, const vector<float>& H)
		{
			float cosThetaH = Dot(H, nor);
			return  (_reflectionPower + 1.f) * pow(cosThetaH, _reflectionPower) / (2.f * PI * 4.f * Dot(outDir, H));
		}

	private:
		color _emission;
		color _materialCol;
		float _fresnel;
		float _reflectionPower;
	};

	// http://120.52.51.13/www.cs.cornell.edu/~srm/publications/EGSR07-btdf.pdf
	class GGXMaterial : public Material
	{
	public:
		GGXMaterial(RefType type, color emission, color matCol, float fresnel, float roughness, float reftIndex)
			:Material(type), _emission(emission), _materialCol(matCol), _fresnel(fresnel), _roughness(roughness), _reftIndex(reftIndex){};
		~GGXMaterial() {};

		virtual color getColor() { return _materialCol; }
		virtual float getFresnel() { return _fresnel; }
		virtual color getEmmision() { return _emission; }
		virtual float getRoughness() { return _roughness; }
		virtual vector<float> Simple_BRDF(const vector<float>& inDir, const float rand1, const float rand2, const vector<float>& nor, float *pdf, float *brdf)
		{
			vector<float> w = nor;
			w = Normalize(w);
			vector<float> u = w._z <= 0.5f ? Cross(vector<float>(0.f, 0.f, 1.f), w) : Cross(vector<float>(0.f, 1.f, 0.f), w);
			u = Normalize(u);
			vector<float> v = Cross(w, u);
			v = Normalize(v);

			//use GGX distribution
			float tan2ThetaH = _roughness * _roughness*rand1 / (1.f - rand1);
			float cosThetaH = 1.f / sqrt(1 + tan2ThetaH);
			float sinThetaH = sqrt(1 - cosThetaH * cosThetaH);
			float phi = 2.f * PI*rand2;
			vector<float> H = u * sinThetaH*cos(phi) + v * sinThetaH* sin(phi) + w * cosThetaH;

			vector<float> dir = inDir - H * Dot(inDir, H) * 2.f;

			float roughness2 = _roughness * _roughness;
			float D = roughness2 / (PI* pow(cosThetaH, 4.f)*(roughness2 + tan2ThetaH));

			float cosThetaIn2 = Dot(dir, nor)*Dot(dir, nor);
			float cosThetaOut2 = Dot(-inDir, nor)*Dot(-inDir, nor);
			float tan2ThetaIn = (1.f - cosThetaIn2) / cosThetaIn2;
			float tan2ThetaOut = (1.f - cosThetaOut2) / cosThetaOut2;

			float GIn = 2.f / (1.f + sqrt(1.f + roughness2 * tan2ThetaIn));
			float GOut = 2.f / (1.f + sqrt(1.f + roughness2 * tan2ThetaOut));
			float G = GIn * GOut;

			float fresnel = _fresnel + (1 - _fresnel)*(pow(1 + Dot(inDir, nor), 5));		//反射项所占比例
			*brdf = D * G*fresnel / (4.f*Dot(-inDir, nor)*Dot(nor, dir));

			*pdf = D * cosThetaH / (4.f* Dot(dir, H));

			if (Dot(dir, nor) < 0.f) *pdf = -1.f;

			return dir;
		}
		virtual float CalculateBRDF(const vector<float>& nor, const vector<float>& inDir, const vector<float>& outDir, const vector<float>& H)
		{
			float roughness2 = _roughness * _roughness;
			float cosThetaH = Dot(H, nor);
			float cosThetaH2 = cosThetaH * cosThetaH;
			float tan2ThetaH = (1.f - cosThetaH2) / cosThetaH2;
			float D = roughness2 / (PI* pow(cosThetaH, 4.f)*(roughness2 + tan2ThetaH));	

			float cosThetaIn2 = Dot(outDir, nor)*Dot(outDir, nor);
			float cosThetaOut2 = Dot(-inDir, nor)*Dot(-inDir, nor);
			float tan2ThetaIn = (1.f - cosThetaIn2) / cosThetaIn2;
			float tan2ThetaOut = (1.f - cosThetaOut2) / cosThetaOut2;

			float GIn = 2.f / (1.f + sqrt(1.f + roughness2 * tan2ThetaIn));
			float GOut = 2.f / (1.f + sqrt(1.f + roughness2 * tan2ThetaOut));
			float G = GIn * GOut;

			float fresnel = _fresnel + (1 - _fresnel)*(pow(1 + Dot(inDir, nor), 5));		//反射项所占比例
			return D * G*fresnel / (4.f*Dot(-inDir, nor)*Dot(nor, outDir));
		}
		virtual float CalculatePDF(const vector<float>& nor, const vector<float>& inDir, const vector<float>& outDir, const vector<float>& H)
		{
			float roughness2 = _roughness * _roughness;
			float cosThetaH = Dot(H, nor);
			float cosThetaH2 = cosThetaH * cosThetaH;
			float tan2ThetaH = (1.f - cosThetaH2) / cosThetaH2;
			float D = roughness2 / (PI* pow(cosThetaH, 4.f)*(roughness2 + tan2ThetaH));

			return  D * cosThetaH / (4.f* Dot(outDir, H));
		}

		virtual void Simple_BSDF(const vector<float>& inDir, const float rand1, const float rand2, const vector<float>& nor,
			float *pdf_r, float *brdf, float *pdf_t, float *btdf, float *fresnel,
			vector<float> *rdir, vector<float> *tdir)
		{
			vector<float> w = nor;
			w = Normalize(w);
			vector<float> u = w._z <= 0.5f ? Cross(vector<float>(0.f, 0.f, 1.f), w) : Cross(vector<float>(0.f, 1.f, 0.f), w);
			u = Normalize(u);
			vector<float> v = Cross(w, u);
			v = Normalize(v);

			//use GGX distribution
			//reflection
			float tan2ThetaH = _roughness * _roughness*rand1/(1.f-rand1);
			float cosThetaH = 1.f / sqrt(1 + tan2ThetaH);
			float sinThetaH = sqrt(1 - cosThetaH * cosThetaH);
			float phi = 2.f * PI*rand2;
			vector<float> H = u * sinThetaH*cos(phi) + v * sinThetaH* sin(phi) + w * cosThetaH;

			*rdir = inDir - H * Dot(inDir, H) * 2.f;

			float roughness2 = _roughness * _roughness;
			float D = roughness2 / (PI* pow(cosThetaH, 4.f)*(roughness2 + tan2ThetaH));

			float cosThetaIn2 = Dot(*rdir, nor)*Dot(*rdir, nor);
			float cosThetaOut2 = Dot(-inDir, nor)*Dot(-inDir, nor);
			float tan2ThetaIn = (1.f - cosThetaIn2) / cosThetaIn2;
			float tan2ThetaOut = (1.f - cosThetaOut2) / cosThetaOut2;

			float GIn = 2.f / (1.f + sqrt(1.f + roughness2 * tan2ThetaIn));
			float GOut = 2.f / (1.f + sqrt(1.f + roughness2 * tan2ThetaOut));
			float G = GIn * GOut;

			*fresnel = _fresnel + (1 - _fresnel)*(pow(1 + Dot(inDir, nor), 5));		//反射项所占比例
			*brdf = D * G*(*fresnel) / (4.f*Dot(-inDir, nor)*Dot(nor, *rdir));

			*pdf_r = D * cosThetaH / (4.f* Dot(*rdir, H));

			//refraction
			float co = Dot(inDir, H);
			float si2 = 1.f - co * co;
			if(Dot(H, inDir) <= 0)
			{

				*tdir = (inDir - H * Dot(inDir, H)) / _reftIndex - H *pow(1 - si2 / (_reftIndex*_reftIndex), 0.5f);
			}
			else
			{
				if (si2 > 1.f / (_reftIndex*_reftIndex)) ;
				*tdir = (inDir - H * Dot(inDir, H)) * _reftIndex - H*pow(1 - (si2 *_reftIndex*_reftIndex), 0.5f);
			}

			float factor = abs(Dot(*tdir, H)*Dot(-inDir, H) / (Dot(*tdir, nor)*Dot(-inDir, nor)));
			*btdf = factor * (1.f - (*fresnel))*G*D / pow(_reftIndex*Dot(*tdir, H) + Dot(-inDir, H), 2);

			float eta = Dot(H, inDir) <= 0 ? 1.f/_reftIndex : _reftIndex;
			float sqrtDenom = Dot(*rdir, H) + eta * Dot(-inDir, H);
			float dwh_dwi =
				std::abs((eta * eta * Dot(-inDir, H)) / (sqrtDenom * sqrtDenom));
			*pdf_t = D * cosThetaH * dwh_dwi;
		}

	private:
		color _emission;
		color _materialCol;
		float _fresnel;
		float _roughness;
		float _reftIndex;
	};

	//=========================================
	//model、light

	//前置声明
	class Model;

	//Intersection class
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

	//Model and light base class
	class Model
	{
	public:
		Model(){};
		virtual ~Model(){};

		//virtual bool intersectAndSetT(Ray ray) = 0;
		virtual bool inModel(point po) = 0;
		virtual Intersection getIntersection(Ray ray) = 0;
		//rand范围：0~1
		virtual color Sample_L(float rand1, float rand2, point rayOrigin, vector<float> *lightRay, float* pdf) { return color(); };
		virtual float CalculatePdf(const point& lightOrigin, float t, const vector<float>& lightDir) { return -1.f; };
		virtual float GetArea() { return 0; };
		virtual bool IsLight()
		{
			return false;
		}
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
				_t = pro + pow(_redius*_redius - inLenSquare, 0.5f);
				po = ray._origin + ray._direction*_t;
				_normal = Normalize(_center - po);
			}
			else
			{
				_t = pro - pow(_redius*_redius - inLenSquare, 0.5f);
				po = ray._origin + ray._direction*_t;
				_normal = Normalize(po - _center);
			}

			if (_t <= 0.f){
				return Intersection::_empty;
			}
			else
			{
				return Intersection(ray._origin + ray._direction*_t, _normal, _t-DELTA*0.5f, _material->_reftype, _material, this);
			}
		}

	protected:
		point _center;
		float _redius;
		Material* _material;
	};

	//Disk_light
	class Disk_Light : public Model
	{
	public:
		Disk_Light(point center, float redius, vector<float> direction, bool islight, Material*  material)
			:_center(center), _redius(redius), _isLight(islight), _material(material) 
		{
			_direction = Normalize(direction);
		}
		~Disk_Light() {};

		bool inModel(point po) {
			return false;
		}

		Intersection getIntersection(Ray ray)
		{
			if (Dot(ray._direction, _direction) > DELTA)
				return Intersection::_empty;

			float d = -Dot(_center, _direction);

			float t = -(d + Dot(_direction, ray._origin)) / Dot(_direction, ray._direction);
			if (t <= 0.f)
				return Intersection::_empty;

			point interPoint = ray._origin + ray._direction*t;
			return Length(interPoint - _center) > _redius + DELTA ? Intersection::_empty : Intersection(interPoint, _direction, t - DELTA*0.5f, _material->_reftype, _material, this);
		}

		//light

		bool IsLight()
		{
			return _isLight;
		}

		virtual color Sample_L(float rand1, float rand2, point rayOrigin, vector<float> *lightRay, float* pdf)
		{
			vector<float> u, v;
			if (_direction._x == 0) u = Normalize(Cross(vector<float>(1, 0, 0), _direction));
			else u = Normalize(Cross(vector<float>(0, 1, 1), _direction));
			v = Normalize(Cross(_direction, u));

			float phi = rand2* 2 * PI;		//phi进行分层采样
			float r = std::sqrt(rand1)*_redius;
			point pos =  _center + u * r*cos(phi) + v * r*sin(phi);

			*lightRay = pos - rayOrigin;
			 vector<float> dir = Normalize(*lightRay);
			 if (Dot(-dir, _direction) <= 0)		//面光源具有方向性
				 *pdf = -1.f;
			 else
				 *pdf = Dot(*lightRay, *lightRay) /( GetArea()*Dot(-dir, _direction));

			 return _material->getEmmision();
		}

		virtual float CalculatePdf(const point& lightOrigin, float t, const vector<float>& lightDir)
		{
				return  t*t / (GetArea()*Dot(-lightDir, _direction));
		}

		float GetArea()
		{
			return PI * _redius*_redius;
		}

	private:
		point _center;
		float _redius;
		vector<float> _direction;
		Material* _material;

		bool _isLight;
	};

	//Sphere_light
	class Sphere_Light : public Sphere
	{
	public:
		Sphere_Light(point center, float redius,  bool islight, Material*  material)
			: Sphere(center, redius, material), _isLight(islight)
		{
		}
		~Sphere_Light() {};

		bool IsLight()
		{
			return _isLight;
		}

		virtual color Sample_L(float rand1, float rand2, point rayOrigin, vector<float> *lightRay, float* pdf)
		{
			vector<float> u, v, w;
			w = Normalize(_center - rayOrigin);
			u = w._z <= 0.5f ? Cross(vector<float>(0.f, 0.f, 1.f), w) : Cross(vector<float>(0.f, 1.f, 0.f), w);
			u = Normalize(u);
			v = Cross(w, u);
			v = Normalize(v);

			//针对球体对应的体积角进行采样
			float phi = rand2 * 2.f * PI;		//phi进行分层采样
			float sinThetaMax = _redius / Distance(rayOrigin, _center);
			float cosThetaMin = sqrt(1.f - sinThetaMax * sinThetaMax);
			float z = 1.f - (1.f - cosThetaMin)*rand1;
			float x = std::cos(phi)*std::sqrt(1.f - z * z);
			float y = std::sin(phi)*std::sqrt(1.f - z * z);

			*lightRay = u * x + v * y + w * z;
			*pdf = 1.f / (2.f * PI * (1.f - cosThetaMin));

			return _material->getEmmision();
		}

		virtual float CalculatePdf(const point& lightOrigin, float t, const vector<float>& lightDir) 
		{
			float sinThetaMax = _redius / Distance(lightOrigin, _center);
			float cosThetaMin = sqrt(1.f - sinThetaMax * sinThetaMax);
			return 1.f / (2.f * PI * (1.f - cosThetaMin));
		}

		float GetArea()
		{
			return PI * _redius*_redius * 4.f;
		}

	private:
		bool _isLight;
	};

	//axis align Box
	class Box : public Model
	{
	public:
		Box(point minPo, point maxPo, Material* mate) :_min(minPo), _max(maxPo), _material(mate){}

		bool inModel(point po)
		{
			return abs(po._x - _min._x + _max._x - po._x) + abs(po._y - _min._y + _max._y - po._y) + abs(po._z - _min._z + _max._z - po._z)
				< _max._x - _min._x + _max._y - _min._y + _max._z - _min._z + DELTA*0.5f;
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

			if (t <= 0.f){
				return Intersection::_empty;
			}
			else
			{
				return Intersection(ray._origin + ray._direction*t, nor, t - DELTA*0.5f, _material->_reftype, _material, this);
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
			float t = -(_d + Dot(_normal, ray._origin)) / Dot(_normal, ray._direction);
			return t <= 0.f ? Intersection::_empty : Intersection(ray._origin + ray._direction*t, _normal, t-DELTA*0.5f, _material->_reftype, _material, this);
		}

	private:
		float _d;
		normal _normal;
		Material* _material;
	};

	//以逆时针顺序初始化 法线由右手准则确定 仅支持凸四边形
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
			return Dot(_normal, _po1 - po) >= -DELTA;
		}

		Intersection getIntersection(Ray ray)
		{
			normal nor =  Dot(ray._direction, _normal) > 0 ? -_normal :_normal;
			float t = -(_d + Dot(_normal, ray._origin)) / Dot(_normal, ray._direction);
			if (t <= DELTA * 0.25f || inModel(ray._origin))
				return Intersection::_empty;
			else
			{
				point po = ray._origin + ray._direction*t;
				vector<float> v1 = Normalize(_po1 - po);
				vector<float> v2 = Normalize(_po2 - po);
				vector<float> v3 = Normalize(_po3 - po);
				vector<float> v4 =Normalize( _po4 - po);
				if (acos(Dot(v1,v2)) + acos(Dot(v2,v3))+acos(Dot(v3,v4))+acos(Dot(v4,v1)) > 2*PI- DELTA)			//根据夹角相加是否等于2PI来判断是否在多边形内部
					return Intersection(ray._origin + ray._direction*t, nor, t - DELTA*0.5f, _material->_reftype, _material, this);
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
			return Dot(_normal, _po1 - po) > -DELTA;
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
				return t <= DELTA * 0.25f ? Intersection::_empty : Intersection(ray._origin + ray._direction*t, nor, t - DELTA*0.5f, _material->_reftype, _material, this);
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

	//=========================================
	//camera

	//Camera; 原点为左下角
	class Camera
	{
	public:
		Camera(){};
		virtual~Camera(){};

		//x,y~(0,1)	图片左下角为原点
		virtual Ray getRay(float screenX, float screenY, float randx, float randy) = 0;
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
		Ray getRay(float screenX, float screenY, float randx, float randy)
		{
			vector<float> dir = _up * std::tan(_fov / 2) *(screenY*2.0f-1.0f) + _right * std::tan(_fov / 2)*_aspect*(screenX*2.0f-1.0f) + _front;		//射线从origin发射
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
		OrthoCamera(float height, float aspect, point origin, vector<float> front, vector<float> up)
			: _height(height), _aspect(aspect),  _origin(origin)
		{
			_right = Normalize(Cross(front, up));
			_up = Normalize(Cross(_right, front));
			_front = Normalize(Cross(_up, _right));
		}
		~OrthoCamera() {}

		//归一化屏幕坐标，范围：0~1
		Ray getRay(float screenX, float screenY, float randx, float randy)
		{
			point origin = _origin + _right * (_aspect*(screenX*2.f - 1.f)) + _up * (screenY*2.f - 1.f);		//射线从像素发射
			return Ray(origin, _front);
		}

	private:
		float _height;
		float _aspect;		//width/height
		point _origin;
		vector<float> _front;
		vector<float> _up;
		vector<float> _right;
	};

	//Pinhole Camera		//针对pinhole相机，不在从屏幕像素中采样，而直接从透镜中采样。
	class PinholeCamera :public Camera
	{
	public:
		PinholeCamera(float apertureRedius, float focusDistance, float imageDistance, float height, float aspect, point origin, vector<float> front, vector<float> up)
			: _apertureRedius(apertureRedius), _focusDistance(focusDistance), _imageDistance(imageDistance), _height(height), _aspect(aspect), _origin(origin)
		{
			_right = Normalize(Cross(front, up));
			_up = Normalize(Cross(_right, front));
			_front = Normalize(Cross(_up, _right));
		}
		~PinholeCamera() {}

		//针孔成像后为倒立图像，此处并没有对倒立过程进行纠正
		Ray getRay(float screenX, float screenY, float randx, float randy)
		{
			point screenPoint = _origin + _right * (_aspect*(screenX*2.f - 1.f)*_height) + _up * ((screenY*2.f - 1.f)*_height);
			point lensPoint = _origin + _front * _imageDistance;
			point focusPoint = lensPoint + (lensPoint - screenPoint) * (_focusDistance / _imageDistance);

			point origin = lensPoint + _right * (std::sqrt(randx) * std::cos(randy*2.f*PI)) + _up * (std::sqrt(randx) * std::sin(randy*2.f*PI)) *_apertureRedius;	//射线从透镜发射
			vector<float> dir = Normalize(focusPoint - origin);
			return Ray(origin, dir);
		}

	private:
		//透镜描述
		float _apertureRedius;	//光圈半径
		float _focusDistance;		//焦距
		float _imageDistance;		//像距

		//画布描述
		float _height;
		float _aspect;

		point _origin;
		vector<float> _front;
		vector<float> _up;
		vector<float> _right;
	};

	//Enviroment Camera		//环境相机， aspect将不起作用
	class EnviromentCamera :public Camera
	{
	public:
		EnviromentCamera(point origin, vector<float> front, vector<float> up)
			: _origin(origin)
		{
			_right = Normalize(Cross(front, up));
			_up = Normalize(Cross(_right, front));
			_front = Normalize(Cross(_up, _right));
		}
		~EnviromentCamera() {}

		Ray getRay(float screenX, float screenY, float randx, float randy)
		{
			float theta = PI*(1.f - screenY);
			float phi = 2*PI * screenX;
			vector<float> dir = _right * (std::sin(theta)*std::cos(phi)) + _front * (std::sin(theta)*std::sin(phi)) + _up * std::cos(theta);
			return Ray(_origin, dir);
		}

	private:

		point _origin;
		vector<float> _front;
		vector<float> _up;
		vector<float> _right;
	};

}

#endif