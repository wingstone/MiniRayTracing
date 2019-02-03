#pragma once

//#define TRANDITION_TRACING_METHOD

#include <vector>
#include <sstream>
#include <iostream>

#include "Core.h"
#include "bmp.h"

using namespace mini;

namespace mini{


	//Scene
	class Scene
	{
	public:
		Scene(int width, int height, int spp, int depth, color bgColor):_samplePerPixel(spp), _depth(depth),_bgColor(bgColor){
			_materialList = new std::vector<Material*>();
			_modelList = new std::vector<Model*>();
			_lightList = new std::vector<Model*>();
			_camera = NULL;
			_buffer = new UCHAR[width*height * 3];
			_bmp = NULL;

			_width = width; _height = height;
		}
		~Scene(){
			for (auto i = _materialList->begin(); i < _materialList->end(); i++)
			{
				delete *i;
			}
			for (auto i = _modelList->begin(); i < _modelList->end(); i++)
			{
				delete *i;
			}
			delete _materialList;
			delete _modelList;
			delete _lightList;
			if (_buffer != NULL) delete[] _buffer;
			if (_bmp != NULL) delete _bmp;
		}

		//向场景添加模型，内存由Scene对象释放
		void addModel(Model* model){
			_modelList->push_back(model);
		}

		//设置相机
		void setCamera(Camera* camera){
			delete _camera;
			_camera = camera;
		}

		//从文件中初始化场景
		void initFromFile(const char* fileName)
		{
			int stack = 0;
			std::fstream filestream(fileName, std::fstream::in);
			if (!filestream.is_open())
			{
				std::cout << "Read " << fileName << " failed~" << std::endl;
				throw std::exception("Can't find rt file");
			}

			char c;
			char cstr[256];	//注释缓存
			std::string str;

			while (filestream.get(c))
			{
				if (isAlpha(c))
				{
					filestream.putback(c);
					str = getNextWord(filestream);
					if (str.compare("Material") == 0)
					{
						stack = 0;
						str = std::string();
						str.reserve(2048);
						while (filestream.get(c))
						{
							str.push_back(c);
							if (c == '{') stack++;
							if (c == '}') stack--;
							if (stack == 0) break;
						}
						str.shrink_to_fit();
						getMaterial(std::stringstream(str));
					}
					else if (str.compare("Model") == 0)
					{
						stack = 0;
						str = std::string();
						str.reserve(2048);
						while (filestream.get(c))
						{
							str.push_back(c);
							if (c == '{') stack++;
							if (c == '}') stack--;
							if (stack == 0) break;
						}
						str.shrink_to_fit();
						getModel(std::stringstream(str));
					}
					else if (str.compare("Camera") == 0)
					{
						str = getNextWord(filestream);
						if (str.compare("PerspectiveCamera") == 0)
						{
							float asp;
							float fov;
							vector<float> up, fron;
							point po;

							fov = getNextNumber(filestream);
							asp = getNextNumber(filestream);
							po = getNextVector(filestream);
							fron = getNextVector(filestream);
							up = getNextVector(filestream);
							_camera = new PerspectiveCamera(fov, asp, po, fron, up );
						}
						else if (str.compare("OrthoCamera") == 0)
						{
							float height = getNextNumber(filestream);
							float aspect = getNextNumber(filestream);
							point origin = getNextVector(filestream);
							vector<float> front = getNextVector(filestream);
							vector<float> up = getNextVector(filestream);
							_camera = new OrthoCamera(height, aspect, origin, front, up);
						}
						else if (str.compare("PinholeCamera") == 0)
						{
							float apertureRedius = getNextNumber(filestream);
							float focusDistance = getNextNumber(filestream);	
							float imageDistance = getNextNumber(filestream);
							float height = getNextNumber(filestream);
							float aspect = getNextNumber(filestream);
							point origin = getNextVector(filestream);
							vector<float> front = getNextVector(filestream);
							vector<float> up = getNextVector(filestream);
							_camera = new PinholeCamera(apertureRedius, focusDistance, imageDistance, height, aspect, origin, front, up);
						}
						else if (str.compare("EnviromentCamera") == 0)
						{
							point origin = getNextVector(filestream);
							vector<float> front = getNextVector(filestream);
							vector<float> up = getNextVector(filestream);
							_camera = new EnviromentCamera(origin, front, up);
						}
					}
					else if (str.compare("Light") == 0)
					{
						stack = 0;
						str = std::string();
						str.reserve(2048);
						while (filestream.get(c))
						{
							str.push_back(c);
							if (c == '{') stack++;
							if (c == '}') stack--;
							if (stack == 0) break;
						}
						str.shrink_to_fit();
						getLight(std::stringstream(str));
					}
				}
				else if (c == '/')
				{
					filestream.get(c);
					if (c == '/')
						filestream.getline(cstr, 256);
				}
			}

			filestream.close();
		}

		//绘制和写入文件RayTrace.bmp
		void renderAndWrite(const char* filename)
		{
			if (_bmp != NULL) delete _bmp;
			_bmp = new BmpClass(filename);
			float spp_1 = float(1)/_samplePerPixel;
			int sum = 0;
			std::cout << "\nRendering prograss: ";
			std::cout.width(4);
			std::cout << sum;
			std::cout << "%";

			for (int i = 0; i < 10; i++)
			{		//从左到右
				for (int j = 0; j < 10; j++)
				{		//从下到上
					renderRegion(i, j, spp_1);
					sum++;
					std::cout << "\b\b\b\b\b";
					std::cout.width(4);
					std::cout << sum;
					std::cout << "%";
					std::cout.flush();
				}
			}

			_bmp->writeFile(_buffer, _width, _height, _width*_height * 3);
			std::cout << "\nRendering and White success~" << std::endl;
		}

	private:
		std::vector<Material*>* _materialList;
		std::vector<Model*>*  _modelList;
		std::vector<Model*>* _lightList;		//光源同时存在于模型列表与光源列表中
		Camera* _camera;
		UCHAR* _buffer;
		BmpClass* _bmp;

		color _bgColor;
		int _width;
		int _height;
		int _samplePerPixel;
		int _depth;

		//file function
		inline bool isAlpha(char c)
		{
			return c >= 'a'&&c <= 'z' || c >= 'A'&&c <= 'Z' || c == '_';
		}
		inline bool isNumber(char c)
		{
			return c >= '0' && c <= '9' || c == '-' || c == '.';
		}
		std::string getNextWord(std::iostream& sstream)
		{
			bool haveWordFlag = false;
			std::string str;
			str.reserve(20);
			char c;
			char cstr[256];

			while (sstream.get(c))
			{
				if (isAlpha(c))
				{
					str.push_back(c);
					haveWordFlag = true;
				}
				else if (c == '/')
				{
					sstream.get(c);
					if (c == '/')
						sstream.getline(cstr, 256);
				}
				else if (haveWordFlag)
				{
					sstream.putback(c);
					break;
				}
				else
				{
					continue;
				}
			}
			str.shrink_to_fit();
			return str;
		}
		float getNextNumber(std::iostream& sstream)
		{
			float num;
			char c;
			char cstr[256];

			while (sstream.get(c))
			{
				if (isNumber(c))
				{
					sstream.putback(c);
					sstream >> num;
					return num;
				}
				else if (c == '/')
				{
					sstream.get(c);
					if (c == '/')
						sstream.getline(cstr, 256);
				}
				else
				{
					continue;
				}
			}
		}
		vector<float> getNextVector(std::iostream& sstream)
		{
			vector<float> vec;
			vec._x = getNextNumber(sstream);
			vec._y = getNextNumber(sstream);
			vec._z = getNextNumber(sstream);
			return vec;
		}
		void getMaterial(std::iostream& sstream)
		{
			std::string str;
			while (true)
			{
				str = getNextWord(sstream);

				if (str.compare("DiffuseMaterial") == 0)
				{
					color matCol, emmision;
					emmision = getNextVector(sstream);
					matCol = getNextVector(sstream);

					_materialList->push_back(new DiffuseMaterial(DIFFUSE, emmision, matCol));
				}
				else if (str.compare("ReflectMaterial") == 0)
				{
					color matCol;
					matCol = getNextVector(sstream);
					_materialList->push_back(new MirrorMaterial(MIRROR, matCol));
				}
				else if (str.compare("TranspMaterial") == 0)
				{
					float fresnel, reftIndex;
					color matCol;
					fresnel = getNextNumber(sstream);
					reftIndex = getNextNumber(sstream);
					matCol = getNextVector(sstream);
					_materialList->push_back(new TramsparentMaterial(TRANSP, fresnel, reftIndex, matCol));
				}
				else if (str.compare("CookMaterial") == 0)
				{
					color emission, matcol;
					float fresnel, roughness;
					emission = getNextVector(sstream);
					matcol = getNextVector(sstream);
					fresnel = getNextNumber(sstream);
					roughness = getNextNumber(sstream);
					_materialList->push_back(new CookTorranceMaterial(COOK, emission, matcol, fresnel, roughness));
				}
				else if (str.compare("BlinMaterial") == 0)
				{
					color emission, matcol;
					float fresnel, reflectionPower;
					emission = getNextVector(sstream);
					matcol = getNextVector(sstream);
					fresnel = getNextNumber(sstream);
					reflectionPower = getNextNumber(sstream);
					_materialList->push_back(new BlinMaterial(BLIN, emission, matcol, fresnel, reflectionPower));
				}
				else
					break;
			}


		}
		void getModel(std::iostream& sstream)
		{
			std::string str;
			while (true)
			{
				str = getNextWord(sstream);
				if (str.compare("Sphere") == 0)
				{
					point center;
					float redius;
					int matId;
					center = getNextVector(sstream);
					redius = getNextNumber(sstream);
					matId = (int)getNextNumber(sstream);
					_modelList->push_back(new Sphere(center, redius, _materialList->at(matId)));
				}
				else if (str.compare("Box") == 0)
				{
					point minPo;
					point maxPo;
					int matId;
					minPo = getNextVector(sstream);
					maxPo = getNextVector(sstream);
					matId = (int)getNextNumber(sstream);
					_modelList->push_back(new Box(minPo, maxPo, _materialList->at(matId)));
				}
				else if (str.compare("Plane") == 0)
				{
					point center;
					normal nor;
					int matId;
					center = getNextVector(sstream);
					nor = getNextVector(sstream);
					matId = (int)getNextNumber(sstream);
					_modelList->push_back(new Plane(center, nor, _materialList->at(matId)));
				}
				else if (str.compare("Triangle") == 0)
				{
					point po1;
					point po2;
					point po3;
					int matId;
					po1 = getNextVector(sstream);
					po2 = getNextVector(sstream);
					po3 = getNextVector(sstream);
					matId = (int)getNextNumber(sstream);
					_modelList->push_back(new Triangle(po1, po2, po3, _materialList->at(matId)));
				}
				else if (str.compare("Polygon") == 0)
				{
					point po1;
					point po2;
					point po3;
					point po4;
					int matId;
					po1 = getNextVector(sstream);
					po2 = getNextVector(sstream);
					po3 = getNextVector(sstream);
					po4 = getNextVector(sstream);
					matId = (int)getNextNumber(sstream);
					_modelList->push_back(new Polygon(po1, po2, po3, po4, _materialList->at(matId)));
				}
				else
				{
					break;
				}
			}

		}
		void getLight(std::iostream& sstream)
		{
			std::string str;
			while (true)
			{
				str = getNextWord(sstream);
				if (str.compare("Disk_Light") == 0)
				{
					point center;
					float redius;
					vector<float> nor;
					int matId;
					center = getNextVector(sstream);
					redius = getNextNumber(sstream);
					nor = getNextVector(sstream);
					matId = (int)getNextNumber(sstream);
					Disk_Light* disklight = new Disk_Light(center, redius, nor, true, _materialList->at(matId));
					_lightList->push_back(disklight);
					_modelList->push_back(disklight);
				}
				else if (str.compare("Sphere_Light") == 0)
				{
					point center;
					float redius;
					int matId;
					center = getNextVector(sstream);
					redius = getNextNumber(sstream);
					matId = (int)getNextNumber(sstream);
					Sphere_Light* sphLight = new Sphere_Light(center, redius, true, _materialList->at(matId));
					_lightList->push_back(sphLight);
					_modelList->push_back(sphLight);
				}
				else
				{
					break;
				}
			}

		}
	
		//render function
		void renderRegion(int x, int y, float spp_1)
		{

			int sub_wid = _width / 10, sub_hei = _height / 10;
			float r1 = .0f, r2 = .0f;
			for (int i = 0; i < sub_wid; i++)
			{
				for (int j = 0; j < sub_hei; j++)
				{
					color col = color();

					//采用分层采样进行采样,针对于每个像素
					int sqrt_spp = (int)std::round(std::sqrt(_samplePerPixel));
					for (int s1 = 0; s1 < sqrt_spp; s1++)
						for (int s2 = 0; s2 < sqrt_spp; s2++) {
							r1 = (s1 + RandNumber()) / float(sqrt_spp);
							r2 = (s2 + RandNumber())/ float(sqrt_spp);
							float screenX = ( (s1 + r1) / float(sqrt_spp) + i + x*sub_wid) / _width;
							float screenY = ( (s1 + r2) / float(sqrt_spp) + j + y*sub_hei) / _height;

							Ray ray = _camera->getRay(screenX, screenY, r1, r2);
							col = col + rayTrace(ray, _depth);
						}

					//Gamma 校正
					col = col / (float)_samplePerPixel;
					col = Power(col, 1.f/2.2f);

					//文件写入顺序：以bgr顺序写入三个字节
					_buffer[3 * (x * sub_wid + i + (y * sub_hei + j)*_width)] = int(Saturate(col._z) * 255);
					_buffer[3 * (x * sub_wid + i + (y * sub_hei + j)*_width) + 1] = int(Saturate(col._y) * 255);
					_buffer[3 * (x * sub_wid + i + (y * sub_hei + j)*_width) + 2] = int(Saturate(col._x) * 255);
				}
			}

		}
		
		//跟踪射线
		color rayTrace(Ray ray, int depth)
		{
			color col = color();
			vector<float> mask(1.f, 1.f, 1.f);
			bool hitSpeculer = false;

#ifdef TRANDITION_TRACING_METHOD

			if (depth == 0) return _bgColor;

			//select closet model;
			Intersection inter = Intersection::_empty;
			for (auto i = _modelList->begin(); i != _modelList->end(); i++) {

				Intersection tempInter = (*i)->getIntersection(ray);
				if (inter._type != EMPTY && tempInter._type != EMPTY && tempInter._t < inter._t) {
					inter = tempInter;
				}
				else if (inter._type == EMPTY && tempInter._type != EMPTY) {
					inter = tempInter;
				}
			}

			if (inter._type == EMPTY)
			{
				return _bgColor;
			}

			//俄罗斯轮盘赌
			color mat_col = inter._material->getColor();
			float p = Max(mat_col);
			if(_depth -  depth >4)
				if ( RandNumber() < p)
				{
					mat_col = mat_col / p;
				}
				else
				{
					return inter._material->getEmmision();
				}

			if (inter._type == DIFFUSE || inter._type == BLIN || inter._type == COOK) {		//积分漫反射材质

				ray._origin = inter._pos;
				vector<float> z = inter._nor;
				z = Normalize(z);
				vector<float> x = z._z <= 0.5f ? Cross(vector<float>(0.f, 0.f, 1.f), z) : Cross(vector<float>(0.f, 1.f, 0.f), z);
				x = Normalize(x);
				vector<float> y = Cross(z, x);
				y = Normalize(y);

				float rand1 = RandNumber();
				float rand2 = RandNumber();
				float pdf = 0.f;
				float brdf = 0.f;
				vector<float> reflDir = inter._material->Simple_BRDF(ray._direction, rand1, rand2, inter._nor, &pdf, &brdf);
				ray._direction = reflDir;

				float cosTheta = Dot(inter._nor, ray._direction);
				col =  inter._material->getEmmision() + mat_col *rayTrace(ray, --depth)*brdf*cosTheta / pdf;
			}
			else if (inter._type == MIRROR) {			//积分镜面反射材质
				vector<float> dir = ray._direction - inter._nor * 2 * Dot(ray._direction, inter._nor);
				ray._origin = inter._pos;
				ray._direction = Normalize(dir);
				col = mat_col * rayTrace(ray, --depth);
			}
			else if (inter._type == TRANSP) {			//积分透射材质
				float reftIndex = inter._material->getReftrectRate();
				float fresnel0 = inter._material->getFresnel();
				float fresnel = fresnel0 + (1 - fresnel0)*(pow(1 + Dot(ray._direction, inter._nor), 5));		//反射项所占比例
				vector<float> reflDir = ray._direction - inter._nor * 2 * Dot(ray._direction, inter._nor);
				float reftRate = inter._material->getReftrectRate();
				vector<float> reftDir;
				float co = Dot(ray._direction, inter._nor);
				float si2 = 1.f - co * co;
				if (inter._model->inModel(ray._origin))
				{
					if (si2 > 1.f / (reftIndex*reftIndex)) return _bgColor;
					reftDir = (ray._direction - inter._nor * Dot(ray._direction, inter._nor)) * reftRate - inter._nor*pow(1 - (si2 *reftRate*reftRate), 0.5f);
				}
				else
				{
					reftDir = (ray._direction - inter._nor * Dot(ray._direction, inter._nor)) / reftRate - inter._nor*pow(1 - si2 / (reftRate*reftRate), 0.5f);
				}
				ray._origin = inter._pos;

				Ray refRay(inter._pos, reflDir);
				Ray reftRay(inter._pos, reftDir);
				depth--;
				col = mat_col * rayTrace(refRay, depth)*fresnel + mat_col * rayTrace(reftRay, depth)*(1.f- fresnel);
			}
			else
			{
				col = _bgColor;
			}

#else

			//光线每一步进行迭代
			for (int i = 0; i < depth; i++)
			{
				//select closet model;
				Intersection inter = Intersection::_empty;
				for (auto i = _modelList->begin(); i != _modelList->end(); i++) {

					Intersection tempInter = (*i)->getIntersection(ray);
					if (inter._type != EMPTY && tempInter._type != EMPTY&& tempInter._t < inter._t) {
						inter = tempInter;
					}
					else if (inter._type == EMPTY&& tempInter._type != EMPTY){
						inter = tempInter;
					}
				}

				if (inter._type == EMPTY)
				{
					col = _bgColor * mask + col;
					break;
				}

				//return normal for debug
				//if(i == 1)
				//	return (inter._nor + vector<float>(1.f,1.f,1.f)) / 2.f;	

				//第一次击中光源or经过镜面
				if (inter._model->IsLight() && i == 0 || hitSpeculer)
				{
					col = col + mask*inter._material->getEmmision();
					hitSpeculer = false;
				}

				if (inter._type == MIRROR){			//积分镜面反射材质
					hitSpeculer = true;

					vector<float> dir = ray._direction - inter._nor * 2 * Dot(ray._direction, inter._nor);
					ray._origin = inter._pos;
					ray._direction = Normalize(dir);
					mask = mask * inter._material->getColor();
				}
				else if (inter._type == TRANSP){			//积分透射材质
					hitSpeculer = true;

					float reftIndex = inter._material->getReftrectRate();
					float fresnel0 = inter._material->getFresnel();
					float fresnel = fresnel0 + (1 - fresnel0)*(pow(1 + Dot(ray._direction, inter._nor), 5));		//反射项所占比例
					vector<float> reflDir = ray._direction - inter._nor * 2 * Dot(ray._direction, inter._nor);
					float reftRate = inter._material->getReftrectRate();
					vector<float> reftDir;
					float co = Dot(ray._direction, inter._nor);
					float si2 = 1.f - co*co;
					if (inter._model->inModel(ray._origin))
					{
						if (si2 > 1.f / (reftIndex*reftIndex)) break;
						reftDir = (ray._direction - inter._nor * Dot(ray._direction, inter._nor)) * reftRate - inter._nor*pow(1 - (si2 *reftRate*reftRate), 0.5f);
					}
					else
					{
						reftDir = (ray._direction - inter._nor * Dot(ray._direction, inter._nor)) / reftRate - inter._nor*pow(1 - si2 / (reftRate*reftRate), 0.5f);
					}


					ray._origin = inter._pos;
					float pdf = fresnel;
					mask = mask * inter._material->getColor();
					if (RandNumber() <= pdf)
					{
						ray._direction = Normalize(reflDir);
						//mask = mask / pdf * fresnel;			//除以是因为重要性采样，乘以是因为基于菲涅尔现象，刚好抵消
					}
					else
					{
						ray._direction = Normalize(reftDir);
						//mask = mask / (1.f - pdf) * （1.f - fresnel);
					}
				}
				else if (inter._type == BLIN || inter._type == COOK)
				{
					hitSpeculer = false;

					//direct lighting
					if (!inter._model->IsLight())
					{
						size_t light_samples = _lightList->size();		//全光源采样
						color directColor;
						for (size_t i = 0; i < light_samples; i++)
						{
							vector<float> lightDir;
							float brdf;
							float light_pdf;
							float brdf_pdf;
							float weight;
							color lightCol;

							//MIS: light sampling
							lightCol = _lightList->at(i)->Sample_L(RandNumber(), RandNumber(), inter._pos, &lightDir, &light_pdf);
							lightDir = Normalize(lightDir);
							if (light_pdf < 0.f) continue;

							//test visibility
							Intersection shadowInter = Intersection::_empty;
							Ray  shadowRay(inter._pos, lightDir);
							for (auto i = _modelList->begin(); i != _modelList->end(); i++) {
								Intersection tempInter = (*i)->getIntersection(shadowRay);
								if (shadowInter._type != EMPTY && tempInter._type != EMPTY && tempInter._t < shadowInter._t) {
									shadowInter = tempInter;
								}
								else if (shadowInter._type == EMPTY && tempInter._type != EMPTY) {
									shadowInter = tempInter;
								}
							}
							if (shadowInter._model == _lightList->at(i))
							{
								vector<float> nor = inter._nor;
								vector<float> inDir = ray._direction;
								vector<float> outDir = lightDir;
								vector<float> H = Normalize(lightDir - ray._direction);
								brdf = inter._material->CalculateBRDF(nor, inDir, outDir, H);

								brdf_pdf = inter._material->CalculatePDF(nor, inDir, inDir, H);
								weight = light_pdf * light_pdf / (light_pdf*light_pdf + brdf_pdf * brdf_pdf);
								directColor = directColor + mask * weight * lightCol * inter._material->getColor()* Dot(lightDir, inter._nor) * brdf / light_pdf;
							}

							//MIS: brdf sampling
							float rand1 = RandNumber();
							float rand2 = RandNumber();
							lightDir = inter._material->Simple_BRDF(ray._direction, rand1, rand2, inter._nor, &brdf_pdf, &brdf);

							//test visibility
							shadowInter = Intersection::_empty;
							shadowRay._direction = lightDir;
							for (auto i = _modelList->begin(); i != _modelList->end(); i++) {

								Intersection tempInter = (*i)->getIntersection(shadowRay);
								if (shadowInter._type != EMPTY && tempInter._type != EMPTY && tempInter._t < shadowInter._t) {
									shadowInter = tempInter;
								}
								else if (shadowInter._type == EMPTY && tempInter._type != EMPTY) {
									shadowInter = tempInter;
								}
							}
							if (shadowInter._model == _lightList->at(i))
							{
								light_pdf = _lightList->at(i)->CalculatePdf(inter._pos, inter._t, lightDir);
								weight = brdf_pdf * brdf_pdf / (light_pdf*light_pdf + brdf_pdf * brdf_pdf);
								directColor = directColor + mask * weight * lightCol * inter._material->getColor()* Dot(lightDir, inter._nor) * brdf / brdf_pdf;
							}

						}
						col = col + directColor / (float)light_samples;
					}

					//undirect lighting
					ray._origin = inter._pos;
					float rand1 = RandNumber();
					float rand2 = RandNumber();
					float pdf = 0.f;
					float brdf = 0.f;
					vector<float> reflDir = inter._material->Simple_BRDF(ray._direction, rand1, rand2, inter._nor, &pdf, &brdf);
					ray._direction = reflDir;

					float cosTheta = Dot(inter._nor, ray._direction);
					mask = mask * inter._material->getColor() * cosTheta * brdf / pdf;
				}
				else if (inter._type == DIFFUSE) {		//积分漫反射材质
					hitSpeculer = false;

					//direct lighting
					if (!inter._model->IsLight())
					{
						size_t light_samples = _lightList->size();		//全光源采样
						color directColor;
						for (size_t i = 0; i < light_samples; i++)
						{
							vector<float> lightDir;
							float pdf;
							color lightCol = _lightList->at(i)->Sample_L(RandNumber(), RandNumber(), inter._pos, &lightDir, &pdf);
							lightDir = Normalize(lightDir);
							if (pdf < 0.f) continue;

							//test visibility
							Intersection shadowInter = Intersection::_empty;
							Ray  shadowRay(inter._pos, lightDir);
							for (auto i = _modelList->begin(); i != _modelList->end(); i++) {

								Intersection tempInter = (*i)->getIntersection(shadowRay);
								if (shadowInter._type != EMPTY && tempInter._type != EMPTY && tempInter._t < shadowInter._t) {
									shadowInter = tempInter;
								}
								else if (shadowInter._type == EMPTY && tempInter._type != EMPTY) {
									shadowInter = tempInter;
								}
							}
							if (shadowInter._model != _lightList->at(i))
								continue;

							vector<float> nor = inter._nor;
							vector<float> inDir = ray._direction;
							vector<float> outDir = lightDir;
							vector<float> H = Normalize(lightDir - ray._direction);
							float brdf = inter._material->CalculateBRDF(nor, inDir, outDir, H);
							directColor = directColor + mask * lightCol * inter._material->getColor()* Dot(lightDir, inter._nor) * brdf / pdf;
						}
						col = col + directColor / (float)light_samples;
					}

					//undirect lighting
					ray._origin = inter._pos;
					float rand1 = RandNumber();
					float rand2 = RandNumber();
					float pdf = 0.f;
					float brdf = 0.f;
					vector<float> reflDir = inter._material->Simple_BRDF(ray._direction, rand1, rand2, inter._nor, &pdf, &brdf);
					ray._direction = reflDir;

					float cosTheta = Dot(inter._nor, ray._direction);
					mask = mask * inter._material->getColor() * cosTheta * brdf / pdf;
				}
				else
				{
					col = _bgColor * mask + col;
					break;
				}

				//俄罗斯轮盘赌
				color mat_col = inter._material->getColor();
				float p = Max(mat_col);
				if (_depth - depth >4)
					if (RandNumber() < p)
					{
						mask = mask / p;
					}
					else
					{
						break;
					}
			}
			
#endif
			return col;
		};
	};
}