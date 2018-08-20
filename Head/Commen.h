#ifndef  MINI_COMMEN
#define  MINI_COMMEN

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

//采样OpenMP并行处理		//如果采用并行处理，字符输出可能顺序混乱
//#pragma omp parallel for
			for (int i = 0; i < 10; i++)
			{		//从左到右
				for (int j = 0; j < 10; j++)
				{		//从下到上
					renderRegion(i, j, spp_1);
//#pragma omp atomic
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
			return c >= 'a'&&c <= 'z' || c >= 'A'&&c <= 'Z';
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
					int sqrt_spp = std::round(std::sqrt(_samplePerPixel));
					for (int s1 = 0; s1 < sqrt_spp; s1++)
						for (int s2 = 0; s2 < sqrt_spp; s2++) {
							r1 = (s1 + RandNumber()) / float(sqrt_spp);
							r2 = (s2 + RandNumber())/ float(sqrt_spp);
							float screenX = ( (s1 + r1) / float(sqrt_spp) + i + x*sub_wid) / _width;
							float screenY = ( (s1 + r2) / float(sqrt_spp) + j + y*sub_hei) / _height;

							Ray ray = _camera->getRay(screenX, screenY, r1, r2);
							col = col + rayTrace(ray, 0)*spp_1;
						}

					//Gamma 校正
					col = Power(col);

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
			if (depth >= _depth) return color(0,0,0);
			
			//select closest model;
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
				return _bgColor;

			//俄罗斯轮盘赌
			vector<float> mat_col = inter._material->getColor();
			float russian_p = Max(mat_col);
			if (depth > 5)
				if(RandNumber() < russian_p)
					mat_col = mat_col / russian_p;
				else
					return inter._material->getEmmision();

			if (inter._type == DIFFUSE) {		//积分漫反射材质

				vector<float> u, v;
				if (inter._nor._x == 0) u = Normalize(Cross(vector<float>(1, 0, 0), inter._nor));
				else u = Normalize(Cross(vector<float>(0, 1, 1), inter._nor));
				v = Normalize(Cross(inter._nor, u));

				//采用cos加权方法采样	：单位圆采样投至单位球
				size_t brdf_samples = 4;
				float brdf = 1.0f / PI;
				vector<float> color;
				for (size_t i = 0; i < brdf_samples; i++)
				{
					float rd = RandNumber();
					float phi = (float(i) + RandNumber() ) / float(brdf_samples) * 2 * PI;		//phi进行分层采样
					float r = pow(rd, 0.5);
					vector<float> dir = u*r*cos(phi) + v*r*sin(phi) + inter._nor * std::sqrt(1 - r*r);

					float cos_theta = Dot(inter._nor, dir);
					float pdf = cos_theta / PI;
					color =color + inter._material->getEmmision() + mat_col *rayTrace(Ray(inter._pos, dir), depth + 1) * cos_theta * brdf / pdf;
				}

				return color / float(brdf_samples);
			}
			else if (inter._type == MIRROR){			//积分镜面反射材质	针对Alpha分布，不采用MC积分
				vector<float> dir = ray._direction - inter._nor * 2 * Dot(ray._direction, inter._nor);
				return mat_col * rayTrace(Ray(inter._pos, dir), depth + 1);
			}
			else if (inter._type == TRANSP){			//积分透射材质		针对Alpha分布，不采用MC积分
				float fresnel0 = inter._material->getFresnel();
				float fresnel = fresnel0 + (1 - fresnel0)*(pow(1 + Dot(ray._direction, inter._nor), 5));		//反射项所占比例
				vector<float> reflDir = ray._direction - inter._nor * 2 * Dot(ray._direction, inter._nor);
				float reftRate = inter._material->getReftrectRate();
				vector<float> reftDir;
				float co = Dot(ray._direction, inter._nor);
				float si2 = 1 - co*co;
				if (inter._model->inModel(ray._origin))
				{
					if (si2 > 1 / (1.5*1.5))return rayTrace(Ray(inter._pos, reflDir), depth + 1);
					reftDir = (ray._direction - inter._nor * Dot(ray._direction, inter._nor)) * reftRate - inter._nor*pow(1 - (si2 *reftRate*reftRate), 0.5);
				}
				else
				{
					reftDir = (ray._direction - inter._nor * Dot(ray._direction, inter._nor)) / reftRate - inter._nor*pow(1 - si2 / (reftRate*reftRate), 0.5);
				}
				return mat_col * (rayTrace(Ray(inter._pos, reflDir), depth + 1)*fresnel + rayTrace(Ray(inter._pos, reftDir), depth + 1)*(1 - fresnel) );
			}
			//还需要改进
			else if (inter._type == COOK){			//积分Cook材质	

				vector<float> u, v;
				if (inter._nor._x == 0) u = Normalize(Cross(vector<float>(1, 0, 0), inter._nor));
				else u = Normalize(Cross(vector<float>(0, 1, 1), inter._nor));
				v = Normalize(Cross(inter._nor, u));

				float roughness = inter._material->getRoughness();
				vector<float> color;
				size_t brdf_samples = 4;
				for (size_t i = 0; i < brdf_samples; i++)
				{
					//重要性采样  根据D的分布来进行采样		参考：PBRT第三版808页		//个人感觉有错误，该pdf并没有在半球范围内进行归一化
					float rd = RandNumber();
					float phi = (float(i) + RandNumber()) / float(brdf_samples) * 2.f * PI;		//分层采样
					float tan2_theta = -roughness* roughness*std::log(rd);
					float cos_theta_h = 1.f/std::sqrt(1.f + tan2_theta);
					float sin_theta_h = std::sqrt(max(0.f, 1.0f - cos_theta_h * cos_theta_h));
					vector<float> H = Normalize(u*sin_theta_h*std::cos(phi) + v * sin_theta_h*std::sin(phi) + inter._nor * cos_theta_h);
					vector<float> refldir = ray._direction + H * Dot(-ray._direction, H) * 2.f;

					//计算BRDF
					float g1 = 2 * Dot(H, inter._nor)*Dot(inter._nor, -ray._direction) / Dot(-ray._direction, H);
					float g2 = 2 * Dot(H, inter._nor)*Dot(inter._nor, refldir) / Dot(-ray._direction, H);
					float G = g1 < g2 ? g1 : g2;
					G = 1 < g1 ? 1 : g1;
					G = G < 0 ? 0 : G;

					float fresnel0 = inter._material->getFresnel();
					float fresnel = fresnel0 + (1 - fresnel0)*(pow(1 + Dot(ray._direction, inter._nor), 5));		//反射项所占比例

					float dot_N_H = Dot(inter._nor, H);
					float D = std::exp(-tan2_theta / (roughness*roughness)) / (roughness*roughness * std::pow(cos_theta_h, 4)* PI);

					float brdf = D * G*fresnel / (4.f*Dot(inter._nor, -ray._direction)*Dot(inter._nor, refldir));
					if (brdf < 0.f) brdf = 0.f;

					//计算 PDF
					float pdf_h = D * g1*Dot(H, -ray._direction) / Dot(inter._nor, -ray._direction);
					float pdf = pdf_h / (4.f * Dot(H, -ray._direction));

					vector<float> incolor = rayTrace(Ray(inter._pos, refldir), depth + 1);
					color = color + inter._material->getEmmision() +
						mat_col * incolor * brdf *Dot(inter._nor, refldir) / pdf;
				}
				return color / float(brdf_samples);
			}
			else if (inter._type == BLIN) {			//积分Blin微面元材质	

				vector<float> u, v;
				if (inter._nor._x == 0) u = Normalize(Cross(vector<float>(1, 0, 0), inter._nor));
				else u = Normalize(Cross(vector<float>(0, 1, 1), inter._nor));
				v = Normalize(Cross(inter._nor, u));

				float power = inter._material->getReflectionPower();
				vector<float> color;
				size_t brdf_samples = 4;
				for (size_t i = 0; i < brdf_samples; i++)
				{
					//重要性采样  根据D的分布来进行采样
					float rd = RandNumber();
					float phi = (float(i) + RandNumber() ) / float(brdf_samples) * 2.f * PI;
					float cos_theta_h = std::pow(rd, 1.0f/ (power + 1.f));
					float sin_theta_h = std::sqrt(max(0.f, 1.0f - cos_theta_h * cos_theta_h));
					vector<float> H = Normalize(u*sin_theta_h*std::cos(phi) + v * sin_theta_h*std::sin(phi) + inter._nor * cos_theta_h);
					vector<float> refldir = ray._direction + H * Dot(-ray._direction, H) * 2.f;

					//计算 PDF
					float pdf = (power + 1.f)*std::pow(cos_theta_h, power) / (4.f * Dot(-ray._direction, H) * 2.f * PI);

					//计算BRDF
					float g1 = 2 * Dot(H, inter._nor)*Dot(inter._nor, -ray._direction) / Dot(-ray._direction, H);
					float g2 = 2 * Dot(H, inter._nor)*Dot(inter._nor, refldir) / Dot(-ray._direction, H);
					float G = g1 < g2 ? g1 : g2;
					G = 1 < g1 ? 1 : g1;
					G = G < 0 ? 0 : G;

					float fresnel0 = inter._material->getFresnel();
					float fresnel = fresnel0 + (1 - fresnel0)*(pow(1 + Dot(ray._direction, inter._nor), 5));		//反射项所占比例

					float dot_N_H = Dot(inter._nor, H);
					float D = (power + 2.f)*std::pow(dot_N_H, power) / (2.f*PI);

					float brdf = D * G*fresnel / (4.f*Dot(inter._nor, -ray._direction)*Dot(inter._nor, refldir));
					if (brdf < 0.f) brdf = 0.f;
					vector<float> incolor = rayTrace(Ray(inter._pos, refldir), depth + 1);
					color = color + inter._material->getEmmision() +
						mat_col * incolor * brdf *Dot(inter._nor, refldir) / pdf;
				}
				return color / float(brdf_samples);
			}
			else
			{
				return _bgColor;
			}
		};
	};
}

#endif