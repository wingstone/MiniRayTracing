#ifndef  MINI_RTFILE
#define  MINI_RTFILE

#include <sstream>

#include "Core.h"

namespace mini
{
	class RTFileClass
	{
	public:
		RTFileClass();
		~RTFileClass();

		void initFromFile(const char* fileName)
		{
			int stack = 0;
			std::fstream filestream(fileName, std::fstream::in);
			if (filestream.is_open())
			{
				std::cout << "Read " << fileName << " success" << std::endl;
			}
			else
			{
				return;
			}
			char c;
			char cstr[256];	//×¢ÊÍ»º´æ
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
						str.reserve(1024);
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
						str.reserve(1024);
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

							asp = getNextNumber(filestream);
							fov = getNextNumber(filestream);
							po = getNextVector(filestream);
							up = getNextVector(filestream);
							fron = getNextVector(filestream);
							_camera = new PerspectiveCamera(asp, fov, po, up, fron);
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

		void getMaterialList(std::vector<Material*>* materialList)
		{
			for (auto i = _materialList->begin(); i < _materialList->end();i++)
			{
				materialList->push_back(*i);
			}
		}
		void getModelList(std::vector<Model*>* modelList)
		{
			for (auto i = _modelList->begin(); i < _modelList->end(); i++)
			{
				modelList->push_back(*i);
			}
		}
		void getCamera(Camera** camera)
		{
			*camera = _camera;
		}

	private:
		std::vector<Material*>* _materialList;
		std::vector<Model*>*  _modelList;
		Camera* _camera;

		inline bool isAlpha(char c)
		{
			return c >= 'a'&&c <= 'z' || c >= 'A'&&c <= 'Z';
		}
		inline bool isNumber(char c)
		{
			return c >= '0' && c <= '9' || c == '-' ||c == '.';
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
				else if ( haveWordFlag)
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
			char cstr[256];
			std::string str;
			while (true)
			{
				str = getNextWord(sstream);

				if (str.compare("DiffuseMaterial") == 0)
				{
					color matCol, emmision;
					matCol = getNextVector(sstream);
					emmision = getNextVector(sstream);

					_materialList->push_back(new DiffuseMaterial(emmision, matCol, DIFFUSE));
				}
				else if (str.compare("ReflectMaterial") == 0)
				{
					_materialList->push_back(new MirrorMaterial(MIRROR));
					sstream.getline(cstr, 256);
					sstream.getline(cstr, 256);
				}
				else if (str.compare("TranspMaterial") == 0)
				{
					float fresnel, reftIndex;
					fresnel = getNextNumber(sstream);
					reftIndex = getNextNumber(sstream);
					_materialList->push_back(new TramsparentMaterial(fresnel, reftIndex, TRANSP));
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
					matId = getNextNumber(sstream);
					_modelList->push_back(new Sphere(center, redius, _materialList->at(matId)));
				}
				else
				{
					break;
				}
			}

		}
	};

	RTFileClass::RTFileClass()
	{
		_materialList = new std::vector<Material*>();
		_modelList = new std::vector<Model*>();
		_camera = NULL;
	}

	RTFileClass::~RTFileClass()
	{
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
		delete _camera;
	}
}

#endif