#ifndef  MINI_BMP
#define  MINI_BMP

#include <fstream>
#include <string>
#include <iostream>

#include "windowDecl.h"

using namespace mini;

//bmp文件读取类
class BmpClass
{
public:
	BmpClass(const char* bmpName);
	~BmpClass();

	//从文件中读取数据
	bool readFile()
	{
		char temp[100];
		std::fstream infile(_bmpName, std::ios_base::in | std::ios_base::binary);
		if (!infile.is_open())
		{
			std::cout << "read file " << _bmpName << " filed!" << std::endl;
			return false;
		}
		BITMAPINFOHEADER head;
		
		infile.read(temp, sizeof(BITMAPFILEHEADER));
		infile.read((char*)&head, sizeof(BITMAPINFOHEADER));

		_bmpWidth = head.biWidth;
		_bmpHeight = head.biHeight;
		_biBitCount = head.biBitCount;

		int lineByte = (_bmpWidth * _biBitCount / 8 );
		_pBmpBuf = new unsigned char[lineByte * _bmpHeight];
		infile.read((char*)_pBmpBuf, lineByte * _bmpHeight);
		infile.close();

		return true;
	}

	//以24位真彩色写入文件				//从左下到右上读取，以行为优先
	bool writeFile(UCHAR* frameBuffer, UINT width, UINT height, UINT len)
	{
		BITMAPFILEHEADER fileHead;
		fileHead.bfType = 0x4D42;
		int a = sizeof(BITMAPFILEHEADER);
		fileHead.bfSize = sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER)+len;
		fileHead.bfReserved1 = fileHead.bfReserved2 = 0;
		fileHead.bfOffBits = sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER);
		BITMAPINFOHEADER infoHead;
		infoHead.biSize = sizeof(BITMAPINFOHEADER);
		infoHead.biWidth = width;
		infoHead.biHeight = height;
		infoHead.biPlanes = 1;
		infoHead.biBitCount = 24;
		infoHead.biCompression = BI_RGB;
		infoHead.biSizeImage = 0;
		infoHead.biXPelsPerMeter = 0;
		infoHead.biYPelsPerMeter = 0;
		infoHead.biClrUsed = 0;
		infoHead.biClrImportant = 0;

		std::fstream outfile(_bmpName, std::ios_base::out | std::fstream::binary);
		if (!outfile.is_open())
		{
			std::cout << "\nWritefile " << _bmpName << " filed!" << std::endl;
			return false;
		}
		outfile.write((char*)&fileHead, sizeof(BITMAPFILEHEADER));
		outfile.write((char*)&infoHead, sizeof(BITMAPINFOHEADER));
		outfile.write((char*)frameBuffer, len);
		outfile.close();

		return true;
	}

	//根据坐标获取图片颜色,左下角(0,0)-(1,1)
	UINT GetRGB(float x, float y)
	{
		//Clamp to 0-1
		x = x - (int)x;
		if (x < 0) x++;
		y = y - (int)y;
		if (y < 0)y++;

		UINT co;
		int index = (int(x*_bmpWidth + 0.5f) + int(y * _bmpHeight + 0.5)* _bmpWidth);

		if (_biBitCount == 24)
		{
			co = _pBmpBuf[index * 3] << 16 | _pBmpBuf[index * 3+1] << 8 | _pBmpBuf[index * 3+2];
		}
		if (_biBitCount == 32)
		{
			UINT* tmp = (UINT*)_pBmpBuf;
			co = tmp[index] >> 8;
		}
		return co;
	}
private:
	unsigned char *_pBmpBuf;
	std::string _bmpName;
	int _bmpWidth;
	int _bmpHeight;
	int _biBitCount;//图像类型，每像素位数  

};

BmpClass::BmpClass(const char* bmpName)
{
	_bmpName = std::string(bmpName);
	_pBmpBuf = NULL;
	_bmpWidth = 0;
	_bmpHeight = 0;
	_biBitCount = 0;
}

BmpClass::~BmpClass()
{
	if(_pBmpBuf != NULL)delete[] _pBmpBuf;
}

#endif