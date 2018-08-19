#ifndef MINI_WINDOWDELC
#define MINI_WINDOWDELC

//移值linux系统使用，用于bmp类的一些宏声明
#ifndef _WIN32
namespace mini
{
#define NULL 0
	typedef int INT;
	typedef unsigned int UINT;
	typedef unsigned char UCHAR;
	typedef unsigned short WORD;
	typedef unsigned int DWORD;
	typedef int LONG;

#pragma pack (2)	//修改默认的字节对齐
	typedef struct tagBITMAPFILEHEADER
	{
		WORD    bfType;
		DWORD   bfSize;
		WORD    bfReserved1;
		WORD    bfReserved2;
		DWORD   bfOffBits;
	} BITMAPFILEHEADER;
#pragma pack ()

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L
#define BI_JPEG       4L
#define BI_PNG        5L

	typedef struct tagBITMAPINFOHEADER
	{
		DWORD  biSize;
		LONG   biWidth;
		LONG   biHeight;
		WORD   biPlanes;
		WORD   biBitCount;
		DWORD  biCompression;
		DWORD  biSizeImage;
		LONG   biXPelsPerMeter;
		LONG   biYPelsPerMeter;
		DWORD  biClrUsed;
		DWORD  biClrImportant;
	} BITMAPINFOHEADER;
}
#endif

#endif