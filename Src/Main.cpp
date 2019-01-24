#include "../Head/Commen.h"
#include <thread>
#include <chrono> 
#include <exception>

using namespace mini;

/*
使用方法：
	前提条件：将程序miniRT.exe与场景文件放置同一文件夹；
	在当前文件夹打开CMD或Power Shell，打开方法：在文件夹空白处“shift+鼠标右键”，然后选择“在此处打开CMD或Power Shell”；
	在CMD或Power Shell中，输入“./miniRT.exe [Scene文件，如：Scene1.rt]”
	然后根据提示输入单像素采样数，采样深度。
*/

//真个程序在右手坐标系下构建

//主函数
int main(int argc, char* argv[])
{
	int spp = 1, depth = 1;
#ifndef _DEBUG
	if (argc != 2)
	{
		std::cout << "Can't read scene file, please input scene file as program parameter!\n"
			<< "For Example :\n" << "[in cmd.exe] mimiRT Scene.rt" << std::endl;
		std::this_thread::sleep_for(std::chrono::seconds(3));
		return 0;
	}
	std::cout << "Render parmater:\n"
		<< "Please input SamplePerPixel( Squre number is better):";
	std::cin >> spp;
	std::cout << "Please input TraceDepth:";
	std::cin >> depth;

#endif
	try
	{
		Scene* scene = new Scene(600, 600, spp, depth, color(0, 0, 0));
#ifndef _DEBUG
		scene->initFromFile(argv[1]);
#else
		scene->initFromFile("MyScene.rt");
#endif // !_DEBUG
		scene->renderAndWrite("RayTracing.bmp");
		delete scene;
	}
	catch( std::exception e){
		std::cout << e.what() << std::endl;
	}
	
	system("pause");
	return 0;
}