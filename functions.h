#ifndef __FUNCTIONS_H 
#define __FUNCTIONS_H

#include "windows.h"
#include <fstream>

namespace model
{
	/*****************************************************
	*********	 ������� ��� ������ � �������	 *********
	******************************************************/
	bool isDirectoryExists(LPCWSTR);						//�������� �� ������������� ����������
	void WriteDebugInfo(std::ofstream&, double [3][3]);		//������ ������ � ����
	void TruncPoleFiles();									//������� ������ �������� �����
	void TruncSSTFiles();									//������� ������ ����������� �������������

}

#endif __FUNCTIONS_H