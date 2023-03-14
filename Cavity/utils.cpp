#include <cmath>
#include <iomanip>
#include "utils.h"

using namespace std;

char* Get_Time()
{
	struct tm* ptr;
	time_t lt;

	lt = time(NULL);
	ptr = localtime(&lt);

	return asctime(ptr);
}

void Memory_Allocation_Error()
{
	printf("\n");
	printf("Memory allocation error!\n");
	printf("Maybe there is not enough memory to run the program.\n");
	printf("So vote for another person as the President.\n");
	printf("Hope he would increase the budget for education & science.\n");
	exit(1);
}

void Openning_File_Error(const char* filename)
{
	cout << "\nError: cannot open the file " << filename << endl;
	cout << "Please make sure it exists." << endl;
	exit(1);
}

void PDB_Format_Error(const char* filename)
{
	printf("\n");
	printf("Error: %s lacks necessary information.\n", filename);
	printf("1.This file may not be in PDB format.\n");
	printf("2.Maybe you used a strict filter ruler in the index file.\n");
	printf("Please check it and try again.\n");
	exit(1);
}

void Mol2_Format_Error(const char* filename)
{
	printf("\n");
	printf("Error: %s lacks necessary information.\n", filename);
	printf("This file may not be in Mol2 format.\n");
	printf("Please check it and try again.\n");
	exit(1);
}

void Missing_Parameter_Error(const char* name)
{
	cout << "\nError: you have probably forgot telling me " << name << "." << endl;
	return;
}

void Check_Directory(string& name, char sep)
{
	if (name[name.length() - 1] != sep) name = name + sep;
	return;
}

void Create_Ball_Grid(int radius, float cutoff, int***& fullball, int& border_min_x, int& border_max_x, int*& border_min_y, int*& border_max_y, int**& border_min_z, int**& border_max_z)
{
	int length = radius * 2 + 1;
	
	fullball = new int** [length];
	for (int i = 0; i < length; i++)
	{
		fullball[i] = new int* [length];
		for (int j = 0; j < length; j++)
		{
			fullball[i][j] = new int[length];
			for (int k = 0; k < length; k++) fullball[i][j][k] = 0;
		}
	}
	// store the border of the ball
	border_min_y = new int[length];
	border_max_y = new int[length];
	border_min_z = new int* [length];
	border_max_z = new int* [length];
	for (int i = 0; i < length; i++)
	{
		border_min_y[i] = -1;
		border_max_y[i] = -1;
		border_min_z[i] = new int[length];
		border_max_z[i] = new int[length];
		for (int j = 0; j < length; j++)
		{
			border_min_z[i][j] = -1;
			border_max_z[i][j] = -1;
		}
	}
	int id, jd, kd;
	double r;
	for (id = 0; id < length; id++)
	{
		for (jd = 0; jd < length; jd++)
		{
			for (kd = 0; kd < length; kd++)
			{
				r = sqrt((id - radius) * (id - radius) + (jd - radius) * (jd - radius) + (kd - radius) * (kd - radius));
				if (r <= cutoff)
				{
					fullball[id][jd][kd] = 1;
					if (border_min_z[id][jd] == -1)
					{
						border_min_z[id][jd] = kd;
						border_max_z[id][jd] = kd + 1;
					}
					else
					{
						border_max_z[id][jd] = kd + 1;
					}
				}
			}
			if (border_min_z[id][jd] != -1)
			{
				if (border_min_y[id] == -1)
				{
					border_min_y[id] = jd;
					border_max_y[id] = jd + 1;
				}
				else
				{
					border_max_y[id] = jd + 1;
				}
			}
		}
		if (border_min_y[id] != -1)
		{
			if (border_min_x == -1)
			{
				border_min_x = id;
				border_max_x = id + 1;
			}
			else
			{
				border_max_x = id + 1;
			}
		}
	}
}

void Clear_Ball_Grid(int radius, int***& fullball, int& border_min_x, int& border_max_x, int*& border_min_y, int*& border_max_y, int**& border_min_z, int**& border_max_z)
{
	int length = radius * 2 + 1;
	for (int i = 0; i < length; i++)
	{
		for (int j = 0; j < length; j++)
		{
			delete[] fullball[i][j];
		}
		delete[] fullball[i];
	}
	delete fullball;
	delete border_max_y, border_min_y;
	for (int i = 0; i < length; i++)
	{
		delete[] border_min_z[i];
		delete[] border_max_z[i];
	}
	delete border_min_z, border_max_z;
}

void Test()
{
	return;
}