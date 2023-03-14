#pragma once
#include <time.h>
#include <iostream>
#include <string>
#include <sstream>
#include <set>

using std::string;

class Clock
{
public:
	Clock() : _begin(0), _end(0), _start_time(nullptr), _stop_time(nullptr) {};
	inline void begin() { _begin = clock();  time_t now = time(NULL); _start_time = localtime(&now); }
	inline void end() { _end = clock(); time_t now = time(NULL); _stop_time = localtime(&now); }
	inline double duration() { return (double)(_end - _begin) / CLOCKS_PER_SEC; }
	inline string strstart() { return string(asctime(_start_time)); }
	inline string strstop() { return string(asctime(_stop_time)); }
	inline string strcurrent() { time_t now = time(NULL); struct tm* curr = localtime(&now); return string(asctime(curr)); }
private:
	clock_t _begin; // begin clock time, for calculating duration
	clock_t _end; // end clock time, for calculating duration
	struct tm* _start_time; // begin date and time, for display
	struct tm* _stop_time; // end date and time for display
};

char* Get_Time();

void Memory_Allocation_Error();

void Openning_File_Error(const char* filename);

void PDB_Format_Error(const char* filename);

void Mol2_Format_Error(const char* filename);

void Missing_Parameter_Error(const char* name);

void Check_Directory(string& name, char sep='/');

template <typename T>
void string_to(string str, T& value)
{
	if (str.size() == 0) return;
	std::stringstream s;
	s << str;
	s >> value;
}

void Create_Ball_Grid(int radius, float cutoff, int***& fullball, int& border_min_x, int& border_max_x, int*& border_min_y, int*& border_max_y, int**& border_min_z, int**& border_max_z);

void Clear_Ball_Grid(int radius, int***& fullball, int& border_min_x, int& border_max_x, int*& border_min_y, int*& border_max_y, int**& border_min_z, int**& border_max_z);

void Test();

