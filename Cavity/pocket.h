#pragma once
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <set>
#include "grid.h"
#include "utils.h"

# define HB_D_R 1.86
# define HB_A_R 1.67 
# define HYDROPHOBIC_R 1.94
# define OLD_HB_ANGLE_CUTOFF 60.0
# define VB  -0.168
# define SHB  0.593
# define MHB  0.216
# define WHB  0.141
# define SWH  0.291
# define MWH -0.708
# define WWH  0.327
# define MB   0.916
# define HM   0.473
# define SB_AWARD 0.150

using namespace std;

class Cavity
{
public:
	int valid;
	int id;
	int re;
	int surfacenum;
	int vacant_num;
	int index;
	float rank;
	//float vacant1, vacant2, vacant3, vacant4, vacant5;
	float surface;
	float edge;
	int area;
	float dep1, dep2;
	float min_x, min_y, min_z;
	float max_x, max_y, max_z;
	float atom_min_x, atom_min_y, atom_min_z;
	float atom_max_x, atom_max_y, atom_max_z;

	Position3D area_core;
	Position3D surface_core;
	Point3D ab_vector;
	Point3D core_vector;
	float ab_vector_norm;
	float core_vector_norm;

	float hydrophobic_surface, acceptor_surface, donor_surface, bind_surface;
	float hydrophobic_s, acceptor_s, donor_s;
	int hydrophobic_vacant, acceptor_vacant, donor_vacant, bind_vacant;
	int hydrophobic_v, acceptor_v, donor_v;

	int surface_dep, vacant_dep;
	int step_vacant[1000];

	float drugscore;

	// for ligand mode
	int percent; // 配体原子所占据的口袋的'V'/'S'格点的数量
	int ligmatch; // 配体所占据的口袋'V'格点的数量
	int ligall; // 配体所占据的全部格点的数量

	vector<string> slines; // surface file
	vector<string> vlines; // vacant file
	vector<string> clines; // cavity file
	string s_str;
	string v_str;
	string c_str;
public:
	Cavity();
};

bool comp(Cavity& c1, Cavity& c2);

class bGrid
{
public:
	int i, j, k;
	int cavity;
public:
	bGrid(int ii, int jj, int kk, int c);
};

class Pocket
{
public:
	GridBox grid;
	FGridBox fgrid;
	int shift_x;
	int shift_y;
	int shift_z;
	vector<bGrid> bgrid;

	int cavity_valid;
	int totalcavity;
	int new_total;

	// define the ball
	int fullradius;
	int length;
	int*** fullball;
	int border_min_x, border_max_x;
	int* border_min_y;
	int* border_max_y;
	int** border_min_z;
	int** border_max_z;
	map<int, vector<Position3D>> shiftgrid;

	map<int, vector<Position3D>> vacantset;
	vector<Position3D> surfacegrid;
	map<int, vector<Position3D>> surfaceset;
	vector<Position3D> excludedset;
	vector<Position3D> lpointgrid; // for ligand mode
	Cavity* cavity;
	
public:
	Pocket();
	void Cavity_Search();
	void Initialize();
	void Recover(int maxdep);
	void Erase();
	void Erase_Cumulative(Ppair& p, vector<Ppair>& list);
	void Define_Surface();
	void Define_Vacant();
	void Reedge();
	void Get_Real_Depth(int n);
	int Check(int n, int seed);
	int Detect_V_Cluster(int n, int seed, int* vnum, vector<int>& tmpvnum);
	int Get_Depth_Predict(int n);
	int Get_Predict_V(int n, int shrink);
	void Backup_Adjust(int n);
	void Adjust(int n, int edge, int shrink, int* re);
	void Restore_Adjust(int n);
	void Core_Adjust(int n, Position3D& core_s, Position3D& core_v);
	void Backup(int ruler);
	void Restore(int n, int ruler);
	int Choose(int seed, int ruler);
	void Cavity_Info(int seed, int ruler);
	int Get_Surface_Dep(int n);
	int Get_Vacant_Dep(int n);
	void Prepare_Output(int n, int k);
	void Output_All(int n, int k);
	void Prepare_Out_Surface(int seed, int k);
	void Prepare_Out_Vacant(int seed, int k);
	void Prepare_Out_Cavity(int seed, int k);
	void Output_Vacant(string fname, vector<Position3D>& vset);
	void Output_Middle_Vacant(int start, int total);
	
	// for ligand mode
	void Define_Ligand_Grid();
	void Calculate_Ligand_Percent(int seed);
};