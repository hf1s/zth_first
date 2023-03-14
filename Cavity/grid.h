#pragma once
//#include <set>
#include "atom.h"

class Grid
{
public:
	Point3D coord;
	
	char type;
	char type_bak;
	char type_final;
	int neib;
	
	float grid_surface;
	int cavity;
	int cavity_final;
	int inner;
	int inner_bak;
	int inner2;

	int re;
	int re_bak;
	int predict;
	int predict2;

	int dep;
	float rdep;
	float donor_score, acceptor_score, hydrophobic_score;
	float donor_score_bak, acceptor_score_bak, hydrophobic_score_bak;
	float acceptor_overlap, donor_overlap;
	int acceptor_wtype, donor_wtype;
	float logp;
	char adh;
	char adh_bak;

	int mark;

	int lpoint; //=1: ligand atom occupies this grid
	
public:
	Grid();
};

class FGrid
{
public:
	Point3D coord;
	char type;
	char type_bak;
	int neib;
	FGrid();
};

class GridBox
{
public:
	int extend_length;
	int min_x, min_y, min_z;
	int max_x, max_y, max_z;
	int dx, dy, dz;
	int num_grid;
	float space;
	Grid*** grid;
public:
	GridBox();
	void Initialize(int extend, float space, int min_xx, int max_xx, int min_yy, int max_yy, int min_zz, int max_zz);
	int OutOfBoundary(int id, int kd, int jd);

	Grid& operator()(Point3D& p);
	Grid& operator()(Position3D& p);
	Grid& operator()(const Position3D& p);

	void Clear();
};

class FGridBox
{
public:
	int extend_length;
	int min_x, min_y, min_z;
	int max_x, max_y, max_z;
	int dx, dy, dz;
	int num_grid;
	float space;
	FGrid*** grid;

	FGridBox();
	void Initialize(int extend, float space, int min_xx, int max_xx, int min_yy, int max_yy, int min_zz, int max_zz);
	int OutOfBoundary(int id, int kd, int jd);
	FGrid& operator()(Position3D& p);
	void Clear();
};