#include "grid.h"
#include "parameter.h"
#include "protein.h"

using namespace std;

Grid::Grid()
{
	type = 'V';
	type_bak = 'V';
	type_final = 'O';
	neib = 0;
	grid_surface = 0;
	cavity = 0;
	cavity_final = 0;
	inner = -1;
	inner_bak = -1;
	inner2 = -1;

	re = -1;
	re_bak = -1;
	predict = -1;
	predict2 = -1;

	dep = 0;
	rdep = -1;
	donor_score = acceptor_score = hydrophobic_score = 0.0;
	donor_score_bak = acceptor_score_bak = hydrophobic_score_bak = 0.0;

	donor_overlap = 10;
	acceptor_overlap = 10;
	donor_wtype = 0;
	acceptor_wtype = 0;
	logp = 0;
	adh = 'N';
	adh_bak = 'N';

	mark = 0;
	
	lpoint = 0;
}

FGrid::FGrid()
{
	type = 'V';
	type_bak = 'V';
	neib = 0;
}

GridBox::GridBox()
{
	extend_length = 0;
	min_x = min_y = min_z = max_x = max_y = max_z = 0;
	dx = dy = dz = 0;
	num_grid = 0;
	space = 1;
	grid = nullptr;
}

void GridBox::Initialize(int extend, float space, int min_xx, int max_xx, int min_yy, int max_yy, int min_zz, int max_zz)
{
	extend_length = extend;
	min_x = min_xx - extend;
	max_x = max_xx + extend;
	min_y = min_yy - extend;
	max_y = max_yy + extend;
	min_z = min_zz - extend;
	max_z = max_zz + extend;
	
	dx = (max_x - min_x) * 2 + 1;
	dy = (max_y - min_y) * 2 + 1;
	dz = (max_z - min_z) * 2 + 1;
	
	grid = new Grid * *[dx];
	for (int i = 0; i < dx; i++)
	{
		grid[i] = new Grid * [dy];
		for (int j = 0; j < dy; j++)
		{
			grid[i][j] = new Grid[dz];
			for (int k = 0; k < dz; k++)
			{
				grid[i][j][k].coord.x = min_x + i * space;
				grid[i][j][k].coord.y = min_y + j * space;
				grid[i][j][k].coord.z = min_z + k * space;
				grid[i][j][k].type = 'V';
			}
		}
	}
	num_grid = dx * dy * dz;
}

int GridBox::OutOfBoundary(int id, int jd, int kd)
{
	if (id >= 0 && id < dx && jd >= 0 && jd < dy && kd >= 0 && kd < dz)
	{
		return 0;
	}
	else return 1;
}

Grid& GridBox::operator()(Point3D& p)
{
	int i, j, k;
	i = (int)p.x;
	j = (int)p.y;
	k = (int)p.z;
	return grid[i][j][k];
}

Grid& GridBox::operator()(Position3D& p)
{
	return grid[p.i][p.j][p.k];
}

Grid& GridBox::operator()(const Position3D& p)
{
	return grid[p.i][p.j][p.k];
}

void GridBox::Clear()
{
	for (int i = 0; i < dx; i++)
	{
		for (int j = 0; j < dy; j++)
		{
			delete[] grid[i][j];
		}
		delete[] grid[i];
	}
	delete[] grid;
}


FGridBox::FGridBox()
{
	extend_length = 0;
	min_x = min_y = min_z = max_x = max_y = max_z = 0;
	dx = dy = dz = 0;
	num_grid = 0;
	space = 1;
	grid = nullptr;
}

void FGridBox::Initialize(int extend, float space, int min_xx, int max_xx, int min_yy, int max_yy, int min_zz, int max_zz)
{
	extend_length = extend;
	min_x = min_xx - extend;
	max_x = max_xx + extend;
	min_y = min_yy - extend;
	max_y = max_yy + extend;
	min_z = min_zz - extend;
	max_z = max_zz + extend;

	dx = (max_x - min_x) * 2 + 1;
	dy = (max_y - min_y) * 2 + 1;
	dz = (max_z - min_z) * 2 + 1;

	grid = new FGrid * *[dx];
	for (int i = 0; i < dx; i++)
	{
		grid[i] = new FGrid * [dy];
		for (int j = 0; j < dy; j++)
		{
			grid[i][j] = new FGrid[dz];
			for (int k = 0; k < dz; k++)
			{
				grid[i][j][k].coord.x = min_x + i * space;
				grid[i][j][k].coord.y = min_y + j * space;
				grid[i][j][k].coord.z = min_z + k * space;
				grid[i][j][k].type = 'V';
			}
		}
	}
	num_grid = dx * dy * dz;
}

int FGridBox::OutOfBoundary(int id, int jd, int kd)
{
	if (id >= 0 && id < dx && jd >= 0 && jd < dy && kd >= 0 && kd < dz)
	{
		return 0;
	}
	else return 1;
}

FGrid& FGridBox::operator()(Position3D& p)
{
	return grid[p.i][p.j][p.k];
}


void FGridBox::Clear()
{
	for (int i = 0; i < dx; i++)
	{
		for (int j = 0; j < dy; j++)
		{
			delete[] grid[i][j];
		}
		delete[] grid[i];
	}
	delete[] grid;
}
