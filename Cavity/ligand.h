#pragma once
#include <vector>
#include "parameter.h"
#include "atom.h"

using namespace std;

vector<string> get_atom_def();

struct ATOM_DEF
{
	string type;
	float r;
	float eps;
	int weight;
};

class Ligand
{
public:
	int num_atom;
	int num_bond;
	int num_heavy;
	vector<Latom> atoms;
	vector<Bond> bonds;
	vector<ATOM_DEF> atomdef;
	float volume;

	int min_x, min_y, min_z;
	int max_x, max_y, max_z;

public:
	Ligand();
	void Read_ATOM_DEF(const char* dirname);
	void Read_ATOM_DEF();
	void Read_From_Mol2(const char* filename);
	float Calc_Volume();
	void Ligand_Box();
};