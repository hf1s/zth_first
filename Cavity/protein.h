#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include "atom.h"

using std::string;
using std::vector;
using std::map;
using std::set;

vector<string> get_residue_def();

class Residue
{
public:
	string name;
	map<string, Patom> atoms;
	int num_atom;
public:
	Residue();
};

class Protein
{
public:
	int num_resi_def;
	map<string, Residue> resi_defs;
	
	int num_atom;
	vector<Patom> atoms;
	map<char, vector<int>> chain_record;

	Point3D center_coord;
	int min_x, min_y, min_z;
	int max_x, max_y, max_z;

public:
	Protein();
	void Read_RESIDUE(const char* dirname);
	void Read_RESIDUE();
	void Show_RESIDUE();
	void Read_PDB(const char* filename);
	void Protein_Box();

	void Value_Atom();
	void Cal_HB_Root();

	void cal_geometric_center();

};

