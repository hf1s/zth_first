#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include "cmdline.h"

using std::string;

class Parameter
{
public:
	string pdb_file;
	string lig_file;
	string autoname;
	string out_dir;

	// output file
	string _v_surface_file;
	string _v_vacant_file;
	string _v_atom_file;

	// detect method
	int detect_mode; // 0: whole protein mode 1: ligand detection mode 2: area detection mode
	int judge; // 0=surface 1=vacant 2=vacant-surface

	// for erase ball 
	float radius_length; //radius or half lenth [step:1A]
	float radius_mis; // bigger to flat and thick and slow, smaller to opposite. Use bigger value when radius is bigger
	float SA_radius; // solvent accessible radius

	// reading pdb file
	int hetmetal; // 1: allow 0: forbidden 
	int hetwater; // 1: allow 0: forbidden 

	// for ligand mode
	int single_ligand; // 1:only process the largest compoents of ligand
	
	// Parameter for vacant / vacant - surface method
	// some criterion
	int chip_v; // skip cavity separate step if vacant < limit
	int chip_depth; // skip cavity fill step if depth < limit
	int max_depth_vacant; 
	int max_abstract_depth;
	int min_depth;
	int max_abstract_limit_v;
	int max_limit_v;
	int min_abstract_depth;
	int rigid_limit;
	bool soft_separate;
	int softsep_limit;
	bool rescue;
	int rescue_limit;

	int edge_adjust;
	int vacant_adjust;

	// output filter
	int ruler_1;
	float output_rank;

	int out_slop;
	int in_slop;

	int info;

	// for area mode or ligand mode
	float min_x, min_y, min_z;
	float max_x, max_y, max_z;

	// output 
	int visual_output; //0:close 1: only output surface results 2:output surface and vacant results 
	int pdb_output; // 0: no output 1:output all atoms 2:output cavity atoms
	int v_num; // larger for smaller visual pdb size, smaller for better quality(min = 1)
	int v_mod; // 0:random point 1:chips
	float distance; // cutoff for cavity atoms

	string version;
	bool debug;

	Parameter();
	void Read_Index(const char* filename);
	void Initialize(cmdline::parser& cmdargs);
	void Process();
};