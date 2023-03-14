#include <fstream>
#include <vector>
#include "ligand.h"
#include "utils.h"
#include "atom.h"

using namespace std;

vector<string> get_atom_def()
{
	vector<string> atom_def;
	atom_def.push_back(string("# Data cited from the Tripos force field as implemented in SYBYL 6.3"));
	atom_def.push_back(string("# The 1st column: numeric ID"));
	atom_def.push_back(string("# The 2nd column: atom type symbol"));
	atom_def.push_back(string("# The 3rd column: van der Waals radius"));
	atom_def.push_back(string("# The 4th column: van der Walls epsilon value"));
	atom_def.push_back(string("# The 5th column: atomic weight"));
	atom_def.push_back(string("# The 6th column: description"));
	atom_def.push_back(string("#"));
	atom_def.push_back(string("1	C.3	1.70	0.107	12	0.000   carbon sp3"));
	atom_def.push_back(string("2	C.2	1.70	0.107	12	0.000   carbon sp2"));
	atom_def.push_back(string("3	C.1	1.70	0.107	12	0.000   carbon sp"));
	atom_def.push_back(string("4	C.cat	1.70	0.107	12	1.000   carbocation (C+) used in guadinium"));
	atom_def.push_back(string("5	C.ar	1.70	0.107	12	0.000   carbon aromatic"));
	atom_def.push_back(string("6	H	1.20	0.042	1	0.000   hydrogen"));
	atom_def.push_back(string("7	N.4	1.55	0.095	14	1.000   nitrogen sp3 positively charged"));
	atom_def.push_back(string("8	N.3	1.55	0.095	14	0.000   nitrogen sp3"));
	atom_def.push_back(string("9	N.2	1.55	0.095	14	0.000   nitrogen sp2"));
	atom_def.push_back(string("10	N.1	1.55	0.095	14	0.000   nitrogen sp"));
	atom_def.push_back(string("11	N.ar	1.55	0.095	14	0.000   nitrogen aromatic"));
	atom_def.push_back(string("12	N.pl3	1.55	0.095	14	0.000   nitrogen trigonal planar"));
	atom_def.push_back(string("13	N.am	1.55	0.095	14	0.000   nitrogen amide"));
	atom_def.push_back(string("14	O.3	1.52	0.116	16	0.000   oxygen sp3	"));
	atom_def.push_back(string("15	O.2	1.52	0.116	16	0.000   oxygen sp2"));
	atom_def.push_back(string("16	O.co2	1.52	0.116	16     -0.500   oxygen in carboxylate and phosphate"));
	atom_def.push_back(string("17	P.3	1.80	0.314	31	0.000   phosphorous sp3"));
	atom_def.push_back(string("18	S.3	1.80	0.314	32	0.000   sulfur sp3"));
	atom_def.push_back(string("19	S.2	1.80	0.314	32	0.000   sulfur sp2"));
	atom_def.push_back(string("20	S.o	1.70	0.314	32	0.000   sulfoxide sulfur"));
	atom_def.push_back(string("21	S.o2	1.70	0.314	32	0.000   sulfone sulfur"));
	atom_def.push_back(string("22	F	1.47	0.109	19	0.000   fluorine"));
	atom_def.push_back(string("23	Cl	1.75	0.314	35	0.000   chlorine"));
	atom_def.push_back(string("24	Br	1.85	0.434	80	0.000   bromine"));
	atom_def.push_back(string("25	I	1.98	0.623	127	0.000   iodine"));
	atom_def.push_back(string("26	Du	0	0.000	0	0.000   dummy atom"));
	atom_def.push_back(string("27	Un	0	0.000	0	0.000   unknown atom"));
	atom_def.push_back(string("28	H.spc	1.20	0.042	1	0.000   seed hydrogen"));
	atom_def.push_back(string("#"));
	atom_def.push_back(string("# Notice: the VDW radius of hydrogen has been changed from original 1.2 to "));
	atom_def.push_back(string("# current 1.0 to reproduce reasonable conformations"));
	atom_def.push_back(string("# now 1.2 applied"));
	atom_def.push_back(string("#"));
	return atom_def;
}

Ligand::Ligand()
{
	num_atom = 0;
	num_bond = 0;
	num_heavy = 0;
	min_x = min_y = min_z = max_x = max_y = max_z = 0;
	volume = 0;
}

void Ligand::Read_ATOM_DEF(const char* dirname)
{
	string filename(dirname);
	filename = filename + "ATOM_DEF";
	fstream infile;
	infile.open(filename, ios::in);
	if (!infile.is_open())
	{
		Openning_File_Error(filename.c_str());
	}
	string line, head, value;
	stringstream sline;
	while (getline(infile, line))
	{
		if (line.size() == 0) { continue; }// blank line
		else if (line[0] == '#') { continue; } // comment line
		else
		{
			sline << line;
			ATOM_DEF atmdef;
			sline >> head >> atmdef.type >> atmdef.r >> atmdef.eps >> atmdef.weight;
			sline.str(std::string());
			sline.clear();
			atomdef.push_back(atmdef);
		}
	}
	infile.close();
	return;
}

void Ligand::Read_ATOM_DEF()
{
	vector<string> atom_def_lines = get_atom_def();
	string head, value;
	stringstream sline;
	for(auto & line : atom_def_lines)
	{
		if (line.size() == 0) { continue; }// blank line
		else if (line[0] == '#') { continue; } // comment line
		else
		{
			sline << line;
			ATOM_DEF atmdef;
			sline >> head >> atmdef.type >> atmdef.r >> atmdef.eps >> atmdef.weight;
			sline.str(std::string());
			sline.clear();
			atomdef.push_back(atmdef);
		}
	}
}

void Ligand::Read_From_Mol2(const char* filename)
{
	extern Parameter* parm;
	ifstream infile;
	infile.open(filename);
	if (!infile.is_open()) { Openning_File_Error(filename); }

	string line;
	stringstream sline;
	int success = 0;
	while (getline(infile, line))
	{
		if (line == "@<TRIPOS>MOLECULE") 
		{
			success = 1; break;
		}
	}
	if (success == 0) Mol2_Format_Error(filename);
	getline(infile, line);
	getline(infile, line);
	sline << line;
	sline >> num_atom >> num_bond;
	sline.clear(); sline.str("");
	success = 0;
	while (getline(infile, line))
	{
		if (line == "@<TRIPOS>ATOM") { success = 1; break; }
	}
	if (success == 0) Mol2_Format_Error(filename);

	for (int i = 0; i < num_atom; i++)
	{
		getline(infile, line);
		Latom atm;
		atm.Read_MOL2_line(line);
		atoms.push_back(atm);
	}
	success = 0;
	while (getline(infile, line))
	{
		if (line == "@<TRIPOS>BOND") { success = 1; break; }
	}
	if (success == 0) Mol2_Format_Error(filename);
	for (int i = 0; i < num_bond; i++)
	{
		getline(infile, line);
		Bond bon;
		bon.Read_MOL2_line(line);
		bonds.push_back(bon);
	}
	infile.close();
	if (parm->single_ligand == 1)
	{
		vector<int> component;
		int cid = 0;
		int maxnum = 0, maxid = 1;
		int curr_id = 0;
		for (int i=0; i< num_atom;i++)
		{
			if (atoms[i].belong != 0) continue;
			cid++;
			atoms[i].belong = cid;
			component.push_back(i);
			for (int j = 0; j < component.size(); j++)
			{
				curr_id = atoms[j].atom_id;
				for (int k = 0; k < num_bond; k++)
				{
					if (bonds[k].belong != 0) continue;
					if (curr_id == bonds[k].atom_1 && atoms[bonds[k].atom_2 - 1].belong == 0)
					{
						atoms[bonds[k].atom_2 - 1].belong = cid;
						component.push_back(bonds[k].atom_2 - 1);
					}
					else if (curr_id == bonds[k].atom_2 && atoms[bonds[k].atom_1 - 1].belong == 0)
					{
						atoms[bonds[k].atom_1 - 1].belong = cid;
						component.push_back(bonds[k].atom_1 - 1);
					}
				}
			}
			if (cid > 1) cout << "The ligand has more than one components, only the biggest part ..." << endl;
			if (component.size() > maxnum) 
			{
				maxid = cid;
				maxnum = component.size();
			}
			for (int k = 0; k < num_atom; k++)
			{
				if (atoms[k].belong != maxid) atoms[k].valid = 0;
			}
		}
		volume = Calc_Volume();
	}
	for (int i = 0; i < num_atom; i++)
	{
		if (atoms[i].valid == 1 && atoms[i].tripos_type != "H") num_heavy++;
	}
	return;
}

float Ligand::Calc_Volume()
{
	int mark = 0;
	int lmin_x, lmin_y, lmin_z, lmax_x, lmax_y, lmax_z;
	lmin_x = lmin_y = lmin_z = lmax_x = lmax_y = lmax_z = 0;
	for (auto& atm : atoms)
	{
		if (atm.valid == 0) continue;
		for (auto& adef : atomdef)
		{
			if (atm.tripos_type == adef.type)
			{
				atm.vdw_r = adef.r;
				break;
			}
		}
		if (atm.vdw_r == 0) { cout << "Unknown Ligand Atom Type " << atm.tripos_type << ", Skiped..." << endl; }
		if (mark == 0)
		{
			lmin_x = (int)(atm.coord.x); lmax_x = (int)(atm.coord.x);
			lmin_y = (int)(atm.coord.y); lmax_y = (int)(atm.coord.y);
			lmin_z = (int)(atm.coord.z); lmax_z = (int)(atm.coord.z);
			mark = 1;
		}
		else
		{
			lmin_x = atm.coord.x < lmin_x ? (int)(atm.coord.x) : lmin_x;
			lmax_x = atm.coord.x > lmax_x ? (int)(atm.coord.x) : lmax_x;
			lmin_y = atm.coord.y < lmin_y ? (int)(atm.coord.y) : lmin_y;
			lmax_y = atm.coord.y > lmax_y ? (int)(atm.coord.y) : lmax_y;
			lmin_z = atm.coord.z < lmin_z ? (int)(atm.coord.z) : lmin_z;
			lmax_z = atm.coord.z > lmax_z ? (int)(atm.coord.z) : lmax_z;
		}
	}
	lmax_x += 3; lmax_y += 3; lmax_z += 3;
	lmin_x -= 3; lmin_y -= 3; lmin_z -= 3;
	int ldx, ldy, ldz;
	ldx = (lmax_x - lmin_x) * 10 + 1;
	ldy = (lmax_y - lmin_y) * 10 + 1;
	ldz = (lmax_z - lmin_z) * 10 + 1;
	int*** lgrid = new int * *[ldx];
	for (int i = 0; i < ldx; i++)
	{
		lgrid[i] = new int * [ldy];
		for (int j = 0; j < ldy; j++)
		{
			lgrid[i][j] = new int[ldz];
			for (int k = 0; k < ldz; k++)
			{
				lgrid[i][j][k] = 0;
			}
		}
	}
	float it, jt, kt;
	int id, jd, kd;
	Point3D v1, v2;
	for (auto& atm : atoms)
	{
		it = ((float)((int)((atm.coord.x) * 10))) / 10;
		jt = ((float)((int)((atm.coord.y) * 10))) / 10;
		kt = ((float)((int)((atm.coord.z) * 10))) / 10;

		id = (int)((it - lmin_x) * 10);
		jd = (int)((jt - lmin_y) * 10);
		kd = (int)((kt - lmin_z) * 10);
		v1.set_xyz((float)(id) / 10, (float)(jd) / 10, (float)(kd) / 10);
		for (int i = id - 20; i < id + 21; i++)
			for (int j = jd - 20; j < jd + 21; j++)
				for (int k = kd - 20; k < kd + 21; k++)
				{
					if (i >= 0 && i < ldx && j >= 0 && j < ldy && k >= 0 && k < ldz)
					{
						v2.set_xyz((float)(i) / 10, (float)(j) / 10, (float)(k) / 10);
						if (Distance(v1, v2) <= atm.vdw_r)lgrid[i][j][k] = 1;
					}
				}
	}
	float total = 0;
	for (int i = 0; i < ldx; i++)
		for (int j = 0; j < ldy; j++)
			for (int k = 0; k < ldz; k++)
			{
				if (lgrid[i][j][k] == 1)total++;
			}
	total = total / 1000;
	for (int i = 0; i < ldx; i++)
	{
		for (int j = 0; j < ldy; j++)
		{
			delete[] lgrid[i][j];
		}
		delete[] lgrid[i];
	}
	delete[] lgrid;
	return total;
}

void Ligand::Ligand_Box()
{
	extern Parameter* parm;
	int mark = 0;
	for (auto& it : atoms)
	{
		if (it.valid == 0) { continue; }
		if (mark == 0)
		{
			min_x = (int)(it.coord.x); max_x = (int)(it.coord.x);
			min_y = (int)(it.coord.y); max_y = (int)(it.coord.y);
			min_z = (int)(it.coord.z); max_z = (int)(it.coord.z);
			mark = 1;
		}
		else
		{
			min_x = it.coord.x < min_x ? (int)(it.coord.x) : min_x;
			max_x = it.coord.x > max_x ? (int)(it.coord.x) : max_x;
			min_y = it.coord.y < min_y ? (int)(it.coord.y) : min_y;
			max_y = it.coord.y > max_y ? (int)(it.coord.y) : max_y;
			min_z = it.coord.z < min_z ? (int)(it.coord.z) : min_z;
			max_z = it.coord.z > max_z ? (int)(it.coord.z) : max_z;
		}
	}
	if (parm->info == 1)
	{
		printf("**INFO: The area of ligand is:\n");
		printf("Max_x = %d\tMin_x = %d\n", max_x, min_x);
		printf("Max_y = %d\tMin_y = %d\n", max_y, min_y);
		printf("Max_z = %d\tMin_z = %d\n", max_z, min_z);
	}
	// add 12.0 A margin to the box
	max_x += 12; max_y += 12; max_z += 12;
	min_x -= 12; min_y -= 12; min_z -= 12;
	parm->min_x = min_x;
	parm->max_x = max_x;
	parm->min_y = min_y;
	parm->max_y = max_y;
	parm->min_z = min_z;
	parm->max_z = max_z;
	return;
}
