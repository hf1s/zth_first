#pragma once
#include <iostream>
#include <sstream>
#include <iomanip>

using std::string;

class Position3D 
{
public:
	int i;
	int j;
	int k;
	Position3D() { i = 0; j = 0; k = 0; }
	Position3D(int a, int b, int c) { i = a; j = b; k = c; }
	Position3D(const Position3D& p);
	void set(int a, int b, int c);
	friend bool operator < (Position3D p1, Position3D p2);
};

typedef std::pair<Position3D, Position3D> Ppair;

class Point3D
{
public:
	float x, y, z;
public:
	Point3D();
	Point3D(float a, float b, float c);
	Point3D(const Point3D& p);

	void set_xyz(float a, float b, float c); 
	void set_xyz(int a, int b, int c);
	void reset();
	double norm2();
	float norm();
	
	friend Point3D operator + (Point3D& p1, Point3D& p2);
	friend Point3D operator - (Point3D& p1, Point3D& p2);
	friend float operator * (Point3D& p1, Point3D& p2);
	Point3D& operator +=(const Point3D& p);
	template <typename T>
	Point3D& operator/=(T a)
	{
		if (a == 0)
		{
			std::cout << "Warning: can not divide 0. Nothing will be done." << std::endl;
			return *this;
		}
		this->x = this->x / a;
		this->y = this->y / a;
		this->z = this->z / a;
		return *this;
	}
	friend double Angle_Of_Two_Points(Point3D& p1, Point3D& p2);
	friend float Distance(Point3D& p1, Point3D& p2);
	friend double Distance2(Point3D& p1, Point3D& p2);
};

class Bond
{
public:
	int bond_id;
	int atom_1;
	int atom_2;
	string type;
	int valid;
	int belong;
	Bond();
	void Read_MOL2_line(string line);
};

class Atom
{
public:
	int atom_id;
	string name;
	string tripos_type; // old tripos_type
	string amber_type; // old hettype
	string element; // element type

	Point3D coord;

	// property
	float logp; // atomic hydrophobic scale
	float q; // atomic charge
	float vdw_r; // VDW radius, old r

	// hydrogen bond
	string hb_type; // hydrogen bond type, old hb
	int max_anum; // max acceptor num, old anum
	int max_dnum; // max donor num, old dnum
	int hblevel;
	float hangle;

	int valid;

public:
	Atom();
};

class Patom :public Atom
{
public:
	string head;
	int resi_id;
	string resi_name; // old residue
	char altloc; // Alternate location indicator
	char iCode; // Code for insertion of residues.
	string charge; // charge in PDB line

	float occupancy;
	float bfactor;

	char chain;

	string pdbline;

	int neib[20];
	int num_neib;
	Point3D root;
	Point3D plane;

public:
	Patom();
	void Read_RESIDUEDEF_line(string line);
	void Read_PDB_line(string line);
	string Get_PDB_line();
	bool is_longpair();
	bool is_hydrogen();
};

class Latom : public Atom
{
public:
	int belong;
public:
	Latom();
	void Read_MOL2_line(string line);
};
