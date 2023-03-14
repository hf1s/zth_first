#include <cmath>
#include "atom.h"
#include "utils.h"

using namespace std;

Position3D::Position3D(const Position3D& p)
{
	i = p.i; j = p.j; k = p.k;
}

void Position3D::set(int a, int b, int c)
{
	i = a; j = b; k = c;
}

bool operator<(Position3D p1, Position3D p2)
{
	if (p1.i < p2.i) return true;
	else if (p1.i == p2.i && p1.j < p2.j) return true;
	else if (p1.i == p2.i && p1.j == p2.j && p1.k < p2.k)return true;
	else return false;
}

Point3D::Point3D()
{
	x = y = z = 0;
}

Point3D::Point3D(float a, float b, float c)
{
	x = a; y = b; z = c;
}

Point3D::Point3D(const Point3D& p)
{
	x = p.x;
	y = p.y;
	z = p.z;
}

void Point3D::set_xyz(float a, float b, float c)
{
	x = a; y = b; z = c;
}

void Point3D::set_xyz(int a, int b, int c)
{
	x = (float)a;
	y = (float)b;
	z = (float)c;
}

void Point3D::reset()
{
	x = 0; y = 0; z = 0;
}

double Point3D::norm2()
{
	return (double) (x * x + y * y + z * z);
}

float Point3D::norm()
{
	return (float)(sqrt(norm2()));
}

Point3D operator +(Point3D& p1, Point3D& p2)
{
	return Point3D(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
}

Point3D operator -(Point3D& p1, Point3D& p2)
{
	return Point3D(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
}

float operator * (Point3D& p1, Point3D& p2)
{
	return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}
Point3D& Point3D::operator+=(const Point3D& p)
{
	this->x += p.x;
	this->y += p.y;
	this->z += p.z;
	return *this;
}

double Angle_Of_Two_Points(Point3D& p1, Point3D& p2)
{
	double acosvalue;
	float angle;
	acosvalue = (double) ((p1 * p2) / (p1.norm() * p2.norm()));
	if (acosvalue > 1)acosvalue = 1;
	else if (acosvalue < -1)acosvalue = -1;
	angle = acos(acosvalue);
	angle = angle / 3.14159265 * 180.0;
	return angle;
}

float Distance(Point3D& p1, Point3D& p2)
{
	float x, y, z;
	x = p1.x - p2.x;
	y = p1.y - p2.y;
	z = p1.z - p2.z;
	return sqrt(x * x + y * y + z * z);
}

double Distance2(Point3D& p1, Point3D& p2)
{
	float x, y, z;
	x = p1.x - p2.x;
	y = p1.y - p2.y;
	z = p1.z - p2.z;
	return double(x * x + y * y + z * z);
}

Bond::Bond()
{
	bond_id = 0;
	atom_1 = 0;
	atom_2 = 0;
	type = "";
	valid = 0;
	belong = 0;
}

void Bond::Read_MOL2_line(string line)
{
	stringstream sline;
	sline << line;
	sline >> bond_id
		>> atom_1
		>> atom_2
		>> type;
	valid = 1;
}

Atom::Atom()
{
	atom_id = 0;

	logp = 0;
	q = 0;
	vdw_r = 0;
	max_anum = 0;
	max_dnum = 0;
	hblevel = 0;
	hangle = 0;
	valid = 0;
}

Patom::Patom()
{
	resi_id = 0;
	occupancy = 0;
	bfactor = 0;
	chain = ' ';
	num_neib = 0;
	altloc = ' ';
	iCode = ' ';
	for (int i = 0; i < 20; i++) neib[i] = 0;
}

void Patom::Read_RESIDUEDEF_line(string line)
{
	stringstream sline;
	sline << line;
	string skip;
	sline >> skip
		>> name
		>> amber_type
		>> logp
		>> hb_type
		>> q
		>> vdw_r
		>> tripos_type
		>> max_anum
		>> max_dnum
		>> hblevel;
	if (hb_type != "M")
	{
		if (max_anum > 0 || max_dnum > 0) hb_type = "P";
		else hb_type = "N";
	}
}

void Patom::Read_PDB_line(string line)
{
	float xx, yy, zz;

	pdbline = line;

	head = line.substr(0, 6);
	string_to(line.substr(6, 5), atom_id);
	string_to(line.substr(12, 4), name);
	altloc = line[16];
	string_to(line.substr(17, 3), resi_name);
	chain = line[21];
	string_to(line.substr(22, 4), resi_id);
	iCode = line[26];
	string_to(line.substr(30, 8), xx);
	string_to(line.substr(38, 8), yy);
	string_to(line.substr(46, 8), zz);
	string_to(line.substr(54, 6), occupancy);
	string_to(line.substr(60, 6), bfactor);
	string_to(line.substr(76, 2), element);
	string_to(line.substr(78, 2), charge);

	coord.set_xyz(xx, yy, zz);
}

string Patom::Get_PDB_line()
{
	stringstream s;
	s << setw(6) << left << head  // 1-6
		<< setw(5) << atom_id  // 7-11
		<< " "                // 12
		<< setw(4) << name    // 13-16
		<< altloc             // 17
		<< setw(3) << resi_name // 18-20
		<< " "                  // 21
		<< chain                // 22
		<< setw(4) << resi_id   // 23-26
		<< iCode                // 27
		<< "   "                // 28-30
		<< setw(8) << setprecision(3) << coord.x  // 31-38
		<< setw(8) << setprecision(3) << coord.y  // 39-46
		<< setw(8) << setprecision(3) << coord.z  // 47-54
		<< setw(6) << setprecision(2) << occupancy // 55-60
		<< setw(6) << setprecision(2) << bfactor // 61-66
		<< "          "                           // 67-76
		<< setw(2) << element     // 77-78
		<< setw(2) << charge  // 79-80
		<< endl;
	string str;
	getline(s, str);
	return str;
}

bool Patom::is_longpair()
{
	if (name[0] == 'L' && name[1] == 'P') return true;
	else return false;
}

bool Patom::is_hydrogen()
{
	if (name[0] == 'H') return true;
	else if ((name[0] == '1' || name[0] == '2' || name[0] == '3') && name[1] == 'H')return true;
	else return false;
}

Latom::Latom()
{
	belong = 0;
}
void Latom::Read_MOL2_line(string line)
{
	float xx, yy, zz;
	stringstream sline;
	sline << line;
	sline >> atom_id
		>> name
		>> xx
		>> yy
		>> zz
		>> tripos_type;
	valid = 1;
	belong = 0;
	coord.set_xyz(xx, yy, zz);
}