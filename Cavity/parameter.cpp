# include "parameter.h"
# include "utils.h"
#ifdef __GNUC__
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#else
#include <io.h>
# include <direct.h>
#endif
#ifdef linux

#endif

using namespace std;

Parameter::Parameter()
{
	pdb_file = "";
	lig_file = "";
	out_dir = "";
	detect_mode = -1;
	judge = 2;
	radius_length = 10.0;
	radius_mis = 0.5;
	SA_radius = 1.5;

	min_x = min_y = min_z = 9999;
	max_x = max_y = max_z = -9999;
	hetmetal = 1;
	hetwater = 0;
	single_ligand = 1;
	chip_v = 300;
	chip_depth = 1;
	max_depth_vacant = 100;
	max_abstract_depth = 20;
	min_depth = 8;
	max_abstract_limit_v = 1500;
	max_limit_v = 6000;
	min_abstract_depth = 2;
	rigid_limit = 0;
	soft_separate = false;
	softsep_limit = 3000;
	rescue = false;
	rescue_limit = 1000;

	edge_adjust = 0;
	vacant_adjust = 0;

	ruler_1 = 100;
	output_rank = 1.5;

	out_slop = 3;
	in_slop = 3;

	info = 0;

	visual_output = 1;
	pdb_output = 2;

	v_num = 1;
	v_mod = 0;
	distance = 5;

	version = "No version information";
	debug = false;
}

void Parameter::Read_Index(const char* filename)
{
	fstream infile;
	infile.open(filename, ios::in);
	if (!infile.is_open()) { Openning_File_Error(filename); }

	string line, key, value;
	stringstream sline;
	while (getline(infile, line))
	{
		if (line.size() == 0) { continue; } // blank line
		else if (line[0] == '#') { continue; } // comment line
		
		sline << line;
		sline >> key;
		switch (key[0])
		{
		case 'A':
			if (key == "ATOM_RADIUS_ADJUST") { sline >> SA_radius; }
			break;
		case 'B':
			//if (key == "BOUNDARY") { sline >> boundary; }
			break;
		case 'D':
			if (key == "DETECT_MODE") { sline >> detect_mode; }
			else if (key == "DISTANCE") { sline >> distance; }
			break;
		case 'E':
			if (key == "EDGE_ADJUST") { sline >> edge_adjust; }
			break;
		case 'H':
			if (key == "HETMETAL") { sline >> hetmetal; }
			else if (key == "HETWATER") { sline >> hetwater; }
			break;
		case 'I':
			if (key == "INCLUDE") {
				sline >> value;
				Read_Index(value.c_str());
			}
			else if (key == "IN_SLOP") { sline >> in_slop; }
			break;
		case 'J':
			if (key == "JUDGE") { sline >> judge; }
			break;
		case 'L': 
			if (key == "LIGAND_FILE") { sline >> lig_file; }
			break;
		case 'M':
			if (key == "MIN_X") { sline >> min_x; }
			else if (key == "MIN_Y") { sline >> min_y; }
			else if (key == "MIN_Z") { sline >> min_z; }
			else if (key == "MAX_X") { sline >> max_x; }
			else if (key == "MAX_Y") { sline >> max_y; }
			else if (key == "MAX_Z") { sline >> max_z; }
			else if (key == "MAX_DEPTH_VACANT") { sline >> max_depth_vacant; }
			else if (key == "MAX_ABSTRACT_DEPTH") { sline >> max_abstract_depth; }
			else if (key == "MAX_ABSTRACT_LIMIT_V") { sline >> max_abstract_limit_v; }
			else if (key == "MIN_ABSTRACT_DEPTH") { sline >> min_abstract_depth; }
			break;
		case 'O':
			if (key == "OUT_SLOP") { sline >> out_slop; }
			else if (key == "OUTPUT_RANK") { sline >> output_rank; }
			break;
		case 'P':
			if (key == "PDB_OUTPUT") { sline >> pdb_output; }
			break;
		case 'R':
			if (key == "RECEPTOR_FILE") {
				sline >> pdb_file;
				autoname = pdb_file.substr(0, pdb_file.rfind('.'));
				_v_surface_file = autoname + "_surface.pdb";
				_v_vacant_file = autoname + "_vacant.pdb";
				_v_atom_file = autoname + "_cavity.pdb";
			}
			else if (key == "RADIUS_LENGTH") { sline >> radius_length; }
			else if (key == "RADIUS_MIS") { sline >> radius_mis; }
			else if (key == "RIGID_LIMIT") { sline >> rigid_limit; }
			else if (key == "RULER_1") { sline >> ruler_1; }
			else if (key == "RUN_INFO") { sline >> info; }
			break;
		case 'S':
			if (key == "SEPARATE_CHIP_LIMIT_V") { sline >> chip_v; }
			else if (key == "SEPARATE_MIN_DEPTH") { sline >> min_depth; }
			else if (key == "SEPARATE_CHIP_DEPTH") { sline >> chip_depth; }
			else if (key == "SEPARATE_MAX_LIMIT_V") { sline >> max_limit_v; }
			else if (key == "SINGLE_LIGAND") { sline >> single_ligand; }
			break;
		case 'V':
			if (key == "VACANT_ADJUST") { sline >> vacant_adjust; }
			else if (key == "VISUAL_OUTPUT") { sline >> visual_output; }
			else if (key == "V_POINT_NUM") { sline >> v_num; }
			else if (key == "V_OUTPUT_MOD") { sline >> v_mod; }
			break;
		default: break;
		}
		// clear string stream
		sline.str(std::string());
		sline.clear();
	}
	infile.close();
}

void Parameter::Initialize(cmdline::parser& cmdargs)
{
	debug = cmdargs.exist("debug");
	string in_dir, protname, ext, tmppdb_file;
	visual_output = cmdargs.get<int>("visual_output");
	if (cmdargs.exist("verbose")) info = 1;
	pdb_file = cmdargs.get<string>("receptor");
	lig_file = cmdargs.get<string>("ligand");
	//detect_mode = cmdargs.get<int>("mode");
	if (lig_file != "") detect_mode = 1;
	else detect_mode = 0;

	radius_length = cmdargs.get<float>("radius_length");
	min_depth = cmdargs.get<int>("min_depth");
	max_abstract_limit_v = cmdargs.get<int>("min_abstract_limit_v");
	max_limit_v = cmdargs.get<int>("max_limit_v");
	min_abstract_depth = cmdargs.get<int>("min_abstract_depth");
	ruler_1 = cmdargs.get<int>("minVnum");
	output_rank = cmdargs.get<float>("rank_score");

	chip_v = cmdargs.get<int>("chip_v");
	soft_separate = cmdargs.exist("softsep");
	softsep_limit = cmdargs.get<int>("softsep_limit");
	rescue = cmdargs.exist("rescue");
	rescue_limit = cmdargs.get<int>("rescue_limit");
	out_dir = cmdargs.get<string>("outdir");
	tmppdb_file = pdb_file;
	ext = pdb_file.substr(pdb_file.find_last_of('.'));
	if (ext == ".gz")
		pdb_file = pdb_file.erase(pdb_file.find_last_of('.'));
	int flag = 0;
#ifdef __GNUC__
	in_dir = pdb_file.substr(0, pdb_file.find_last_of('/'));
	protname = pdb_file.substr(pdb_file.find_last_of('/') + 1);
	protname = protname.erase(protname.find_last_of('.'));
	if (out_dir == "")
	{
		autoname = in_dir + "/" + protname;
	}
	else
	{
		if (access(out_dir.c_str(), 0) == -1)
		{
			cout << out_dir << " is not exist. Now create it ..." << endl;
			flag = mkdir(out_dir.c_str(), 0775);
			if (flag != 0)
			{
				cout << "Cannot create this directory: " << out_dir << ". Use protein's directory." << endl;
				autoname = in_dir + "/" + protname;
			}
			else
			{
				Check_Directory(out_dir);
				autoname = out_dir + protname;
			}
		}
		else
		{
			Check_Directory(out_dir);
			autoname = out_dir + protname;
		}
	}
#elif defined(_MSC_VER)
	in_dir = pdb_file.substr(0, pdb_file.find_last_of('\\'));
	protname = pdb_file.substr(pdb_file.find_last_of('\\') + 1);
	protname = protname.erase(protname.find_last_of('.'));
	if (out_dir == "")
	{
		autoname = in_dir + "\\" + protname;
	}
	else 
	{
		if (_access(out_dir.c_str(), 0) == -1)
		{
			cout << out_dir << " is not exist. Now create it ..." << endl;
			flag = _mkdir(out_dir.c_str());
			if (flag != 0)
			{
				cout << "Cannot create this directory: " << out_dir << ". Use protein's directory." << endl;
				autoname = in_dir + "\\" + protname;
			}
			else
			{
				Check_Directory(out_dir, '\\');
				autoname = out_dir + protname;
			}
		}
		else
		{
			Check_Directory(out_dir, '\\');
			autoname = out_dir + protname;
		}
	}
#else
		cout << "Cannot create this directory: " << out_dir << ". Use protein's directory." << endl;
		autoname = pdb_file.substr(0, pdb_file.rfind('.'));
#endif
	_v_surface_file = autoname + "_surface.pdb";
	_v_vacant_file = autoname + "_vacant.pdb";
	_v_atom_file = autoname + "_cavity.pdb";
	
	pdb_file = tmppdb_file;
	
}

void Parameter::Process()
{
	int num_error = 0;
	if (pdb_file == "") {
		Missing_Parameter_Error("RECEPTOR_FILE");
		num_error++;
	}
	if (lig_file != "" && detect_mode != 1)
	{
		cout << "\n***Warning: Ligand file is given but not used for detect_mode 0 \n\n";
	}
	if ((detect_mode != 0) && (detect_mode != 1) && (detect_mode != 2))
	{
		Missing_Parameter_Error("DETECT_MODE");
		num_error++;
	}
	if (detect_mode == 1 && lig_file == "")
	{
		Missing_Parameter_Error("LIGAND_FILE");
		num_error++;
	}
	//Check_Directory(parameter_dir);
	if (num_error != 0)
	{
		cout << "\n" << num_error << " errors have been detected in the parameter file." << endl;
		cout << "Please correct them and try again." << endl;
		exit(1);
	}
	return;
}

