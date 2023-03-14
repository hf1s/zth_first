#include <iostream>
#include <string>

#include "utils.h"
#include "parameter.h"
#include "protein.h"
#include "ligand.h"
#include "pocket.h"
#include "cmdline.h"

Parameter* parm;
Protein* prot;
Ligand* lig;
using namespace std;
int hetres = -1;

int main(int argc, char* argv[])
{
	Clock myclock;
	myclock.begin();

	// --receptor example\2ikg\2ikgA.pdb -v
	// --receptor example\2nvd\2nvdA.pdb -v
	// --receptor example\3dn5\3dn5A.pdb -v
	string version = "Version 2022.09";
	cmdline::parser cmdargs(version);
	cmdargs.add("help", 'h', "show this message");
	cmdargs.add<string>("receptor", 'r', "receptor file", true, "");
	cmdargs.add<string>("ligand", 'l', "ligand file (necessary when mode=1)", false, "");
	//cmdargs.add<int>("mode", 'm', "detect mode: 0: whole protein mode, 1: ligand detection mode, 2: area detection mode", false, 0, cmdline::oneof<int>(0, 1, 2));
	cmdargs.add<string>("outdir", 'o', "output directory(default is the protein's directory", false, "");
	cmdargs.add<int>("visual_output", '\0', "control output files. 0: close 1: only output surface grid 2: output surface and vacant grid", false, 2, cmdline::oneof<int>(0, 1, 2));
	
	cmdargs.add<float>("radius_length", '\0', "the radius length of eraser ball", false, 10);
	cmdargs.add<int>("chip_v", '\0', "control the minimal size of shrinked cavities. Cavity size which is less than chip_v will be ignored", false, 300);
	cmdargs.add<int>("min_depth", '\0', "separate minimal depth", false, 8);
	cmdargs.add<int>("min_abstract_limit_v", '\0', "max_abstract_limit_v", false, 1500);
	cmdargs.add<int>("max_limit_v", '\0', "separate_max_limit_v", false, 6000);
	cmdargs.add<int>("min_abstract_depth", '\0', "min_abstract_depth", false, 2);
	cmdargs.add<int>("minVnum", '\0', "control the minimal size of final cavities. Cavity size which is less than minVnum will be skipped.", false, 100);
	cmdargs.add<float>("rank_score", '\0', "cutoff of rank score. Only output cavities whose rank score is larger then cutoff.", false, 1.5);

	cmdargs.add("softsep", '\0', "Do not separate small cavities");
	cmdargs.add<int>("softsep_limit", '\0', "conrol the softsep mode. The maximal size of shrinked cavities which are going to apply softsep", false, 3000);
	cmdargs.add("rescue", '\0', "rescue large cavities which was accidentally excluded");
	cmdargs.add<int>("rescue_limit", '\0', "control the rescue mode. The maximal size of shrinked cavities which are going to apply rescue.", false, 1000);
	cmdargs.add("verbose", 'v', "show more info");
	cmdargs.add("version", 'V', "show the version information.");
	cmdargs.add("debug", '\0', "debug mode");
	stringstream oss;
	oss << "\n    example 1: " << argv[0] << " --receptor=protein.pdb" << endl;
	oss << "    example 1: " << argv[0] << " --receptor=protein.pdb --ligand=ligand.mol2 --mode=1" << endl;
	string fstr = oss.str();
	cmdargs.footer(fstr);

	cmdargs.parse_check(argc, argv);

	cout << version << endl;
	cout << myclock.strcurrent();
	cout << "\nCavity detection starts ..." << endl;

	cout << "*** Now parsing parameters ..." << endl;
	Parameter parameter; parm = &parameter;
	parameter.version = version;
	parameter.Initialize(cmdargs);
	parameter.Process();

	Protein protein; prot = &protein;
	protein.Read_RESIDUE();
	Ligand ligand; lig = &ligand;
	
	if (parm->judge == 0) { cout << "Search Rule: Surface Separate ..." << endl; }
	else if (parm->judge == 1) { cout << "Search Rule: Depth Limit ..." << endl; }
	else if (parm->judge == 2)
	{ 
		cout << "Search Rule: Depth-Volume Limit ..." << endl;
	}
	else
	{
		cout << "Wrong judge rule!" << endl;
		exit(1);
	}
	
	Pocket pocket;
	if (parm->detect_mode == 0)
	{
		cout << "Search Mode: Whole Protein Mode ..." << endl;
		cout << "Now reading the protein from '" << parm->pdb_file << "' ..." << endl;
		protein.Read_PDB(parm->pdb_file.c_str());
		cout << "Now defining the box size of protein ..." << endl;
		protein.Protein_Box();
	}
	else if (parm->detect_mode == 1)
	{
		ligand.Read_ATOM_DEF();
		cout << "Search Mode : Ligand Mode ..." << endl;
		cout << "Now reading the ligand from '" << parm->lig_file << "' ..." << endl;
		ligand.Read_From_Mol2(parm->lig_file.c_str());
		cout << "Now define the box by ligand ..." << endl;
		ligand.Ligand_Box();
		cout << "Now reading the protein from '" << parm->pdb_file << "' ..." << endl;
		protein.Read_PDB(parm->pdb_file.c_str());
	}

	cout << "Now analyzing the protein ..." << endl;
	protein.Value_Atom();
	
	cout << "Now detecting cavity in the box ..." << endl;
	if (parm->detect_mode == 0) { cout << "Whole protein search need more time, please wait ..." << endl; }
	pocket.Cavity_Search();

	cout << "Cavity has done the job successfully!\n" << endl;
	myclock.end();
	cout << myclock.strcurrent();
	cout << "Running time: " << myclock.duration() << " seconds." << endl;
	return 0;
}
