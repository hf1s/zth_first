#include "protein.h"
#include "utils.h"
#include "parameter.h"
#include <set>
#include <vector>
#include <cmath>
#ifdef __GNUC__
#include <zlib.h>
#endif

using namespace std;

vector<string> get_residue_def()
{
	vector<string> residue_def;
	residue_def.push_back(string("#"));
	residue_def.push_back(string("# This table contains the parameters applied to the target protein."));
	residue_def.push_back(string("# Twenty amino acids as well as several specific residues are presented."));
	residue_def.push_back(string("# The lines headed by \"ATOM\" are organized as following:"));
	residue_def.push_back(string("# The 1st column: heading"));
	residue_def.push_back(string("# The 2nd column: name of the atom used in PDB file"));
	residue_def.push_back(string("# The 3rd column: atom type used in AMBER force field"));
	residue_def.push_back(string("# The 4th column: atomic hydrophobic scale, derived from XLOGP2"));
	residue_def.push_back(string("# The 5th column: hydrogen bonding character"));
	residue_def.push_back(string("# The 6th column: atomic charge"));
	residue_def.push_back(string("# The 7th column: VDW radius from AMBER force field"));
	residue_def.push_back(string("# The 8th column: atom type used in Tripos force field "));
	residue_def.push_back(string("# The 9th column: max acceptor num of atom"));
	residue_def.push_back(string("# The 10th column: max donor num of atom"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("#----------------------------------------------------------------------------"));
	residue_def.push_back(string("#	 	 XlogP2	HB  Charge  r  	 Type   Acceptor Donor"));
	residue_def.push_back(string("#----------------------------------------------------------------------------"));
	residue_def.push_back(string("RESI ACE  11       0.00  !  Acetyl  "));
	residue_def.push_back(string("ATOM CA   CT   	 0.267	N -0.3662  1.94  C.3	0 0 0  ! acetylated N-terminus"));
	residue_def.push_back(string("ATOM CH3  CT     0.267  N -0.3662  1.94  C.3	0 0 0  !"));
	residue_def.push_back(string("ATOM CAY  CT   	 0.267	N -0.3662  1.94  C.3	0 0 0  ! HY1 HY2 HY3              "));
	residue_def.push_back(string("ATOM HA1  HC   	 0.000	N  0.1123  0.00  H	0 0 0  !    \\ | /    "));
	residue_def.push_back(string("ATOM HA2  HC   	 0.000	N  0.1123  0.00  H	0 0 0  !     CAY     "));
	residue_def.push_back(string("ATOM HA3  HC   	 0.000	N  0.1123  0.00  H	0 0 0  !      |      "));
	residue_def.push_back(string("ATOM HY1  HC   	 0.000	N  0.1123  0.00  H	0 0 0  !      C=O  "));
	residue_def.push_back(string("ATOM HY2  HC   	 0.000	N  0.1123  0.00  H	0 0 0  !      |      "));
	residue_def.push_back(string("ATOM HY3  HC   	 0.000	N  0.1123  0.00  H	0 0 0  "));
	residue_def.push_back(string("ATOM C    C    	-0.030	N  0.5972  1.90  C.2	0 0 0  "));
	residue_def.push_back(string("ATOM O    O    	-0.399	A -0.5679  1.66  O.2	2 0 2	"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI ALA  11      0.00  !  Alanine"));
	residue_def.push_back(string("ATOM N    N    	-0.096	D -0.4157  1.83  N.am	0 1 1  !     |"));
	residue_def.push_back(string("ATOM H    H    	 0.000	N  0.2719  0.00  H	0 0 0  !  HN-N"));
	residue_def.push_back(string("ATOM CA   CT   	-0.305	N  0.0337  1.94  C.3	0 0 0  !     |     HB1"));
	residue_def.push_back(string("ATOM HA   H1   	 0.000	N  0.0823  0.00  H	0 0 0  !     |    /"));
	residue_def.push_back(string("ATOM CB   CT   	 0.528	N -0.1825  1.94  C.3	0 0 0  !  HA-CA--CB-HB2            ALANINE"));
	residue_def.push_back(string("ATOM HB1  HC   	 0.000	N  0.0603  0.00  H	0 0 0  !     |    \\      "));
	residue_def.push_back(string("ATOM HB2  HC   	 0.000	N  0.0603  0.00  H	0 0 0  !     |     HB2 "));
	residue_def.push_back(string("ATOM HB3  HC   	 0.000	N  0.0603  0.00  H	0 0 0  !   O=C   "));
	residue_def.push_back(string("ATOM C    C    	-0.030	N  0.5973  1.90  C.2	0 0 0  !     |"));
	residue_def.push_back(string("ATOM O    O    	-0.399	A -0.5679  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5679  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI ARG  25      1.00  !  Arginine"));
	residue_def.push_back(string("ATOM N    N    	-0.096	D -0.3479  1.83  N.am	0 1 1   !     |                      HH11"));
	residue_def.push_back(string("ATOM H    H    	 0.000	N  0.2747  0.00  H	0 0 0   !  HN-N                       |"));
	residue_def.push_back(string("ATOM CA   CT   	-0.305	N -0.2637  1.94  C.3	0 0 0   !     |   HB1 HG1 HD1 HE     NH1-HH12"));
	residue_def.push_back(string("ATOM HA   H1   	 0.000	N  0.1560  0.00  H	0 0 0   !     |   |   |   |   |    //(+)  "));
	residue_def.push_back(string("ATOM CB   CT   	 0.358 	N -0.0007  1.94  C.3	0 0 0   !  HA-CA--CB--CG--CD--NE--CZ          "));
	residue_def.push_back(string("ATOM HB1  HC   	 0.000	N  0.0327  0.00  H	0 0 0   !     |   |   |   |         \\         "));
	residue_def.push_back(string("ATOM HB2  HC   	 0.000	N  0.0327  0.00  H	0 0 0   !     |   HB2 HG2 HD2        NH2-HH22 "));
	residue_def.push_back(string("ATOM CG   CT   	 0.358	N  0.0390  1.94  C.3	0 0 0   !   O=C                       |       "));
	residue_def.push_back(string("ATOM HG1  HC   	 0.000	N  0.0285  0.00  H	0 0 0   !     |                      HH21     "));
	residue_def.push_back(string("ATOM HG2  HC   	 0.000	N  0.0285  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CD   CT   	-0.137  N  0.0486  1.94  C.3	0 0 0"));
	residue_def.push_back(string("ATOM HD1  H1   	 0.000	N  0.0687  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HD2  H1   	 0.000	N  0.0687  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM NE   N2   	-0.096	D -0.5295  1.83  N.pl3  0 1 1"));
	residue_def.push_back(string("ATOM HE   H    	 0.000	N  0.3456  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CZ   CA   	 0.005	N  0.8076  1.90  C.cat	0 0 0"));
	residue_def.push_back(string("ATOM NH1  N2   	-0.646	D -0.8627  1.83  N.pl3  0 2 2"));
	residue_def.push_back(string("ATOM HH11 H    	 0.000	N  0.4478  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HH12 H    	 0.000	N  0.4478  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM NH2  N2   	-0.646	D -0.8627  1.83  N.pl3	0 2 2"));
	residue_def.push_back(string("ATOM HH21 H    	 0.000	N  0.4478  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HH22 H    	 0.000	N  0.4478  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM C    C    	-0.030	N  0.7341  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM O    O    	-0.399	A -0.5894  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5894  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI ASN  17      0.00  !  Asparagine"));
	residue_def.push_back(string("ATOM N    N    	-0.096	D -0.4157  1.83  N.am	0 1 1 !     |                                   "));
	residue_def.push_back(string("ATOM H    H    	 0.000	N  0.2719  0.00  H	0 0 0 !  HN-N                                   "));
	residue_def.push_back(string("ATOM CA   CT   	-0.305	N  0.0143  1.94  C.3	0 0 0 !     |   HB1 OD1    HD21 (cis to OD1)    "));
	residue_def.push_back(string("ATOM HA   H1   	 0.000	N  0.1048  0.00  H	0 0 0 !     |   |   ||    /                     "));
	residue_def.push_back(string("ATOM CB   CT   	-0.008	N -0.2041  1.94  C.3	0 0 0 !  HA-CA--CB--CG--ND2                     "));
	residue_def.push_back(string("ATOM HB1  HC   	 0.000	N  0.0797  0.00  H	0 0 0 !     |   |         \\                     "));
	residue_def.push_back(string("ATOM HB2  HC   	 0.000	N  0.0797  0.00  H	0 0 0 !     |   HB2        HD22 (trans to OD1)  "));
	residue_def.push_back(string("ATOM CG   C    	-0.030	N  0.7130  1.90  C.2	0 0 0 !   O=C                                   "));
	residue_def.push_back(string("ATOM OD1  O    	-0.399	A -0.5931  1.66  O.2	2 0 2 !     |                                   "));
	residue_def.push_back(string("ATOM AD1  O     -0.399  A -0.5931  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM ND2  N    	-0.646	D -0.9196  1.83  N.am	0 2 2"));
	residue_def.push_back(string("ATOM AD2  N     -0.646  D -0.9196  1.83  N.am	0 2 2"));
	residue_def.push_back(string("ATOM HD21 H    	 0.000	N  0.4196  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HD22 H    	 0.000	N  0.4196  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM C    C    	-0.030	N  0.5973  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM O    O    	-0.399	A -0.5679  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5679  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI ASP  14     -1.00  !  Aspartic acid"));
	residue_def.push_back(string("ATOM N    N    	-0.096	D -0.5163  1.83  N.am	0 1 1  !     |                 "));
	residue_def.push_back(string("ATOM H    H    	 0.000	N  0.2936  0.00  H	0 0 0  !  HN-N                 "));
	residue_def.push_back(string("ATOM CA   CT   	-0.305	N  0.0381  1.94  C.3	0 0 0  !     |   HB1   OD1     "));
	residue_def.push_back(string("ATOM HA   H1   	 0.000	N  0.0880  0.00  H	0 0 0  !     |   |    //       "));
	residue_def.push_back(string("ATOM CB   CT   	 0.008	N -0.0303  1.94  C.3	0 0 0  !  HA-CA--CB--CG        "));
	residue_def.push_back(string("ATOM HB1  HC   	 0.000	N -0.0122  0.00  H	0 0 0  !     |   |    \\        "));
	residue_def.push_back(string("ATOM HB2  HC   	 0.000	N -0.0122  0.00  H	0 0 0  !     |   HB2   OD2(-)  HD2"));
	residue_def.push_back(string("ATOM CG   C    	-0.030	N  0.7994  1.90  C.2	0 0 0  !   O=C                 "));
	residue_def.push_back(string("ATOM OD1  O2   	-0.880	A -0.8014  1.66  O.co2	2 0 2  !     |                 "));
	residue_def.push_back(string("ATOM OD2  O2   	-0.880	A -0.8014  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("ATOM HD2  HO	 0.000  N  0.4641  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM C    C    	-0.030	N  0.5366  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM O    O    	-0.399	A -0.5819  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5819  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI CYS  12      0.00  !  Cysteine"));
	residue_def.push_back(string("ATOM N    N    	-0.096	D -0.4157  1.83  N.am	0 1 1  !     |              "));
	residue_def.push_back(string("ATOM H    H    	 0.000	N  0.2719  0.00  H	0 0 0  !  HN-N              "));
	residue_def.push_back(string("ATOM CA   CT   	-0.305	N  0.0213  1.94  C.3	0 0 0  !     |   HB1        "));
	residue_def.push_back(string("ATOM HA   H1   	 0.000	N  0.1124  0.00  H	0 0 0  !     |   |          "));
	residue_def.push_back(string("ATOM CB   CT   	 0.358 	N -0.1231  1.94  C.3	0 0 0  !  HA-CA--CB--SG     or SG--SG"));
	residue_def.push_back(string("ATOM HB1  H1   	 0.000	N  0.1112  0.00  H	0 0 0  !     |   |     \\    "));
	residue_def.push_back(string("ATOM HB2  H1   	 0.000	N  0.1112  0.00  H	0 0 0  !     |   HB2    HG1 "));
	residue_def.push_back(string("ATOM SG   SH   	 0.419	N -0.3119  2.09  S.3	0 1 3  !   O=C              "));
	residue_def.push_back(string("ATOM HG   HS   	 0.000	N  0.1933  0.00  H	0 0 0  !     |              "));
	residue_def.push_back(string("ATOM C    C    	-0.030	N  0.5973  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM O    O    	-0.399	A -0.5679  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5679  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI GLN  20      0.00  !  Glutamine"));
	residue_def.push_back(string("ATOM N    N    	-0.096	D -0.4157  1.83  N.am	0 1 1    !     |                                     "));
	residue_def.push_back(string("ATOM H    H    	 0.000	N  0.2719  0.00  H	0 0 0    !  HN-N                                     "));
	residue_def.push_back(string("ATOM CA   CT   	-0.305	N -0.0031  1.94  C.3	0 0 0    !     |   HB1 HG1 OE1   HE21 (cis to OE1)   "));
	residue_def.push_back(string("ATOM HA   H1   	 0.000	N  0.0850  0.00  H	0 0 0    !     |   |   |   ||    /                   "));
	residue_def.push_back(string("ATOM CB   CT   	 0.358	N -0.0036  1.94  C.3	0 0 0    !  HA-CA--CB--CG--CD--NE2                   "));
	residue_def.push_back(string("ATOM HB1  HC   	 0.000	N  0.0171  0.00  H	0 0 0    !     |   |   |         \\                   "));
	residue_def.push_back(string("ATOM HB2  HC   	 0.000	N  0.0171  0.00  H	0 0 0    !     |   HB2 HG2       HE22 (trans to OE1) "));
	residue_def.push_back(string("ATOM CG   CT   	-0.008	N -0.0645  1.94  C.3	0 0 0    !   O=C                                     "));
	residue_def.push_back(string("ATOM HG1  HC   	 0.000	N  0.0352  0.00  H	0 0 0    !     |                                     "));
	residue_def.push_back(string("ATOM HG2  HC   	 0.000	N  0.0352  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CD   C    	-0.030	N  0.6951  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM OE1  O    	-0.399	A -0.6086  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM AE1  O     -0.399  A -0.6086  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM NE2  N    	-0.646	D -0.9407  1.83  N.am	0 2 2"));
	residue_def.push_back(string("ATOM AE2  N     -0.646  D -0.9407  1.83  N.am	0 2 2"));
	residue_def.push_back(string("ATOM HE21 H    	 0.000	N  0.4251  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HE22 H    	 0.000	N  0.4251  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM C    C    	-0.030	N  0.5973  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM O    O    	-0.399	A -0.5679  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5679  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI GLU  17     -1.00  !  Glutamic acid"));
	residue_def.push_back(string("ATOM N    N    	-0.096	D -0.5163  1.83  N.am	0 1 1  !     |                       "));
	residue_def.push_back(string("ATOM H    H    	 0.000	N  0.2936  0.00  H	0 0 0  !  HN-N                       "));
	residue_def.push_back(string("ATOM CA   CT   	-0.305	N  0.0397  1.94  C.3	0 0 0  !     |   HB1 HG1   OE1       "));
	residue_def.push_back(string("ATOM HA   H1   	 0.000	N  0.1105  0.00  H	0 0 0  !     |   |   |    //         "));
	residue_def.push_back(string("ATOM CB   CT   	 0.358	N  0.0560  1.94  C.3	0 0 0  !  HA-CA--CB--CG--CD          "));
	residue_def.push_back(string("ATOM HB1  HC   	 0.000	N -0.0173  0.00  H	0 0 0  !     |   |   |    \\          "));
	residue_def.push_back(string("ATOM HB2  HC   	 0.000	N -0.0173  0.00  H	0 0 0  !     |   HB2 HG2   OE2(-)   HE2 "));
	residue_def.push_back(string("ATOM CG   CT   	-0.008	N  0.0136  1.94  C.3	0 0 0  !   O=C                       "));
	residue_def.push_back(string("ATOM HG1  HC   	 0.000	N -0.0425  0.00  H	0 0 0  !     |                       "));
	residue_def.push_back(string("ATOM HG2  HC   	 0.000	N -0.0425  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CD   C    	-0.030	N  0.8054  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM OE1  O2   	-0.880	A -0.8188  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("ATOM OE2  O2   	-0.880	A -0.8188  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("ATOM HE2  HO	 0.000  N  0.4641  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM C    C    	-0.030	N  0.5366  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM O    O    	-0.399	A -0.5819  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5819  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI GLY  8       0.00  !  Glycine		            !     |       "));
	residue_def.push_back(string("ATOM N    N    	-0.096	D -0.4157  1.83  N.am	0 1 1       !     N-H     "));
	residue_def.push_back(string("ATOM H    H    	 0.000	N  0.2719  0.00  H	0 0 0       !     |       "));
	residue_def.push_back(string("ATOM CA   CT   	-0.303	N -0.0252  1.94  C.3	0 0 0       !     |       "));
	residue_def.push_back(string("ATOM HA1  H1   	 0.000	N  0.0698  0.00  H	0 0 0       ! HA1-CA-HA2  "));
	residue_def.push_back(string("ATOM HA2  H1   	 0.000	N  0.0698  0.00  H	0 0 0       !     |       "));
	residue_def.push_back(string("ATOM C    C    	-0.030	N  0.5973  1.90  C.2	0 0 0       !     |       "));
	residue_def.push_back(string("ATOM O    O    	-0.399	A -0.5679  1.66  O.2	2 0 2       !     C=O     "));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5679  1.66  O.co2	2 0 2       !     |       "));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI HIS  19      0.90  !  Histidine (uncharged wit h proton on nd1)"));
	residue_def.push_back(string("ATOM N    N    	-0.096	D -0.4157  1.83  N.am	0 1 1  !     |          HD1    HE1     !     |                 HE1    !     |          HD1    HE1 "));
	residue_def.push_back(string("ATOM H    H    	 0.000	N  0.2719  0.00  H	0 0 0  !  HN-N           |     /       !  HN-N             __  /      !  HN-N           |     /   "));
	residue_def.push_back(string("ATOM CA   CT   	-0.305	N  0.0188  1.94  C.3	0 0 0  !     |   HB1    ND1--CE1       !     |   HB1    ND1--CE1      !     |   HB1    ND1--CE1   "));
	residue_def.push_back(string("ATOM HA   H1   	 0.000	N  0.0881  0.00  H	0 0 0  !     |   |     /      ||       !     |   |     /      |       !     |   |     /      ||   "));
	residue_def.push_back(string("ATOM CB   CT   	-0.008	N -0.0462  1.94  C.3	0 0 0  !  HA-CA--CB--CG       ||       !  HA-CA--CB--CG       |       !  HA-CA--CB--CG       ||   "));
	residue_def.push_back(string("ATOM HB1  HC   	 0.000	N  0.0402  0.00  H	0 0 0  !     |   |     \\     ||       !     |   |     \\     |       !     |   |     \\     ||   "));
	residue_def.push_back(string("ATOM HB2  HC   	 0.000	N  0.0402  0.00  H	0 0 0  !     |   HB2    CD2--NE2       !     |   HB2    CD2--NE2      !     |   HB2    CD2--NE2(+)"));
	residue_def.push_back(string("ATOM CG   CC   	-0.027	N -0.0266  1.90  C.ar	0 0 0  !   O=C           |             !   O=C           |     \\      !   O=C           |     \\   "));
	residue_def.push_back(string("ATOM ND1  NA   	 0.135	D -0.3811  1.87  N.ar	1 1 1  !     |          HD2            !     |          HD2    HE2    !     |          HD2    HE2 "));
	residue_def.push_back(string("ATOM HD1  H    	 0.000	N  0.3649  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CD2  CV   	-0.310	N  0.1292  1.90  C.ar	0 0 0"));
	residue_def.push_back(string("ATOM HD2  H4   	 0.000	N  0.1147  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CE1  CR   	-0.310	N  0.2057  1.90  C.ar	0 0 0"));
	residue_def.push_back(string("ATOM HE1  H5   	 0.000	N  0.1392  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM NE2  NB   	 0.135	D -0.5727  1.87  N.ar	1 1 1"));
	residue_def.push_back(string("ATOM HE2  H5   	 0.000	N  0.3339  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM C    C    	-0.030	N  0.5973  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM O    O    	-0.399	A -0.5679  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5679  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI ILE  20      0.00  !  Isoleucine"));
	residue_def.push_back(string("ATOM N    N    	-0.096	D -0.4157  1.83  N.am	0 1 1  !     |    HG21 HG22     "));
	residue_def.push_back(string("ATOM H    H    	 0.000	N  0.2719  0.00  H	0 0 0  !  HN-N      | /         "));
	residue_def.push_back(string("ATOM CA   CT   	-0.305	N -0.0597  1.94  C.3	0 0 0  !     |     CG2--HG23    "));
	residue_def.push_back(string("ATOM HA   H1   	 0.000	N  0.0869  0.00  H	0 0 0  !     |    /             "));
	residue_def.push_back(string("ATOM CB   CT   	 0.127	N  0.1303  1.94  C.3	0 0 0  !  HA-CA--CB-HB    HD1   "));
	residue_def.push_back(string("ATOM HB   HC   	 0.000	N  0.0187  0.00  H	0 0 0  !     |    \\       /     "));
	residue_def.push_back(string("ATOM CG2  CT   	 0.528	N -0.3204  1.94  C.3	0 0 0  !     |     CG1--CD--HD2 "));
	residue_def.push_back(string("ATOM HG21 HC   	 0.000	N  0.0882  0.00  H	0 0 0  !   O=C    / \\     \\     "));
	residue_def.push_back(string("ATOM HG22 HC   	 0.000	N  0.0882  0.00  H	0 0 0  !     | HG11 HG12  HD3   "));
	residue_def.push_back(string("ATOM HG23 HC   	 0.000	N  0.0882  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CG1  CT   	 0.358	N -0.0430  1.94  C.3	0 0 0"));
	residue_def.push_back(string("ATOM HG11 HC   	 0.000	N  0.0236  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HG12 HC   	 0.000	N  0.0236  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CD1  CT   	 0.528	N -0.0660  1.94  C.3	0 0 0"));
	residue_def.push_back(string("ATOM HD11 HC   	 0.000	N  0.0186  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HD12 HC   	 0.000	N  0.0186  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HD13 HC   	 0.000	N  0.0186  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM C    C    	-0.030	N  0.5973  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM O    O    	-0.399	A -0.5679  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5679  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI LEU  20      0.00  !  Leucine"));
	residue_def.push_back(string("ATOM N    N    	-0.096	D -0.4157  1.83  N.am	0 1 1    !     |        HD11 HD12  "));
	residue_def.push_back(string("ATOM H    H    	 0.000	N  0.2719  0.00  H	0 0 0    !  HN-N          | /      "));
	residue_def.push_back(string("ATOM CA   CT   	-0.305	N -0.0518  1.94  C.3	0 0 0    !     |   HB1   CD1--HD13 "));
	residue_def.push_back(string("ATOM HA   H1   	 0.000	N  0.0922  0.00  H	0 0 0    !     |   |    /          "));
	residue_def.push_back(string("ATOM CB   CT   	 0.358	N -0.1102  1.94  C.3	0 0 0    !  HA-CA--CB--CG-HG       "));
	residue_def.push_back(string("ATOM HB1  HC   	 0.000	N  0.0457  0.00  H	0 0 0    !     |   |    \\          "));
	residue_def.push_back(string("ATOM HB2  HC   	 0.000	N  0.0457  0.00  H	0 0 0    !     |   HB2   CD2--HD23 "));
	residue_def.push_back(string("ATOM CG   CT   	 0.127	N  0.3531  1.94  C.3	0 0 0    !   O=C          | \\      "));
	residue_def.push_back(string("ATOM HG   HC   	 0.000	N -0.0361  0.00  H	0 0 0    !     |        HD21 HD22  "));
	residue_def.push_back(string("ATOM CD1  CT   	 0.528	N -0.4121  1.94  C.3	0 0 0"));
	residue_def.push_back(string("ATOM HD11 HC   	 0.000	N  0.1000  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HD12 HC   	 0.000	N  0.1000  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HD13 HC   	 0.000	N  0.1000  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CD2  CT   	 0.528	N -0.4121  1.94  C.3	0 0 0"));
	residue_def.push_back(string("ATOM HD21 HC   	 0.000	N  0.1000  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HD22 HC   	 0.000	N  0.1000  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HD23 HC   	 0.000	N  0.1000  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM C    C    	-0.030	N  0.5973  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM O    O    	-0.399	A -0.5679  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5679  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI LYS  23      1.00  !  Lysine"));
	residue_def.push_back(string("ATOM N    N    	-0.096	D -0.3479  1.83  N.am	0 1 1   !     |                             "));
	residue_def.push_back(string("ATOM H    H    	 0.000	N  0.2747  0.00  H	0 0 0   !  HN-N                             "));
	residue_def.push_back(string("ATOM CA   CT   	-0.305	N -0.2400  1.94  C.3	0 0 0   !     |   HB1 HG1 HD1 HE1    HZ1    "));
	residue_def.push_back(string("ATOM HA   HC   	 0.000	N  0.1426  0.00  H	0 0 0   !     |   |   |   |   |     /       "));
	residue_def.push_back(string("ATOM CB   CT   	 0.358	N -0.0094  1.94  C.3	0 0 0   !  HA-CA--CB--CG--CD--CE--NZ--HZ2   "));
	residue_def.push_back(string("ATOM HB1  HC   	 0.000	N  0.0362  0.00  H	0 0 0   !     |   |   |   |   |     \\       "));
	residue_def.push_back(string("ATOM HB2  HC   	 0.000	N  0.0362  0.00  H	0 0 0   !     |   HB2 HG2 HD2 HE2    HZ3    "));
	residue_def.push_back(string("ATOM CG   CT   	 0.358	N  0.0187  1.94  C.3	0 0 0   !   O=C                             "));
	residue_def.push_back(string("ATOM HG1  HC   	 0.000	N  0.0103  0.00  H	0 0 0   !     |                             "));
	residue_def.push_back(string("ATOM HG2  HC   	 0.000	N  0.0103  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CD   CT   	 0.358	N -0.0479  1.94  C.3	0 0 0"));
	residue_def.push_back(string("ATOM HD1  HC   	 0.000	N  0.0621  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HD2  HC   	 0.000	N  0.0621  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CE   CT   	-0.137	N -0.0143  1.94  C.3	0 0 0"));
	residue_def.push_back(string("ATOM HE1  HP   	 0.000	N  0.1135  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HE2  HP   	 0.000	N  0.1135  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM NZ   N3   	-2.000	D -0.3854  1.86  N.4	0 3 3 "));
	residue_def.push_back(string("ATOM HZ1  H    	 0.000	N  0.3400  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HZ2  H    	 0.000	N  0.3400  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HZ3  H    	 0.000	N  0.3400  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM C    C   	-0.030	N  0.7341  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM O    O    	-0.399	A -0.5894  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5894  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI MET  18      0.00  !  Methionine"));
	residue_def.push_back(string("ATOM N    N    	-0.096	D -0.4157  1.83  N.am	0 1 1   !     |                        "));
	residue_def.push_back(string("ATOM H    H    	 0.000	N  0.2719  0.00  H	0 0 0   !  HN-N                        "));
	residue_def.push_back(string("ATOM CA   CT   	-0.305	N -0.0237  1.94  C.3	0 0 0   !     |   HB1 HG1     HE1      "));
	residue_def.push_back(string("ATOM HA   H1   	 0.000	N  0.0880  0.00  H	0 0 0   !     |   |   |       |        "));
	residue_def.push_back(string("ATOM CB   CT   	 0.358	N  0.0342  1.94  C.3	0 0 0   !  HA-CA--CB--CG--SD--CE--HE3  "));
	residue_def.push_back(string("ATOM HB1  HC   	 0.000	N  0.0241  0.00  H	0 0 0   !     |   |   |       |        "));
	residue_def.push_back(string("ATOM HB2  HC   	 0.000	N  0.0241  0.00  H	0 0 0   !     |   HB2 HG2     HE2      "));
	residue_def.push_back(string("ATOM CG   CT   	 0.358	N  0.0018  1.94  C.3	0 0 0   !   O=C                        "));
	residue_def.push_back(string("ATOM HG1  H1   	 0.000	N  0.0440  0.00  H	0 0 0   !     |                        "));
	residue_def.push_back(string("ATOM HG2  H1   	 0.000	N  0.0440  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM SD   S    	 0.255	N -0.2737  2.09  S.3	0 0 0"));
	residue_def.push_back(string("ATOM CE   CT   	 0.528	N -0.0536  1.94  C.3	0 0 0"));
	residue_def.push_back(string("ATOM HE1  H1   	 0.000	N  0.0684  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HE2  H1   	 0.000	N  0.0684  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HE3  H1   	 0.000	N  0.0684  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM C    C    	-0.030	N  0.5973  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM O    O    	-0.399	A -0.5679  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5679  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI NME  6       0.00  !  N-Methylamide C-terminus 0"));
	residue_def.push_back(string("ATOM N    N    	-0.096	D -0.4157  1.83  N.am	0 1 1    !      |        "));
	residue_def.push_back(string("ATOM H    H    	 0.000	N  0.2719  0.00  H	0 0 0    !      N-HN"));
	residue_def.push_back(string("ATOM CA   CT   	-0.032	N -0.1490  1.94  C.3	0 0 0    !      |        "));
	residue_def.push_back(string("ATOM HA1  H1   	 0.000	N  0.0976  0.00  H	0 0 0    ! HT1-CAT-HT3   "));
	residue_def.push_back(string("ATOM HA2  H1   	 0.000	N  0.0976  0.00  H	0 0 0    !      |        "));
	residue_def.push_back(string("ATOM HA3  H1   	 0.000	N  0.0976  0.00  H	0 0 0    !     HT2       "));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI PHE  21      0.00  !  Phenylalanine"));
	residue_def.push_back(string("ATOM N    N   	-0.096	D -0.4157  1.83  N.am	0 1 1   !     |        HD1  HE1       "));
	residue_def.push_back(string("ATOM H    H   	 0.000	N  0.2719  0.00  H	0 0 0   !  HN-N         |    |        "));
	residue_def.push_back(string("ATOM CA   CT    -0.305	N -0.0024  1.94  C.3	0 0 0   !     |   HB1  CD1--CE1       "));
	residue_def.push_back(string("ATOM HA   H1     0.000	N  0.0978  0.00  H	0 0 0   !     |   |    //     \\      "));
	residue_def.push_back(string("ATOM CB   CT    -0.008	N -0.0343  1.94  C.3	0 0 0   !  HA-CA--CB--CG      CZ--HZ  "));
	residue_def.push_back(string("ATOM HB1  HC     0.000	N  0.0295  0.00  H	0 0 0   !     |   |    \\  __  /       "));
	residue_def.push_back(string("ATOM HB2  HC     0.000	N  0.0295  0.00  H	0 0 0   !     |   HB2  CD2--CE2       "));
	residue_def.push_back(string("ATOM CG   CA     0.296	N  0.0118  1.85  C.ar	0 0 0   !   O=C         |    |        "));
	residue_def.push_back(string("ATOM CD1  CA     0.337	N -0.1256  1.85  C.ar	0 0 0   !     |        HD2  HE2       "));
	residue_def.push_back(string("ATOM HD1  HA     0.000	N  0.1330  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CD2  CA     0.337	N -0.1256  1.85  C.ar	0 0 0"));
	residue_def.push_back(string("ATOM HD2  HA     0.000	N  0.1330  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CE1  CA     0.337	N -0.1704  1.85  C.ar	0 0 0"));
	residue_def.push_back(string("ATOM HE1  HA     0.000	N  0.1430  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CE2  CA     0.337	N -0.1704  1.85  C.ar	0 0 0"));
	residue_def.push_back(string("ATOM HE2  HA     0.000	N  0.1430  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CZ   CA     0.337	N -0.1072  1.85  C.ar	0 0 0"));
	residue_def.push_back(string("ATOM HZ   HA     0.000	N  0.1297  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM C    C    	-0.030	N  0.5973  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM O    O     -0.399	A -0.5679  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5679  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI PRO  15      0.00  !  Proline"));
	residue_def.push_back(string("ATOM N    N      0.078	N -0.2548  1.83  N.am	0 0 0   !       HD1 HD2       "));
	residue_def.push_back(string("ATOM CD   CT    -0.137	N  0.0192  1.94  C.3	0 0 0   !     |   \\ /         "));
	residue_def.push_back(string("ATOM HD1  H1     0.000	N  0.0391  0.00  H	0 0 0   !     N---CD   HG1    "));
	residue_def.push_back(string("ATOM HD2  H1     0.000	N  0.0391  0.00  H	0 0 0   !     |     \\  /      "));
	residue_def.push_back(string("ATOM CA   CT    -0.305	N -0.0266  1.94  C.3	0 0 0   !     |      CG       "));
	residue_def.push_back(string("ATOM HA   H1     0.000	N  0.0641  0.00  H	0 0 0   !     |     /  \\      "));
	residue_def.push_back(string("ATOM CB   CT     0.358	N -0.0070  1.94  C.3	0 0 0   !  HA-CA--CB   HG2    "));
	residue_def.push_back(string("ATOM HB1  HC     0.000	N  0.0253  0.00  H	0 0 0   !     |   / \\         "));
	residue_def.push_back(string("ATOM HB2  HC     0.000	N  0.0253  0.00  H	0 0 0   !     | HB1 HB2       "));
	residue_def.push_back(string("ATOM CG   CT     0.358	N  0.0189  1.94  C.3	0 0 0   !   O=C               "));
	residue_def.push_back(string("ATOM HG1  HC     0.000	N  0.0213  0.00  H	0 0 0   !     |                     "));
	residue_def.push_back(string("ATOM HG2  HC     0.000	N  0.0213  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM C    C     -0.030	N  0.5896  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM O    O     -0.399	A -0.5748  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5748  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI SER  12      0.00  !  Serine"));
	residue_def.push_back(string("ATOM N    N     -0.096	D -0.4157  1.83  N.am	0 1 1   !     |               "));
	residue_def.push_back(string("ATOM H    H      0.000	N  0.2719  0.00  H	0 0 0   !  HN-N               "));
	residue_def.push_back(string("ATOM CA   CT    -0.305	N -0.0249  1.94  C.3	0 0 0   !     |   HB1         "));
	residue_def.push_back(string("ATOM HA   H1     0.000	N  0.0843  0.00  H	0 0 0   !     |   |           "));
	residue_def.push_back(string("ATOM CB   CT    -0.137	N  0.2117  1.94  C.3	0 0 0   !  HA-CA--CB--OG      "));
	residue_def.push_back(string("ATOM HB1  H1     0.000	N  0.0352  0.00  H	0 0 0   !     |   |     \\     "));
	residue_def.push_back(string("ATOM HB2  H1     0.000	N  0.0352  0.00  H	0 0 0   !     |   HB2    HG1  "));
	residue_def.push_back(string("ATOM OG   OH    -0.467	N -0.6546  1.74  O.3	2 1 3   !   O=C               "));
	residue_def.push_back(string("ATOM HG   HO     0.000	N  0.4275  0.00  H	0 0 0   !     |               "));
	residue_def.push_back(string("ATOM C    C     -0.030	N  0.5973  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM O    O     -0.399	A -0.5679  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5679  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI THR  15      0.00  !  Threonine"));
	residue_def.push_back(string("ATOM N    N     -0.096	D -0.4157  1.83  N.am	0 1 1  !     |               "));
	residue_def.push_back(string("ATOM H    H      0.000	N  0.2719  0.00  H	0 0 0  !  HN-N               "));
	residue_def.push_back(string("ATOM CA   CT    -0.305	N -0.0389  1.94  C.3	0 0 0  !     |     OG1--HG1  "));
	residue_def.push_back(string("ATOM HA   H1     0.000	N  0.1007  0.00  H	0 0 0  !     |    /          "));
	residue_def.push_back(string("ATOM CB   CT    -0.205	N  0.3654  1.94  C.3	0 0 0  !  HA-CA--CB-HB       "));
	residue_def.push_back(string("ATOM HB   H1     0.000	N  0.0043  0.00  H	0 0 0  !     |    \\          "));
	residue_def.push_back(string("ATOM OG1  OH    -0.467	N -0.6761  1.74  O.3	2 1 3  !     |     CG2--HG21 "));
	residue_def.push_back(string("ATOM HG1  HO     0.000	N  0.4102  0.00  H	0 0 0  !   O=C    / \\        "));
	residue_def.push_back(string("ATOM CG2  CT     0.528	N -0.2438  1.94  C.3	0 0 0  !     | HG21 HG22     "));
	residue_def.push_back(string("ATOM HG21 HC     0.000	N  0.0642  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HG22 HC     0.000	N  0.0642  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HG23 HC     0.000	N  0.0642  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM C    C     -0.030	N  0.5973  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM O    O     -0.399	A -0.5679  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5679  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI TRP  25      0.00  !  Tryptophan"));
	residue_def.push_back(string("ATOM N    N     -0.096	D -0.4157  1.83  N.am	0 1 1  !     |                  HE3         "));
	residue_def.push_back(string("ATOM H    H      0.000	N  0.2719  0.00  H	0 0 0  !  HN-N                   |          "));
	residue_def.push_back(string("ATOM CA   CT    -0.305	N -0.0275  1.94  C.3	0 0 0  !     |   HB1            CE3         "));
	residue_def.push_back(string("ATOM HA   H1     0.000	N  0.1123  0.00  H	0 0 0  !     |   |             /  \\\\        "));
	residue_def.push_back(string("ATOM CB   CT    -0.008	N -0.0050  1.94  C.3	0 0 0  !  HA-CA--CB---CG-----CD2   CZ3-HZ3  "));
	residue_def.push_back(string("ATOM HB1  HC     0.000	N  0.0339  0.00  H	0 0 0  !     |   |    ||     ||     |       "));
	residue_def.push_back(string("ATOM HB2  HC     0.000	N  0.0339  0.00  H	0 0 0  !     |   HB2  CD1    CE2   CH2-HH2  "));
	residue_def.push_back(string("ATOM CG   C*     0.013	N -0.1415  1.90  C.ar	0 0 0  !   O=C       /   \\   / \\  //        "));
	residue_def.push_back(string("ATOM CD2  CB     0.296	N  0.1243  1.85  C.ar	0 0 0  !     |     HD1    NE1   CZ2         "));
	residue_def.push_back(string("ATOM CE2  CN    -0.151	N  0.1380  1.85  C.ar	0 0 0  !                   |     |          "));
	residue_def.push_back(string("ATOM CE3  CA     0.337	N -0.2387  1.85  C.ar	0 0 0  !                  HE1   HZ2         "));
	residue_def.push_back(string("ATOM HE3  HA     0.000	N  0.1700  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CD1  CW    -0.310	N -0.1638  1.90  C.ar	0 0 0"));
	residue_def.push_back(string("ATOM HD1  H4     0.000	N  0.2062  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM NE1  NA     0.545	D -0.3418  1.87  N.ar	0 1 1"));
	residue_def.push_back(string("ATOM HE1  H      0.000	N  0.3412  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CZ2  CA     0.337	N -0.2601  1.85  C.ar	0 0 0"));
	residue_def.push_back(string("ATOM HZ2  HA     0.000	N  0.1572  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CZ3  CA     0.337	N -0.1972  1.85  C.ar	0 0 0"));
	residue_def.push_back(string("ATOM HZ3  HA     0.000	N  0.1447  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CH2  CA     0.337	N -0.1134  1.85  C.ar	0 0 0"));
	residue_def.push_back(string("ATOM HH2  HA     0.000	N  0.1417  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM C    C     -0.030	N  0.5973  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM O    O     -0.399	A -0.5679  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5679  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI TYR  22      0.00  !  Tyrosine"));
	residue_def.push_back(string("ATOM N    N     -0.096	D -0.4157  1.83  N.am	0 1 1  !     |        HD1  HE1         "));
	residue_def.push_back(string("ATOM H    H      0.000	N  0.2719  0.00  H	0 0 0  !  HN-N         |    |          "));
	residue_def.push_back(string("ATOM CA   CT    -0.305	N -0.0014  1.94  C.3	0 0 0  !     |   HB1  CD1--CE1         "));
	residue_def.push_back(string("ATOM HA   H1     0.000	N  0.0876  0.00  H	0 0 0  !     |   |   //      \\\\        "));
	residue_def.push_back(string("ATOM CB   CT    -0.008	N -0.0152  1.94  C.3	0 0 0  !  HA-CA--CB--CG      CZ--OH    "));
	residue_def.push_back(string("ATOM HB1  HC     0.000	N  0.0295  0.00  H	0 0 0  !     |   |    \\  __  /     \\   "));
	residue_def.push_back(string("ATOM HB2  HC     0.000	N  0.0295  0.00  H	0 0 0  !     |   HB2  CD2--CE2     HH  "));
	residue_def.push_back(string("ATOM CG   CA     0.296	N -0.0011  1.85  C.ar	0 0 0  !   O=C         |    |          "));
	residue_def.push_back(string("ATOM CD1  CA     0.337	N -0.1906  1.85  C.ar	0 0 0  !     |        HD2  HE2         "));
	residue_def.push_back(string("ATOM HD1  HA     0.000	N  0.1699  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CD2  CA     0.337	N -0.1906  1.85  C.ar	0 0 0"));
	residue_def.push_back(string("ATOM HD2  HA     0.000	N  0.1699  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CE1  CA     0.337	N -0.2341  1.85  C.ar	0 0 0"));
	residue_def.push_back(string("ATOM HE1  HA     0.000	N  0.1656  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CE2  CA     0.337	N -0.2341  1.85  C.ar	0 0 0"));
	residue_def.push_back(string("ATOM HE2  HA     0.000	N  0.1656  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CZ   C     -0.151	N  0.3226  1.85  C.ar	0 0 0"));
	residue_def.push_back(string("ATOM OH   OH     0.082	N -0.5579  1.74  O.3	2 1 3"));
	residue_def.push_back(string("ATOM HH   HO     0.000	N  0.3992  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM C    C     -0.030	N  0.5973  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM O    O     -0.399	A -0.5679  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5679  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI VAL  17      0.00  !  Valine"));
	residue_def.push_back(string("ATOM N    N     -0.096	D -0.4157  1.83  N.am	0 1 1  !     |    HG11 HG12   "));
	residue_def.push_back(string("ATOM H    H      0.000	N  0.2719  0.00  H	0 0 0  !  HN-N      | /       "));
	residue_def.push_back(string("ATOM CA   CT    -0.305	N -0.0875  1.94  C.3	0 0 0  !     |     CG1--HG13  "));
	residue_def.push_back(string("ATOM HA   H1     0.000	N  0.0969  0.00  H	0 0 0  !     |    /           "));
	residue_def.push_back(string("ATOM CB   CT     0.127	N  0.2985  1.94  C.3	0 0 0  !  HA-CA--CB-HB        "));
	residue_def.push_back(string("ATOM HB   HC     0.000	N -0.0297  0.00  H	0 0 0  !     |    \\           "));
	residue_def.push_back(string("ATOM CG1  CT     0.528	N -0.3192  1.94  C.3	0 0 0  !     |     CG2--HG21  "));
	residue_def.push_back(string("ATOM HG11 HC     0.000	N  0.0791  0.00  H	0 0 0  !   O=C    / \\         "));
	residue_def.push_back(string("ATOM HG12 HC     0.000	N  0.0791  0.00  H	0 0 0  !     | HG21 HG22      "));
	residue_def.push_back(string("ATOM HG13 HC     0.000	N  0.0791  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM CG2  CT     0.528	N -0.3192  1.94  C.3	0 0 0"));
	residue_def.push_back(string("ATOM HG21 HC     0.000	N  0.0791  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HG22 HC     0.000	N  0.0791  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM HG23 HC     0.000	N  0.0791  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM C    C     -0.030	N  0.5973  1.90  C.2	0 0 0"));
	residue_def.push_back(string("ATOM O    O     -0.399	A -0.5679  1.66  O.2	2 0 2"));
	residue_def.push_back(string("ATOM OXT  O2    -0.399  A -0.5679  1.66  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("RESI HET  22"));
	residue_def.push_back(string("ATOM O    OW    -1.000	N -0.8340  1.77  O.w	2 2 3	"));
	residue_def.push_back(string("ATOM OH2  OW    -1.000	N -0.8340  1.77  O.w	2 2 3	"));
	residue_def.push_back(string("ATOM H1   H  	 0.000	N  0.4170  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM H2   H      0.000	N  0.4170  0.00  H	0 0 0"));
	residue_def.push_back(string("ATOM LI   LI    -2.000	M  1.0000  1.20  Li  	0 0 0"));
	residue_def.push_back(string("ATOM NA   NA    -2.000	M  1.0000  1.20  Na	0 0 0"));
	residue_def.push_back(string("ATOM K    K     -2.000	M  1.0000  1.20  K  	0 0 0"));
	residue_def.push_back(string("ATOM MG   MG    -2.000	M  2.0000  1.90  Mg     6 6 0"));
	residue_def.push_back(string("ATOM CA   CA    -2.000	M  2.0000  1.20  Ca	6 6 0"));
	residue_def.push_back(string("ATOM AL   AL    -2.000	M  3.0000  1.20  Al	6 6 0"));
	residue_def.push_back(string("ATOM MN   MN	-2.000	M  2.0000  1.90  Mn  	6 6 0"));
	residue_def.push_back(string("ATOM FE   FE    -2.000	M  2.0000  1.84  Fe     6 6 0"));
	residue_def.push_back(string("ATOM CO   CO    -2.000  M  2.0000  1.73  Co	6 6 0"));
	residue_def.push_back(string("ATOM NI   NI    -2.000  M  2.0000  1.73  Ni     6 6 0"));
	residue_def.push_back(string("ATOM CU   CU	-2.000  M  2.0000  1.73  Cu     6 6 0"));
	residue_def.push_back(string("ATOM ZN   ZN    -2.000	M  2.0000  1.68  Zn	6 6 0"));
	residue_def.push_back(string("ATOM P	  P	-0.477	N  0.0000  1.80  P.3    0 0 0"));
	residue_def.push_back(string("ATOM S	  SO4	-0.168	N  0.0000  1.80  S.o2   0 0 0"));
	residue_def.push_back(string("ATOM O1	  SO4	-0.399	A -0.5000  1.52  O.co2  2 0 2"));
	residue_def.push_back(string("ATOM O2   SO4   -0.399	A -0.5000  1.52  O.co2  2 0 2"));
	residue_def.push_back(string("ATOM O3   SO4   -0.399	A -0.5000  1.52  O.co2	2 0 2"));
	residue_def.push_back(string("ATOM O4	  SO4   -0.399	A -0.5000  1.52  O.co2	2 0 2"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string("END"));
	residue_def.push_back(string("#"));
	residue_def.push_back(string(""));
	return residue_def;
}


Residue::Residue()
{
	num_atom = 0;
}

Protein::Protein()
{
	num_resi_def = 0;
	
	num_atom = 0;

	min_x = min_y = min_z = max_x = max_y = max_z = 0;
}

void Protein::Read_RESIDUE(const char* dirname)
{
	string filename(dirname);
	filename = filename + "RESIDUE_DEF";
	fstream infile;
	infile.open(filename, ios::in);
	if (!infile.is_open()) 
	{
		Openning_File_Error(filename.c_str());
	}
		
	string line, head, value;
	stringstream sline;
	int num_atom = 0;
	while (getline(infile, line))
	{
		if (line.size() == 0) { continue; }// blank line
		else if (line[0] == '#') { continue; } // comment line
		sline << line;
		sline >> head;
		if (head == "END") break;
		else if (head == "RESI")
		{
			if (!resi_defs.empty())
			{
				int i = resi_defs[value].num_atom;
				if (num_atom != i)
				{
					printf("\nError when reading file: %s. Only %d of %d atoms were successfully read.\n", filename.c_str(), i, num_atom);
					printf("Please check the parameter file.\n");
					exit(1);
				}
			}
			Residue tmpres;
			sline >> value >> num_atom;
			tmpres.name = value;
			resi_defs.insert(map<string, Residue>::value_type(value, tmpres));
		}
		else if (head == "ATOM")
		{
			Patom tmpatom;
			tmpatom.Read_RESIDUEDEF_line(line);
			resi_defs[value].atoms.insert(map<string, Patom>::value_type(tmpatom.name, tmpatom));
			resi_defs[value].num_atom++;
		}
		sline.str(std::string());
		sline.clear();
	}
	infile.close();
	// skip read nucleotide acid parameters
	num_resi_def = int(resi_defs.size());
}

void Protein::Read_RESIDUE()
{
	string head, value;
	stringstream sline;
	int num_atom = 0;
	vector<string> residue_def_lines = get_residue_def();
	for (auto& line : residue_def_lines)
	{
		if (line.size() == 0) { continue; }// blank line
		else if (line[0] == '#') { continue; } // comment line
		sline << line;
		sline >> head;
		if (head == "END") break;
		else if (head == "RESI")
		{
			if (!resi_defs.empty())
			{
				int i = resi_defs[value].num_atom;
				if (num_atom != i)
				{
					cout << "\nError when reading predefined residue parameters.\n";
					exit(1);
				}
			}
			Residue tmpres;
			sline >> value >> num_atom;
			tmpres.name = value;
			resi_defs.insert(map<string, Residue>::value_type(value, tmpres));
		}
		else if (head == "ATOM")
		{
			Patom tmpatom;
			tmpatom.Read_RESIDUEDEF_line(line);
			resi_defs[value].atoms.insert(map<string, Patom>::value_type(tmpatom.name, tmpatom));
			resi_defs[value].num_atom++;
		}
		sline.str(std::string());
		sline.clear();
	}
	// skip read nucleotide acid parameters
	num_resi_def = int(resi_defs.size());
}

void Protein::Show_RESIDUE()
{
	cout << "The number of residue templates is: " << num_resi_def << endl;
	for (auto& iter : resi_defs)
	{
		printf("Residue\t%3s\t%d\n", iter.first.c_str(), iter.second.num_atom);
		for (auto& it : iter.second.atoms)
		{
			printf("ATOM %8s %6.3f %4s %6.3f %6.3f %8s\n",
				it.second.name.c_str(),
				it.second.logp,
				it.second.hb_type.c_str(),
				it.second.q,
				it.second.vdw_r,
				it.second.tripos_type.c_str());
		}
		printf("\n");
	}
}

void Protein::Read_PDB(const char* filename)
{
	extern Parameter* parm;
	string fname(filename);
	string ext = fname.substr(fname.find_last_of('.'));
	vector<string> pdblines;
#ifdef _MSC_VER
	if (ext == ".gz")
	{
		cout << "Reader for .gz file does not implement for this complier. Please use .pdb format." << endl;
		exit(1);
	}
#elif defined(__GNUC__)
	if (ext == ".gz")
	{
		gzFile file = gzopen(fname.c_str(), "r");
		char bufline[100];
		memset(bufline, 0, sizeof(bufline));
		string line;
		while (gzgets(file, bufline, sizeof(bufline)))
		{
			line = bufline;
			pdblines.push_back(line);
		}
		gzclose(file);
	}
#else
	if (ext == ".gz")
	{
		cout << "Reader for .gz file does not implement for this complier. Please use .pdb format." << endl;
		exit(1);
	}
#endif
	// Now read the information of receptor atoms from the PDB file
	else
	{
		if (ext != ".pdb")
		{
			cout << "Warning: this protein file is not a PDB format file. Try parsing it as PDB format ..." << endl;
		}
		ifstream infile;
		infile.open(filename);
		if (!infile.is_open()) { Openning_File_Error(filename); }
		string line;
		while (getline(infile, line))
		{
			if (line.size() == 0) { continue; }// blank line
			pdblines.push_back(line);
		}
		infile.close();
	}
	
	if (pdblines.size() == 0)
	{
		cout << "Error. Cant not read pdb file." << endl;
		exit(1);
	}

	string head;
	int num_model = 0;
	vector<string> tmplines;
	int num_stdatm = 0;
	int valid_model = 0;
	for (auto& line : pdblines)
	{
		if (line.size() == 0) { continue; }
		string head = line.substr(0, 6);
		head.erase(head.find_last_not_of(' ') + 1);
		if (head == "MODEL")
		{
			num_model += 1;
			if (valid_model == 0)
			{
				tmplines.clear();
				num_stdatm = 0;
				tmplines.push_back(line);
			}
		}
		else
		{
			if (valid_model == 0)
			{
				tmplines.push_back(line);
				if (head == "ATOM") { num_stdatm += 1; }
				else if (head == "ENDMDL")
				{
					if (num_stdatm > 0) { valid_model = 1; }
				}
			}
		}
	}
	if (num_model > 1)
	{
		pdblines = tmplines;
		cout << "Warning: this protein file contains more than one model. Only the first model is used." << endl;
	}
	
	cout << "Detecting chain: ";
	
	int terflag = 0, ter = 0;
	int num_hetwater = 0, num_hetmetal = 0;
	int atom_id = 0;
	int unknown = 0;
	set<string> unknown_AA;
	for (auto & line : pdblines)
	{
		if (line.size() == 0) { continue; }// blank line
		head = line.substr(0, 6);
		head.erase(head.find_last_not_of(' ') + 1);
		if (head == "END") break;
		else if (head == "ENDMDL") break;
		else if (head == "TER") { terflag = 1; ter++; }
		else if (head == "ATOM")
		{
			terflag = 0;
			Patom tmpatm;
			tmpatm.Read_PDB_line(line);
			// skip long pairs and hydrogens
			if (tmpatm.is_longpair() || tmpatm.is_hydrogen()) continue;
			if (parm->detect_mode != 0)
			{
				if (tmpatm.coord.x < parm->min_x || tmpatm.coord.x > parm->max_x)continue;
				if (tmpatm.coord.y < parm->min_y || tmpatm.coord.y > parm->max_y)continue;
				if (tmpatm.coord.z < parm->min_z || tmpatm.coord.z > parm->max_z)continue;
			}
			tmpatm.valid = true;
			atoms.push_back(tmpatm);
			if (chain_record.find(tmpatm.chain) != chain_record.end())
			{
				chain_record[tmpatm.chain].push_back(atom_id);
			}
			else
			{
				vector<int> idlist; 
				idlist.push_back(atom_id);
				chain_record.insert(map<char, vector<int>>::value_type(tmpatm.chain, idlist));
				cout << tmpatm.chain;
			}
			atom_id++;
		}
		else if (head == "HETATM") // read water molecules or metal atoms
		{
			if (terflag != 1) 
			{
				unknown = 1;
				Patom tmpatm;
				tmpatm.Read_PDB_line(line);
				unknown_AA.insert(tmpatm.resi_name + "_" + to_string(tmpatm.resi_id) + "_" + to_string(tmpatm.chain));
			}
			if (parm->hetmetal == 0 && parm->hetwater == 0) { continue; }
			
			Patom tmpatm;
			tmpatm.Read_PDB_line(line);
			tmpatm.amber_type = tmpatm.resi_name;
			if (parm->hetwater == 0) // remove all water molecules
			{
				if (tmpatm.name[0] == 'O') { continue; }
			}
			if (parm->hetmetal == 0) // remove all non-water atoms
			{
				if (tmpatm.name[0] != 'O') { continue; }
			}
			// skip Hydrogens
			if (tmpatm.name[0] == 'H') { continue; }
			else if (tmpatm.name[0] == '1' || tmpatm.name[0] == '2' || tmpatm.name[0] == '3') { continue; }
			// boundary check
			if (parm->detect_mode != 0)
			{
				if (tmpatm.coord.x < parm->min_x || tmpatm.coord.x > parm->max_x)continue;
				if (tmpatm.coord.y < parm->min_y || tmpatm.coord.y > parm->max_y)continue;
				if (tmpatm.coord.z < parm->min_z || tmpatm.coord.z > parm->max_z)continue;
			}

			if (parm->hetwater == 1 && (tmpatm.resi_name == "OW" || tmpatm.resi_name == "HOH")) { num_hetwater++; }
			else if (parm->hetmetal == 1)
			{
				if (resi_defs["HET"].atoms.find(tmpatm.name) != resi_defs["HET"].atoms.end())
				{
					if (resi_defs["HET"].atoms[tmpatm.name].amber_type == tmpatm.amber_type &&
						resi_defs["HET"].atoms[tmpatm.name].hb_type == "M")
						num_hetmetal++;
					else { continue; }
				}
				else { continue; }
			}
			else continue;
			tmpatm.resi_name = "HET";
			tmpatm.valid = true;
			atoms.push_back(tmpatm);
			if (chain_record.find(tmpatm.chain) != chain_record.end())
			{
				chain_record[tmpatm.chain].push_back(atom_id);
			}
			else
			{
				vector<int> idlist;
				idlist.push_back(atom_id);
				chain_record.insert(map<char, vector<int>>::value_type(tmpatm.chain, idlist));
			}
			atom_id++;
		}
	}

	num_atom = int(atoms.size());
	if (num_atom == 0) PDB_Format_Error(filename);

	cal_geometric_center();

	if (num_hetmetal > 0) cout << " + HETMETAL";
	if (num_hetwater > 0) cout << " + HETWATER";
	cout << endl;

	// check information
	// 1. check chain length
	int peptide = 0;
	for (map<char, vector<int>>::iterator it = chain_record.begin(); it != chain_record.end(); it++)
	{
		if (it->second.size() < 100) { peptide = 1; break; }
	}
	int occupancyflag = 0;
	for (auto& iter : atoms)
	{
		if (iter.occupancy < 1) { occupancyflag += 1; }
	}
	int warning = peptide + unknown;
	if (occupancyflag > 0) warning += 1;
	if (warning == 0) { cout << "Protein checking: Pass" << endl; }
	else if (warning == 1) { cout << "Protein checking: 1 warning" << endl; }
	else { cout << "Protein checking: " << warning << " warnings." << endl; }
	if (peptide == 1) { cout << "--> Detected small peptide (number of heavy atoms/chain less than 100)" << endl; }
	if (unknown == 1)
	{
		cout << "--> Detected undefined residue: ";
		for (auto& iter : unknown_AA)
		{
			cout << iter << " ";
		}
		cout <<"Skipped." << endl;
	}
	if (occupancyflag >= 1) 
	{
		cout << "--> Detected "<<occupancyflag<<" uncertain atom(s) (occupancy < 1.0)." << endl;
	}

}

void Protein::Protein_Box()
{
	extern Parameter* parm;
	int mark = 0;
	for (auto& it : atoms)
	{
		if (it.valid == 0) { continue; }
		if (mark == 0)
		{
			min_x = (int)floor(it.coord.x); max_x = (int)ceil(it.coord.x);
			min_y = (int)floor(it.coord.y); max_y = (int)ceil(it.coord.y);
			min_z = (int)floor(it.coord.z); max_z = (int)ceil(it.coord.z);
			mark = 1;
		}
		else
		{
			min_x = it.coord.x < min_x ? (int)floor(it.coord.x) : min_x;
			max_x = it.coord.x > max_x ? (int)ceil(it.coord.x) : max_x;
			min_y = it.coord.y < min_y ? (int)floor(it.coord.y) : min_y;
			max_y = it.coord.y > max_y ? (int)ceil(it.coord.y) : max_y;
			min_z = it.coord.z < min_z ? (int)floor(it.coord.z) : min_z;
			max_z = it.coord.z > max_z ? (int)ceil(it.coord.z) : max_z;
		}
	}
	if (parm->info == 1)
	{
		printf("**INFO: The area of protein is:\n");
		printf("Max_x = %d\tMin_x = %d\n", max_x, min_x);
		printf("Max_y = %d\tMin_y = %d\n", max_y, min_y);
		printf("Max_z = %d\tMin_z = %d\n", max_z, min_z);
	}
	parm->min_x = min_x;
	parm->max_x = max_x;
	parm->min_y = min_y;
	parm->max_y = max_y;
	parm->min_z = min_z;
	parm->max_z = max_z;
	return;
}

void Protein::Value_Atom()
{
	for (vector<Patom>::iterator it = atoms.begin(); it != atoms.end(); it++)
	{
		// First, check the atom to see whether it has a legal residue
		if (resi_defs.find(it->resi_name) == resi_defs.end())
		{
			cout << "Warning: cannot find parameter for the protein atom " << it->atom_id
				<< " " << it->name << " " << it->resi_name << endl;
			cout << "This atom is skipped." << endl;
			it->valid = false; 
			continue;
		}
		// Second, check the atom in the specific residue
		Residue& res = resi_defs[it->resi_name];
		if (res.atoms.find(it->name) != res.atoms.end())
		{
			Patom& resatm = res.atoms[it->name];
			it->tripos_type = resatm.tripos_type;
			it->vdw_r = resatm.vdw_r;
			it->q = resatm.q;
			it->logp = resatm.logp;
			it->max_anum = resatm.max_anum;
			it->max_dnum = resatm.max_dnum;
			it->hblevel = resatm.hblevel;
			it->hb_type = resatm.hb_type;
			it->valid = true;
		}
		else
		{
			cout << "Warning: cannot find parameter for the protein atom " << it->atom_id << " " << it->name << " " << it->resi_name << endl;
			cout << "This atom is skipped." << endl;
			it->valid = false;
			continue;
		}
	}
}

void Protein::Cal_HB_Root()
{
	// check CYS-CYS
	vector<int> tmp_index;
	double d;
	for (int i = 0; i < num_atom; i++) 
	{
		atoms[i].num_neib = 0;
		if (atoms[i].valid <= 0) continue;
		if (atoms[i].amber_type != "S.3") continue;
		else tmp_index.push_back(i);
	}
	for (int i = 0; i < tmp_index.size(); i++)
	{
		for (int j = i + 1; j < tmp_index.size(); j++)
		{
			d = Distance2(atoms[i].coord, atoms[j].coord);
			if (d < 6.25)
			{
				atoms[i].max_dnum = 0;
				atoms[j].max_dnum = 0;
				atoms[i].hb_type = "N";
				atoms[j].hb_type = "N";
				break;
			}
		}
	}
	// get neighbor
	Patom* atm_i, * atm_j;
	for (int i = 0; i < num_atom; i++)
	{
		atm_i = &atoms[i];
		if (atm_i->valid <= 0) continue;
		atm_i->root.reset();
		atm_i->plane.reset();
		if (atm_i->hb_type == "M") continue; // metal
		else if (atm_i->amber_type == "O.w") continue; // water
		
		int stop_i = i + 51;
		if (stop_i > num_atom) stop_i = num_atom;
		for (int j = i + 1; j < stop_i; j++)
		{
			atm_j = &atoms[j];
			if (atm_j->valid <= 0)continue;
			if (atm_j->hb_type == "M") continue; // metal
			d = Distance2(atm_i->coord, atm_j->coord);
			if (d > 4) continue;
			else
			{
				atm_i->neib[atm_i->num_neib] = j; atm_i->num_neib++;
				atm_j->neib[atm_j->num_neib] = i; atm_j->num_neib++;
			}
		}
	}
	// check ADpair 
}

void Protein::cal_geometric_center()
{
	if (num_atom == 0) return;
	for (int i = 0; i < num_atom; i++)
	{
		center_coord += atoms[i].coord;
	}
	center_coord /= num_atom;
}

