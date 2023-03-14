#include <cmath>
#include <algorithm>
#include "pocket.h"
#include "grid.h"
#include "protein.h"
#include "ligand.h"
#include "parameter.h"

using namespace std;

Cavity::Cavity()
{
	valid = 0;
	id = 0;
	re = 0;
	surfacenum = 0;
	vacant_num = 0;
	index = 0;
	rank = 0;
	edge = 0;
	area = 0;
	dep1 = 0;
	dep2 = 0;
	min_x = min_y = min_z = 99999;
	max_x = max_y = max_z = -99999;
	atom_min_x = atom_min_y = atom_min_z = 99999;
	atom_max_x = atom_max_y = atom_max_z = -99999;
	ab_vector_norm = 0;
	core_vector_norm = 0;
	hydrophobic_surface = acceptor_surface = donor_surface = bind_surface = 0;
	hydrophobic_vacant = acceptor_vacant = donor_vacant = bind_vacant = 0;
	hydrophobic_s = acceptor_s = donor_s = 0;
	hydrophobic_v = acceptor_v = donor_v = 0;
	surface_dep = vacant_dep = 0;
	drugscore = 0;
	surface = 0;
	for (int i = 0; i < 1000; i++)step_vacant[i] = 0;
	// for ligand mod
	percent = 0;
	ligall = 0;
	ligmatch = 0;
}

bool comp(Cavity& c1, Cavity& c2)
{
	return c1.rank > c2.rank;
}

bGrid::bGrid(int ii, int jj, int kk, int c)
{
	i = ii, j = jj, k = kk;
	cavity = c;
}

Pocket::Pocket()
{
	shift_x = shift_y = shift_z = 0;
	cavity_valid = 0;
	totalcavity = 0;
	new_total = 0;
	cavity = nullptr;
	fullradius = 0;
	length = 0;
	fullball = nullptr;
	border_min_x = -1;
	border_max_x = -1;
	border_min_y = nullptr;
	border_max_y = nullptr;
	border_min_z = nullptr;
	border_max_z = nullptr;
}

void Pocket::Cavity_Search() // main function for cavity search
{
	extern Parameter* parm;
	// initialize grid box
	grid.Initialize(3, 0.5, parm->min_x, parm->max_x, parm->min_y, parm->max_y, parm->min_z, parm->max_z);
	fgrid.Initialize((int)(parm->radius_length + parm->radius_mis + 3) + 1, 0.5, parm->min_x, parm->max_x, parm->min_y, parm->max_y, parm->min_z, parm->max_z);
	shift_x = (fgrid.dx - grid.dx) / 2;
	shift_y = (fgrid.dy - grid.dy) / 2;
	shift_z = (fgrid.dz - grid.dz) / 2;

	// first, check if the grid bumps to any atom and initialize the grids
	Initialize();

	// global eraser
	Erase();

	if (parm->info == 1) cout << "**INFO: Now detecting valid matrix ...\n";

	Define_Surface(); // find the surface grid of protein surface
	Define_Vacant(); // find the vacant grid
	
	if (parm->detect_mode == 1)
	{
		Define_Ligand_Grid();
	}
	// Get Real Depth
	if (parm->info == 1) cout << "**INFO: Now detecting real depth ...\n";
	Reedge();
	Get_Real_Depth(0);

	// Check and GO
	int vnum_size = (grid.dx * grid.dy * grid.dz) / parm->chip_v + 1;
	int* vnum = new int[vnum_size];
	for (int i = 0; i < vnum_size; i++) vnum[i] = 0;
	vector<int> tmpvnum;
	if (parm->judge == 0) totalcavity = Check(0, 0);
	else
	{
		totalcavity = Detect_V_Cluster(0, 0, vnum, tmpvnum);
	}

	if (parm->debug) Output_Middle_Vacant(1, totalcavity+1);

	new_total = totalcavity * 300;
	if (new_total > 300000) new_total = 300000;
	int* re = new int[new_total * 10]; // re[i]用来记录该口袋缩了几层
	re[0] = vnum[0] = new_total * 10;
	for (int i = 1; i < re[0]; i++) { re[i] = -1; }

	cavity = new Cavity[new_total];
	if (cavity == NULL) Memory_Allocation_Error();
	// initialize all the cavity
	for (int i = 0; i < new_total; i++)
	{
		cavity[i].id = i;
	}

	// Separate and Choose
	int nowtotal;
	nowtotal = totalcavity;
	new_total = totalcavity;
	int tmp;

	if (parm->judge == 0)
	{
		printf("Please use depth or d-v mode\n");
		exit(0);
	}
	else
	{
		int predict_depth = 0, predict_vacant = 0;
		vector<Position3D> tmpvacantset;
		int softsep_mark = 0, rescue_mark = 0;
		re[0] = 1;
		// 该循环中，对Check_V()检测出的口袋不停的收缩，然后重新检测口袋。
		// 其中，判断口袋是否要继续收缩的首要条件是：
		// 1. predict_depth>parm->max_abstract_depth(default=20),即口袋中层格点数>max_depth_vacant的层数；
		// 2. re[i] < parm->min_depth(default=8), 即该口袋还没有收缩够8次；
		// 3. parm->judge=2且vnum[i]>parm->max_abstract_limit_v (1500), 即口袋中格点数>max_abstract_limit_v;
		// 4. parm->judge=2且predict_vacant > parm->max_limit_v (6000)
		// 在以上条件有一个满足的情况下，需要：
		// 1. predict_depth > min_abstract_depth (2)
		// 2. 收缩后口袋格点数不少于chip_v。这样的口袋
		// 
		for (int i = 1; i <= nowtotal; i++)
		{
			if (re[i] == -1) re[i] = 0;
			predict_depth = Get_Depth_Predict(i);
			if (predict_depth > parm->max_abstract_depth || re[i] < parm->min_depth) predict_vacant = -1;
			else if (parm->judge == 2)
			{
				predict_vacant = Get_Predict_V(i, re[i]);
			}
			if (predict_depth > parm->max_abstract_depth
				|| re[i] < parm->min_depth
				|| (parm->judge == 2 && parm->max_abstract_limit_v < vnum[i])
				|| (parm->judge == 2 && predict_vacant > parm->max_limit_v))
			{
				/*if (predict_depth > parm->max_abstract_depth || re[i] < parm->min_depth)
				{
					if (predict_depth == 0)continue;
				}*/
				if (parm->judge == 2)
				{
					if (predict_depth <= parm->min_abstract_depth)continue;
				}
				if (re[i] != 0 && re[i] % 10 == 0)
				{
					Get_Real_Depth(i);
				}
				vnum[0] = vnum[i];
				softsep_mark = 0; rescue_mark = 0;
				if (parm->soft_separate || parm->rescue)
				{
					if (re[i] < parm->min_depth) { predict_vacant = Get_Predict_V(i, re[i]); }
					if (vnum[i] < parm->softsep_limit)
					{
						softsep_mark = 1;
					}
					if (vnum[i] < parm->rescue_limit && predict_vacant > parm->max_limit_v)
					{
						rescue_mark = 1;
					}
				}
				if (softsep_mark == 1 || rescue_mark == 1)
				{
					tmpvacantset.clear();
					tmpvacantset = vacantset[i];
					Backup_Adjust(i);
				}

				Adjust(i, 0, -1, re); // 将口袋i的V表面收缩一层
				vnum[i] = (int)vacantset[i].size();
				if (vnum[i] == 0)
				{
					re[i] = -1;
					continue;
				}
				else if (vnum[i] == vnum[0])
				{
					re[i] = 0;
					continue;
				}
				tmp = Detect_V_Cluster(nowtotal, i, vnum, tmpvnum); // 重新检查原来口袋i中现在的口袋情况，新口袋从nowtotal开始编号
				if (softsep_mark == 1) // shrink和separation过程会持续进行，但当口袋已经比较小时，separation已经没有必要，可能会导致最终的口袋缺失一小块，这里会检测这种情况并提前停止收缩
				{
					if ((tmp - nowtotal) >= 1)
					{
						if (tmpvnum.size() != 1)
						{
							int countnum = 0;
							for (auto& a : tmpvnum) { if (a > parm->max_depth_vacant) countnum++; }
							if (countnum > 1)
							{
								if (parm->info) { cout << "Soft separation is specified and ignore the separation of cavity " << i << endl; }
								vacantset[i] = tmpvacantset;
								vnum[i] = (int)vacantset[i].size();
								for (int t = nowtotal + 1; t <= tmp; t++) { vacantset.erase(t); }
								Restore_Adjust(i);
								continue;
							}
						}
					}
					else if (tmp == nowtotal)
					{
						if (predict_vacant / re[i] > 500 && re[i] > 2) // 较平坦的口袋比较浅，容易在收缩后被舍弃，尝试召回这类口袋。同时，平坦口袋受擦除球影响较大，扩展时多扩展一次
						{
							if (parm->info) { cout << "Soft separation is specified and ignore the separation of cavity " << i << endl; }
							vacantset[i] = tmpvacantset;
							vnum[i] = (int)vacantset[i].size();
							re[i] = re[i] + 1;
							Restore_Adjust(i);
							continue;
						}
					}
				}
				if (rescue_mark == 1) 
				// 一些大口袋，收缩到最后，可能会出现一些特殊情况，例如:
				// 1. predict_depth>2，但再收缩一次后，vnum<300，这样会导致该口袋被舍弃；
				// 2. 相反，也会出现vnum>300但predict_depth=0的情况，predict_depth=0会导致口袋在之后仍然被丢弃，该类口袋会在之后召回，这里不处理;
				// 3. 当收缩面与蛋白表面平行时，最后几层收缩后口袋可能被分成碎块而遭到舍弃；
				// 这里将尝试召回这类会被舍弃的口袋
				{
					if (tmp == nowtotal)
					{
						if (parm->info) { cout << "Rescue cavity " << i << " whose potential size is larger than parm->max_limit_v (=" << parm->max_limit_v << ")" << endl; }
						vacantset[i] = tmpvacantset;
						vnum[i] = (int)vacantset[i].size();
						Restore_Adjust(i);
						continue;
					}
				}
				if (parm->debug) Output_Middle_Vacant(nowtotal + 1, tmp + 1);
				for (int t = nowtotal + 1; t <= tmp; t++)
				{
					re[t] = re[i] + 1;
				}
				re[i] = -1;
				nowtotal = tmp;
			}
		}
		int isvalid = true;
		for (int i = 1; i <= nowtotal; i++)
		{
			isvalid = true;
			if (re[i] < 0 || vnum[i] < parm->chip_v || i >= vnum_size) { isvalid = false; }
			else
			{
				predict_depth = Get_Depth_Predict(i);
				if (parm->rescue && predict_depth < parm->chip_depth) // 之前提到，大口袋缩小到最后时，会出现predict_depth=0但vnum>300的情况，这里尝试召回该类口袋
				{
					predict_vacant = Get_Predict_V(i, re[i]);
					if (predict_vacant > parm->max_limit_v)
					{
						predict_depth = parm->chip_depth;
					}
				}
				if (predict_depth >= parm->chip_depth)
				{
					new_total++;
					cavity[new_total].id = i;
					cavity[new_total].re = re[i];
					Get_Real_Depth(i);
					for (auto& vit : vacantset[i])
					{
						if (grid(vit).re == -1) grid(vit).re = re[i] + grid(vit).dep;
					}
					isvalid = true;
				}
				else { isvalid = false; }
			}
			if (!isvalid)
			{
				if (vacantset.count(i) != 0)
				{
					for (auto& vit : vacantset[i])
					{
						grid(vit).cavity = -1;
						grid(vit).type = 'O';
					}
				}
				for (auto& sit : surfacegrid)
				{
					if (grid(sit).cavity == i) grid(sit).cavity = -1;
				}
			}
		}
		printf("**INFO: Now processing cavities ...\n");
		Backup(1);  // backup 'V' grid
		int flag = 0;
		int k = 1;

		for (int i = totalcavity + 1; i <= new_total; i++)
		{
			Restore(cavity[i].id, 1); // 所有格点的inner恢复为inner_bak的值，且令口袋cavity[i].id格点以外的所有格点的cavity=-1, 所有'V'格点的type='O'.
			re[0] = 1;
			for (int j = 0; j < re[cavity[i].id]; j++)
			{
				Adjust(cavity[i].id, 0, 1, re); // 将口袋cavity[i].id的表面往外扩
				if (parm->judge == 2 && parm->rigid_limit == 1)
				{
					vnum[cavity[i].id] = (int)vacantset[cavity[i].id].size();
					if (vnum[cavity[i].id] > parm->max_limit_v) break;
				}
			}
			re[0] = 0;
			if (parm->edge_adjust != 0 || parm->vacant_adjust != 0)
			{
				for (int j = 0; j < abs(parm->edge_adjust); j++)
				{
					if (parm->edge_adjust > 0) Adjust(cavity[i].id, 1, 0, re);
					else Adjust(cavity[i].id, -1, 0, re);
				}
				for (int j = 0; j < abs(parm->vacant_adjust); j++)
				{
					if (parm->vacant_adjust > 0) Adjust(cavity[i].id, 0, 1, re);
					else Adjust(cavity[i].id, 0, -1, re);
				}
			}
			vnum[cavity[i].id] = (int)vacantset[cavity[i].id].size();
			cavity[i].surfacenum = 0;
			cavity[i].vacant_num = vnum[cavity[i].id];
			cavity[i].valid = Choose(i, 2);
			if (cavity[i].valid == 1) cavity_valid++;
			else continue;

			// get surface
			Cavity_Info(i, 9); // 将口袋cavity[i].id的surface格点的cavity值标记为cavity[i].id
			Cavity_Info(i, 0); // count surface， 修改口袋的totalnum的值
			Cavity_Info(i, 10);  // 找到口袋的盒子的边界，并将totalnum和vacant_num的值加到cavity[0]上
			Reedge(); // 重新定义边界
			Cavity_Info(i, 3); // get dep1 = edge/surface dep2 = area/vacant
			Cavity_Info(i, 7); // get cavity core 把口袋缩了一圈，然后找到内表面和外表面的中心点，然后算ab_vector, core_vector
			Cavity_Info(i, 5); // Mark bind surface and vacant 统计hydrophobic,donor, acceptor格点个数
			Cavity_Info(i, 6); // Get Surface & Vacant Dep

			cavity[i].index = 0;
			if (cavity[i].vacant_num - 2 * cavity[i].area + cavity[i].hydrophobic_v <= 0
				|| cavity[i].surface + cavity[i].edge - (cavity[i].donor_s + cavity[i].acceptor_s) <= 0)
			{
				cavity[i].rank = 0;
			}
			else
			{
				cavity[i].rank = ((float)(cavity[i].vacant_num - 2 * cavity[i].area + cavity[i].hydrophobic_v)) / ((float)(cavity[i].surface + cavity[i].edge - (cavity[i].donor_s + cavity[i].acceptor_s)));
			}
			if (cavity[i].rank < parm->output_rank)
			{
				cavity[i].valid = 0; 
				flag++;
			}
			if (cavity[i].valid == 0)continue;
			cavity[i].drugscore = (float)(cavity[i].bind_vacant * cavity[i].hydrophobic_vacant / cavity[i].vacant_num - 6000 * cavity[i].area / cavity[i].vacant_num);
			
			if (parm->detect_mode == 1)
			{
				Calculate_Ligand_Percent(i);
			}
			Prepare_Output(i, k);
			k++;
		}
		cavity_valid -= flag;
		if (flag != 0)
		{
			printf("Rank score limit skipped %d cavities...\n", flag);
		}
		printf("Total valid cavity num is %d\n", cavity_valid);
		if (cavity_valid == 0) printf("Warning: No cavity detected in protein, please change search ruler and try again!\n");
		else printf("Cavity will output % d cavity file(s)\n", cavity_valid);
		vector<Cavity> validcavity;
		for (int i = totalcavity + 1; i <= new_total; i++)
		{
			if (cavity[i].valid != 1) continue;
			validcavity.push_back(cavity[i]);
		}
		sort(validcavity.begin(), validcavity.end(), comp);
		int j = 0;
		for (int i = 0; i < validcavity.size(); i++)
		{
			validcavity[i].index = i + 1;
			for (j = totalcavity + 1; j <= new_total; j++)
			{
				if (cavity[j].id == validcavity[i].id) break;
			}
			Output_All(j, i + 1);
		}
	}
	delete[] re;
	delete[] vnum;
	delete[] cavity;
	return;
}

void Pocket::Initialize()
// 该函数的任务是：
// 1. 将蛋白原子占据的格点的type标记为'A'；
// 2. 将每个蛋白原子周围的空格点，根据到蛋白原子的距离，标记其neib
// 3. 计算每个原子周围grid空格点的donor_score, acceptor_score和hydrophobic_score
// 4. 统计grid格点的各项分数，并标记adh
// 5. 根据neib对每个grid的type进行分类
// 6. 根据neib，标记fgrid中的'X'
// 
// 经过该函数后，
// 1. 找到了fgrid中的所有'X'格点，即between (atom vdw_radius + probe radius) ~ (atom vdw_radius + probe radius + probe mis)的格点。'X'格点将用于global eraser过程
// 2. 定义了grid格点的三类分数，该分数会被用于最后计算分数；
// 3. 标记了grid中的格点类型：
//	E: 原子占据的格点,以及溶剂可及表面空间中经过Recover()后剩余的格点
//  V: 原子周围(atom vdw_radius + probe radius + probe mis)以内的空格点
//  O: 原子周围(atom vdw_radius + probe radius + probe mis)以外的空格点
// 
{
	extern Parameter* parm;
	extern Protein* prot;

	float id, jd, kd;
	int it = 0, jt = 0, kt = 0;
	int fit = 0, fjt = 0, fkt = 0;
	int neib_scale = 0;
	int neib_min_x, neib_max_x, neib_min_y, neib_max_y, neib_min_z, neib_max_z;
	neib_min_x = neib_max_x = neib_min_y = neib_max_y = neib_min_z = neib_max_z = 0;
	float d = 0, d_tmp = 0, d0 = 0, d1 = 0, d2 = 0;

	float r_N = HB_D_R; // the probe atom is N.4 (donor)
	float r_O = HB_A_R; // the probe atom is O.co2 (acceptor)
	float r_C = HYDROPHOBIC_R; // the probe atom is C.3 (hydrophobic)

	Point3D v1, v2;
	float angle;

	//printf("**WARNING: pocket.cpp:Pocket:Initialize(): When calculating donor_score and acceptor_score for each grid, the atom.pos was not correctly defined.\n");
	float x, y, z, cutoff;
	float c2_2, c2_3, c2_4;

	// define the global size 
	cutoff = (float)fgrid.extend_length * 2;
	neib_scale = fgrid.extend_length * 2;
	if (cutoff > neib_scale) { cout << "Error." << endl; }
	int*** ball = nullptr;
	int border_x_min = -1, border_x_max = -1;
	int* border_y_min = nullptr, * border_y_max = nullptr;
	int** border_z_min = nullptr, ** border_z_max = nullptr;
	Create_Ball_Grid(neib_scale, cutoff, ball, border_x_min, border_x_max, border_y_min, border_y_max, border_z_min, border_z_max);

	int count_0, count_1, count_2, count_3, count_4;
	count_0 = count_1 = count_2 = count_3 = count_4 = 0;
	int start_x, start_y, start_z, end_x, end_y, end_z;
	for (auto& iter : prot->atoms)
	{
		if (iter.valid == false) { continue; }
		// grid coordinates of the current atom
		id = (float)((int)(iter.coord.x * 2)) / 2;
		jd = (float)((int)(iter.coord.y * 2)) / 2;
		kd = (float)((int)(iter.coord.z * 2)) / 2;
		it = (int)((id - grid.min_x) * 2);
		jt = (int)((jd - grid.min_y) * 2);
		kt = (int)((kd - grid.min_z) * 2);
		grid.grid[it][jt][kt].type = 'A';
		grid.grid[it][jt][kt].neib = 10000;
		fit = it + shift_x; fjt = jt + shift_y; fkt = kt + shift_z;
		fgrid.grid[fit][fjt][fkt].type = 'A';
		fgrid.grid[fit][fjt][fkt].neib = 10000;

		int small_cutoff = 6;
		neib_min_x = max(it - small_cutoff * 2 + 1, 0);
		neib_max_x = min(it + small_cutoff * 2 + 1, grid.dx - 1);
		neib_min_y = max(jt - small_cutoff * 2 + 1, 0);
		neib_max_y = min(jt + small_cutoff * 2 + 1, grid.dy - 1);
		neib_min_z = max(kt - small_cutoff * 2 + 1, 0);
		neib_max_z = min(kt + small_cutoff * 2 + 1, grid.dz - 1);
		for (int ii = neib_min_x; ii <= neib_max_x; ii++)
		{
			x = iter.coord.x - grid.grid[ii][0][0].coord.x;
			if (x > small_cutoff) continue;
			else if (x < -small_cutoff) continue;
			for (int jj = neib_min_y; jj <= neib_max_y; jj++)
			{
				y = iter.coord.y - grid.grid[ii][jj][0].coord.y;
				if (y > small_cutoff) continue;
				else if (y < -small_cutoff) continue;
				for (int kk = neib_min_z; kk <= neib_max_z; kk++)
				{
					// first, find out vacant grids and excluded grids
					z = iter.coord.z - grid.grid[ii][jj][kk].coord.z;
					if (z > small_cutoff) continue;
					else if (z < -small_cutoff) continue;
					d = (float)(sqrt(x * x + y * y + z * z));
					//if (d < iter.vdw_r) { grid.grid[ii][jj][kk].neib = 10000; }
					//else if (grid.grid[ii][jj][kk].neib < 9999 && d < iter.vdw_r + parm->SA_radius) { grid.grid[ii][jj][kk].neib = 1000; }
					if (d < iter.vdw_r + parm->SA_radius) { grid.grid[ii][jj][kk].neib = 1000; }
					else if (grid.grid[ii][jj][kk].neib < 999 && d < iter.vdw_r + parm->radius_length) { grid.grid[ii][jj][kk].neib = 100; }
					else if (grid.grid[ii][jj][kk].neib < 99 && d < iter.vdw_r + parm->radius_length + parm->radius_mis) { grid.grid[ii][jj][kk].neib = 10; }
					fgrid.grid[ii + shift_x][jj + shift_y][kk + shift_z].neib = grid.grid[ii][jj][kk].neib;
					if (grid.grid[ii][jj][kk].type == 'A')continue;
					if (d > small_cutoff) continue;
					//second, calculate donor score for each grid
					d0 = r_N + iter.vdw_r;
					if ((d - d0) < -0.6) { grid.grid[ii][jj][kk].donor_score += (VB); }
					if (d < 3.5 && iter.q < 0) { grid.grid[ii][jj][kk].donor_score += (SB_AWARD); }
					if (iter.max_anum > 0 && iter.hb_type != "M")
					{
						if (d < d0) // maximum h-bond length
						{
							v1 = iter.root - iter.coord;
							v2 = iter.coord - grid.grid[ii][jj][kk].coord;
							angle = 180 - Angle_Of_Two_Points(v1, v2);
							angle = fabs(angle - iter.hangle);
							if (angle <= OLD_HB_ANGLE_CUTOFF) // h-bond angle cutoff
							{
								// the probe may form HB pair with nearby pocket atoms, only the shortest HBs are taken into account
								if (d - d0 < grid.grid[ii][jj][kk].donor_overlap)
								{
									if (iter.tripos_type == "O.w") grid.grid[ii][jj][kk].donor_wtype = 2;
									else grid.grid[ii][jj][kk].donor_wtype = 1;
									grid.grid[ii][jj][kk].donor_overlap = d - d0;
								}
							}
						}
					}
					// third, calculate the acceptor score for each grid
					d1 = r_O + iter.vdw_r;
					if ((d - d1) < -0.6) { grid.grid[ii][jj][kk].acceptor_score += (VB); }
					if (d < 3.5 && iter.q >0) { grid.grid[ii][jj][kk].acceptor_score += (SB_AWARD); }
					if (iter.max_dnum > 0)
					{
						if (iter.hb_type != "M")
						{
							if (d < d1)
							{
								v1 = iter.root - iter.coord;
								v2 = iter.coord - grid.grid[ii][jj][kk].coord;
								angle = 180.0 - Angle_Of_Two_Points(v1, v2);
								angle = fabs(angle - iter.hangle);
								if (angle <= OLD_HB_ANGLE_CUTOFF)
								{
									if (d - d1 < grid.grid[ii][jj][kk].acceptor_overlap)
									{
										if (iter.tripos_type == "O.w")grid.grid[ii][jj][kk].acceptor_wtype = 2;
										else grid.grid[ii][jj][kk].acceptor_wtype = 1;
										grid.grid[ii][jj][kk].acceptor_overlap = d - d1;
									}
								}
							}
						}
						else if (iter.hb_type == "M")
						{
							if (d < 2.0)grid.grid[ii][jj][kk].acceptor_score += (MB);
							else if (d < 3.0) grid.grid[ii][jj][kk].acceptor_score += ((3.0 - d) * (MB));
						}
					}
					// fourth, calculate the hydrophobic score for each grid
					d2 = r_C + iter.vdw_r;
					if ((d - d2) < -0.6)grid.grid[ii][jj][kk].hydrophobic_score += (VB);
					if (d < 5.0)grid.grid[ii][jj][kk].logp += iter.logp;
				}
			}
		}

		// traverse every grid around the atom 
		cutoff = iter.vdw_r + parm->radius_length + parm->radius_mis;
		//c2_1 = iter.vdw_r * iter.vdw_r;
		c2_2 = (iter.vdw_r + parm->SA_radius) * (iter.vdw_r + parm->SA_radius);
		c2_3 = (iter.vdw_r + parm->radius_length) * (iter.vdw_r + parm->radius_length);
		c2_4 = cutoff * cutoff;
		start_x = fit + border_x_min - neib_scale;
		end_x = fit + border_x_max - neib_scale;
		for (int ii = start_x; ii < end_x; ii++)
		{
			x = iter.coord.x - fgrid.grid[ii][0][0].coord.x;
			if (x > cutoff) { continue; }
			else if (x < -cutoff) { continue; }
			start_y = fjt + border_y_min[ii + neib_scale - fit] - neib_scale;
			end_y = fjt + border_y_max[ii + neib_scale - fit] - neib_scale;
			for (int jj = start_y; jj < end_y; jj++)
			{
				y = iter.coord.y - fgrid.grid[ii][jj][0].coord.y;
				if (y > cutoff) { continue; }
				else if (y < -cutoff) { continue; }
				d_tmp = x * x + y * y;
				if (d_tmp > c2_4) { continue; }
				start_z = fkt + border_z_min[ii + neib_scale - fit][jj + neib_scale - fjt] - neib_scale;
				end_z = fkt + border_z_max[ii + neib_scale - fit][jj + neib_scale - fjt] - neib_scale;
				for (int kk = start_z; kk < end_z; kk++)
				{
					//if (fgrid.grid[ii][jj][kk].neib > 100) { continue; }
					z = iter.coord.z - fgrid.grid[ii][jj][kk].coord.z;
					if (z > cutoff) { continue; }
					else if (z < -cutoff) continue;
					d = d_tmp + z * z;
					//if (d < c2_1) { fgrid.grid[ii][jj][kk].neib = 10000; }
					//else if (d < c2_2 && fgrid.grid[ii][jj][kk].neib < 9999) { fgrid.grid[ii][jj][kk].neib = 1000; }
					if (d < c2_2) { fgrid.grid[ii][jj][kk].neib = 1000; }
					else if (d < c2_3 && fgrid.grid[ii][jj][kk].neib < 999) { fgrid.grid[ii][jj][kk].neib = 100; }
					else if (d < c2_4 && fgrid.grid[ii][jj][kk].neib < 99) { fgrid.grid[ii][jj][kk].neib = 10; }
					/* neib/type:
							10000/A: it is the atom
							10000/V: in atom vdw_radius
							1000/V:  between atom vdw_radius ~ (atom vdw_radius + water_access)
							100/V:   between (atom vdw_radius + water_access) ~ (atom vdw_radius + probe radius)
							10/V:    between (atom vdw_radius + probe radius) ~ (atom vdw_radius + probe radius + probe mis)
							0/V:	 others
					*/
				}
			}
		}
	}

	for (int i = 0; i < grid.dx; i++)
		for (int j = 0; j < grid.dy; j++)
			for (int k = 0; k < grid.dz; k++)
			{
				if (grid.grid[i][j][k].donor_wtype == 1)
				{
					if (grid.grid[i][j][k].donor_overlap < -0.6) grid.grid[i][j][k].donor_score += SHB;
					else if (grid.grid[i][j][k].donor_overlap < -0.3) grid.grid[i][j][k].donor_score += MHB;
					else if (grid.grid[i][j][k].donor_overlap < 0.0) grid.grid[i][j][k].donor_score += WHB;
				}
				else if (grid.grid[i][j][k].donor_wtype == 2)
				{
					if (grid.grid[i][j][k].donor_overlap < -0.6) grid.grid[i][j][k].donor_score += SWH;
					else if (grid.grid[i][j][k].donor_overlap < -0.3) grid.grid[i][j][k].donor_score += MWH;
					else if (grid.grid[i][j][k].donor_overlap < 0.0) grid.grid[i][j][k].donor_score += WWH;
				}
				if (grid.grid[i][j][k].acceptor_wtype == 1)
				{
					if (grid.grid[i][j][k].acceptor_overlap < -0.6) grid.grid[i][j][k].acceptor_score += SHB;
					else if (grid.grid[i][j][k].acceptor_overlap < -0.3) grid.grid[i][j][k].acceptor_score += MHB;
					else if (grid.grid[i][j][k].acceptor_overlap < 0.0) grid.grid[i][j][k].acceptor_score += WHB;
				}
				else if (grid.grid[i][j][k].acceptor_wtype == 2)
				{
					if (grid.grid[i][j][k].acceptor_overlap < -0.6) grid.grid[i][j][k].acceptor_score += SWH;
					else if (grid.grid[i][j][k].acceptor_overlap < -0.3) grid.grid[i][j][k].acceptor_score += MWH;
					else if (grid.grid[i][j][k].acceptor_overlap < 0.0) grid.grid[i][j][k].acceptor_score += WWH;
				}
				if (grid.grid[i][j][k].logp > 0.1)grid.grid[i][j][k].hydrophobic_score += HM;
			}
	// classify the grid according to its scores
	for (int i = 0; i < grid.dx; i++)
		for (int j = 0; j < grid.dy; j++)
			for (int k = 0; k < grid.dz; k++)
			{
				if ((grid.grid[i][j][k].donor_score >= grid.grid[i][j][k].acceptor_score)
					&& (grid.grid[i][j][k].donor_score >= grid.grid[i][j][k].hydrophobic_score)
					&& (grid.grid[i][j][k].donor_score > 0.001))
				{
					grid.grid[i][j][k].adh = 'D';
				}
				else if ((grid.grid[i][j][k].acceptor_score >= grid.grid[i][j][k].donor_score)
					&& (grid.grid[i][j][k].acceptor_score >= grid.grid[i][j][k].hydrophobic_score)
					&& (grid.grid[i][j][k].acceptor_score > 0.001))
				{
					grid.grid[i][j][k].adh = 'A';
				}
				else if ((grid.grid[i][j][k].hydrophobic_score >= grid.grid[i][j][k].donor_score)
					&& (grid.grid[i][j][k].hydrophobic_score >= grid.grid[i][j][k].acceptor_score)
					&& (grid.grid[i][j][k].hydrophobic_score > 0.001))
				{
					grid.grid[i][j][k].adh = 'H';
				}
				else
				{
					grid.grid[i][j][k].adh = 'N';
				}
			}
	for (int i = 0; i < grid.dx; i++)
		for (int j = 0; j < grid.dy; j++)
			for (int k = 0; k < grid.dz; k++)
			{
				grid.grid[i][j][k].neib = fgrid.grid[i + shift_x][j + shift_y][k + shift_z].neib;
			}

	for (int i = 0; i < grid.dx; i++)
		for (int j = 0; j < grid.dy; j++)
			for (int k = 0; k < grid.dz; k++)
			{
				if (grid.grid[i][j][k].type == 'A' || grid.grid[i][j][k].neib > 900) grid.grid[i][j][k].type = 'E';
				else if (grid.grid[i][j][k].neib > 9) grid.grid[i][j][k].type = 'V';
				else grid.grid[i][j][k].type = 'O';
			}
	Recover(int(parm->SA_radius / 0.5));

	int count = 0;
	for (int i = 0; i < fgrid.dx; i++)
		for (int j = 0; j < fgrid.dy; j++)
			for (int k = 0; k < fgrid.dz; k++)
			{
				if (fgrid.grid[i][j][k].neib == 10) { fgrid.grid[i][j][k].type = 'X'; count++; }
			}
	Clear_Ball_Grid(neib_scale, ball, border_x_min, border_x_max, border_y_min, border_y_max, border_z_min, border_z_max);
	return;
}

//void Pocket::Initialize()
//// 该函数的任务是：
//// 1. 将蛋白原子占据的格点的type标记为'A'；
//// 2. 将每个蛋白原子周围的空格点，根据到蛋白原子的距离，标记其neib
//// 3. 计算每个原子周围grid空格点的donor_score, acceptor_score和hydrophobic_score
//// 4. 统计grid格点的各项分数，并标记adh
//// 5. 根据neib对每个grid的type进行分类
//// 6. 根据neib，标记fgrid中的'X'
//// 
//// 经过该函数后，
//// 1. 找到了fgrid中的所有'X'格点，即between (atom vdw_radius + probe radius) ~ (atom vdw_radius + probe radius + probe mis)的格点。'X'格点将用于global eraser过程
//// 2. 定义了grid格点的三类分数，该分数会被用于最后计算分数；
//// 3. 标记了grid中的格点类型：
////	E: 原子占据的格点,以及溶剂可及表面空间中经过Recover()后剩余的格点
////  V: 原子周围(atom vdw_radius + probe radius + probe mis)以内的空格点
////  O: 原子周围(atom vdw_radius + probe radius + probe mis)以外的空格点
//// 
//{
//	extern Parameter* parm;
//	extern Protein* prot;
//
//	float id, jd, kd;
//	int it = 0, jt = 0, kt = 0;
//	int fit = 0, fjt = 0, fkt = 0;
//	int neib_scale = 0;
//	int neib_min_x, neib_max_x, neib_min_y, neib_max_y, neib_min_z, neib_max_z;
//	neib_min_x = neib_max_x = neib_min_y = neib_max_y = neib_min_z = neib_max_z = 0;
//	float d = 0, d_tmp = 0, d0 = 0, d1 = 0, d2 = 0;
//
//	float r_N = HB_D_R; // the probe atom is N.4 (donor)
//	float r_O = HB_A_R; // the probe atom is O.co2 (acceptor)
//	float r_C = HYDROPHOBIC_R; // the probe atom is C.3 (hydrophobic)
//
//	Point3D v1, v2;
//	float angle;
//
//	//printf("**WARNING: pocket.cpp:Pocket:Initialize(): When calculating donor_score and acceptor_score for each grid, the atom.pos was not correctly defined.\n");
//	float x, y, z, cutoff;
//	float c2_2, c2_3, c2_4;
//
//	// define the global size 
//	cutoff = (float)fgrid.extend_length * 2;
//	neib_scale = fgrid.extend_length * 2;
//	if (cutoff > neib_scale) { cout << "Error." << endl; }
//	int*** ball = nullptr;
//	int border_x_min = -1, border_x_max = -1;
//	int* border_y_min = nullptr, * border_y_max = nullptr;
//	int** border_z_min = nullptr, ** border_z_max = nullptr;
//	Create_Ball_Grid(neib_scale, cutoff, ball, border_x_min, border_x_max, border_y_min, border_y_max, border_z_min, border_z_max);
//
//	int count_0, count_1, count_2, count_3, count_4;
//	count_0=count_1 = count_2 = count_3 = count_4 = 0;
//	int start_x, start_y, start_z, end_x, end_y, end_z;
//	for (auto& iter : prot->atoms)
//	{
//		if (iter.valid == false) { continue; }
//		// grid coordinates of the current atom
//		id = (float)((int)(iter.coord.x * 2)) / 2;
//		jd = (float)((int)(iter.coord.y * 2)) / 2;
//		kd = (float)((int)(iter.coord.z * 2)) / 2;
//		it = (int)((id - grid.min_x) * 2);
//		jt = (int)((jd - grid.min_y) * 2);
//		kt = (int)((kd - grid.min_z) * 2);
//		grid.grid[it][jt][kt].type = 'A';
//		grid.grid[it][jt][kt].neib = 10000;
//		fit = it + shift_x; fjt = jt + shift_y; fkt = kt + shift_z;
//		fgrid.grid[fit][fjt][fkt].type = 'A';
//		fgrid.grid[fit][fjt][fkt].neib = 10000;
//
//		int small_cutoff = 6;
//		neib_min_x = max(it - small_cutoff * 2 + 1, 0);
//		neib_max_x = min(it + small_cutoff * 2 + 1, grid.dx - 1);
//		neib_min_y = max(jt - small_cutoff * 2 + 1, 0);
//		neib_max_y = min(jt + small_cutoff * 2 + 1, grid.dy - 1);
//		neib_min_z = max(kt - small_cutoff * 2 + 1, 0);
//		neib_max_z = min(kt + small_cutoff * 2 + 1, grid.dz - 1);
//		for (int ii = neib_min_x; ii <= neib_max_x; ii++) 
//		{
//			x = iter.coord.x - grid.grid[ii][0][0].coord.x;
//			if (x > small_cutoff) continue;
//			else if (x < -small_cutoff) continue;
//			for (int jj = neib_min_y; jj <= neib_max_y; jj++) 
//			{
//				y = iter.coord.y - grid.grid[ii][jj][0].coord.y;
//				if (y > small_cutoff) continue;
//				else if (y < -small_cutoff) continue;
//				for (int kk = neib_min_z; kk <= neib_max_z; kk++)
//				{
//					// first, find out vacant grids and excluded grids
//					z = iter.coord.z - grid.grid[ii][jj][kk].coord.z;
//					if (z > small_cutoff) continue;
//					else if (z < -small_cutoff) continue;
//					d = (float)(sqrt(x * x + y * y + z * z));
//					//if (d < iter.vdw_r) { grid.grid[ii][jj][kk].neib = 10000; }
//					//else if (grid.grid[ii][jj][kk].neib < 9999 && d < iter.vdw_r + parm->SA_radius) { grid.grid[ii][jj][kk].neib = 1000; }
//					if (d < iter.vdw_r + parm->SA_radius) { grid.grid[ii][jj][kk].neib = 1000; }
//					else if (grid.grid[ii][jj][kk].neib < 999 && d < iter.vdw_r + parm->radius_length) { grid.grid[ii][jj][kk].neib = 100; }
//					else if (grid.grid[ii][jj][kk].neib < 99 && d < iter.vdw_r + parm->radius_length + parm->radius_mis) { grid.grid[ii][jj][kk].neib = 10; }
//					fgrid.grid[ii + shift_x][jj + shift_y][kk + shift_z].neib = grid.grid[ii][jj][kk].neib;
//					if (grid.grid[ii][jj][kk].type == 'A')continue;
//					if (d > small_cutoff) continue;
//					//second, calculate donor score for each grid
//					d0 = r_N + iter.vdw_r;
//					if ((d - d0) < -0.6) { grid.grid[ii][jj][kk].donor_score += (VB); }
//					if (d < 3.5 && iter.q < 0) { grid.grid[ii][jj][kk].donor_score += (SB_AWARD); }
//					if (iter.max_anum > 0 && iter.hb_type != "M")
//					{
//						if (d < d0) // maximum h-bond length
//						{
//							v1 = iter.root - iter.coord;
//							v2 = iter.coord - grid.grid[ii][jj][kk].coord;
//							angle = 180 - Angle_Of_Two_Points(v1, v2);
//							angle = fabs(angle - iter.hangle);
//							if (angle <= OLD_HB_ANGLE_CUTOFF) // h-bond angle cutoff
//							{
//								// the probe may form HB pair with nearby pocket atoms, only the shortest HBs are taken into account
//								if (d - d0 < grid.grid[ii][jj][kk].donor_overlap)
//								{
//									if (iter.tripos_type == "O.w") grid.grid[ii][jj][kk].donor_wtype = 2;
//									else grid.grid[ii][jj][kk].donor_wtype = 1;
//									grid.grid[ii][jj][kk].donor_overlap = d - d0;
//								}
//							}
//						}
//					}
//					// third, calculate the acceptor score for each grid
//					d1 = r_O + iter.vdw_r;
//					if ((d - d1) < -0.6) { grid.grid[ii][jj][kk].acceptor_score += (VB); }
//					if (d < 3.5 && iter.q >0) { grid.grid[ii][jj][kk].acceptor_score += (SB_AWARD); }
//					if (iter.max_dnum > 0)
//					{
//						if (iter.hb_type != "M")
//						{
//							if (d < d1)
//							{
//								v1 = iter.root - iter.coord;
//								v2 = iter.coord - grid.grid[ii][jj][kk].coord;
//								angle = 180.0 - Angle_Of_Two_Points(v1, v2);
//								angle = fabs(angle - iter.hangle);
//								if (angle <= OLD_HB_ANGLE_CUTOFF)
//								{
//									if (d - d1 < grid.grid[ii][jj][kk].acceptor_overlap)
//									{
//										if (iter.tripos_type == "O.w")grid.grid[ii][jj][kk].acceptor_wtype = 2;
//										else grid.grid[ii][jj][kk].acceptor_wtype = 1;
//										grid.grid[ii][jj][kk].acceptor_overlap = d - d1;
//									}
//								}
//							}
//						}
//						else if (iter.hb_type == "M")
//						{
//							if (d < 2.0)grid.grid[ii][jj][kk].acceptor_score += (MB);
//							else if (d < 3.0) grid.grid[ii][jj][kk].acceptor_score += ((3.0 - d) * (MB));
//						}
//					}
//					// fourth, calculate the hydrophobic score for each grid
//					d2 = r_C + iter.vdw_r;
//					if ((d - d2) < -0.6)grid.grid[ii][jj][kk].hydrophobic_score += (VB);
//					if (d < 5.0)grid.grid[ii][jj][kk].logp += iter.logp;
//				}
//			}
//		}
//
//		// traverse every grid around the atom 
//		cutoff = iter.vdw_r + parm->radius_length + parm->radius_mis;
//		//c2_1 = iter.vdw_r * iter.vdw_r;
//		c2_2 = (iter.vdw_r + parm->SA_radius) * (iter.vdw_r + parm->SA_radius);
//		c2_3 = (iter.vdw_r + parm->radius_length) * (iter.vdw_r + parm->radius_length);
//		c2_4 = cutoff * cutoff;
//		start_x = fit + border_x_min - neib_scale;
//		end_x = fit + border_x_max - neib_scale;
//		for (int ii = start_x; ii < end_x; ii++) 
//		{
//			x = iter.coord.x - fgrid.grid[ii][0][0].coord.x;
//			if (x > cutoff) { continue; }
//			else if (x < -cutoff) { continue; }
//			start_y = fjt + border_y_min[ii + neib_scale - fit] - neib_scale;
//			end_y = fjt + border_y_max[ii + neib_scale - fit] - neib_scale;
//			for (int jj = start_y; jj < end_y; jj++) 
//			{
//				y = iter.coord.y - fgrid.grid[ii][jj][0].coord.y;
//				if (y > cutoff) { continue; }
//				else if (y < -cutoff) { continue; }
//				d_tmp = x * x + y * y;
//				if (d_tmp > c2_4) { continue; }
//				start_z = fkt + border_z_min[ii + neib_scale - fit][jj + neib_scale - fjt] - neib_scale;
//				end_z = fkt + border_z_max[ii + neib_scale - fit][jj + neib_scale - fjt] - neib_scale;
//				for (int kk = start_z; kk < end_z; kk++)
//				{
//					//if (fgrid.grid[ii][jj][kk].neib > 100) { continue; }
//					z = iter.coord.z - fgrid.grid[ii][jj][kk].coord.z;
//					if (z > cutoff) { continue; }
//					else if (z < -cutoff) continue;
//					d = d_tmp + z * z;
//					//if (d < c2_1) { fgrid.grid[ii][jj][kk].neib = 10000; }
//					//else if (d < c2_2 && fgrid.grid[ii][jj][kk].neib < 9999) { fgrid.grid[ii][jj][kk].neib = 1000; }
//					if (d < c2_2) { fgrid.grid[ii][jj][kk].neib = 1000; }
//					else if (d < c2_3 && fgrid.grid[ii][jj][kk].neib < 999) { fgrid.grid[ii][jj][kk].neib = 100; }
//					else if (d < c2_4 && fgrid.grid[ii][jj][kk].neib < 99) { fgrid.grid[ii][jj][kk].neib = 10; }
//					/* neib/type:
//							10000/A: it is the atom
//							10000/V: in atom vdw_radius
//							1000/V:  between atom vdw_radius ~ (atom vdw_radius + water_access)
//							100/V:   between (atom vdw_radius + water_access) ~ (atom vdw_radius + probe radius)
//							10/V:    between (atom vdw_radius + probe radius) ~ (atom vdw_radius + probe radius + probe mis)
//							0/V:	 others
//					*/
//				}
//			}
//		}
//	}
//
//	for (int i = 0; i < grid.dx; i++)
//		for (int j = 0; j < grid.dy; j++)
//			for (int k = 0; k < grid.dz; k++)
//			{
//				if (grid.grid[i][j][k].donor_wtype == 1)
//				{
//					if (grid.grid[i][j][k].donor_overlap < -0.6) grid.grid[i][j][k].donor_score += SHB;
//					else if (grid.grid[i][j][k].donor_overlap < -0.3) grid.grid[i][j][k].donor_score += MHB;
//					else if (grid.grid[i][j][k].donor_overlap < 0.0) grid.grid[i][j][k].donor_score += WHB;
//				}
//				else if (grid.grid[i][j][k].donor_wtype == 2)
//				{
//					if (grid.grid[i][j][k].donor_overlap < -0.6) grid.grid[i][j][k].donor_score += SWH;
//					else if (grid.grid[i][j][k].donor_overlap < -0.3) grid.grid[i][j][k].donor_score += MWH;
//					else if (grid.grid[i][j][k].donor_overlap < 0.0) grid.grid[i][j][k].donor_score += WWH;
//				}
//				if (grid.grid[i][j][k].acceptor_wtype == 1)
//				{
//					if (grid.grid[i][j][k].acceptor_overlap < -0.6) grid.grid[i][j][k].acceptor_score += SHB;
//					else if (grid.grid[i][j][k].acceptor_overlap < -0.3) grid.grid[i][j][k].acceptor_score += MHB;
//					else if (grid.grid[i][j][k].acceptor_overlap < 0.0) grid.grid[i][j][k].acceptor_score += WHB;
//				}
//				else if (grid.grid[i][j][k].acceptor_wtype == 2)
//				{
//					if (grid.grid[i][j][k].acceptor_overlap < -0.6) grid.grid[i][j][k].acceptor_score += SWH;
//					else if (grid.grid[i][j][k].acceptor_overlap < -0.3) grid.grid[i][j][k].acceptor_score += MWH;
//					else if (grid.grid[i][j][k].acceptor_overlap < 0.0) grid.grid[i][j][k].acceptor_score += WWH;
//				}
//				if (grid.grid[i][j][k].logp > 0.1)grid.grid[i][j][k].hydrophobic_score += HM;
//			}
//	// classify the grid according to its scores
//	for (int i = 0; i < grid.dx; i++)
//		for (int j = 0; j < grid.dy; j++)
//			for (int k = 0; k < grid.dz; k++)
//			{
//				if ((grid.grid[i][j][k].donor_score >= grid.grid[i][j][k].acceptor_score)
//					&& (grid.grid[i][j][k].donor_score >= grid.grid[i][j][k].hydrophobic_score)
//					&& (grid.grid[i][j][k].donor_score > 0.001))
//				{
//					grid.grid[i][j][k].adh = 'D'; 
//				}
//				else if ((grid.grid[i][j][k].acceptor_score >= grid.grid[i][j][k].donor_score)
//					&& (grid.grid[i][j][k].acceptor_score >= grid.grid[i][j][k].hydrophobic_score)
//					&& (grid.grid[i][j][k].acceptor_score > 0.001))
//				{
//					grid.grid[i][j][k].adh = 'A'; 
//				}
//				else if ((grid.grid[i][j][k].hydrophobic_score >= grid.grid[i][j][k].donor_score)
//					&& (grid.grid[i][j][k].hydrophobic_score >= grid.grid[i][j][k].acceptor_score)
//					&& (grid.grid[i][j][k].hydrophobic_score > 0.001))
//				{
//					grid.grid[i][j][k].adh = 'H'; 
//				}
//				else
//				{
//					grid.grid[i][j][k].adh = 'N'; 
//				}
//			}
//	for (int i = 0; i < grid.dx; i++)
//		for (int j = 0; j < grid.dy; j++)
//			for (int k = 0; k < grid.dz; k++)
//			{
//				grid.grid[i][j][k].neib = fgrid.grid[i + shift_x][j + shift_y][k + shift_z].neib;
//			}
//
//	for (int i = 0; i < grid.dx; i++)
//		for (int j = 0; j < grid.dy; j++)
//			for (int k = 0; k < grid.dz; k++)
//			{
//				if (grid.grid[i][j][k].type == 'A' || grid.grid[i][j][k].neib > 900) grid.grid[i][j][k].type = 'E';
//				else if (grid.grid[i][j][k].neib > 9) grid.grid[i][j][k].type = 'V';
//				else grid.grid[i][j][k].type = 'O';
//			}
//	Recover(int(parm->SA_radius / 0.5));
//
//	int count = 0;
//	for (int i = 0; i < fgrid.dx; i++)
//		for (int j = 0; j < fgrid.dy; j++)
//			for (int k = 0; k < fgrid.dz; k++)
//			{
//				if (fgrid.grid[i][j][k].neib == 10) { fgrid.grid[i][j][k].type = 'X'; count++; }
//			}
//	Clear_Ball_Grid(neib_scale, ball, border_x_min, border_x_max, border_y_min, border_y_max, border_z_min, border_z_max);
//	return;
//}

void Pocket::Recover(int maxdep)
// The grids between (radius ~ radius + SA_radius) were labeled as 'E' in the Initialize step, now recover them to 'V'
{
	int i, j, k;
	int dep = 0;
	for (dep = 0; dep < maxdep; dep++) // maxdep is the depth of SA_radius (the number of grid layers)
	{
		for (i = 0; i < grid.dx; i++)
			for (j = 0; j < grid.dy; j++)
				for (k = 0; k < grid.dz; k++)
				{
					if (grid.grid[i][j][k].type == 'E')
					{
						if (grid.OutOfBoundary(i - 1, j, k) == 0 && grid.grid[i - 1][j][k].type == 'V') grid.grid[i][j][k].type = 'W';
						else if (grid.OutOfBoundary(i + 1, j, k) == 0 && grid.grid[i + 1][j][k].type == 'V') grid.grid[i][j][k].type = 'W';
						else if (grid.OutOfBoundary(i, j - 1, k) == 0 && grid.grid[i][j - 1][k].type == 'V') grid.grid[i][j][k].type = 'W';
						else if (grid.OutOfBoundary(i, j + 1, k) == 0 && grid.grid[i][j + 1][k].type == 'V') grid.grid[i][j][k].type = 'W';
						else if (grid.OutOfBoundary(i, j, k - 1) == 0 && grid.grid[i][j][k - 1].type == 'V') grid.grid[i][j][k].type = 'W';
						else if (grid.OutOfBoundary(i, j, k + 1) == 0 && grid.grid[i][j][k + 1].type == 'V') grid.grid[i][j][k].type = 'W';
						else continue;
					}
				}
		for (i = 0; i < grid.dx; i++)
			for (j = 0; j < grid.dy; j++)
				for (k = 0; k < grid.dz; k++)
				{
					if (grid.grid[i][j][k].type == 'W')grid.grid[i][j][k].type = 'V';
				}
	}
	for (i = 0; i < grid.dx; i++)
		for (j = 0; j < grid.dy; j++)
			for (k = 0; k < grid.dz; k++)
			{
				if (grid.grid[i][j][k].type == 'E')
				{
					excludedset.push_back(Position3D(i, j, k)); // store the 'E' grid
					grid.grid[i][j][k].predict = -2;
					grid.grid[i][j][k].predict2 = -2;
				}
			}
	return;
}

void Pocket::Erase()
{
	extern Parameter* parm;
	if (parm->info == 1)	cout << "**INFO: Create full eraser ball, radius is " << setw(5) << setprecision(2) << fixed << parm->radius_length << "A" << endl;
	// define a full ball, radius=?A
	fullradius = (int)(2 * parm->radius_length) + 2;
	float cutoff = (parm->radius_length + parm->radius_mis) * 2;
	Create_Ball_Grid(fullradius, cutoff, fullball, border_min_x, border_max_x, border_min_y, border_max_y, border_min_z, border_max_z);
	length = fullradius * 2 + 1;

	int x, y, z;
	int diff_i, diff_j, diff_k;
	for (int m = 0; m < 6; m++)
	{
		vector<Position3D> shiftset;
		if (m == 0) { diff_i = -1; diff_j = 0; diff_k = 0; }
		else if (m == 1) { diff_i = 1; diff_j = 0; diff_k = 0; }
		else if (m == 2) { diff_i = 0; diff_j = -1; diff_k = 0; }
		else if (m == 3) { diff_i = 0; diff_j = 1; diff_k = 0; }
		else if (m == 4) { diff_i = 0; diff_j = 0; diff_k = -1; }
		else if (m == 5) { diff_i = 0; diff_j = 0; diff_k = 1; }
		for (int i = border_min_x; i < border_max_x; i++)
		{
			for (int j = border_min_y[i]; j < border_max_y[i]; j++)
			{
				for (int k = border_min_z[i][j]; k < border_max_z[i][j]; k++)
				{
					x = i + diff_i;
					y = j + diff_j;
					z = k + diff_k;
					if (x >= border_max_x || x < border_min_x) {}
					else if (y >= border_max_y[x] || y < border_min_y[x]) {}
					else if (z >= border_max_z[x][y] || z < border_min_z[x][y]) {}
					else continue;
					shiftset.push_back(Position3D(i, j, k));
				}
			}
		}
		shiftgrid.insert(pair<int, vector<Position3D>>(m, shiftset));
	}
	
	if (parm->info == 1) cout << "**INFO: Now starting global eraser ..." << endl;
	int min_x, min_y, min_z, max_x, max_y, max_z;
	int tmp_min_x, tmp_min_y, tmp_min_z, tmp_max_x, tmp_max_y, tmp_max_z;
	vector<Ppair> xgrid;
	int i, j, k;
	int ii, jj, kk;
	int id, jd, kd;
	for (i = 0; i < fgrid.dx; i++)
		for (j = 0; j < fgrid.dy; j++)
			for (k = 0; k < fgrid.dz; k++)
			{
				if (fgrid.grid[i][j][k].type == 'X')
				{
					fgrid.grid[i][j][k].type = 'O';
					fgrid.grid[i][j][k].type_bak = 'O';
					min_x = max(fullradius - i, 0);
					max_x = min(length, fgrid.dx + fullradius - i);
					min_y = max(fullradius - j, 0);
					max_y = min(length, fgrid.dy + fullradius - j);
					min_z = max(fullradius - k, 0);
					max_z = min(length, fgrid.dz + fullradius - k);
					tmp_min_x = max(min_x, border_min_x);
					tmp_max_x = min(max_x, border_max_x);
					for (id = tmp_min_x; id < tmp_max_x; id++)
					{
						tmp_min_y = max(min_y, border_min_y[id]);
						tmp_max_y = min(max_y, border_max_y[id]);
						for (jd = tmp_min_y; jd < tmp_max_y; jd++)
						{
							tmp_min_z = max(min_z, border_min_z[id][jd]);
							tmp_max_z = min(max_z, border_max_z[id][jd]);
							for (kd = tmp_min_z; kd < tmp_max_z; kd++)
							{
								ii = i + id - fullradius;
								jj = j + jd - fullradius;
								kk = k + kd - fullradius;
								fgrid.grid[ii][jj][kk].type_bak = 'O';
							}
						}
					}
					if (i > 0 && fgrid.grid[i - 1][j][k].type == 'X') xgrid.push_back(make_pair(Position3D(i, j, k), Position3D(i - 1, j, k)));
					if (i < fgrid.dx - 1 && fgrid.grid[i + 1][j][k].type == 'X')xgrid.push_back(make_pair(Position3D(i, j, k), Position3D(i + 1, j, k)));
					if (j < 0 && fgrid.grid[i][j - 1][k].type == 'X') xgrid.push_back(make_pair(Position3D(i, j, k), Position3D(i, j - 1, k)));
					if (j < fgrid.dy - 1 && fgrid.grid[i][j + 1][k].type == 'X') xgrid.push_back(make_pair(Position3D(i, j, k), Position3D(i, j + 1, k)));
					if (k > 0 && fgrid.grid[i][j][k - 1].type == 'X') xgrid.push_back(make_pair(Position3D(i, j, k), Position3D(i, j, k - 1)));
					if (k < fgrid.dz - 1 && fgrid.grid[i][j][k + 1].type == 'X') xgrid.push_back(make_pair(Position3D(i, j, k), Position3D(i, j, k + 1)));
					for (int it = 0; it < xgrid.size(); it++)
					{
						if (fgrid(xgrid[it].second).type == 'X')
						{
							Erase_Cumulative(xgrid[it], xgrid);
						}
					}
					xgrid.clear();
				}
			}
	for (int i = 0; i < grid.dx; i++)
		for (int j = 0; j < grid.dy; j++)
			for (int k = 0; k < grid.dz; k++)
			{
				if (grid.grid[i][j][k].type == 'V')
				{
					id = i + shift_x;
					jd = j + shift_y;
					kd = k + shift_z;
					if (fgrid.grid[id][jd][kd].type_bak == 'O') grid.grid[i][j][k].type = 'O';
				}
			}

	fgrid.Clear();
	// free memory
	Clear_Ball_Grid(fullradius, fullball, border_min_x, border_max_x, border_min_y, border_max_y, border_min_z, border_max_z);
	return;
}

void Pocket::Erase_Cumulative(Ppair& p, vector<Ppair>& list)
{
	int last_i = p.first.i, last_j = p.first.j, last_k = p.first.k;
	int i = p.second.i, j = p.second.j, k = p.second.k;
	int diff_i = i - last_i, diff_j = j - last_j, diff_k = k - last_k;
	int id, jd, kd;
	int ii, jj, kk;
	fgrid.grid[i][j][k].type = 'O';
	fgrid.grid[i][j][k].type_bak = 'O';
	int min_x = max(fullradius - i, 0);
	int max_x = min(length, fgrid.dx + fullradius - i);
	int min_y = max(fullradius - j, 0);
	int max_y = min(length, fgrid.dy + fullradius - j);
	int min_z = max(fullradius - k, 0);
	int max_z = min(length, fgrid.dz + fullradius - k);
	int m = 0;
	if (diff_i == -1) m = 0;
	else if (diff_i == 1) m = 1;
	else if (diff_j == -1) m = 2;
	else if (diff_j == 1) m = 3;
	else if (diff_k == -1) m = 4;
	else if (diff_k == 1) m = 5;
	for (auto& sit : shiftgrid[m])
	{
		id = sit.i; jd = sit.j; kd = sit.k;
		if (id >= min_x && id < max_x && jd >= min_y && jd < max_y && kd >= min_z && kd < max_z)
		{
			ii = i + id - fullradius;
			jj = j + jd - fullradius;
			kk = k + kd - fullradius;
			fgrid.grid[ii][jj][kk].type_bak = 'O';
		}
	}
	if (i > 0 && fgrid.grid[i - 1][j][k].type == 'X') list.push_back(make_pair(Position3D(i, j, k), Position3D(i - 1, j, k)));
	if (i < fgrid.dx - 1 && fgrid.grid[i + 1][j][k].type == 'X')list.push_back(make_pair(Position3D(i, j, k), Position3D(i + 1, j, k)));
	if (j < 0 && fgrid.grid[i][j - 1][k].type == 'X') list.push_back(make_pair(Position3D(i, j, k), Position3D(i, j - 1, k)));
	if (j < fgrid.dy - 1 && fgrid.grid[i][j + 1][k].type == 'X') list.push_back(make_pair(Position3D(i, j, k), Position3D(i, j + 1, k)));
	if (k > 0 && fgrid.grid[i][j][k - 1].type == 'X') list.push_back(make_pair(Position3D(i, j, k), Position3D(i, j, k - 1)));
	if (k < fgrid.dz - 1 && fgrid.grid[i][j][k + 1].type == 'X') list.push_back(make_pair(Position3D(i, j, k), Position3D(i, j, k + 1)));
	return;
}

void Pocket::Define_Surface()
// find the surface of protein surface
// surface grid分为两类，一类是潜在口袋的surface，其cavity=0, inner=1, 另一类是非口袋的surface，其cavity=-1，inner=0
// 另外，位于口袋交界面上的grid的cavity=0， inner=0。剩余的其他格点，则是默认的cavity=0,inner=-1
// 则在此之后，所有格点的状态应该为4种：
// 1. 口袋表面格点：type='S', cavity=0, inner=1;
// 2. 非口袋表面格点：type='S', cavity=-1, inner=0;
// 3. 交界面的表面格点，type='S', cavity=0, inner=0;
// 4. 其他非表面格点： cavity=0, inner=-1。
// 但inner的值会在reedge里重新初始化为-1，不知道这里给inner赋值有什么意义。
{
	int flag = 0;
	int i, j, k;
	for (auto& eit : excludedset)
	{
		i = eit.i, j = eit.j, k = eit.k;
		flag = 0;

		Position3D neib_point[6];
		neib_point[0].set(i - 1, j, k);
		neib_point[1].set(i + 1, j, k);
		neib_point[2].set(i, j - 1, k);
		neib_point[3].set(i, j + 1, k);
		neib_point[4].set(i, j, k - 1);
		neib_point[5].set(i, j, k + 1);
		for (auto& pi : neib_point)
		{
			if (grid.OutOfBoundary(pi.i, pi.j, pi.k) == 0)
			{
				if (grid(pi).type == 'V')
				{
					grid.grid[i][j][k].type = 'S';
					grid.grid[i][j][k].grid_surface += 1;
					grid.grid[i][j][k].cavity = 0;
					if (flag == 0)grid.grid[i][j][k].inner = 1;
					flag = 1;
				}
				else if (grid(pi).type == 'O')
				{
					grid.grid[i][j][k].type = 'S';
					grid.grid[i][j][k].grid_surface += 1;
					if (flag == 0)grid.grid[i][j][k].cavity = -1;
					grid.grid[i][j][k].inner = 0;
					flag = 1;
				}
			}
		}
		if (grid.grid[i][j][k].type == 'S' && grid.grid[i][j][k].cavity == 0)
		{
			surfacegrid.push_back(eit);
		}
		if (grid.grid[i][j][k].grid_surface > 1) grid.grid[i][j][k].grid_surface = (float)(grid.grid[i][j][k].grid_surface / 2 * 1.414); // EQUAL
	}
	return;
}

void Pocket::Define_Vacant()
// 找到所有cavity=0的'V'格点
{
	vector<Position3D> vlist;
	for (int i = 0; i < grid.dx; i++)
		for (int j = 0; j < grid.dy; j++)
			for (int k = 0; k < grid.dz; k++)
			{
				if (grid.grid[i][j][k].cavity == 0 && grid.grid[i][j][k].type == 'V')
				{
					vlist.push_back(Position3D(i, j, k));
				}
			}
	vacantset.insert(pair<int, vector<Position3D>>(0, vlist));
	return;
}

void Pocket::Reedge()
// 为grid.inner赋值。首先初始化grid.inner=-1 (这样的话，之前给inner赋值就没有意义了，待观察)
// 对所有cavity>=0的格点，根据以下情形分配inner值：
// 1. 当前格点type='S', 周围一圈格点（26个）中也有type='S'的格点，且该邻居格点不是非口袋的表面格点，则inner=1, 若邻居格点是非口袋的表面格点，则inner=0;
// 2. 当前格点的type='V',周围格点（这里是上下左右前后6个点）中有type='O'的格点，或者有type='V'且cavity值与当前格点不一样的格点，则inner=0
// 最终，inner=0 和inner=1的点围成了整个口袋
{
	int it, jt, kt;
	for (auto& sit : surfacegrid)
	{
		grid(sit).inner = -1;
	}
	int i, j, k;
	int ia, ja, ka, ib, jb, kb;
	for (auto& sit : surfacegrid)
	{
		if (grid(sit).cavity >= 0)
		{
			i = sit.i, j = sit.j, k = sit.k;
			ia = max(i - 1, 0), ib = min(i + 2, grid.dx);
			ja = max(j - 1, 0), jb = min(j + 2, grid.dy);
			ka = max(k - 1, 0), kb = min(k + 2, grid.dz);
			for (it = ia; it < ib; it++)
				for (jt = ja; jt < jb; jt++)
					for (kt = ka; kt < kb; kt++)
					{
						if (it == i && jt == j && kt == k)continue;
						if (grid.grid[it][jt][kt].type == 'S' && grid.grid[i][j][k].inner != 0)
						{
							grid.grid[i][j][k].inner = 1;
							if (grid.grid[it][jt][kt].cavity != grid.grid[i][j][k].cavity)
							{
								grid.grid[i][j][k].inner = 0;
							}
							continue;
						}
					}
		}
	}
	for (map<int, vector<Position3D>>::iterator iter = vacantset.begin(); iter != vacantset.end(); iter++)
	{
		for (auto& vit : iter->second)
		{
			i = vit.i, j = vit.j, k = vit.k;
			if (i > 0 && (grid.grid[i - 1][j][k].type == 'O' ||
				(grid.grid[i - 1][j][k].type == 'V' && grid.grid[i - 1][j][k].cavity != grid.grid[i][j][k].cavity)))
			{
				grid.grid[i][j][k].inner = 0;
			}
			else if (i + 1 < grid.dx && (grid.grid[i + 1][j][k].type == 'O' ||
				(grid.grid[i + 1][j][k].type == 'V' && grid.grid[i + 1][j][k].cavity != grid.grid[i][j][k].cavity)))
			{
				grid.grid[i][j][k].inner = 0;
			}
			else if (j > 0 && (grid.grid[i][j - 1][k].type == 'O' ||
				(grid.grid[i][j - 1][k].type == 'V' && grid.grid[i][j - 1][k].cavity != grid.grid[i][j][k].cavity)))
			{
				grid.grid[i][j][k].inner = 0;
			}
			else if (j + 1 < grid.dy && (grid.grid[i][j + 1][k].type == 'O' ||
				(grid.grid[i][j + 1][k].type == 'V' && grid.grid[i][j + 1][k].cavity != grid.grid[i][j][k].cavity)))
			{
				grid.grid[i][j][k].inner = 0;
			}
			else if (k > 0 && (grid.grid[i][j][k - 1].type == 'O' ||
				(grid.grid[i][j][k - 1].type == 'V' && grid.grid[i][j][k - 1].cavity != grid.grid[i][j][k].cavity)))
			{
				grid.grid[i][j][k].inner = 0;
			}
			else if (k + 1 < grid.dz && (grid.grid[i][j][k + 1].type == 'O' ||
				(grid.grid[i][j][k + 1].type == 'V' && grid.grid[i][j][k + 1].cavity != grid.grid[i][j][k].cavity)))
			{
				grid.grid[i][j][k].inner = 0;
			}
			else grid.grid[i][j][k].inner = -1;
		}
	}
	return;
}

void Pocket::Get_Real_Depth(int n)
// n is used to compare with grid.cavity
// used to assign the value of grid.rdep (default=-1) and grid.dep (default=0)
// dep是指空腔内的点到空腔外表面（V/O交界面）之间的最短距离，rdep是实际距离，dep是整数化的距离指标。但这里rdep的赋值方式非常奇怪
{
	if (vacantset.count(n) == 0)
	{
		cout << "No valid vacant grid." << endl;
	}
	vector<Position3D> rgrid;
	map<int, map<int, vector<int>>> newrgrid;
	for (auto& vit : vacantset[n])
	{
		if (grid(vit).inner == 0)
		{
			//rgrid_num++;
			rgrid.push_back(vit);
			grid(vit).rdep = 0;
			grid(vit).dep = 0;
			if (newrgrid.count(vit.i) == 0)
			{
				map<int, vector<int>> m;
				vector<int> tmp;
				tmp.push_back(vit.k);
				m.insert(pair<int, vector<int>>(vit.j, tmp));
				newrgrid.insert(pair<int, map<int, vector<int>>>(vit.i, m));
			}
			else
			{
				if (newrgrid[vit.i].count(vit.j) == 0)
				{
					vector<int> tmp;
					tmp.push_back(vit.k);
					newrgrid[vit.i].insert(pair<int, vector<int>>(vit.j, tmp));
				}
				else
				{
					newrgrid[vit.i][vit.j].push_back(vit.k);
				}
			}
		}
	}
	if (rgrid.size() == 0)
	{
		for (auto& vit : vacantset[n])
		{
			grid(vit).rdep = 0;
			grid(vit).dep = 0;
		}
		return;
	}
	sort(rgrid.begin(), rgrid.end());
	Position3D init;
	int min, now;
	int x, y, z;
	int i, j, k;
	//int rgrid_num = rgrid.size();
	//int m = -1;
	int sum;
	for (auto& vit : vacantset[n])
	{
		if (grid(vit).inner == 0)continue;
		i = vit.i, j = vit.j, k = vit.k;
		min = (i - init.i) * (i - init.i) + (j - init.j) * (j - init.j) + (k - init.k) * (k - init.k);
		//m = -1;
		for (auto & pi : newrgrid)
		{
			x = (pi.first - i) * (pi.first - i);
			if (x >= min) continue;
			for (auto& pj : pi.second)
			{
				y = (pj.first - j) * (pj.first - j);
				if (y >= min) continue;
				sum = x + y;
				if (sum >= min) continue;
				for (auto& pk : pj.second)
				{
					z = (pk - k) * (pk - k);
					if (z >= min) continue;
					now = sum + z;
					if (min > now)
					{
						min = now;
						init = Position3D(pi.first, pj.first, pk);
					}
				}
			}
		}
		grid(vit).rdep = (float)Distance(grid(vit).coord, grid(init).coord) - 1.0;
		grid(vit).dep = (int)(grid(vit).rdep * 2 + 0.5);
	}
	return;
}

int Pocket::Check(int n, int seed)
{
	return 0;
}

int Pocket::Detect_V_Cluster(int n, int seed, int* vnum, vector<int>& tmpvnum)
// 对于当前格点中所有cavity=seed的格点，对于每块连在一起的‘V'格点空间进行查找和格点计数。每块空间从n开始编号
{
	if (vacantset.count(seed) == 0)
	{
		cout << "No valid vacant grid." << endl;
		return 0;
	}
	extern Parameter* parm;
	int ii, jj, kk;
	int n_bak = n;
	tmpvnum.clear();
	for (auto& vit : vacantset[seed])
	{
		if (grid(vit).cavity != seed) continue;
		n++;
		vector<Position3D> vlist;
		grid(vit).cavity = n;
		vlist.push_back(vit);
		for (int it = 0; it < vlist.size(); it++)
		{
			ii = vlist[it].i; jj = vlist[it].j; kk = vlist[it].k;
			if (ii > 0)
			{
				if (grid.grid[ii - 1][jj][kk].type == 'V' && grid.grid[ii - 1][jj][kk].cavity == seed)
				{
					grid.grid[ii - 1][jj][kk].cavity = n;
					vlist.push_back(Position3D(ii - 1, jj, kk));
				}
			}
			if (ii < grid.dx - 1)
			{
				if (grid.grid[ii + 1][jj][kk].type == 'V' && grid.grid[ii + 1][jj][kk].cavity == seed)
				{
					grid.grid[ii + 1][jj][kk].cavity = n;
					vlist.push_back(Position3D(ii + 1, jj, kk));
				}
			}
			if (jj > 0)
			{
				if (grid.grid[ii][jj - 1][kk].type == 'V' && grid.grid[ii][jj - 1][kk].cavity == seed)
				{
					grid.grid[ii][jj - 1][kk].cavity = n;
					vlist.push_back(Position3D(ii, jj - 1, kk));
				}
			}
			if (jj < grid.dy - 1)
			{
				if (grid.grid[ii][jj + 1][kk].type == 'V' && grid.grid[ii][jj + 1][kk].cavity == seed)
				{
					grid.grid[ii][jj + 1][kk].cavity = n;
					vlist.push_back(Position3D(ii, jj + 1, kk));
				}
			}
			if (kk > 0)
			{
				if (grid.grid[ii][jj][kk - 1].type == 'V' && grid.grid[ii][jj][kk - 1].cavity == seed)
				{
					grid.grid[ii][jj][kk - 1].cavity = n;
					vlist.push_back(Position3D(ii, jj, kk - 1));
				}
			}
			if (kk < grid.dz - 1)
			{
				if (grid.grid[ii][jj][kk + 1].type == 'V' && grid.grid[ii][jj][kk + 1].cavity == seed)
				{
					grid.grid[ii][jj][kk + 1].cavity = n;
					vlist.push_back(Position3D(ii, jj, kk + 1));
				}
			}
		}
		vnum[n] = (int)vlist.size();
		tmpvnum.push_back(vnum[n]);
		if (vnum[n] < parm->chip_v) // 口袋n的格点数小于chip_v，舍弃该口袋
		{
			for (auto& it : vlist)
			{
				grid(it).cavity = -1;
				grid(it).type = 'O';
			}
			n--;
			continue;
		}
		else
		{
			vacantset.insert(pair<int, vector<Position3D>>(n, vlist));
		}
		if (seed == 0)
		{
			if (parm->info == 1) printf("**INFO: New matrix detected: No. %d, vnum = %d\n", n, vnum[n]);
		}
		else
		{
			if (parm->info == 1) { printf("**INFO: Sub cavity No. %d is separated from matrix No. %d vnum = %d\n", n, seed, vnum[n]); }
		}
	}
	vacantset.erase(seed);
	if (seed == 0) { printf("Cavity matrix detecting completed.\nTotal matrix num is %d\n", n); }
	else if (n - n_bak > 1 && parm->info == 1) { printf("**INFO: Matrix No. %d has been separated to %d parts\n", seed, n - n_bak); }
	else if (n - n_bak == 1 && parm->info == 1) { printf("**INFO: Just one part left after separate Matrix %d\n", seed); }
	else { if (parm->info == 1) printf("**INFO: Matrix %d has nothing left after operation\n", seed); }
	/*printf("  erased grid: ");
	for (auto& a : tmpvnum) { printf("%d ", a); }
	printf("\n");*/
	return n;
}

int Pocket::Get_Depth_Predict(int n)
//  返回的是口袋n中,层格点数大于max_depth_vacant的层数，层在这里由dep定义
//
{
	extern Parameter* parm;
	int* depth=new int[300];
	for (int i = 0; i < 300; i++) depth[i] = 0;
	int max = 0, min = 300, limit = 1;
	int dep = 0, mark = 0;
	for (auto& it : vacantset[n])
	{
		dep = grid(it).dep;
		if (dep < 300) depth[dep]++;
		else mark = 1;
		if (dep > max) max = dep;
		if (dep < min) min = dep;
	}
	if (mark == 1) cout << "**WARNING: grid.dep >= 300" << endl;
	if (min < 0)min = 0;
	limit = min;
	for (int i = min; i <= max; i++)
	{
		if (i >= 300) break;
		if (depth[i] > parm->max_depth_vacant) limit++;
	}
	delete[] depth;
	return (limit - min);
}

int Pocket::Get_Predict_V(int n, int shrink)
// 口袋n收缩了shrink次，将口袋n的shrink次收缩倒回去，看看收缩前的格点数
{
	extern Parameter* parm;
	int i, j, k;
	vector<Position3D> vset0;
	vector<Position3D> vset1;
	for (auto& vit : vacantset[n])
	{
		if (grid(vit).inner == 0)
		{
			grid(vit).predict = 0;
			vset0.push_back(vit);
		}
		else if (grid(vit).inner == 1)
		{
			grid(vit).predict = 1;
			vset1.push_back(vit);
		}
	}

	int flag;
	int it, jt, kt;
	for (flag = 0; flag < shrink; flag++)
	{
		vector<Position3D> vset2;
		for (auto& iter : vset0)
		{
			i = iter.i; j = iter.j; k = iter.k;
			grid(iter).predict = 1;
			for (it = i - 1; it < i + 2; it++)
				for (jt = j - 1; jt < j + 2; jt++)
					for (kt = k - 1; kt < k + 2; kt++)
					{
						if (it == i && jt == j && kt == k) continue;
						if (grid.OutOfBoundary(it, jt, kt) == 0)
						{
							if (grid.grid[it][jt][kt].predict == -1)
							{
								if (grid(iter).re - grid.grid[it][jt][kt].re > 1) continue;
								grid.grid[it][jt][kt].predict = 2;
								vset2.push_back(Position3D(it, jt, kt));
							}
						}
					}
		}
		vset1.insert(vset1.end(), vset0.begin(), vset0.end());
		vset0.clear();
		for (auto& iter : vset2)
		{
			grid(iter).predict = 0;
			vset0.push_back(iter);
		}
	}
	int total = 0;
	total = (int)(vset1.size() + vset0.size());
	for (auto& vit : vset1)
	{
		grid(vit).predict = grid(vit).predict2;
	}
	for (auto& vit : vset0)
	{
		grid(vit).predict = grid(vit).predict2;
	}
	if (parm->debug)
	{
		vector<Position3D> vset;
		vset.insert(vset.end(), vset1.begin(), vset1.end());
		vset.insert(vset.end(), vset0.begin(), vset0.end());
		string fname = parm->_v_vacant_file;
		fname.insert(parm->_v_vacant_file.find_last_of('.'), "_middle_" + to_string(n) + "_predicted");
		Output_Vacant(fname, vset);
	}
	return total;
}

void Pocket::Backup_Adjust(int n)
{
	int i, j, k;
	for (auto& vit : vacantset[n])
	{
		i = vit.i, j = vit.j, k = vit.k;
		grid.grid[i][j][k].inner_bak = grid.grid[i][j][k].inner;
		grid.grid[i][j][k].re_bak = grid.grid[i][j][k].re;
		grid.grid[i][j][k].type_bak = grid.grid[i][j][k].type;
	}
	return;
}

void Pocket::Adjust(int n, int edge, int shrink, int* re)
// 边界调整：目标是口袋中type='S',inner=0的格点，edge>0,则将边界往外扩，edge<0则往里缩
// 收缩调整：目标是口袋中type='V',inner=0的格点，shrink<0, 则往里收缩边界，shrink>0则往外扩边界
// re[0]=0: 仅处理上下左右前后的邻居格点； re[0]=1: 处理3*3空间内的邻居格点
{
	if (vacantset.count(n) == 0)
	{
		cout << "No valid vacant grid." << endl;
		return;
	}
	extern Parameter* parm;
	// Edge Adjust 0.5A * edge
	int flag;
	int i, j, k;
	int it, jt, kt;
	int neib;

	for (flag = 0; flag < abs(edge); flag++)
	{
		for (i = 0; i < grid.dx; i++)
			for (j = 0; j < grid.dy; j++)
				for (k = 0; k < grid.dz; k++)
				{
					if (grid.grid[i][j][k].type == 'S' && grid.grid[i][j][k].cavity == n && grid.grid[i][j][k].inner == 0)
					{
						if (edge > 0)
						{
							if (re[0] == 1 && grid.grid[i][j][k].re < (edge - flag))continue;
							grid.grid[i][j][k].inner = 1;
						}
						else
						{
							grid.grid[i][j][k].inner = -1;
							grid.grid[i][j][k].cavity = -1;
							grid.grid[i][j][k].re = re[n] + 1;
						}
						for (it = i - 1; it < i + 2; it++)
							for (jt = j - 1; jt < j + 2; jt++)
								for (kt = k - 1; kt < k + 2; kt++)
								{
									if (it == i && jt == j && kt == k) continue;
									neib = abs(it - i) + abs(jt - j) + abs(kt - k);
									if (grid.OutOfBoundary(it, jt, kt) == 0)
									{
										if (edge > 0 && grid.grid[it][jt][kt].cavity == -1 && grid.grid[it][jt][kt].type == 'S')
										{
											if (neib < parm->out_slop)
											{
												grid.grid[it][jt][kt].inner = 2;
												grid.grid[it][jt][kt].cavity = n;
											}
										}
										if (edge < 0 && grid.grid[it][jt][kt].cavity == n && grid.grid[it][jt][kt].type == 'S' && grid.grid[it][jt][kt].inner != 0)
										{
											if (neib < parm->in_slop)
											{
												grid.grid[it][jt][kt].inner = 2;
												grid.grid[it][jt][kt].re = re[n] + 2;
											}
										}
									}
								}
					}
				}

		for (i = 0; i < grid.dx; i++)
for (j = 0; j < grid.dy; j++)
	for (k = 0; k < grid.dz; k++)
	{
		if (grid.grid[i][j][k].inner == 2)grid.grid[i][j][k].inner = 0;
	}
	}
	// Shrink Adjust 0.5A * shrink
	vector<Position3D> addset;
	for (flag = 0; flag < abs(shrink); flag++)
	{
		addset.clear();
		for (auto& vit : vacantset[n])
		{
			i = vit.i, j = vit.j, k = vit.k;
			if (grid.grid[i][j][k].inner == 0)
			{
				if (shrink > 0)
				{
					grid.grid[i][j][k].inner = -1;
				}
				else
				{
					grid.grid[i][j][k].inner = -1;
					grid.grid[i][j][k].cavity = -1;
					grid.grid[i][j][k].re = re[n] + 1;
					grid.grid[i][j][k].type = 'O';
				}
				for (it = i - 1; it < i + 2; it++)
					for (jt = j - 1; jt < j + 2; jt++)
						for (kt = k - 1; kt < k + 2; kt++)
						{
							neib = abs(it - i) + abs(jt - j) + abs(kt - k);
							if (re[0] != 1 && neib != 1)continue;
							if (grid.OutOfBoundary(it, jt, kt) == 0)
							{
								if (shrink > 0 && grid.grid[it][jt][kt].type == 'O')
								{
									if (re[0] == 1 && grid.grid[i][j][k].re > grid.grid[it][jt][kt].re + 1)continue;
									grid.grid[it][jt][kt].type = 'V';
									grid.grid[it][jt][kt].inner = 2;
									grid.grid[it][jt][kt].cavity = n;
									addset.push_back(Position3D(it, jt, kt));
								}
								if (shrink < 0 && grid.grid[it][jt][kt].type == 'V' && grid.grid[it][jt][kt].cavity == n && grid.grid[it][jt][kt].inner != 0)
								{
									if (re[0] == 1 && grid.grid[i][j][k].dep + 1 < grid.grid[it][jt][kt].dep)continue;
									grid.grid[it][jt][kt].inner = 2;
									grid.grid[it][jt][kt].re = re[n] + 2;
								}
							}
						}
			}
		}
		if (shrink > 0)
		{
			for (auto& ait : addset)
			{
				grid(ait).inner = 0;
				vacantset[n].push_back(ait);
			}
		}
		else if (shrink < 0)
		{
			for (auto& vit : vacantset[n])
			{
				if (grid(vit).inner == 2) 
					grid(vit).inner = 0;
			}
		}
	}
	vector<Position3D> vlist;
	if (shrink < 0)
	{
		for (auto& vit : vacantset[n])
		{
			if (grid(vit).cavity == n) 
				vlist.push_back(vit);
		}
		vacantset[n] = vlist;
	}
	
	return;
}

void Pocket::Restore_Adjust(int n)
{
	int i, j, k;
	for (auto& vit : vacantset[n])
	{
		i = vit.i, j = vit.j, k = vit.k;
		grid.grid[i][j][k].cavity = n;
		grid.grid[i][j][k].inner = grid.grid[i][j][k].inner_bak;
		grid.grid[i][j][k].re = grid.grid[i][j][k].re_bak;
		grid.grid[i][j][k].type = grid.grid[i][j][k].type_bak;
	}
	return;
}

void Pocket::Core_Adjust(int n, Position3D& core_s, Position3D& core_v)
// 把口袋的边界整体缩了一圈,一直缩到最后，获取剩下的最后一个点
{
	extern Parameter* parm;
	int i, j, k, it, jt, kt;
	int neib;
	int s_valid = 0, v_valid = 0;
	vector<Position3D> sset0;
	for (auto& sit : surfaceset[n])
	{
		grid(sit).inner2 = grid(sit).inner;
		if (grid(sit).inner2 == 0) sset0.push_back(sit);
	}
	for (auto& vit : vacantset[n])
	{
		if (grid(vit).inner == 0) grid(vit).inner2 = 1;
	}
	vector<Position3D> sset2; // protein surface edge
	vector<Position3D> vset0; // vacant surface edge
	if (sset0.size() != 0) s_valid = 1;
	for (auto& sit : sset0)
	{
		i = sit.i, j = sit.j, k = sit.k;
		grid(sit).inner2 = -1;
		core_s = sit;
		for (it = i - 1; it < i + 2; it++)
			for (jt = j - 1; jt < j + 2; jt++)
				for (kt = k - 1; kt < k + 2; kt++)
				{
					neib = abs(it - i) + abs(jt - j) + abs(kt - k);
					if (grid.OutOfBoundary(it, jt, kt) == 0 && neib != 0)
					{
						if (grid.grid[it][jt][kt].cavity == n && grid.grid[it][jt][kt].inner2 == 1)
						{
							if (grid.grid[it][jt][kt].type == 'V')
							{
								grid.grid[it][jt][kt].inner2 = 0;
								vset0.push_back(Position3D(it, jt, kt));
							}
							else if (grid.grid[it][jt][kt].type == 'S' && neib < parm->in_slop)
							{
								grid.grid[it][jt][kt].inner2 = 2;
								sset2.push_back(Position3D(it, jt, kt));
							}
						}
					}
				}
	}
	// find protein surface core 
	while (sset2.size() != 0)
	{
		sset0.clear();
		for (auto& sit : sset2)
		{
			grid(sit).inner2 = 0;
			sset0.push_back(sit);
		}
		sset2.clear();
		for (auto& sit : sset0)
		{
			i = sit.i, j = sit.j, k = sit.k;
			grid(sit).inner2 = -1;
			core_s = sit;
			for (it = i - 1; it < i + 2; it++)
				for (jt = j - 1; jt < j + 2; jt++)
					for (kt = k - 1; kt < k + 2; kt++)
					{
						neib = abs(it - i) + abs(jt - j) + abs(kt - k);
						if (grid.OutOfBoundary(it, jt, kt) == 0 && neib != 0)
						{
							if (grid.grid[it][jt][kt].cavity == n && grid.grid[it][jt][kt].inner2 == 1 && grid.grid[it][jt][kt].type == 'S')
							{
								if (neib < parm->in_slop)
								{
									grid.grid[it][jt][kt].inner2 = 2;
									sset2.push_back(Position3D(it, jt, kt));
								}
							}
						}
					}
		}
	}
	// find vacant surface core
	vector<Position3D> newvset0;
	if (vset0.size() != 0) v_valid = 1;
	while (vset0.size() != 0)
	{
		newvset0.clear();
		newvset0.insert(newvset0.end(), vset0.begin(), vset0.end());
		vset0.clear();
		for (auto& vit : newvset0)
		{
			i = vit.i, j = vit.j, k = vit.k;
			grid(vit).inner = -1;
			core_v = vit;
			for (it = i - 1; it < i + 2; it++)
				for (jt = j - 1; jt < j + 2; jt++)
					for (kt = k - 1; kt < k + 2; kt++)
					{
						neib = abs(it - i) + abs(jt - j) + abs(kt - k);
						if (grid.OutOfBoundary(it, jt, kt) == 0 && neib != 0)
						{
							if (grid.grid[it][jt][kt].cavity == n && grid.grid[it][jt][kt].inner2 == 1 && grid.grid[it][jt][kt].type == 'V')
							{
								grid.grid[it][jt][kt].inner2 = 0;
								vset0.push_back(Position3D(it, jt, kt));
							}
						}
					}
		}
	}
	if (s_valid == 1) grid(core_s).inner2 = 8;
	if (v_valid == 1) grid(core_v).inner2 = 8;
	return;
}

void Pocket::Backup(int ruler)
// ruler = 0: backup 'S' grid
// ruler = others: backup 'V' grid
{
	if (ruler == 0)
	{
		for (int i = 0; i < grid.dx; i++)
			for (int j = 0; j < grid.dy; j++)
				for (int k = 0; k < grid.dz; k++)
				{
					if (grid.grid[i][j][k].type == 'S')
					{
						bgrid.push_back(bGrid(i, j, k, grid.grid[i][j][k].cavity));
					}
					grid.grid[i][j][k].inner_bak = grid.grid[i][j][k].inner;
				}
	}
	else
	{
		for (int i = 0; i < grid.dx; i++)
			for (int j = 0; j < grid.dy; j++)
				for (int k = 0; k < grid.dz; k++)
				{
					grid.grid[i][j][k].inner_bak = grid.grid[i][j][k].inner;
				}
	}
	return;
}

void Pocket::Restore(int n, int ruler)
// ruler = 0  bGrid备份的是'S' grid
// ruler = other bGrid备份的是'V' grid
// 该函数将所有格点的cavity初始化为-1,然后将口袋n的格点的cavity恢复为n，
// 对于ruler!=0, 还会将除了口袋n内'V'格点以外的'V'格点变为'O'格点。
{
	if (ruler == 0)
	{
		for (int i = 0; i < grid.dx; i++)
			for (int j = 0; j < grid.dy; j++)
				for (int k = 0; k < grid.dz; k++)
				{
					grid.grid[i][j][k].inner = grid.grid[i][j][k].inner_bak;
					grid.grid[i][j][k].cavity = -1;
				}
		for (auto& it : bgrid)
		{
			if (it.cavity != n)continue;
			else
			{
				grid.grid[it.i][it.j][it.k].cavity = n;
			}
		}
	}
	else
	{
		for (int i = 0; i < grid.dx; i++)
			for (int j = 0; j < grid.dy; j++)
				for (int k = 0; k < grid.dz; k++)
				{
					grid.grid[i][j][k].inner = grid.grid[i][j][k].inner_bak;
					grid.grid[i][j][k].cavity = -1;
					if (grid.grid[i][j][k].type == 'V') grid.grid[i][j][k].type = 'O';
				}
		for (auto& vit : vacantset[n])
		{
			grid(vit).cavity = n;
			grid(vit).type = 'V';
		}
	}
	return;
}

int Pocket::Choose(int seed, int ruler)
// 选择条件：
// ruler = 0: totalnum > ruler_1
// ruler = 1: total != 0
// ruler = 2: vacant_num > ruler_1
{
	extern Parameter* parm;
	if (ruler == 0)
	{
		if (cavity[seed].surfacenum > parm->ruler_1) return 1;
	}
	if (ruler == 1)
	{
		if (cavity[seed].surfacenum != 0) return 1;
	}
	if (ruler == 2)
	{
		if (cavity[seed].vacant_num > parm->ruler_1) return 1;
	}
	return 0;
}

void Pocket::Cavity_Info(int seed, int ruler)
// ruler = 9: 找surface
// ruler = 0: 统计surface，并赋值给cavity[seed].totalnum 
// ruler = 10: 找到口袋cavity[seed].id的盒子边界，并将当前口袋的totalnum和vacant_sum加到cavity[0]上
{
	if (ruler == 0) // count cavity point num for prechoose
	{
		if (surfaceset.count(cavity[seed].id) == 0) return;
		cavity[seed].surface = 0;
		for (auto& sit : surfaceset[cavity[seed].id])
		{
			cavity[seed].surface += grid(sit).grid_surface;
		}
		cavity[seed].surfacenum = (int)(cavity[seed].surface);
		return;
	}
	if (ruler == 3) // get dep1 = edge/surface dep2 = area/vacant
	{
		if (cavity[seed].valid != 1) return;
		for (auto& vit : vacantset[cavity[seed].id])
		{
			if (grid(vit).inner == 0) cavity[seed].area++;
		}
		for (auto& sit : surfaceset[cavity[seed].id])
		{
			if (grid(sit).inner == 0) { cavity[seed].edge += grid(sit).grid_surface; }
		}
		cavity[seed].dep1 = cavity[seed].edge / (float)(cavity[seed].surfacenum);
		cavity[seed].dep2 = cavity[seed].area / (float)(cavity[seed].vacant_num);
		return;
	}
	if (ruler == 4) // Mark cavity atoms
	{
		extern Protein* prot;
		extern Parameter* parm;
		float r = parm->distance * parm->distance;
		for (auto& atm : prot->atoms)
		{
			if (atm.valid == 2) atm.valid = 1;
		}
		float x, y, z, sum, tmp;
		int i, j, k;
		vector<Position3D> cset;
		
		for (auto& vit : vacantset[cavity[seed].id])
		{
			if (grid(vit).inner == 0) { cset.push_back(vit); }
		}
		cset.insert(cset.end(), surfaceset[cavity[seed].id].begin(), surfaceset[cavity[seed].id].end());
		map<int, map<int, vector<int>>> newcset;
		for (auto& vit : cset)
		{
			if (newcset.count(vit.i) == 0)
			{
				map<int, vector<int>> m;
				vector<int> tmp;
				tmp.push_back(vit.k);
				m.insert(pair<int, vector<int>>(vit.j, tmp));
				newcset.insert(pair<int, map<int, vector<int>>>(vit.i, m));
			}
			else
			{
				if (newcset[vit.i].count(vit.j) == 0)
				{
					vector<int> tmp;
					tmp.push_back(vit.k);
					newcset[vit.i].insert(pair<int, vector<int>>(vit.j, tmp));
				}
				else
				{
					newcset[vit.i][vit.j].push_back(vit.k);
				}
			}
		}
		vector<Patom> validatm;
		float min_x, min_y, min_z, max_x, max_y, max_z;
		min_x = cavity[seed].min_x - parm->distance;
		max_x = cavity[seed].max_x + parm->distance;
		min_y = cavity[seed].min_y - parm->distance;
		max_y = cavity[seed].max_y + parm->distance;
		min_z = cavity[seed].min_z - parm->distance;
		max_z = cavity[seed].max_z + parm->distance;
		for (auto& atm : prot->atoms)
		{
			x = atm.coord.x, y = atm.coord.y, z = atm.coord.z;
			if ((x < min_x || x > max_x) && (y <  min_y || y >max_y) && (z < min_z || z >max_z)) atm.valid = 0;
			else atm.valid = 1;
		}
		int mark = 0;
		for (auto& atm : prot->atoms)
		{
			if (atm.valid == 0) continue;
			x = atm.coord.x, y = atm.coord.y, z = atm.coord.z;
			if ((x < min_x || x > max_x) && (y <  min_y || y >max_y) && (z < min_z || z >max_z)) continue;
			mark = 0;
			for (auto& pi : newcset)
			{
				i = pi.first;
				x = atm.coord.x - grid.grid[i][0][0].coord.x;
				if (x > parm->distance) { continue; }
				else if (x < -parm->distance) { continue; }
				for (auto& pj : pi.second)
				{
					j = pj.first;
					y = atm.coord.y - grid.grid[i][j][0].coord.y;
					if (y > parm->distance) { continue; }
					else if (y < -parm->distance) { continue; }
					tmp = x * x + y * y;
					if (tmp > r) continue;
					for (auto& pk : pj.second)
					{
						k = pk;
						z = atm.coord.z - grid.grid[i][j][k].coord.z;
						if (z > parm->distance) continue;
						else if (z < -parm->distance) continue;
						sum = tmp + z * z;
						if (sum <= r)
						{
							atm.valid = 2;
							if (cavity[seed].atom_max_x < atm.coord.x) cavity[seed].atom_max_x = atm.coord.x;
							if (cavity[seed].atom_min_x > atm.coord.x) cavity[seed].atom_min_x = atm.coord.x;
							if (cavity[seed].atom_max_y < atm.coord.y) cavity[seed].atom_max_y = atm.coord.y;
							if (cavity[seed].atom_min_y > atm.coord.y) cavity[seed].atom_min_y = atm.coord.y;
							if (cavity[seed].atom_max_z < atm.coord.z) cavity[seed].atom_max_z = atm.coord.z;
							if (cavity[seed].atom_min_z > atm.coord.z) cavity[seed].atom_min_z = atm.coord.z;
							mark = 1;
							break;
						}
					}
					if (mark == 1) break;
				}
				if (mark == 1) break;
			}
		}
		return;
	}
	if (ruler == 5) // Mark bind surface and vacant
	{
		if (cavity[seed].valid != 1) return;
		int flag;
		for (auto& vit : vacantset[cavity[seed].id])
		{
			flag = 0;
			if (grid(vit).hydrophobic_score > 0.001) { cavity[seed].hydrophobic_vacant++; flag = 1; }
			if (grid(vit).acceptor_score > 0.001) { cavity[seed].acceptor_vacant++; flag = 1; }
			if (grid(vit).donor_score > 0.001) { cavity[seed].donor_vacant++; flag = 1; }
			if (flag == 1) cavity[seed].bind_vacant++;
			if (grid(vit).adh == 'D') cavity[seed].donor_v++;
			else if (grid(vit).adh == 'A') cavity[seed].acceptor_v++;
			else if (grid(vit).adh == 'H') cavity[seed].hydrophobic_v++;
		}
		for (auto& sit : surfaceset[cavity[seed].id])
		{
			flag = 0;
			if (grid(sit).hydrophobic_score > 0.001) { cavity[seed].hydrophobic_surface += grid(sit).grid_surface; flag = 1; }
			if (grid(sit).acceptor_score > 0.001) { cavity[seed].acceptor_surface += grid(sit).grid_surface; flag = 1; }
			if (grid(sit).donor_score > 0.001) { cavity[seed].donor_surface += grid(sit).grid_surface; flag = 1; }
			if (flag == 1) cavity[seed].bind_surface += grid(sit).grid_surface;
			if (grid(sit).adh == 'D') cavity[seed].donor_s += grid(sit).grid_surface;
			else if (grid(sit).adh == 'A') { cavity[seed].acceptor_s += grid(sit).grid_surface; }
			else if (grid(sit).adh == 'H') cavity[seed].hydrophobic_s+= grid(sit).grid_surface;
		}
		return;
	}
	if (ruler == 6)
	{
		if (cavity[seed].valid != 1)return;
		cavity[seed].surface_dep = Get_Surface_Dep(cavity[seed].id);
		cavity[seed].vacant_dep = Get_Vacant_Dep(seed);
		return;
	}
	if (ruler == 7) // get cavity core
	{
		extern Protein* prot;
		Position3D core_s, core_v;
		int area_mark = 0;
		if (cavity[seed].valid != 1)return;
		Core_Adjust(cavity[seed].id, core_s, core_v); // 寻找口袋内表面和外表面的中心点
		cavity[seed].area_core = core_v;
		cavity[seed].surface_core = core_s;
		
		if (grid(core_v).inner2 != 8) core_v = core_s;
		cavity[seed].ab_vector = grid(core_v).coord - grid(core_s).coord; // 内表面指向外表面的向量
		cavity[seed].ab_vector_norm = (float)cavity[seed].ab_vector.norm();
		cavity[seed].core_vector = grid(core_s).coord - prot->center_coord; // 蛋白中心指向内表面的向量
		cavity[seed].core_vector_norm = (float)cavity[seed].core_vector.norm();
		return;
	}
	if (ruler == 9) // get surface
	{
		int i, j, k;
		int flag;
		if (cavity[seed].valid != 1) return;
		vector<Position3D> slist;
		for (auto& sit : surfacegrid)
		{
			flag = 0;
			i = sit.i, j = sit.j, k = sit.k;
			if (i > 0 && grid.grid[i - 1][j][k].type == 'V' && cavity[seed].id == grid.grid[i - 1][j][k].cavity) flag = 1;
			else if (i < grid.dx - 1 && grid.grid[i + 1][j][k].type == 'V' && cavity[seed].id == grid.grid[i + 1][j][k].cavity) flag = 1;
			else if (j > 0 && grid.grid[i][j - 1][k].type == 'V' && cavity[seed].id == grid.grid[i][j - 1][k].cavity) flag = 1;
			else if (j < grid.dy - 1 && grid.grid[i][j + 1][k].type == 'V' && cavity[seed].id == grid.grid[i][j + 1][k].cavity) flag = 1;
			else if (k > 0 && grid.grid[i][j][k - 1].type == 'V' && cavity[seed].id == grid.grid[i][j][k - 1].cavity) flag = 1;
			else if (k < grid.dz - 1 && grid.grid[i][j][k + 1].type == 'V' && cavity[seed].id == grid.grid[i][j][k + 1].cavity) flag = 1;
			else continue;
			if (flag == 1)
			{
				grid.grid[i][j][k].cavity = cavity[seed].id;
				slist.push_back(sit);
			}
		}
		surfaceset.insert(pair<int, vector<Position3D>>(cavity[seed].id, slist));
		return;
	}
	if (ruler == 10) // Get cavity box
	{
		if (cavity[seed].valid != 1) return;
		for (auto& vit : vacantset[cavity[seed].id])
		{
			if (grid(vit).coord.x > cavity[seed].max_x) cavity[seed].max_x = grid(vit).coord.x;
			if (grid(vit).coord.x < cavity[seed].min_x) cavity[seed].min_x = grid(vit).coord.x;
			if (grid(vit).coord.y > cavity[seed].max_y) cavity[seed].max_y = grid(vit).coord.y;
			if (grid(vit).coord.y < cavity[seed].min_y) cavity[seed].min_y = grid(vit).coord.y;
			if (grid(vit).coord.z > cavity[seed].max_z) cavity[seed].max_z = grid(vit).coord.z;
			if (grid(vit).coord.z < cavity[seed].min_z) cavity[seed].min_z = grid(vit).coord.z;
		}
		for (auto& sit : surfaceset[cavity[seed].id])
		{
			if (grid(sit).coord.x > cavity[seed].max_x) cavity[seed].max_x = grid(sit).coord.x;
			if (grid(sit).coord.x < cavity[seed].min_x) cavity[seed].min_x = grid(sit).coord.x;
			if (grid(sit).coord.y > cavity[seed].max_y) cavity[seed].max_y = grid(sit).coord.y;
			if (grid(sit).coord.y < cavity[seed].min_y) cavity[seed].min_y = grid(sit).coord.y;
			if (grid(sit).coord.z > cavity[seed].max_z) cavity[seed].max_z = grid(sit).coord.z;
			if (grid(sit).coord.z < cavity[seed].min_z) cavity[seed].min_z = grid(sit).coord.z;
		}

		cavity[0].vacant_num += cavity[seed].vacant_num;
		cavity[0].surface += cavity[seed].surfacenum;
		return;
	}
}

int Pocket::Get_Surface_Dep(int n)
// 获取surface表面的深度
{
	extern Parameter* parm;
	int dep, neib;
	int i, j, k, it, jt, kt;
	// initialize 
	vector<Position3D> validset;
	for (auto& sit : surfaceset[n])
	{
		grid(sit).inner2 = grid(sit).inner;
		if (grid(sit).inner == 0) validset.push_back(sit);
	}
	vector<Position3D> validset2;
	for (dep = 1;; dep++)
	{
		for (auto& vit : validset)
		{
			i = vit.i; j = vit.j; k = vit.k;
			grid.grid[i][j][k].inner = -1;
			for (it = i - 1; it < i + 2; it++)
				for (jt = j - 1; jt < j + 2; jt++)
					for (kt = k - 1; kt < k + 2; kt++)
					{
						neib = abs(it - i) + abs(jt - j) + abs(kt - k);
						if (grid.OutOfBoundary(it, jt, kt) == 0)
						{
							if (grid.grid[it][jt][kt].type == 'S' && grid.grid[it][jt][kt].cavity == n && grid.grid[it][jt][kt].inner == 1)
							{
								if (neib < parm->in_slop)
								{
									grid.grid[it][jt][kt].inner = 2;
									validset2.push_back(Position3D(it, jt, kt));
								}
							}
						}
					}
		}
		validset.clear();
		for (auto& vit : validset2)
		{
			grid.grid[vit.i][vit.j][vit.k].inner = 0;
			validset.push_back(vit);
		}
		validset2.clear();
		if (validset.size() == 0)break;
	}
	//Reedge(); // 恢复边界
	for (auto& sit : surfaceset[n])
	{
		grid(sit).inner = grid(sit).inner2;
	}
	return dep;
}

int Pocket::Get_Vacant_Dep(int n)
{
	extern Parameter* parm;
	int max = 0, min = 300;
	int depth[300] = { 0 };
	for (auto& vit : vacantset[cavity[n].id])
	{
		if (grid(vit).re == 0) grid(vit).dep = 0;
		else grid(vit).dep = grid(vit).re - 1;
		depth[grid(vit).dep]++;
		if (grid(vit).dep > max) max = grid(vit).dep;
		if (grid(vit).dep < min) min = grid(vit).dep;
	}
	int limit = min;
	//float deptmp, flagtmp;
	for (int i = min; i <= max; i++)
	{
		if (i < 300)
		{
			cavity[n].step_vacant[i] = depth[i];
			if (depth[i] > 0)
			{
				limit++;
			}
		}
	}
	return (limit - min);
}

void Pocket::Prepare_Output(int seed, int k)
{
	extern Parameter* parm;
	int tmpsurface = 0, tmpvacant = 0;
	if (parm->visual_output > 0)
	{
		Prepare_Out_Surface(seed, k);
		Prepare_Out_Vacant(seed, k);
	}
	if (parm->pdb_output == 2) 	Prepare_Out_Cavity(seed, k);
}

void Pocket::Prepare_Out_Surface(int seed, int k)
// seed: 当前口袋的id
// k: 输出文件中，POK字段之后的数值，暂时意义不明，输入值为口袋的编号
// mod：控制输出格点的密度，mod=1会输出每个格点，否则，每mod个格点中取1个格点输出
// vmod： 控制输出格点，vmod=1，每隔固定距离输出一个格点，vmod=0, 随机确定当前格点是否输出。需与mod一起使用
{
	extern Parameter* parm;
	extern Ligand* lig;
	int mod = parm->v_num;
	int vmod = parm->v_mod;
	stringstream outfp;
	int n = cavity[seed].id;
	outfp << "HEADER    User defined\n";
	outfp << "COMPND    Visual Surface\n";
	outfp << "AUTHOR    Generated by " << parm->version << "\n";
	outfp << "REMARK   1\n";
	outfp << "REMARK   1 Creation time " << Get_Time();
	outfp << "REMARK   1\n";
	outfp << "REMARK   2 Carbon atoms represent surface" << "\n";
	outfp << "REMARK   2" << "\n";
	
	outfp << "REMARK   3 Area :" << "\n";
	outfp << "REMARK   3 MIN\tX\tY\tZ" << "\n";
	outfp << "REMARK   3    \t" << fixed << setprecision(2) << setw(4) << cavity[seed].min_x << "\t"
		<< setw(4) << setprecision(2) << cavity[seed].min_y << "\t"
		<< setw(4) << setprecision(2) << cavity[seed].min_z << "\n";
	outfp << "REMARK   3 MAX\tX\tY\tZ" << "\n";
	outfp << "REMARK   3    \t" << setw(4) << setprecision(2) << cavity[seed].max_x << "\t"
		<< setw(4) << setprecision(2) << cavity[seed].max_y << "\t"
		<< setw(4) << setprecision(2) << cavity[seed].max_z << "\n";
	outfp << "REMARK   3 " << "\n";
	outfp << "REMARK   4 " << "\n";
	outfp << "REMARK   4 Total Surface num is " << cavity[seed].surfacenum << "\n";
	outfp << "REMARK   4 Total Surface area is " << fixed << (float)(cavity[seed].surfacenum * 0.25) << " A^2" << "\n";
	outfp << "REMARK   4 " << "\n";
	outfp << "REMARK   4 Edge  = " << (int)cavity[seed].edge << "\n";
	outfp << "REMARK   4 Hydrophobic  = " << (int)cavity[seed].hydrophobic_surface << "\n";
	outfp << "REMARK   4 Acceptor  = " << (int)cavity[seed].acceptor_surface << "\n";
	outfp << "REMARK   4 Donor  = " << (int)cavity[seed].donor_surface << "\n";
	outfp << "REMARK   4 Bind  = " << (int)cavity[seed].bind_surface << "\n";
	outfp << "REMARK   4 Hydrophobic2  = " << (int)cavity[seed].hydrophobic_s << "\n";
	outfp << "REMARK   4 Acceptor2  = " << (int)cavity[seed].acceptor_s << "\n";
	outfp << "REMARK   4 Donor2  = " << (int)cavity[seed].donor_s << "\n";
	outfp << "REMARK   4 Dep1 = " << cavity[seed].surface_dep << "\n";
	outfp << "REMARK   4 " << "\n";
	if (parm->detect_mode == 1)
	{
		outfp << "REMARK   5 Percent " << setprecision(8) << (float)(cavity[seed].percent / lig->num_heavy) << "\n";
		outfp << "REMARK   5 Lig_Match " << setprecision(6) << (float)cavity[seed].ligmatch * 0.125 << " / " << (float)cavity[seed].ligall * 0.125
			<< " " << (float)cavity[seed].ligmatch / (float)cavity[seed].vacant_num << "\n";
		outfp << "REMARK   5 Pecision " << (float)cavity[seed].ligmatch / (float)(cavity[seed].ligall + cavity[seed].vacant_num - cavity[seed].ligmatch) << "\n";
		outfp << "REMARK   5 Lig_Num " << lig->num_atom << "\n";
		if (parm->single_ligand == 1) outfp << "REMARK   5 Lig_Volume " << lig->volume << "\n";
	}
	outfp << "REMARK   5 Round " << cavity[seed].re << "\n";
	outfp << "REMARK   5 Step "<< cavity[seed].step_vacant[1];
	for (int i = 2; i < cavity[seed].vacant_dep; i++)
	{
		outfp << "_" << cavity[seed].step_vacant[i];
	}
	outfp << "\n";
	outfp << "REMARK   5 RankScore " << setprecision(6) << cavity[seed].rank << "\n";
	if (cavity[seed].area <= 50)
	{
		outfp << "REMARK   5 This cavity is closed, which may make large deviation in following prediction!" << "\n";
	}
	if (cavity[seed].rank <= 5.18)
	{
		outfp << "REMARK   5 Predict Maximal pKd: " << setw(5) << setprecision(2) << (cavity[seed].rank * 1.796 + 2.719) << "\n";
	}
	else
	{
		outfp << "REMARK   5 Predict Maximal pKd: " << setw(5) << setprecision(2) << (cavity[seed].rank * -1.312 + 18.73) << "\n";
	}
	if (cavity[seed].rank <= 5.6)
	{
		outfp << "REMARK   5 Predict Average pKd: " << setw(5) << setprecision(2) << (cavity[seed].rank * 0.615 + 3.553) << "\n";
	}
	else
	{
		outfp << "REMARK   5 Predict Average pKd: " << setw(5) << setprecision(2) << (cavity[seed].rank * -0.086 + 7.478) << " (Large deviation)" << "\n";
	}
	outfp << "REMARK   6 DrugScore: " << fixed << setw(5) << setprecision(2) << cavity[seed].drugscore << "\n";
	outfp << "REMARK   6 Druggability: ";
	if (cavity[seed].drugscore >= 600) outfp << "Strong" << "\n";
	else if (cavity[seed].drugscore >= -180) outfp << "Medium" << "\n";
	else outfp << "Weak" << "\n";
	outfp << "REMARK   6" << "\n";
	outfp << "REMARK   7" << "\n";

	int count = 1;
	int tmp = -1;
	char buff[200];
	for (auto& sit : surfaceset[n])
	{
		tmp++;
		if (mod != 1)
		{
			if (vmod == 1)
			{
				if ((tmp % mod) != 0)continue;
			}
			if (vmod == 0)
			{
				if ((tmp % mod) != (rand() % mod)) continue;
			}
		}
		snprintf(buff, 200, "HETATM%6d %-3c POK   %-5d   %7.3f %7.3f %7.3f  1.00\n", count, 'C', k, grid(sit).coord.x, grid(sit).coord.y, grid(sit).coord.z);
		outfp << buff;
		count++;
	}
	cavity[seed].slines.clear();
	string s;
	while (getline(outfp, s))
	{
		cavity[seed].slines.push_back(s);
	}
}

void Pocket::Prepare_Out_Vacant(int seed, int k)
{
	extern Parameter* parm;
	extern Ligand* lig;
	int mod = parm->v_num;
	int vmod = parm->v_mod;
	stringstream outfp;
	int n = cavity[seed].id;
	outfp << "HEADER    User defined" << "\n";
	outfp << "COMPND    Visual Vacant" << "\n";
	outfp << "AUTHOR    Generated by " << parm->version << "\n";
	outfp << "REMARK   1" << "\n";
	outfp << "REMARK   1 Creation time " << Get_Time();
	outfp << "REMARK   1" << "\n";
	outfp << "REMARK   2 Oxygen atoms represent vacant" << "\n";
	outfp << "REMARK   2" << "\n";
	outfp << "REMARK   3 Area :" << "\n";
	outfp << "REMARK   3 MIN\tX\tY\tZ" << "\n";
	outfp << "REMARK   3    \t" << setw(4) << fixed << setprecision(2) << cavity[seed].min_x << "\t"
		<< setw(4) << setprecision(2) << cavity[seed].min_y << "\t"
		<< setw(4) << setprecision(2) << cavity[seed].min_z << "\n";
	outfp << "REMARK   3 MAX\tX\tY\tZ" << "\n";
	outfp << "REMARK   3    \t" << setw(4) << setprecision(2) << cavity[seed].max_x << "\t"
		<< setw(4) << setprecision(2) << cavity[seed].max_y << "\t"
		<< setw(4) << setprecision(2) << cavity[seed].max_z << "\n";
	outfp << "REMARK   3 " << "\n";
	outfp << "REMARK   4 " << "\n";
	outfp << "REMARK   4 Total Vacant num is " << cavity[seed].vacant_num << "\n";
	outfp << "REMARK   4 Total Vacant volume is " << fixed << (float)(cavity[seed].vacant_num * 0.125) << " A^3" << "\n";
	outfp << "REMARK   4 " << "\n";
	outfp << "REMARK   4 Area  = " << (int)cavity[seed].area << "\n";
	outfp << "REMARK   4 Hydrophobic  = " << (int)cavity[seed].hydrophobic_vacant << "\n";
	outfp << "REMARK   4 Acceptor  = " << (int)cavity[seed].acceptor_vacant << "\n";
	outfp << "REMARK   4 Donor  = " << (int)cavity[seed].donor_vacant << "\n";
	outfp << "REMARK   4 Bind  = " << (int)cavity[seed].bind_vacant << "\n";
	outfp << "REMARK   4 Hydrophobic2  = " << (int)cavity[seed].hydrophobic_v << "\n";
	outfp << "REMARK   4 Acceptor2  = " << (int)cavity[seed].acceptor_v << "\n";
	outfp << "REMARK   4 Donor2  = " << (int)cavity[seed].donor_v << "\n";
	outfp << "REMARK   4 Dep2 = " << cavity[seed].vacant_dep << "\n";
	outfp << "REMARK   4 " << "\n";
	if (parm->detect_mode == 1)
	{
		outfp << "REMARK   5 Percent " << setprecision(6) << (float)(cavity[seed].percent / lig->num_heavy) << "\n";
		outfp << "REMARK   5 Lig_Match " << setprecision(6) << (float)cavity[seed].ligmatch * 0.125 << " / " << (float)cavity[seed].ligall * 0.125
			<< " " << (float)cavity[seed].ligmatch / (float)cavity[seed].vacant_num << "\n";
		outfp << "REMARK   5 Pecision " << (float)cavity[seed].ligmatch / (float)(cavity[seed].ligall + cavity[seed].vacant_num - cavity[seed].ligmatch) << "\n";
		outfp << "REMARK   5 Lig_Num " << lig->num_atom << "\n";
		if (parm->single_ligand == 1) outfp << "REMARK   5 Lig_Volume " << lig->volume << "\n";
	}
	outfp << "REMARK   5 Round " << cavity[seed].re << "\n";
	outfp << "REMARK   5 Step " << cavity[seed].step_vacant[1];
	for (int i = 2; i < cavity[seed].vacant_dep; i++)
	{
		outfp << "_" << cavity[seed].step_vacant[i];
	}
	outfp << "\n";
	outfp << "REMARK   5 RankScore " << setprecision(6) << cavity[seed].rank << "\n";
	if (cavity[seed].area <= 50)
	{
		outfp << "REMARK   5 This cavity is closed, which may make large deviation in following prediction!" << "\n";
	}
	if (cavity[seed].rank <= 5.18)
	{
		outfp << "REMARK   5 Predict Maximal pKd: " << setw(5) << setprecision(2) << (cavity[seed].rank * 1.76 + 2.719) << "\n";
	}
	else
	{
		outfp << "REMARK   5 Predict Maximal pKd: " << setw(5) << setprecision(2) << (cavity[seed].rank * -1.312 + 18.73) << "\n";
	}
	if (cavity[seed].rank <= 5.6)
	{
		outfp << "REMARK   5 Predict Average pKd: " << setw(5) << setprecision(2) << (cavity[seed].rank * 0.615 + 3.553) << "\n";
	}
	else
	{
		outfp << "REMARK   5 Predict Average pKd: " << setw(5) << setprecision(2) << (cavity[seed].rank * -0.086 + 7.478) << " (Large deviation)" << "\n";
	}
	outfp << "REMARK   6 DrugScore: " << fixed << setw(5) << setprecision(2) << cavity[seed].drugscore << "\n";
	outfp << "REMARK   6 Druggability: ";
	if (cavity[seed].drugscore >= 600) outfp << "Strong" << "\n";
	else if (cavity[seed].drugscore >= -180) outfp << "Medium" << "\n";
	else outfp << "Weak" << "\n";
	outfp << "REMARK   6" << "\n";
	outfp << "REMARK   7" << "\n";

	if (parm->visual_output == 2)
	{
		int count = 1;
		int tmp = -1;
		char buff[200];
		for (auto& vit : vacantset[n])
		{
			tmp++;
			if (mod != 1)
			{
				if (vmod == 1)
				{
					if ((tmp % mod) != 0)continue;
				}
				if (vmod == 0)
				{
					if ((tmp % mod) != (rand() % mod)) continue;
				}
			}
			/*outfp << "HETATM" << right << setw(6) << count << " " << setw(3) << left << "O" << " ";
			outfp << "POK   " << setw(5) << left << k << "   ";
			outfp << right << setw(7) << setprecision(3) << grid(vit).coord.x << " "
				<< setw(7) << setprecision(3) << grid(vit).coord.y << " "
				<< setw(7) << setprecision(3) << grid(vit).coord.z << " " << " 1.00" << "\n";*/
			snprintf(buff, 200, "HETATM%6d %-3c POK   %-5d   %7.3f %7.3f %7.3f  1.00\n", count, 'O', k, grid(vit).coord.x, grid(vit).coord.y, grid(vit).coord.z);
			outfp << buff;
			count++;
		}
	}
	
	cavity[seed].vlines.clear();
	string s;
	while (getline(outfp, s))
	{
		cavity[seed].vlines.push_back(s);
	}
}

void Pocket::Prepare_Out_Cavity(int seed, int k)
{
	extern Protein* prot;
	int n = cavity[seed].id;
	Cavity_Info(seed, 4);
	stringstream outfp;
	for (auto& atm : prot->atoms)
	{
		if (atm.valid == 2) outfp << atm.pdbline << "\n";
	}
	outfp << "ENDMDL\n";
	cavity[seed].clines.clear();
	string s;
	while (getline(outfp, s))
	{
		cavity[seed].clines.push_back(s);
	}
}

void Pocket::Output_All(int seed, int k)
{
	extern Parameter* parm;
	string fname;
	ofstream outfp;
	if (parm->visual_output > 0)
	{
		fname = parm->_v_surface_file;
		fname.insert(parm->_v_surface_file.find_last_of('.'), "_" + to_string(k));
		outfp.open(fname, ios::out);
		for (auto& s : cavity[seed].slines)
		{
			outfp << s << "\n";
		}
		outfp.close();
		fname = parm->_v_vacant_file;
		fname.insert(parm->_v_vacant_file.find_last_of('.'), "_" + to_string(k));
		outfp.open(fname, ios::out);
		for (auto& s : cavity[seed].vlines)
		{
			outfp << s << "\n";
		}
		outfp.close();
	}
	if (parm->pdb_output == 2)
	{
		fname = parm->_v_atom_file;
		fname.insert(parm->_v_atom_file.find_last_of('.'), "_" + to_string(k));
		outfp.open(fname, ios::out);
		for (auto& s : cavity[seed].clines)
		{
			outfp << s << "\n";
		}
		outfp.close();
	}
}

void Pocket::Define_Ligand_Grid()
{
	extern Ligand* lig;
	int it, jt, kt;
	double d;
	for (auto& atm : lig->atoms)
	{
		if (atm.valid != 1 || atm.tripos_type == "H") continue;
		it = (int)((atm.coord.x - grid.min_x) * 2);
		jt = (int)((atm.coord.y - grid.min_y) * 2);
		kt = (int)((atm.coord.z - grid.min_z) * 2);
		for (int i = it - 5; i < it + 6; i++)
			for (int j = jt - 5; j < jt + 6; j++)
				for (int k = kt - 5; k < kt + 6; k++)
				{
					if (grid.OutOfBoundary(i, j, k) == 0)
					{
						d = Distance(grid.grid[i][j][k].coord, atm.coord);
						if (d <= atm.vdw_r * 1.2) { grid.grid[i][j][k].lpoint = 1; }
					}
				}
	}
	for (int i = 0; i < grid.dx; i++)
		for (int j = 0; j < grid.dy; j++)
			for (int k = 0; k < grid.dz; k++)
			{
				if (grid.grid[i][j][k].lpoint == 1) lpointgrid.push_back(Position3D(i, j, k));
			}
}

void Pocket::Calculate_Ligand_Percent(int seed)
{
	extern Ligand* lig;
	float it, jt, kt;
	int id, jd, kd;
	for (auto& atm : lig->atoms)
	{
		if (atm.tripos_type == "H" || atm.valid != 1) continue;
		it = ((float)((int)(atm.coord.x * 2))) / 2;
		jt = ((float)((int)(atm.coord.y * 2))) / 2;
		kt = ((float)((int)(atm.coord.z * 2))) / 2;
		id = (int)((it - grid.min_x) * 2);
		jd = (int)((jt - grid.min_y) * 2);
		kd = (int)((kt - grid.min_z) * 2);
		if (grid.OutOfBoundary(id, jd, kd) == 0)
		{
			if (grid.grid[id][jd][kd].type == 'V' || grid.grid[id][jd][kd].type == 'S')
			{
				if (cavity[seed].id == grid.grid[id][jd][kd].cavity)cavity[seed].percent++;
			}
		}
	}
	for (auto& lit : lpointgrid)
	{
		if (cavity[seed].id == grid(lit).cavity && grid(lit).type == 'V') cavity[seed].ligmatch++;
	}
	cavity[seed].ligall = (int)lpointgrid.size();
}

void Pocket::Output_Vacant(string fname, vector<Position3D>& vset)
{
	stringstream outfp;
	int count = 1;
	int tmp = -1;
	char buff[200];
	for (auto& vit : vset)
	{
		tmp++;
		snprintf(buff, 200, "HETATM%6d %-3c POK   %-5d   %7.3f %7.3f %7.3f  1.00\n", count, 'O', 0, grid(vit).coord.x, grid(vit).coord.y, grid(vit).coord.z);
		outfp << buff;
		count++;
	}
	string s;
	ofstream fp;
	fp.open(fname.c_str(), ios::out);
	while (getline(outfp, s))
	{
		fp << s << endl;
	}
	fp.close();
	outfp.clear();
}

void Pocket::Output_Middle_Vacant(int start, int total)
{
	extern Parameter* parm;
	string fname;
	
	for (int i = start; i < total; i++)
	{
		fname = parm->_v_vacant_file;
		fname.insert(parm->_v_vacant_file.find_last_of('.'), "_middle_" + to_string(i));
		stringstream outfp;
		int count = 1;
		int tmp = -1;
		char buff[200];
		for (auto& vit : vacantset[i])
		{
			tmp++;
			snprintf(buff, 200, "HETATM%6d %-3c POK   %-5d   %7.3f %7.3f %7.3f  1.00\n", count, 'O', i, grid(vit).coord.x, grid(vit).coord.y, grid(vit).coord.z);
			outfp << buff;
			count++;
		}
		string s;
		ofstream fp;
		fp.open(fname.c_str(), ios::out);
		while (getline(outfp, s))
		{
			fp << s << endl;
		}
		fp.close();
		outfp.clear();
	}
}
