#include <iostream>
#include <random>
#include <vector>
#include <utility>
#include <iomanip>
#include <cmath>
#include <array>
#include <ostream>
#include <fstream>

using namespace std;

namespace {
	ofstream take_log_type("data_desig9_type_number0_0.01.log");
	// ofstream take_log_network("data_desig9_network0_0.01.log");
	ofstream take_log_outside("data_desig9_outside0_0.01.log");
	// ofstream take_log_devdev("data_desig9_devdev0_0.01.log");
	ofstream take_log_coef("data_desig9_coef0_0.01.log");
	ofstream take_log_inside("data_desig9_inside0_0.01.log");
	ofstream take_log_inside1("data_desig9_inside1_0.01.log");
	ofstream take_log_inside2("data_desig9_inside2_0.01.log");
	ofstream take_log_inside3("data_desig9_inside3_0.01.log");
	ofstream take_log_inside4("data_desig9_inside4_0.01.log");
	ofstream take_log_inside5("data_desig9_inside5_0.01.log");
	ofstream take_log_inside6("data_desig9_inside6_0.01.log");
	ofstream take_log_inside7("data_desig9_inside7_0.01.log");
	ofstream take_log_inside8("data_desig9_inside8_0.01.log");
	ofstream take_log_inside9("data_desig9_inside9_0.01.log");
	ofstream take_log_inside10("data_desig9_inside10_0.01.log");
	ofstream take_log_growth("data_desig9_growth_0.001.log");
	// ofstream take_log_come("data_desig9_come.log");
	// ofstream take_log_up("data_desig9_up.log");
	// ofstream take_log_sum("data_desig9_sum.log");
	// ofstream take_log_boxcon("data_desig9_boxcon.log");
	// ofstream take_log_grow("data_desig9_grow.log");
	// ofstream take_log_pop("data_desig9_pop.log");
}


#define rep(i, n) for (int i = 0; i < n; i++)

double rd(void)
{
	double get;
	get = (double)(rand() % 1000) * 0.001 + 0.001;
	return get;
}

double rand_normal(double mu, double sigma)
{
	double z = sqrt(- 2.0 * log(rd())) * sin(2.0 * M_PI * rd());
	return mu + sigma * z;
}

double get_rand_normal(double size)
{
	double get;
	get = rand_normal(size, size * 0.25);
	return get;
}

const int cell_max = 1000;
int cell_type;
const int init_cell_type = 1;
// const int time_end = 5000000;
const int time_end = 5000000;
// const int time_end = 200;
const double time_bunkai = 0.001;
const int run_time = 1;

const double box_size = 1.0;
const double outside_size = 10;

const int N = 10;
// const int N = 3;

int cell_number;

double reversible[N][N];
double coef_decrease[N];

double nut_coef;
double nut_reversible;
double aver_nut;
double box_con;

typedef struct
{
	int type;
	double nut;
	int nut_cat;

	double mol[N];

	double nut_zero_coef;
	double coef[N][N];

	double go[N];

	int catalyst[N][N];

	double size;
	double init_last;
} Cell;

Cell cell[cell_max];
Cell def;

double outside_nut;
double outside[N];

double prev_size[1000];


double decide_box_nut(int time)
{
	return aver_nut;
	// return aver_nut * (1.0 + 1.0 * sin(time));
}


double get_size(Cell p)
{
	double sum = p.nut;
	rep(i, N) sum += p.mol[i];
	return sum;
}

array<array<double, N>, 10000> begin_coef;

void sum_init(int i, double max)
{
	double sum = 0;
	array<double, N> keep;
	rep(j, N) {
		keep.at(j) = rand() % rand();
		keep.at(j) *= rand();
		sum += keep.at(j);
	}
	rep(j, N) begin_coef.at(i).at(j) = max * keep.at(j) / sum;
}

void desig_init(int i)
{
	rep(j, N) cin >> begin_coef.at(i).at(j);
}

void all_init(void)
{
	rep(i, N) {
		if (i % 3 != 2) begin_coef.at(0).at(i) = 0.4;
		else begin_coef.at(0).at(i) = 2.4;
	}
}

void zero_init(int i)
{
	rep(j, N) {
		begin_coef.at(i).at(j) = 0;
	}
}

void one_init(int i)
{
	rep(j, N) {
		begin_coef.at(i).at(j) = 1;
	}
}

void ten_init(int i)
{
	rep(j, N) {
		begin_coef.at(i).at(j) = 10;
	}
}

double init_go[N];

void init(void)
{
	//defをいれておいて、typeカウントを防ぐ
	def.type = - 1;
	def.nut = 0;
	rep(i, N) def.mol[i] = 0;
	def.size = 0;
	rep(i, cell_max) cell[i] = def;

	//coefの決め方
	// all_init();
	// sum_init(0, 10);
	// one_init(0);
	ten_init(0);
	// one_init(1);
	// zero_init(0);
	// rep(i, 50) sum_init(i + 1, 1);
	// rep(i, 50) desig_init(i + 1);
	// desig_init(0);
	// desig_init(1);
	// rep(i, 50) {
	// 	rep(j, N) take_log_coef << begin_coef.at(i).at(j) << " ";
	// 	take_log_coef << endl;
	// }

	cell_type = init_cell_type;
	// cell_type = 0;
	cell_number = cell_type;

	nut_coef = 1;
	// nut_coef = 0.001;
	nut_reversible = 0;
	// aver_nut = 0.0001;
	aver_nut = 0.1;

	//outside系はloopの外
	outside_nut = aver_nut;
	rep(i,N) {
		outside[i] = 0;
		coef_decrease[i] = 0;
	}

	rep(i, cell_type) {
		cell[i].nut_zero_coef = begin_coef.at(i).at(0);
		cell[i].type = i;
		cell[i].size = box_size / cell_type;
		cell[i].nut = 1.0 / (double)(N + 1);
		// cell[i].nut = 0.25; //test用
		rep(j, N) {
			cell[i].mol[j] = 1.0 / (double)(N + 1);
			// cell[i].mol[j] = 0.25; //test用
			cell[i].nut_cat = N - 3;
			// cell[i].nut_cat = N - 2;
			// cell[i].nut_cat = 1; //test用
			cell[i].go[j] = (j % 3 / 2) * 0;
			// cell[i].go[j] = init_go[j];
			if (j != N - 1) {
				cell[i].coef[j][j + 1] = begin_coef.at(i).at(j + 1);
				cell[i].catalyst[j][j + 1] = (j - 2 + N) % N;
				// cell[i].catalyst[j][j + 1] = (j - 1 + N) % N;
				reversible[j][j + 1] = 0;
			}
		}
		// cell[i].catalyst[0][1] = 2; //test用
		// cell[i].catalyst[1][2] = 0; //test用
		cell[i].init_last = cell[i].mol[N - 1];
		// cell[i].go[0] = 1;
		// cell[i].nut = 0.5; //test
		// cell[i].mol[0] = 0; //test
		// cell[i].mol[1] = 0.5;   //test
		// cell[i].mol[2] = 0;

		// cell[i].nut_cat = 9;
		// cell[i].catalyst[0][1] = 2;
		// cell[i].catalyst[1][2] = 4;
		// cell[i].catalyst[2][3] = 6;
		// cell[i].catalyst[3][4] = 8;
		// cell[i].catalyst[4][5] = 0;
		// cell[i].catalyst[5][6] = 7;
		// cell[i].catalyst[6][7] = 5;
		// cell[i].catalyst[7][8] = 3;
		// cell[i].catalyst[8][9] = 1;
	}
}

double prev_outside_nut;
double prev_outside[N];

double mmm = 0;

double resize = 0;

Cell internal(Cell p, int t, int num)
{
	//前の値を保存
	double prev[N];
	rep(i, N) prev[i] = p.mol[i];

	double prev_nut = p.nut;
	double growth = 0;

	//細胞内外の栄養の流出入
	p.nut += time_bunkai * nut_coef/* * pow(p.size, - 1.0 / 3.0)*/ * (prev_outside_nut - prev_nut);
	outside_nut -= time_bunkai * nut_coef/* * pow(p.size, - 1.0 / 3.0)*/ * (prev_outside_nut - prev_nut) * prev_size[num] / outside_size;
	growth += time_bunkai * nut_coef * (prev_outside_nut - prev_nut);

	//リアクション
	p.mol[0] += time_bunkai * p.nut_zero_coef * prev_nut * prev[p.nut_cat];
	p.nut -= time_bunkai * p.nut_zero_coef * prev_nut * prev[p.nut_cat];

	rep(i, N - 1) {
		p.mol[i + 1] += time_bunkai * p.coef[i][i + 1] * prev[i] * prev[p.catalyst[i][i + 1]];
		p.mol[i] -= time_bunkai * p.coef[i][i + 1] * prev[i] * prev[p.catalyst[i][i + 1]];
		mmm += time_bunkai * p.coef[i][i + 1] * prev[i] * prev[p.catalyst[i][i + 1]];
		// cout << p.coef[i][i + 1] << " " << prev[i] << " " << prev[p.catalyst[i][i + 1]] << endl;
	}

	//可逆反応
	p.nut += time_bunkai * nut_reversible * prev[0];
	p.mol[0] -= time_bunkai * nut_reversible * prev[0];
	rep(i, N - 1) {
		p.mol[i] += time_bunkai * reversible[i][i + 1] * prev[i + 1];
		p.mol[i + 1] -= time_bunkai * reversible[i][i + 1] * prev[i + 1];
	}
	
	//細胞内外の溶質の流出入
	rep(i, N) {
		p.mol[i] += time_bunkai * p.go[i]/* * pow(p.size, - 1.0 / 3.0)*/ * (prev_outside[i] - prev[i]);
		outside[i] -= time_bunkai * p.go[i]/* * pow(p.size, - 1.0 / 3.0)*/ * (prev_outside[i] - prev[i]) * prev_size[num] / outside_size;
		growth += time_bunkai * p.go[i] * (prev_outside[i] - prev[i]);
	}

	//リサイズ
	double sum = p.nut;
	rep(i, N) sum += p.mol[i];
	// //test
	// // sum = 1.00007;
	// //testend
	// p.size = sum * p.size;
	// p.nut = p.nut / sum;
	// rep(i, N) p.mol[i] = p.mol[i] / sum;
	// take_log_sum << fixed << setprecision(8) << (sum - 1) / time_bunkai << endl;

	//growth調節
	// //test
	// if (num == 0) growth = 0.001;
	// if (num == 1) growth = 0.0001;
	// //test end
	resize += growth * prev_size[num];
	p.size += growth * prev_size[num];
	rep(i, N) p.mol[i] -= prev[i] * growth;
	p.nut -= prev_nut * growth;

	// if (t % 10 == 0) take_log_growth << growth / time_bunkai << " ";

	//test用 take_log_grow
	// take_log_grow << p.nut - prev_nut << " ";
	// double grow = p.nut - prev_nut;
	// rep(i, N) {
	// 	take_log_grow << p.mol[i] - prev[i] << " ";
	// 	grow += p.mol[i] - prev[i];
	// }
	// take_log_grow << grow << endl;

	if (t % 1000 == 0) { //take_log_inside << endl; の解除も忘れない
		take_log_inside << p.nut << " ";
	// rep(i, N) take_log_inside << p.mol[i] << " "; // test用 nodeが一つしかない時
	// }
		take_log_inside1 << p.mol[0] << " ";
		take_log_inside2 << p.mol[1] << " ";
		take_log_inside3 << p.mol[2] << " ";
		take_log_inside4 << p.mol[3] << " ";
		take_log_inside5 << p.mol[4] << " ";
		take_log_inside6 << p.mol[5] << " ";
		take_log_inside7 << p.mol[6] << " ";
		take_log_inside8 << p.mol[7] << " ";
		take_log_inside9 << p.mol[8] << " ";
		take_log_inside10 << p.mol[9] << " ";
	}

	return p;
}

void evolve(void)
{
	// double evolve_size = 1.0 / ((double)cell_type + 1.0);
	double evolve_size = 0.01;
	// double evolve_size = 0.5;
	// rep(i, cell_type) {
	// 	cell[i].size = cell[i].size * (1 - evolve_size);
	// }
	cell[0].size -= evolve_size; // 同時投入用
	int i = cell_type;
	cell[i].nut_zero_coef = begin_coef.at(i).at(0);
	cell[i].type = i;
	cell[i].size = evolve_size;
	cell[i].nut = 1.0 / (double)(N + 1);
	rep(j, N) {
		cell[i].mol[j] = 1.0 / (double)(N + 1);
		cell[i].nut_cat = N - 3;
		cell[i].go[j] = (j % 3 / 2) * 1;
		// cell[i].go[j] = init_go[j];
		if (j != N - 1) {
			cell[i].coef[j][j + 1] = begin_coef.at(i).at(j + 1);
			cell[i].catalyst[j][j + 1] = (j - 2 + N) % N;
			reversible[j][j + 1] = 0;
		}
	}
	cell[i].init_last = cell[i].mol[N - 1];
	// cell[i].go[0] = 1;
	cell_type++;
}

double up = 1;

void process(int t)
{
	resize = 0;

	box_con = decide_box_nut(t);
	//outside_nut = box_con; //test用 outsideを一定にする(box_conからの流入を考慮しない)

	//outside系をprevにいれる
	prev_outside_nut = outside_nut;
	rep(i, N) prev_outside[i] = outside[i];
	rep(i, cell_type) prev_size[i] = cell[i].size;

	//loop回してreaction
	rep(i, cell_type) {
		if (cell[i].size < 10e-8) {
			// cell[i].size = 0;
			// double sum = 0;
			// rep(j, cell_type) sum += cell[j].size;
			// rep(j, cell_type) cell[j].size = cell[j].size * box_size / sum;
		} else cell[i] = internal(cell[i], t, i);
	}
	
	if (t % 1000 == 0) {
		take_log_inside << endl;
	// if (1) take_log_inside << endl;
		take_log_inside1 << endl;
		take_log_inside2 << endl;
		take_log_inside3 << endl;
		take_log_inside4 << endl;
		take_log_inside5 << endl;
		take_log_inside6 << endl;
		take_log_inside7 << endl;
		take_log_inside8 << endl;
		take_log_inside9 << endl;
		take_log_inside10 << endl;
	}
	// if (t % 10 == 0) take_log_growth << endl;
	
	//if (t < 0) {
	//リサイズ
	// double sum = 0;
	// rep(i, cell_type) sum += cell[i].size;
	// rep(i, cell_type) {
	// 	cell[i].size = cell[i].size / sum * box_size;
	// 	cell[i].nut = cell[i].nut / sum * box_size;
	// 	rep(j, N) cell[i].mol[j] = cell[i].mol[j] / sum * box_size;
	// }
	// take_log_up << sum << endl;
	// up = up * sum;
	// take_log_up << up << endl;
	//}

	//リサイズの微分方程式`
	rep(i, cell_type) cell[i].size -= prev_size[i] * resize;

	//outsideの値を更新
	outside_nut += time_bunkai * (box_con - prev_outside_nut);
	//outside_nut = prev_outside_nut; //test用 outsideを一定にする(box_conからの流入を考慮しない)
	// take_log_boxcon << box_con << endl;

	if (t % 1000 == 0) {
	// if (1) {
		take_log_outside << outside_nut << " ";
		rep(i, N) {
			take_log_outside << outside[i] << " ";
		}
		take_log_outside << endl;
	}
}

int main(void)
// void main_(int get_rand)
{
	//randomの種を与える
	// srand(get_rand);
	srand(5);

	//tun_time回走らせる
	rep(l, run_time) {
		init();
		//time_end_秒走らせる
		rep(t, time_end) {
			// if (t == 0) rep(i, 1) evolve();
			process(t);
			// if (t % 50000 == 30000) evolve();
			// if (t == 100000) rep(i, 50) evolve();
			// if (t % 2000 == 1000 && cell_type < 50) rep(i, 1) evolve();
			// if (t == 60000) rep(i, 15) cell[i + 1].go[1] = 0;
			cout << t << " ";
			rep(i, cell_type) cout << fixed << setprecision(8) << cell[i].size << " ";
			cout << endl;
			if (t % 1000 == 0) {
				rep(i, cell_type) take_log_type << cell[i].size << " ";
				take_log_type << endl;
			}
		}
	}
	// cout << "mmm" << mmm << endl;
	
	rep(i, cell_type) {
		rep(j, N) {
			cout << begin_coef.at(i).at(j) << " ";
		}
		cout << endl;
		rep(j, N) cout << cell[i].go[j] << " ";
		cout << endl;
	}

	rep(i, cell_type) {
		cout << cell[i].nut_zero_coef << " ";
		rep(j, N - 1) cout << cell[i].coef[j][j + 1] << " ";
		cout << endl;
	}

	rep(i, cell_type) {
		cout << cell[i].nut_cat << " ";
		rep(j, N - 1) cout << cell[i].catalyst[j][j + 1] << " ";
		cout << endl;
	}

	double max_pop = 0;
	int count_pop = 0;
	rep(i, cell_type) {
		if (i == 0) cout << cell[i].size << " ";
		else {
			if (max_pop < cell[i].size) max_pop = cell[i].size;
			if (cell[i].size > 0.05) count_pop++;
			else if (cell[i].size > 0.01 && prev_size[i] < cell[i].size) count_pop++;
		}
	}
	cout << max_pop << " " << count_pop << endl;
	// take_log_pop << aver_nut << " " << cell[0].size << " " << max_pop << " " << count_pop << endl;

	cout << up << endl;
	cout << "cell_type" << cell_type << endl;

	// return 0;
}
/*
int main(void)
{
	// double give_nut[10] = {0.0001, 0.0005, 0.0007, 0.001, 0.003, 0.005, 0.01, 0.05, 0.1};
	double give_nut[14] = {0.0001, 0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.0015, 0.002, 0.003, 0.005, 0.007, 0.01, 0.015, 0.02};
	
	rep(i, 14) {
		aver_nut = give_nut[i];
		// rep(j, N) init_go[j] = 0;
		// init_go[i] = 1;
		rep(j, 5) main_(j + 1);
	}
	return 0;
}*/
