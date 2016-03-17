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
	ofstream take_log_type("data_test_type_number10_0.01.log");
	ofstream take_log_network("data_test_network10_0.01.log");
	ofstream take_log_outside("data_test_outside10_0.01.log");
	ofstream take_log_devdev("data_test_devdev10_0.01.log");
	ofstream take_log_coef("data_test_coef10_0.01.log");
	ofstream take_log_inside("data_test_inside10_0.01.log");
	ofstream take_log_come("data_test_come.log");
	ofstream take_log_up("data_test_up.log");
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
int cell_type = 1;
const int time_end = 1000000;
const double time_bunkai = 0.1;
const int run_time = 1;

const double box_size = 1;
const double outside_size = 0.1;

const int N = 1;

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

void init(void)
{
	//defをいれておいて、typeカウントを防ぐ
	def.type = - 1;
	rep(i, N) cell[i] = def;

	//coefの決め方
	all_init();
	rep(i, 50) sum_init(i + 1, 1);

	cell_number = cell_type;

	//outside系はloopの外
	outside_nut = aver_nut;
	rep(i,N) {
		outside[i] = 0;
		coef_decrease[i] = 0;
	}

	nut_coef = 0.1;
	nut_reversible = 0;
	aver_nut = 0.001;

	rep(i, cell_type) {
		cell[i].nut_zero_coef = begin_coef.at(i).at(0);
		cell[i].type = i;
		cell[i].size = box_size / cell_type;
		cell[i].nut = 1.0 / (double)(N + 1);
		rep(j, N) {
			cell[i].mol[j] = 1.0 / (double)(N + 1);
			cell[i].nut_cat = 0;
			cell[i].go[j] = (j % 3 / 2) * 0;
			if (j != N - 1) {
				cell[i].coef[j][j + 1] = begin_coef.at(i).at(j + 1);
				cell[i].catalyst[j][j + 1] = (j - 2 + N) % N;
				reversible[j][j + 1] = 0;
			}
		}
		cell[i].size = 1.0;
		cell[i].init_last = cell[i].mol[N - 1];
	}
}

double prev_outside_nut;
double prev_outside[N];

Cell internal(Cell p)
{
	//前の値を保存
	double prev[N];
	rep(i, N) prev[i] = p.mol[i];

	double prev_nut = p.nut;

	//細胞内外の栄養の流出入
	p.nut += time_bunkai * nut_coef/* * pow(p.size, - 1.0 / 3.0)*/ * (prev_outside_nut - prev_nut);
	outside_nut -= time_bunkai * nut_coef/* * pow(p.size, - 1.0 / 3.0)*/ * (prev_outside_nut - prev_nut) * p.size / outside_size;

	//リアクション
	p.mol[0] += time_bunkai * p.nut_zero_coef * prev_nut * prev[p.nut_cat];
	p.nut -= time_bunkai * p.nut_zero_coef * prev_nut * prev[p.nut_cat];

	rep(i, N - 1) {
		p.mol[i + 1] += time_bunkai * p.coef[i][i + 1] * prev[i] * prev[p.catalyst[i][i + 1]];
		p.mol[i] -= time_bunkai * p.coef[i][i + 1] * prev[i] * prev[p.catalyst[i][i + 1]];
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
		outside[i] -= time_bunkai * p.go[i]/* * pow(p.size, - 1.0 / 3.0)*/ * (prev_outside[i] - prev[i]) * p.size / outside_size;
	}

	//リサイズ
	double sum = p.nut;
	rep(i, N) sum += p.mol[i];
	p.size = sum * p.size;
	p.nut = p.nut / sum;
	rep(i, N) p.mol[i] = p.mol[i] / sum;

	take_log_inside << p.nut << " ";

	return p;
}

void evolve(void)
{
	// double evolve_size = 1.0 / ((double)cell_type + 1.0);
	double evolve_size = 0.1;
	rep(i, cell_type) {
		cell[i].size = cell[i].size * (1 - evolve_size);
	}
	int i = cell_type;
	cell[i].nut_zero_coef = begin_coef.at(i).at(0);
	cell[i].type = i;
	cell[i].size = evolve_size;
	cell[i].nut = 1.0 / (double)(N + 1);
	rep(j, N) {
		cell[i].mol[j] = 1.0 / (double)(N + 1);
		cell[i].nut_cat = N - 2;
		cell[i].go[j] = (j % 3 / 2) * 0;
		if (j != N - 1) {
			cell[i].coef[j][j + 1] = begin_coef.at(i).at(j + 1);
			cell[i].catalyst[j][j + 1] = (j - 2 + N) % N;
			reversible[j][j + 1] = 0;
		}
	}
	cell[i].init_last = cell[i].mol[N - 1];
	cell_type++;
}

double up = 1;

void process(int t)
{
	box_con = decide_box_nut(t);

	//outside系をprevにいれる
	prev_outside_nut = outside_nut;
	rep(i, N) prev_outside[i] = outside[i];

	//loop回してreaction
	rep(i, cell_type) {
		if (cell[i].size < 10e-8) {
			cell[i].size = 0;
			double sum = 0;
			rep(j, cell_type) sum += cell[j].size;
			rep(j, cell_type) cell[j].size = cell[j].size * box_size / sum;
		} else cell[i] = internal(cell[i]);
	}
	take_log_inside << endl;

	//リサイズ
	double sum = 0;
	rep(i, cell_type) sum += cell[i].size;
	rep(i, cell_type) {
		cell[i].size = cell[i].size / sum * box_size;
		cell[i].nut = cell[i].nut * box_size / sum;
		rep(j, N) cell[i].mol[j] = cell[i].mol[j] * box_size / sum;
	}
	take_log_up << sum << endl;
	up = up * sum;

	//outsideの値を更新
	outside_nut += time_bunkai * (box_con - prev_outside_nut) * outside_size;

	take_log_outside << outside_nut << " ";
	rep(i, N) {
		take_log_outside << outside[i] << " ";
	}
	take_log_outside << endl;
}

int main(void)
{
	//randomの種を与える
	srand(1);

	//tun_time回走らせる
	rep(l, run_time) {
		init();
		//time_end秒走らせる
		rep(t, time_end) {
			process(t);
			if (t % 50000 == 0) evolve();
			cout << t << " ";
			rep(i, cell_type) cout << fixed << setprecision(8) << cell[i].size << " ";
			cout << endl;
			rep(i, cell_type) take_log_type << cell[i].size << " ";
			take_log_type << endl;
		}
	}
	
	rep(i, cell_type) {
		rep(j, N) {
			cout << begin_coef.at(i).at(j) << " ";
		}
		cout << endl;
	}

	cout << up << endl;

	return 0;
}
