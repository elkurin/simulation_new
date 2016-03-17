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
	ofstream take_log_type("data_same_type_number10_test.log");
	ofstream take_log_network("data_same_network10_test.log");
	ofstream take_log_outside("data_same_outside10_test.log");
	ofstream take_log_devdev("data_same_devdev10_test.log");
	ofstream take_log_coef("data_same_coef10_test.log");
	ofstream take_log_inside("data_same_inside10_test.log");
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

const int cell_max = 100;
int cell_type = 1;
const int time_end = 1000000;
const double time_bunkai = 0.1;
const int run_time = 1;
const int init_box_size = 1200;
const int max_cell_size = 10;
double box_con;

const int N = 10;

int cell_number;
double box_size;

double reversible[N][N];
double coef_decrease[N];

double nut_coef;
double nut_reversible;
double aver_nut;

typedef struct
{
	int type;
	double inside_nut;
	int nut_cat;

	double mol[N];

	double nut_zero_coef;
	double coef[N][N];
	// double go[N];

	int catalyst[N][N];

	double size;
	double init_last;
	double init_last_con;
} Cell;

Cell cell[cell_max];
Cell def;

double outside_nut;
double outside[N];

double go[N];

int _count = 0;
int count_ = 0;

int coexist = 0;

double decide_box_nut(int time)
{
	return aver_nut;
}

double get_size(Cell p)
{
	double sum = p.inside_nut;
	rep(i, N) {
		sum += p.mol[i];
	}
	return sum;
}

double get_box_size(void)
{
	double sum = 0;
	rep(i, cell_number) {
		sum += cell[i].size;
	}
	return init_box_size - sum;
}

array<array<double, N>, 1000> begin_coef;

void sum_ten_init(void)
{
	//coefの合計が一定でその中でrandomをつくる関数
	double sum = 0;
	array<double, N> keep;
	for (int i = 0; i < cell_type; i++) {
		rep(j, N) {
			keep.at(j) = rand() % rand();
			keep.at(j) *= rand();
			sum += keep.at(j);
		}
		rep(j, N) {
			begin_coef.at(i).at(j) = 10 * keep.at(j) / sum;
		}
		sum = 0;
	}
}

void desig_init(void)
{
	rep(i, cell_type) {
		rep(j, N) {
			cin >> begin_coef.at(i).at(j);
		}
	}
}

void all_init(void)
{
	rep(i, N) {
		if (i % 3 != 2) begin_coef.at(0).at(i) = 0.4;
		else begin_coef.at(0).at(i) = 2.4;
	}
}

void network(void)
{
	rep(k, cell_type) {
		rep(i, N) {
			rep(j, N) {
				if (i == j) continue;
				cell[i].coef[i][j] = (double)(rand() % 100) * 0.01;
			}
		}
	}
}

void init(void)
{
	//defをいれておいて、typeカウントを防ぐ
	def.type = - 1;
	rep(i, N) cell[i] = def;

	//coefの定め方
	// sum_ten_init();
	all_init();

	cell_number = cell_type;
	
	//outside系はloopの外
	outside_nut = 20;
	rep(i, N) {
		outside[i] = 0;
		go[i] = (i % 3 / 2) * 1;
		coef_decrease[i] = 0;
	}
	nut_coef = 0.1;
	nut_reversible = 0;
	aver_nut = 0.0001;
	// aver_nut = 1.0;


	//cellごとの詳細設定
	rep(i, cell_type) {
		cell[i].nut_zero_coef = begin_coef.at(i).at(0);
		// cell[i].nut_zero_coef = 0;

		cell[i].type = i;
		rep(j, N) {
			cell[i].mol[j] = 0.2;
			cell[i].inside_nut = 0.2;
			// cell[i].mol[j] = 0;
			// cell[i].inside_nut = 0;

			cell[i].nut_cat = N - 3;
			if (j != N - 1) {
				cell[i].coef[j][j + 1] = begin_coef.at(i).at(j + 1);
				cell[i].catalyst[j][j + 1] = (j - 2 + N) % N;
				reversible[j][j + 1] = 0;
			}
		}
		cell[i].size = get_size(cell[i]);
		cell[i].init_last = cell[i].mol[N - 1];
		cell[i].init_last_con = cell[i].mol[N - 1] / cell[i].size;
	}

	rep(i, 10 * cell_type) {
		cell[i] = cell[0];
	}

	box_size = get_box_size();
	// network();
}


//outside分だけはglobalでここで宣言しておく
double prev_outside_nut_con;
double new_outside_nut_con;
double prev_outside_con[N];
double new_outside_con[N];

int a[1000];

int devdev = 0;

Cell internal(Cell p)
{
	double new_con[N];
	double prev_con[N];

	//数から濃度へ変換
	rep(i, N) {
		prev_con[i] = p.mol[i] / p.size;
		new_con[i] = prev_con[i];
	}
	double prev_nut_con = p.inside_nut / p.size;
	double new_nut_con = prev_nut_con;

	//細胞内外の栄養の流出入
	// cout << pow(p.size, - 1.0 / 3.0) << endl;
	new_nut_con += time_bunkai * nut_coef * pow(p.size, - 1.0 / 3.0) * (prev_outside_nut_con - prev_nut_con);
	new_outside_nut_con -= time_bunkai * nut_coef * pow(p.size, - 1.0 / 3.0) * (prev_outside_nut_con - prev_nut_con) * p.size / box_size;

	//ただのリアクション
	new_con[0] += time_bunkai * p.nut_zero_coef * prev_nut_con * prev_con[p.nut_cat];
	new_nut_con -= time_bunkai * p.nut_zero_coef * prev_nut_con * prev_con[p.nut_cat];
	rep(i, N - 1) {
		new_con[i + 1] += time_bunkai * p.coef[i][i + 1] * prev_con[i] * prev_con[p.catalyst[i][i + 1]];
		prev_con[i] -= time_bunkai * p.coef[i][i + 1] * prev_con[i] * prev_con[p.catalyst[i][i + 1]];
	}
	
	//可逆反応
	new_nut_con += time_bunkai * nut_reversible * prev_con[0];
	new_con[0] -= time_bunkai * nut_reversible * prev_con[0];
	rep(i, N - 1) {
		new_con[i] += time_bunkai * reversible[i][i + 1] * prev_con[i + 1];
		new_con[i + 1] -= time_bunkai * reversible[i][i + 1] * prev_con[i + 1];
	}

	//細胞内外の溶質の流出入
	rep(i, N) {
		new_con[i] += time_bunkai * go[i] * pow(p.size, - 1.0 / 3.0) * (prev_outside_con[i] - prev_con[i]);
		new_outside_con[i] -= time_bunkai * go[i] * pow(p.size, - 1.0 / 3.0) * (prev_outside_con[i] - prev_con[i]) * p.size / box_size;
	}

	//濃度から数に戻して挿入
	rep(i, N) {
		p.mol[i] = new_con[i] * p.size;
	}
	p.inside_nut = new_nut_con * p.size;
	
	//サイズを増やす
	double sum = new_nut_con;
	rep(i, N) {
		sum += new_con[i];
	}
	p.size = sum * p.size;

	return p;
}

void evolve(void)
{
	double sum = 0;
	array<double, N> keep;
	rep(j, N) {
		keep.at(j) = rand();
 		sum += keep.at(j);
	}
	rep(j, N) {
		begin_coef.at(cell_type).at(j) = 1 * keep.at(j) / sum;
	}

	int get = rand() % cell_number;

	// if (cell_number + 1 ==  cell_max) {
	// 	int get;
	// 	// do {
	// 		get = rand() % cell_number;
	// 		// cout << "in while" << endl;
	// 	// } while(get == cell_number - 1);
	//
	// 	cell_number--;
	// 	cell[get] = cell[cell_number];
	// }

	cell[get].nut_zero_coef = begin_coef.at(cell_type).at(1);
	cell[get].type = cell_type;

	rep(j, N) {
		if (j != N - 1) {
			cell[get].coef[j][j + 1] = begin_coef.at(cell_type).at(j + 1);
		}
	}

	cell_type++;
	// cell_number++;
}


pair<Cell, Cell> devide(Cell p)
{
	Cell q, r;
	q = p;
	r = p;
	
	//pが分裂するとき中身は6:4に分裂
	while(1) {
		// cout << "dev while" << endl;
		q.inside_nut = get_rand_normal(0.5) * p.inside_nut;
		if (q.inside_nut >= 0) break;
	}
	r.inside_nut = p.inside_nut - q.inside_nut;

	rep(i, N) {
		while(1) {
			// cout << "dev whil " << i << " " << p.mol[i] << endl;
			q.mol[i] = get_rand_normal(0.5) * p.mol[i];
			if (q.mol[i] >= 0) break;
		}
		r.mol[i] = p.mol[i] - q.mol[i];
	}

	//sizeを入れておく
	q.size = get_size(q);
	r.size = get_size(r);

	//pairで渡す
	pair<Cell, Cell> s = {q, r};
	return s;
}	

double pop_sum = 0;

void process(int t)
{
	box_con = decide_box_nut(t);

	//outside系を濃度に変換
	prev_outside_nut_con = outside_nut / box_size;
	new_outside_nut_con = prev_outside_nut_con;
	rep(i, N) {
		prev_outside_con[i] = outside[i] / box_size;
		new_outside_con[i] = prev_outside_con[i];
	}

	//loop回してreaction
	rep(i, cell_number) {
		if (cell[i].size == 0) {
			cell[i] = cell[cell_number - 1];
			cell_number--;
		} else cell[i] = internal(cell[i]);
	}

	//sizeが2倍以上またはmaxを越えたら分裂、半分になったら消滅
	rep(i, cell_number) {
		if (/* cell[i].size > max_cell_size ||*/ cell[i].mol[N - 1] > 2 * cell[i].init_last) {
			devdev++;
			if (cell[i].mol[N - 1] > 2 * cell[i].init_last) _count++;
			else count_++;
			cell_number++;
			pair<Cell, Cell> dev = devide(cell[i]);
			cell[i] = dev.first;
			if (cell_number == cell_max) {
				int get;
				do {
					// cout << "while for dev" << endl;
					get = rand() % cell_number;
				} while(get == i);

				cell_number--;
				cell[get] = dev.second;
			} else {
				cell[cell_number - 1] = dev.second;
			}
		} else if (cell[i].mol[N] / cell[i].size < 0.5 * cell[i].init_last_con) {
			cell[i] = cell[cell_number - 1];
			cell_number--;
		}
	}

	// cout << "out of dev" << endl;

	//outsideの値を更新
	outside_nut = new_outside_nut_con * box_size + time_bunkai * (box_con - prev_outside_nut_con) * box_size;
	rep(i, N) outside[i] = new_outside_con[i] * box_size;

	// cout << "out of outside" << endl;

	//それぞれのcelltypeの数をカウント、表示
	rep(i, cell_type) a[i] = 0;
	rep(i, cell_number) a[cell[i].type]++;
	// cout << "start counting a " << cell_number << endl;
	cout << t << " ";
	rep(i, cell_type) {
		if (a[i] == 0) {
			take_log_type << a[i] << " ";
			continue;
		}
		cout << setw(4) << a[i] << " ";
		take_log_type << a[i] << " ";
	}
	cout << endl;
	take_log_type << endl;
	//prevの変換が全て終わってからbox_sizeの値を更新
	box_size = get_box_size();
	double x = (double)a[0];
	double y = (double)cell_number;
	if (t >= time_end * 0.2) pop_sum += x / y;

	vector<vector<double>> keep8(cell_type);
	// cout << cell[0].mol[7] << endl;
	rep(i, cell_number) {
		if (cell[i].size != 0) keep8[cell[i].type].push_back(cell[i].mol[7] / cell[i].size);
	}
	rep(i, cell_type) {
		vector<double> sum8(cell_type, 0);
		rep(j, (int)keep8[i].size()) sum8[i] += keep8[i][j];
		if (keep8[i].size() != 0) take_log_inside << sum8[i] / (double)keep8[i].size() << " ";
		else take_log_inside << 0 << " ";
	}
	take_log_inside << endl;
}
	


int main(void)
{
	//randomの種を与える
	srand(10);

	//run_time回だけ走らせる
	rep(l, run_time) {
		init();
		//time_end秒だけ走らせる
		rep(t, time_end) {
			devdev = 0;
			process(t);
			rep(i, N) {
				double sum0 = 0;
				double sum2 = 0;
				int num0 = 0;
				int num2 = 0;
				rep(j, cell_number) {
					if (cell[j].type == 0) {
						sum0 += cell[j].mol[i];
						num0++;
					}
					if (cell[j].type == 16) {
						sum2 += cell[j].mol[i];
						num2++;
					}
				}
				double get0, get2;
				if (num0) get0 = sum0 / (double)num0;
				else get0 = 0;
				if (num2) get2 = sum2 / (double)num2;
				else get2 = 0;
				
				// cout << get0  << " " << get2 << endl;
				// if (i == 8) take_log_network << (double)sum0 / (double)num0  << " " << (double)sum2 / (double)num2 << endl;
				if (get0) take_log_network << (get0 - get2) / get0 << " ";
				else take_log_network << 0 << " ";
			}
			take_log_network << endl;
			// cout << cell[0].size << endl;
			rep(i, N) {
				take_log_outside << outside[i] << " ";
			}
			take_log_outside << endl;
			take_log_devdev << devdev << endl;

			if (t % 20000 == 1000) {
				int count_alive = 0;
				rep(i, cell_type) if (a[i]) count_alive++;
				coexist += count_alive;
				evolve();
			}
		}
	}

	rep(i, cell_type) {
		rep(j, N) {
			cout << begin_coef.at(i).at(j) << " ";
			take_log_coef << begin_coef.at(i).at(j) << " ";
		}
		cout << endl;
		take_log_coef << endl;
	}

	rep(i, cell_type) {
		if (a[i] == 0) {
			continue;
		}
		cout << i << " " << setw(4) << a[i] << endl;
	}
	cout << endl;

	cout << endl << _count << " " << count_ << endl << endl;

	cout << "% = " << pop_sum / ((double)time_end * 0.8) << endl;
	cout << "coexist = " << (double)coexist / (double)(time_end / 20000) << endl << endl;
	return 0;
}

	
