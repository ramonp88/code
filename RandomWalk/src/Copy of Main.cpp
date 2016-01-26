//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <string>
//
//#include <algorithm>
//#include <vector>
//#include <map>
//#include <set>
//
//#include <random>
//#include <numeric>
//
//#include <stdlib.h>
//#include <iomanip>
//
//#include <Eigen/Dense>
//
//#define DELIM ','
//#define V_QUADRATIC_ERROR 1e-20
//#define DAMPING_FACTOR 0.8
//
//#define BETA_alpha 1.0
//#define BETA_beta 1.0
//
//using namespace std;
//using namespace Eigen;
//
//struct Item {
//	Item() :
//		// dummy vote
//				nUps(1), nRatings(1) {
//	}
//
//	unsigned int id;
//	unsigned int nUps;
//	unsigned int nRatings;
//	vector<int> users;
//};
//typedef struct Item Item;
//
//typedef struct {
//	unsigned int itemId;
//	float rating;
//} Vote;
//
//typedef struct {
//	unsigned int id;
//	vector<Vote *> ratings;
//} User;
//
//vector<User *> users;
//vector<Item *> items;
//
//map<int, int> userIds; // Maps a user's ID to an integer
//map<int, int> itemIds; // Maps an item's ID to an integer
//
//map<int, int> rUserIds; // Inverse of the above maps
//map<int, int> rItemIds;
//
//map<int, int>::iterator it;
//
//void scaleRating(User *u) {
//	float avg = 0;
//	for (size_t i = 0; i < u->ratings.size(); ++i) {
//		avg += u->ratings[i]->rating;
//	}
//	avg /= u->ratings.size();
//
//	for (size_t i = 0; i < u->ratings.size(); ++i) {
//		u->ratings[i]->rating = u->ratings[i]->rating >= avg ? 1 : 0;
//	}
//}
//
//void write_csv(MatrixXf & m) {
//	ofstream csv("adjacency.csv");
//
//	cout << m.rows() << "\t" << m.cols() << endl;
//
//	for (int i = 0; i < m.rows(); ++i) {
//		for (int j = 0; j < m.cols(); ++j) {
//			csv << (m(i, j) > 0 ? 1 : 0) << ",";
//		}
//		csv << (m(i, m.cols() - 1) > 0 ? 1 : 0) << endl;
//	}
//
//}
//
//void read_data(const char* filename) {
//	ifstream file("ratings.csv");
//	string line;
//
//	unsigned int userId;
//	unsigned int itemId;
//	float rating;
//
//	getline(file, line); // reading header
//
//	while (getline(file, line)) {
//		stringstream ss(line);
//		string tok;
//
//		getline(ss, tok, DELIM);
//		userId = atoi(tok.c_str());
//
//		getline(ss, tok, DELIM);
//		itemId = atoi(tok.c_str());
//
//		getline(ss, tok, DELIM);
//		rating = atof(tok.c_str());
//
//		//timestamp
//		getline(ss, tok, DELIM);
//
//		//cout << userId << " " << itemId << " " << rating << endl;
//
//		User *u;
//		Item *p;
//
//		it = userIds.find(userId);
//		if (it == userIds.end()) {
//			int id = users.size();
//			rUserIds[id] = userId;
//			userIds[userId] = id;
//			u = new User();
//			u->id = id;
//			users.push_back(u);
//		} else {
//			u = users[it->second];
//		}
//
//		it = itemIds.find(itemId);
//		if (it == itemIds.end()) {
//			int id = items.size();
//			rItemIds[id] = itemId;
//			itemIds[itemId] = id;
//			p = new Item();
//			p->id = id;
//			items.push_back(p);
//		} else {
//			p = items[it->second];
//		}
//		p->users.push_back(u->id);
//
//		Vote *v = new Vote();
//		v->itemId = p->id;
//		v->rating = rating;
//		u->ratings.push_back(v);
//	}
//
//	for (size_t i = 0; i < items.size(); ++i) {
//		Item *p = items[i];
//		sort(p->users.begin(), p->users.end());
//	}
//}
//
//void make_stochastic(MatrixXf & m) {
//	const unsigned int N = items.size();
//	for (unsigned int j = 0; j < N; ++j) {
//		double sum = 0;
//		for (unsigned int i = 0; i < N; ++i) {
//			sum += m(i, j);
//		}
//		for (unsigned int i = 0; i < N; ++i) {
//			m(i, j) /= sum;
//		}
//	}
//}
//
//int intersection_size(vector<int> &x, vector<int> &y) {
//	unsigned int i = 0, j = 0;
//	unsigned int count = 0;
//
//	while (i < x.size() && j < y.size()) {
//		if (x[i] < y[j])
//			i++;
//		else if (y[j] < x[i])
//			j++;
//		else {
//			++count;
//			++i;
//			++j;
//		}
//	}
//	return count;
//}
//
//VectorXf gori_pucci(int userId) {
//
//	User *u = users[userId];
//
//	// creating matrix
//	const unsigned int N = items.size();
//
//	// creating vector d
//	VectorXf d = VectorXf::Zero(N);
//	double sum = 0.0;
//	for (unsigned int i = 0; i < u->ratings.size(); ++i) {
//
//		unsigned int itemId = u->ratings[i]->itemId;
//		d(itemId) = u->ratings[i]->rating;
//		sum += d(itemId);
//	}
//	d /= sum; //making stochastic
//
//	// creating transition matrix
//	MatrixXf m = MatrixXf::Zero(N, N);
//	for (size_t a = 0; a < items.size(); ++a) {
//		Item *i = items[a];
//		for (size_t b = a + 1; b < items.size(); ++b) {
//			Item *j = items[b];
//
//			int size = intersection_size(i->users, j->users);
//
//			m(i->id, j->id) = size;
//			m(j->id, i->id) = size;
//
//		}
//	}
//	make_stochastic(m);
//	if (1) {
//		write_csv(m);
//		return d;
//	}
//
//	// random walk
//	VectorXf last_v;
//	VectorXf v = (1.0 / N) * VectorXf::Ones(N);
//
//	float norm2 = 0.0;
//	int i = 0;
//	do {
//		last_v = v;
//		v = DAMPING_FACTOR * m * v + (1 - DAMPING_FACTOR) * d;
//		norm2 = (v - last_v).squaredNorm();
//		cout << (++i) << "\t" << setprecision(20) << norm2 << endl;
//	} while (norm2 > V_QUADRATIC_ERROR);
//
//	return v;
//}
//
//VectorXf random_walk() {
//
//	// scales ratings according to user's bias and counts up votes
//	for (size_t i = 0; i < users.size(); ++i) {
//		User *u = users[i];
//		scaleRating(u);
//
//		for (size_t j = 0; j < u->ratings.size(); ++j) {
//			Vote * v = u->ratings[j];
//
//			Item *p = items[v->itemId];
//			p->nRatings++;
//			if (v->rating > 0) {
//				p->nUps++;
//			}
//		}
//	}
//
//	// ########### creating matrix
//	const unsigned int N = items.size();
//
//	VectorXf d = VectorXf::Zero(N);
//	double sum = 0.0;
//	for (unsigned int i = 0; i < N; ++i) {
//		Item *p = items[i];
//
//		d(i) = (BETA_alpha + p->nUps) / (BETA_alpha + BETA_beta + p->nRatings); // bayesian trust
//		sum += d(i);
//	}
//	d.normalize();
//
//	MatrixXf m = MatrixXf::Zero(N, N);
//	for (size_t k = 0; k < users.size(); ++k) {
//		for (size_t a = 0; a < users[k]->ratings.size(); ++a) {
//			for (size_t b = a + 1; b < users[k]->ratings.size(); ++b) {
//				int i = users[k]->ratings[a]->itemId;
//				int j = users[k]->ratings[b]->itemId;
//
//				//tansition probability
//				float i_ratio = ((float) items[i]->nUps);
//				float j_ratio = ((float) items[j]->nUps);
//
//				m(i, j) = i_ratio / j_ratio;
//				m(j, i) = 1 / m(i, j);
//			}
//		}
//	}
//
//	make_stochastic(m);
//
//	// ###### RANDOM WALK
//	VectorXf last_v;
//	VectorXf v = (1.0 / N) * VectorXf::Ones(N);
//
//	float norm2 = 0.0;
//	int i = 0;
//	do {
//		last_v = v;
//		v = DAMPING_FACTOR * m * v + (1 - DAMPING_FACTOR) * d;
//		norm2 = (v - last_v).squaredNorm();
//		cout << (++i) << "\t" << setprecision(20) << norm2 << endl;
//	} while (norm2 > V_QUADRATIC_ERROR);
//
//	return v;
//}
//
//int main() {
//	srand(0);
//
//	std::default_random_engine generator;
//	std::uniform_int_distribution<int> distribution(0, 9);
//
//	for (int i = 0; i < 10; ++i)
//		cout << distribution(generator) << endl;
//
//	int numbers[10];
//
//	std::iota(numbers, numbers + 10, 100);
//
//	std::cout << "numbers:";
//	for (int i:numbers) {
//		std::cout << ' ' << i;
//	}
//	std::cout << '\n';
//	cin.get();
//
//	read_data("ratings.csv");
//	cout << users.size() << " " << items.size() << endl;
//
//	if (1) {
//		gori_pucci(0);
//		return 0;
//	}
//	VectorXf v = random_walk();
//
//	cin.get();
//	cout << setprecision(10) << v << endl;
//	return 0;
//}
