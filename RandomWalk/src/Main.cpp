#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <algorithm>
#include <vector>
#include <map>
#include <set>

#include <cmath>
#include <random>
#include <numeric>

#include <stdlib.h>
#include <iomanip>

#include <Eigen/Dense>

#define DELIM ','
#define EPS 1e-8
#define V_QUADRATIC_ERROR 1e-20
#define DAMPING_FACTOR 0.8

#define BETA_alpha 1.0
#define BETA_beta 1.0

#define CLIQUE_STRATEGY 1
#define INTERSECTION_STRATEGY 2
#define GORI_STRATEGY 3

#define BAYESIAN_d 1
#define GORI_d 2

#define RANDOM_ALGORITHM 1
#define GORI_ALGORITHM 2

#define sigmoid(x) (1/(1+exp(-x)))

using namespace std;
using namespace Eigen;

unsigned int K_FOLD;

struct Item {
	Item() :
		// dummy votes
				nUps(1), nRatings(2) {
	}
	~Item() {
		users.clear();
	}

	unsigned int id;
	unsigned int nUps;
	unsigned int nRatings;
	vector<int> users;
	set<int> usersUp;
};
typedef struct Item Item;

typedef struct {
	unsigned int itemId;
	float rating;
	unsigned short int scaled_rating;
} Vote;

struct User {
	unsigned int id;
	set<unsigned int> items;
	vector<Vote *> ratings;
	~User() {
		for (size_t i = 0; i < ratings.size(); ++i)
			delete ratings[i];
		ratings.clear();
	}
};
typedef struct User User;

typedef struct {
	unsigned int userId;
	unsigned int itemId;
	float rating;

	unsigned short int fold;
} Review;

vector<Review*> reviews;
vector<User*> users;
vector<Item*> items;

map<unsigned int, unsigned int> userIds; // Maps real user's ID to id
map<unsigned int, unsigned int> itemIds; // Maps real item's ID to id

map<unsigned int, unsigned int> rUserIds; // Inverse of the above maps
map<unsigned int, unsigned int> rItemIds;

map<unsigned int, unsigned int>::iterator it;

vector<unsigned int> k_count;

struct Gen {
	Gen() :
		n(K_FOLD) {
		v.resize(K_FOLD, 0);
		fill();
	}
	vector<unsigned int> v;
	unsigned int n;

	void fill() {
		n = K_FOLD;
		for (unsigned int i = 0; i < K_FOLD; ++i)
			v[i] = i;
	}

	unsigned int next() {
		if (n < 1)
			fill();
		int i = rand() % n;
		unsigned int number = v[i];
		unsigned int temp = v[n - 1];
		v[n - 1] = number;
		v[i] = temp;
		n--;
		return number;
	}

};
typedef struct Gen Gen;

int cmp(double x, double y, double tol = 1e-19) {
	return (x <= y + tol) ? (x + tol < y) ? -1 : 0 : 1;
}

// returns the intersection elements between two ordered vectors
vector<int> take_intersection(vector<int> &x, vector<int> &y) {
	unsigned int i = 0, j = 0;
	vector<int> intersection;

	while (i < x.size() && j < y.size()) {
		if (x[i] < y[j])
			i++;
		else if (y[j] < x[i])
			j++;
		else {
			intersection.push_back(x[i]);
			++i;
			++j;
		}
	}
	return intersection;
}

void read_data2(const char* filename) {

	Gen generator;

	ifstream file(filename);
	string line;

	unsigned int userId;
	unsigned int itemId;
	float rating;

	getline(file, line); // reading header


	map<unsigned int, vector<int> > m;
	vector<unsigned int> movies;

	while (getline(file, line)) {
		stringstream ss(line);
		string tok;

		getline(ss, tok, DELIM);
		userId = atoi(tok.c_str());

		getline(ss, tok, DELIM);
		itemId = atoi(tok.c_str());

		getline(ss, tok, DELIM);
		rating = atof(tok.c_str());

		m[itemId].push_back(userId);

		movies.push_back(itemId);
	}

	for (size_t i = 0; i < movies.size(); ++i) {
		unsigned int id = movies[i];
		sort(m[id].begin(), m[id].end());
	}

	for (size_t i = 0; i < movies.size(); ++i) {
		unsigned int iId = movies[i];

		for (size_t j = i + 1; j < movies.size(); ++j) {
			unsigned int jId = movies[j];
			if (take_intersection(m[iId], m[jId]).size() > 0) {
				cout << iId << "," << jId << endl;
			}
		}
	}

}

void read_data(const char* filename) {

	Gen generator;

	ifstream file(filename);
	string line;

	unsigned int userId;
	unsigned int itemId;
	float rating;

	getline(file, line); // reading header

	while (getline(file, line)) {
		stringstream ss(line);
		string tok;

		getline(ss, tok, DELIM);
		userId = atoi(tok.c_str());

		getline(ss, tok, DELIM);
		itemId = atoi(tok.c_str());

		getline(ss, tok, DELIM);
		rating = atof(tok.c_str());

		//timestamp
		//getline(ss, tok, DELIM);

		//cout << userId << " " << itemId << " " << rating << endl;
		Review * r = new Review();
		r->userId = userId;
		r->itemId = itemId;
		r->rating = rating;
		r->fold = generator.next();

		k_count[r->fold]++;

		reviews.push_back(r);
	}
}

// makes columns to sum one
void make_stochastic(MatrixXd & m) {
	const unsigned int N = items.size();
	for (unsigned int j = 0; j < N; ++j) {
		double sum = 0;
		for (unsigned int i = 0; i < N; ++i) {
			sum += m(i, j);
		}
		for (unsigned int i = 0; i < N; ++i) {
			m(i, j) /= sum;
		}
	}
}

// scales user rating so that 1 if above user average or 0 otherwise
void scaleRating(User *u) {
	float avg = 0;
	for (size_t i = 0; i < u->ratings.size(); ++i) {
		avg += u->ratings[i]->rating;
	}
	avg /= u->ratings.size();

	for (size_t i = 0; i < u->ratings.size(); ++i) {
		u->ratings[i]->scaled_rating = u->ratings[i]->rating >= avg ? 1 : 0;
	}
}

float doa(const User *u, const vector<Review*> &test, const VectorXd & rank) {
	int nw = 0;
	int count = 0;

	for (size_t i = 0; i < items.size(); ++i) {
		unsigned int k = items[i]->id;

		// checks if user u rated item k
		if (u->items.find(k) != u->items.end())
			continue;

		// since item in test set is not in training set
		nw++;

		for (size_t a = 0; a < test.size(); ++a) {
			if(itemIds.find(test[a]->itemId) == itemIds.end()){
				cout << "error: item test" << test[a]->itemId << " is not in training" << endl;
				exit(0);
			}
			unsigned int j = itemIds[test[a]->itemId];
			if (cmp(rank(j), rank(k), EPS) >= 0) {
				count++;
			}
		}
	}
	return ((float) count) / (test.size() * nw);
}

typedef std::pair<double, int> rankPair;
bool rankComparator(const rankPair & l, const rankPair & r) {
	return (cmp(l.first, r.first, 1e-9) < 0);
}

vector< vector<double> > getPrecisionRecall(const User *u, const vector<Review*> &test,
		const VectorXd & rank, vector<unsigned int> kValues) {

	// precision index - 0
	// recall index - 1
	vector< vector<double> > result(2);
	result[0].resize(kValues.size(), 0);
	result[1].resize(kValues.size(), 0);

	// construindo novo rank
	vector<pair<double, int> > uRank;
	for (size_t i = 0; i < items.size(); ++i) {
		uRank.push_back(make_pair(rank[i], i));
	}
	sort(uRank.begin(), uRank.end(), rankComparator);

	// colocando itens do test no sets
	set<int> uTest;
	for(size_t i=0; i < test.size(); ++i){
		uTest.insert(itemIds[test[i]->itemId]);
	}

	for(size_t i=0; i < kValues.size(); ++i){
		const unsigned int k = kValues[i];

		int count = 0;
		for(size_t j=0; j < k; ++j){
			if(uTest.find(uRank[j].second) != uTest.end()){
				count++;
			}
		}


		result[0][i] = ((double) count)/k; //precision
		result[1][i] = ((double) count)/uTest.size(); //recall
	}
	return result;
}

VectorXd gori_pucci(MatrixXd &m, User *u) {

	//cout << "** gori_pucci \t" << userId << endl;
	const unsigned int N = items.size();

	// creating vector d
	VectorXd d = VectorXd::Zero(N);
	double sum = 0.0;
	for (unsigned int i = 0; i < u->ratings.size(); ++i) {

		unsigned int itemId = u->ratings[i]->itemId;
		d(itemId) = u->ratings[i]->rating;
		sum += d(itemId);
	}
	d /= sum; //making stochastic


	// random walk
	VectorXd last_v;
	VectorXd v = (1.0 / N) * VectorXd::Ones(N);

	float norm2 = 0.0;
	do {
		last_v = v;
		v = DAMPING_FACTOR * m * v + (1 - DAMPING_FACTOR) * d;
		norm2 = (v - last_v).squaredNorm();
	} while (cmp(norm2, 0, V_QUADRATIC_ERROR) > 0);

	return v;
}

void make_matrix_gori(MatrixXd &m, const unsigned int N) {
	// creating transition matrix
	for (size_t a = 0; a < N; ++a) {
		Item *i = items[a];
		for (size_t b = a + 1; b < N; ++b) {
			Item *j = items[b];

			int size = take_intersection(i->users, j->users).size();

			m(i->id, j->id) = size;
			m(j->id, i->id) = size;

		}
	}
}

float run_gori(map<unsigned int, vector<Review*>> &test) {

	// creating matrix
	const unsigned int N = items.size();

	MatrixXd m = MatrixXd::Zero(N, N);

	make_matrix_gori(m, N);
	make_stochastic(m);

	double macroDOA = 0;
	map<unsigned int, vector<Review*>>::iterator it;
	for (it = test.begin(); it != test.end(); ++it) {
		User *u = users[userIds[it->first]];

		VectorXd rank = gori_pucci(m, u);

		macroDOA += doa(u, it->second, rank);
	}
	return macroDOA / users.size();
}

bool check_d(User *u, MatrixXd &m, unsigned itemId) {

	// if item is in user's list
	if (u->items.find(itemId) != u->items.end())
		return false;

	// if item is neighbor
	for (set<unsigned int>::iterator it = u->items.begin(); it
			!= u->items.end(); ++it) {
		if (m(*it, itemId) != 0 || m(itemId, *it) != 0) {
			return true;
		}
	}
	return false;
}

VectorXd random_walk(MatrixXd &m, User* u, int strategy_d) {

	const unsigned int N = items.size();

	VectorXd d = VectorXd::Zero(N);
	double sum = 0.0;

	switch (strategy_d) {
	case GORI_d:
		// GORI
		for (unsigned int i = 0; i < u->ratings.size(); ++i) {

			unsigned int itemId = u->ratings[i]->itemId;
			d(itemId) = u->ratings[i]->rating;
			sum += d(itemId);
		}
		d /= sum; //making stochastic
		break;
	case BAYESIAN_d:
		//BAYESIAN
		for (unsigned int i = 0; i < u->ratings.size(); ++i) {

			unsigned int itemId = u->ratings[i]->itemId;
			Item *p = items[itemId];
			d(p->id) = (BETA_alpha + p->nUps + 1) / (BETA_alpha + BETA_beta
					+ p->nRatings + 1);
			sum += d(p->id);
		}

		d /= sum;
		break;
	default:
		cout << "Error \t strategy-d" << endl;
		exit(0);
	}
	/*for (unsigned int i = 0; i < N; ++i) {
	 Item *p = items[i];

	 if (check_d(u, m, p->id)) {
	 // bayesian trust
	 d(p->id) = (BETA_alpha + p->nUps + 1) / (BETA_alpha + BETA_beta
	 + p->nRatings + 1);
	 sum += d(p->id);
	 }
	 }*/

	// ###### RANDOM WALK
	VectorXd last_v;
	//VectorXd v = (1.0 / N) * VectorXd::Ones(N);

	// teste
	VectorXd v = VectorXd::Zero(N);
	for (unsigned int i = 0; i < u->ratings.size(); ++i) {
		unsigned int itemId = u->ratings[i]->itemId;
		v(itemId) = 1;
	}
	v = (1.0 / u->ratings.size()) * v;

	double norm2 = 0.0;
	do {
		last_v = v;
		v = DAMPING_FACTOR * m * v + (1 - DAMPING_FACTOR) * d;
		norm2 = (v - last_v).squaredNorm();
	} while (cmp(norm2, 0, V_QUADRATIC_ERROR) > 0);

	return v;
}

void make_matrix_intersection(MatrixXd &m, const unsigned int N) {
	for (size_t a = 0; a < N; ++a) {
		Item *i = items[a];
		for (size_t b = a + 1; b < N; ++b) {
			Item *j = items[b];

			//tansition
			float i_ratio = ((float) i->nUps) / (i->nRatings - i->nUps);

			float j_ratio = ((float) j->nUps) / (j->nRatings - j->nUps);

			vector<int> v = take_intersection(i->users, j->users);
			if (v.size() > 0) {
				int countUps = 1;
				for (vector<int>::iterator it = v.begin(); it != v.end(); it++) {
					if (i->usersUp.find(*it) != i->usersUp.end()
							&& j->usersUp.find(*it) != j->usersUp.end()) {
						countUps++;
					}
				}

				float ratio = ((float) countUps) / (v.size() + 1);
				//				cout << i->id << " " << j->id << " " << v.size() << " "
				//						<< ratio << endl;
				m(j->id, i->id) = m(i->id, j->id) = ratio
						* sigmoid(i_ratio / j_ratio);
			}

		}
	}
}

void make_matrix_clique(MatrixXd &m, const unsigned int N) {
	for (size_t k = 0; k < users.size(); ++k) {
		for (size_t a = 0; a < users[k]->ratings.size(); ++a) {
			for (size_t b = a + 1; b < users[k]->ratings.size(); ++b) {
				int i = users[k]->ratings[a]->itemId;
				int j = users[k]->ratings[b]->itemId;

				//tansition
				float i_ratio = ((float) items[i]->nUps) / (items[i]->nRatings
						- items[i]->nUps);

				float j_ratio = ((float) items[j]->nUps) / (items[j]->nRatings
						- items[j]->nUps);

				m(i, j) = i_ratio / j_ratio;
				m(j, i) = 1 / m(i, j);
			}
		}
	}
}

float run_random(map<unsigned int, vector<Review*>> &test, int strategy,
		int strategy_d) {

	// scales ratings according to user's bias and counts up-votes
	for (size_t i = 0; i < users.size(); ++i) {
		User *u = users[i];
		scaleRating(u);

		for (size_t j = 0; j < u->ratings.size(); ++j) {
			Vote * v = u->ratings[j];

			Item *p = items[v->itemId];
			p->nRatings++;
			if (v->scaled_rating > 0) {
				p->nUps++;
				p->usersUp.insert(u->id);
			}
		}
	}

	const unsigned int N = items.size();

	// creating transition matrix
	MatrixXd m = MatrixXd::Zero(N, N);

	switch (strategy) {
	case INTERSECTION_STRATEGY:
		make_matrix_intersection(m, N);
		break;
	case CLIQUE_STRATEGY:
		make_matrix_clique(m, N);
		break;
	case GORI_STRATEGY:
		make_matrix_gori(m, N);
	default:
		cout << "Error \t core strategy" << endl;
		exit(0);
	}
	make_stochastic(m);

	float macroDOA = 0;
	map<unsigned int, vector<Review*>>::iterator it;

	for (it = test.begin(); it != test.end(); ++it) {
		User *u = users[userIds[it->first]];
		VectorXd rank = random_walk(m, u, strategy_d);
		macroDOA += doa(u, it->second, rank);
	}
	return macroDOA / users.size();
}

void clear() {
	userIds.clear();
	rUserIds.clear();

	itemIds.clear();
	rItemIds.clear();

	for (size_t i = 0; i < users.size(); ++i) {
		delete users[i];
	}
	users.clear();

	for (size_t i = 0; i < items.size(); ++i)
		delete items[i];
	items.clear();

}

void kfold(char algorithm, int strategy = 0, int strategy_d = 0) {

	for (unsigned int k = 0; k < K_FOLD; ++k) {
		cout << "####K=" << k << endl;

		map<unsigned int, vector<Review*>> test; // key is  real user id in dataset


		for (size_t z = 0; z < reviews.size(); ++z) {
			Review* r = reviews[z];

			if (K_FOLD > 1 && r->fold == k) {
				test[r->userId].push_back(r);
				continue;
			} else {
				unsigned int userId = r->userId;
				unsigned int itemId = r->itemId;

				User * u;
				Item *p;

				it = userIds.find(userId);
				if (it == userIds.end()) {
					int id = users.size();
					userIds[userId] = id;
					rUserIds[id] = userId;
					u = new User();
					u->id = id;
					users.push_back(u);
				} else {
					u = users[it->second];
				}

				it = itemIds.find(itemId);
				if (it == itemIds.end()) {
					int id = items.size();
					itemIds[itemId] = id;
					rItemIds[id] = itemId;
					p = new Item();
					p->id = id;
					items.push_back(p);
				} else {
					p = items[it->second];
				}
				p->users.push_back(u->id);

				Vote *v = new Vote();
				v->itemId = p->id;
				v->rating = r->rating;

				u->ratings.push_back(v);
				u->items.insert(p->id);
			}
		}

		for (size_t i = 0; i < items.size(); ++i) {
			Item *p = items[i];
			sort(p->users.begin(), p->users.end());
		}

		switch (algorithm) {
		case GORI_ALGORITHM:
			cout << run_gori(test) << endl;
			break;
		case RANDOM_ALGORITHM:
			switch (strategy) {
			case CLIQUE_STRATEGY:
				cout << run_random(test, CLIQUE_STRATEGY, strategy_d) << endl;
				break;
			case INTERSECTION_STRATEGY:
				cout << run_random(test, INTERSECTION_STRATEGY, strategy_d)
						<< endl;
				break;
			case GORI_STRATEGY:
				cout << run_random(test, GORI_STRATEGY, strategy_d) << endl;
				break;
			}
		}

		//#pragma omp parallel sections
		//		{
		//			{
		//				macroDOA[0][k] = run_random(test, CLIQUE_STRATEGY);
		//
		//			}
		//#pragma omp section
		//			{
		//				macroDOA[1][k] = run_random(test, INTERSECTION_STRATEGY);
		//			}
		//
		//#pragma omp section
		//			{
		//				macroDOA[2][k] = run_gori(test);
		//			}
		//		}
		//		cout << macroDOA[0][k] << "\t" << macroDOA[1][k] << "\t"
		//				<< macroDOA[1][k] << endl;
	}
	clear();
}

/*
 * 1- filename
 * 2- number of folds of k-fold
 * 3- algorithm to run
 * 4- strategy chosen to run random-walk
 */

int main(int argc, char **argv) {
	srand(0);

	char* filename = argv[1];
	if (1)
		read_data2(filename);

	K_FOLD = atoi(argv[2]);
	int algorithm = atoi(argv[3]);
	int strategy = atoi(argv[4]);
	int strategy_d = atoi(argv[5]);

	k_count.resize(K_FOLD, 0);

	read_data(filename);

	switch (algorithm) {
	case GORI_ALGORITHM:
		kfold(GORI_ALGORITHM, strategy, strategy_d);
		break;
	case RANDOM_ALGORITHM:
		kfold(RANDOM_ALGORITHM, strategy, strategy_d);
		break;
	default:
		cout << "error" << endl;
	}

	return 0;
}
