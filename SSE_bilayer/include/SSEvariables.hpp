#ifndef _SSEVARIABLES_HPP_DEFINED_
#define _SSEVARIABLES_HPP_DEFINED_

#include <random>
#include <chrono>

#pragma once



class SSEvariables {
	public:
    bool time_series;
    int lx, ly, Ns, Nm, NH, NB, NQ;
    double JH, JQ, JB, Beta, lambda;
    double prob_in, prob_rm, cum_prob[10];
    int iter, isteps;
    long n1, Lc;
    int  nl;
    int **JHsites, **JQsites, **JBsites,  *JBsgnx, *JBsgny, *JHsgnx, *JHsgny, *JQsgnx, *JQsgny;
    double *Hsgn, *QQsgn, *Bsgn;

  	void lattice_sites();
  	void declare_variables();

};


class ran {
  private:
    std::mt19937 mt;
    std::uniform_real_distribution < double > dist;

  public:
    ran(double lower = 0.0, double upper = 1.0): mt(std::chrono::high_resolution_clock::now().time_since_epoch().count()), dist(lower, upper) {}
  double operator()() {
    return dist(mt);
  }
};

//class ran {
//private:
//    std::mt19937_64 mt;
//    std::uniform_real_distribution<double> dist;

//public:
//    ran(double lower = 0.0, double upper = 1.0) : mt(std::random_device{}()), dist(lower, upper) {}

//    double operator()() {
//        return dist(mt);
//    }
//};
#endif
