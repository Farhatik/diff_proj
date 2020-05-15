#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

struct Dividends {
    std::vector<double> payments;
    std::vector<double> years;
};

struct Data {
    double S;  // Option price
    double K;  // Strike price
    double r;   // Risk-free rate (5%)
    double v;    // Volatility of the underlying (20%)
    double T;    // One year until expiry
    Dividends D;
};

// Standard normal probability density function

double norm_pdf(const double &x) {
    return (1.0 / (pow(2 * M_PI, 0.5))) * exp(-0.5 * x * x);
}


/* An approximation to the cumulative distribution function
* for the standard normal distribution */

double norm_cdf(const double &x) {
    double k = 1.0 / (1.0 + 0.2316419 * x);
    double k_sum = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + 1.330274429 * k))));

    if (x >= 0.0) {
        return (1.0 - (1.0 / (pow(2 * M_PI, 0.5))) * exp(-0.5 * x * x) * k_sum);
    } else {
        return 1.0 - norm_cdf(-x);
    }
}


/* This calculates d_j, for j in {1,2}. This term appears in the closed
* form solution for the European call or put price */

double d_j(const int &j, const double &S, const double &K, const double &r, const double &v, const double &T) {
    return (log(S / K) + (r + (pow(-1, j - 1)) * 0.5 * v * v) * T) / (v * (pow(T, 0.5)));
}

double sum(Dividends &dividends, double &r) {
    double sum = 0;

    for (size_t i = 0; i < dividends.payments.size(); ++i) {
        double param = 1 + r;
        sum += 1.0 * dividends.payments[i] / std::pow(param, dividends.years[i]);
    }

    return sum;
}

/* underlying S, strike K, risk-free rate r, volatility of
*  underlying sigma and time to maturity */

double call_price_with_dividents(Data &data) {
    double S = data.S;
    double K = data.K;
    double r = data.r;
    double v = data.v;
    double T = data.T;

    Dividends dividends = data.D;

    double first = S * norm_cdf(d_j(1, S, K, r, v, T)) * exp(-r * T);
    double second = -K * exp(-r * T) * norm_cdf(d_j(2, S, K, r, v, T));
    double third = -norm_cdf(d_j(1, S, K, r, v, T)) * sum(dividends, r) * exp(-r * T);

    return first + second + third;
}

double put_price_with_dividents(Data &data) {
    double S = data.S;
    double K = data.K;
    double r = data.r;
    double v = data.v;
    double T = data.T;

    Dividends dividends = data.D;

    //return -S*norm_cdf(-d_j(1, S, K, r, v, T))+K*exp(-r*T) * norm_cdf(-d_j(2, S, K, r, v, T));

    double first = -S * norm_cdf(-d_j(1, S, K, r, v, T)) * exp(-r * T);
    double second = K * exp(-r * T) * norm_cdf(-d_j(2, S, K, r, v, T));
    double third = norm_cdf(-d_j(1, S, K, r, v, T)) * sum(dividends, r) * exp(-r * T);

    return first + second + third;
}

double call_price(Data &data) {
    double S = data.S;
    double K = data.K;
    double r = data.r;
    double v = data.v;
    double T = data.T;

    Dividends dividends = data.D;

    return S * norm_cdf(d_j(1, S, K, r, v, T)) - K * exp(-r * T) * norm_cdf(d_j(2, S, K, r, v, T));
}


double put_price(Data &data) {
    double S = data.S;
    double K = data.K;
    double r = data.r;
    double v = data.v;
    double T = data.T;

    Dividends dividends = data.D;

    return -S * norm_cdf(-d_j(1, S, K, r, v, T)) + K * exp(-r * T) * norm_cdf(-d_j(2, S, K, r, v, T));
}


int main(int, char **argv) {
    std::ifstream infile(argv[1]);
    std::ofstream outfile(argv[2]);

    Data data;

    // First we create the parameter list
    while (infile >> data.S >> data.K >> data.r >> data.v >> data.T) {
        double current = 0;

        infile >> current;
        while (current > 0) {
            data.D.payments.push_back(current);
            infile >> current;
        }

        infile >> current;
        while (current > 0) {
            data.D.years.push_back(current);
            infile >> current;
        }

        double call = call_price(data);
        double put = put_price(data);
        double call_w_div = call_price_with_dividents(data);
        double put_w_div = put_price_with_dividents(data);

        // Finally we output the parameters and prices
        outfile << "Underlying:      " << data.S << std::endl;
        outfile << "Strike:          " << data.K << std::endl;
        outfile << "Risk-Free Rate:  " << data.r << std::endl;
        outfile << "Volatility:      " << data.v << std::endl;
        outfile << "Maturity:        " << data.T << std::endl;

        outfile << "Dividents:        ";
        for (int i = 0; i < data.D.payments.size(); ++i) {
            outfile << data.D.payments[i] << " ";
        }
        outfile << std::endl;

        if ((call > 0) && (call_w_div > 0)) {
            outfile << "Call Price:                     " << call << std::endl;
            outfile << "Call Price with dividents:      " << call_w_div << std::endl;
        }
        if ((put > 0) && (put_w_div > 0)) {
            outfile << "Put Price:                      " << put << std::endl;
            outfile << "Put Price with dividents:       " << put_w_div << std::endl;
        }

        outfile << std::endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;

        data.D.years.clear();
        data.D.payments.clear();
    }

    return 0;
}
