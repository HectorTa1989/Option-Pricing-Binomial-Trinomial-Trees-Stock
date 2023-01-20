#include <iostream>
#include <math.h>
#include <string>
#include <cmath>

auto u = [](double v, double dt)
{
    return exp(v*sqrt(2.0*dt));
};

auto d = [](double v, double dt)
{
    return 1.0 / u(v, dt);
};

auto p_up = [](double r, double v, double dt)
{
    double top = exp(r*dt/2.0) - exp(-v*sqrt(dt/2.0));
    double bot = exp(v*sqrt(dt/2.0)) - exp(-v*sqrt(dt/2.0));
    return pow(top/bot, 2);
};

auto p_dn = [](double r, double v, double dt)
{
    double top = exp(v*sqrt(dt/2.0)) - exp(r*dt/2.0);
    double bot = exp(v*sqrt(dt/2.0)) - exp(-v*sqrt(dt/2.0));
    return pow(top/bot, 2);
};

auto p_m = [](double r, double v, double dt)
{
    return 1.0 - p_up(r, v, dt) - p_dn(r, v, dt);
};

double ** C(double ** tree, int rows, int cols, int split, double S, double K, double r, double v, double dt, double U, double D, double D_UP, double D_M, double D_DOWN, std::string opType)
{
    int ux = 2;
    int cs = 0;
    while(cs <= cols){
        tree[split][cs] = S;
        for(int i = cs + 1; i < cols; ++i){
            tree[split - (i - cs)*ux][i] = tree[split - (i - cs - 1)*ux][i - 1]*U;
            tree[split - (i - cs - 1)*ux][i] = tree[split - (i - cs - 1)*ux][i - 1];
            tree[split + (i - cs)*ux][i] = tree[split + (i - cs - 1)*ux][i - 1]*D;
            tree[split + (i -  cs - 1)*ux][i] = tree[split + (i - cs - 1)*ux][i - 1];
        }
        cs += 2;
    }

    // Calculate Payoffs
    for(int i = 1; i < rows; ++i){
        if(i % 2 != 0){
            if(opType == "call"){
                tree[i][cols - 1] = std::max(tree[i - 1][cols - 1] - K, 0.0);
            } else {
                tree[i][cols - 1] = std::max(K - tree[i - 1][cols - 1], 0.0);
            }
        }
    }

    // Discount backwards
    int cx = 0;
    int ls = 0;
    int osx = 0;

    for(int i = 1; i < cols; ++i){
        cx = cols - (i + 1);
        ls = 4 + osx;
        while(ls < rows - osx){
            tree[ls - 1][cx] = std::max(exp(-r*dt)*(D_UP*tree[ls - 1 - 2][cx + 1] + D_M*tree[ls - 1][cx + 1] + D_DOWN*tree[ls - 1 + 2][cx + 1]), 0.0);
            ls += 2;
        }
       
        osx += 2;
    }

    return tree;
}

double OptionPrice(double ** tree, int split)
{
    return tree[split + 1][0];
}


int main()
{
    // Declare Inputs
    double S = 100.0; // Stock Price
    double K = 95.0; // Strike Price
    double r = 0.05; // RiskFree Rate
    double v = 0.30; // Volatility
    double t = 30.0/365.0; // Expiry
    std::string opType = "put";

    int nodes = 14;
    double dt = t / (double) nodes;

    bool showtree = true;

    // Declare Formulas
    double U = u(v, dt);
    double D = d(v, dt);

    double D_UP = p_up(r, v, dt);
    double D_DOWN = p_dn(r, v, dt);
    double D_M = p_m(r, v, dt);

    // Declare Tree
    int rows = 4*nodes + 2;
    int cols = nodes + 1;

    double ** tree = new double*[rows];
    double ** treeB = new double*[rows];
    double ** treeC = new double*[rows];
    for(int i = 0; i < rows; ++i){
        tree[i] = new double[cols];
        treeB[i] = new double[cols];
        treeC[i] = new double[cols];
    }

    // Indexes middle of tree array which holds the initial stock price
    int split = rows / 2 - 1;

    // Used to calculate delta and gamma
    double ds = S*0.05;

    tree = C(tree, rows, cols, split, S, K, r, v, dt, U, D, D_UP, D_M, D_DOWN, opType);
    double opPriceA = OptionPrice(tree, split);
    treeB = C(treeB, rows, cols, split, S+ds, K, r, v, dt, U, D, D_UP, D_M, D_DOWN, opType);
    double opPriceB = OptionPrice(treeB, split);
    treeC = C(treeC, rows, cols, split, S-ds, K, r, v, dt, U, D, D_UP, D_M, D_DOWN, opType);
    double opPriceC = OptionPrice(treeC, split);

    // Display the inputs
    std::cout << "Stock Price: " << S << std::endl;
    std::cout << "Strike Price: " << K << std::endl;
    std::cout << "RiskFree Rate: " << r << std::endl;
    std::cout << "Volatility: " << v << std::endl;
    std::cout << "Expiry: " << t << std::endl;
    std::cout << "Option Type: " << opType << std::endl;

    // Option Price
    std::cout << "Option Price: " << opPriceA << std::endl;

    // Calculate Delta
    double delta = (opPriceB - opPriceA)/ds;
    std::cout << "Option Delta: " << delta << std::endl;

    // Calculate Gamma
    double gamma = (opPriceB - 2*opPriceA + opPriceC) / pow(ds, 2);
    std::cout << "Option Gamma: " << gamma << std::endl;

    // Print Tree
    if(showtree == true){
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j){
                if(tree[i][j] == 0){
                    std::cout << "\t";
                } else {
                    std::cout << round(tree[i][j]*100)/100 << "\t";
                }
            }
            std::cout << std::endl;
        }
    }
    return 0;
}

