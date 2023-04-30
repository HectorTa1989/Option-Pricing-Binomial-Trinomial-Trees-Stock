#include <iostream>
#include <math.h>
#include <string>
#include <cmath>

// Calculates the UP factor
auto u = [](double v, double dt)
{
    return exp(v*sqrt(2.0*dt));
};

// Calculates the DOWN factor
auto d = [](double v, double dt)
{
    return 1.0 / u(v, dt);
};

// Calculates probability up
auto p_up = [](double r, double v, double dt)
{
    double top = exp(r*dt/2.0) - exp(-v*sqrt(dt/2.0));
    double bot = exp(v*sqrt(dt/2.0)) - exp(-v*sqrt(dt/2.0));
    return pow(top/bot, 2);
};

// Calculates probability down
auto p_dn = [](double r, double v, double dt)
{
    double top = exp(v*sqrt(dt/2.0)) - exp(r*dt/2.0);
    double bot = exp(v*sqrt(dt/2.0)) - exp(-v*sqrt(dt/2.0));
    return pow(top/bot, 2);
};

// Counts the middle probability
auto p_m = [](double r, double v, double dt)
{
    return 1.0 - p_up(r, v, dt) - p_dn(r, v, dt);
};

// Computes the option price for a call or put
double C(int rows, int cols, int split, double S, double K, double r, double v, double dt, double U, double D, double D_UP, double D_M, double D_DOWN, std::string opType)
{   
    // Declare tree array
    double ** tree = new double*[rows];
    for(int i = 0; i < rows; ++i){
        tree[i] = new double[cols];
    }

    // Set indexors
    int ux = 2;
    int cs = 0;

    // Forward propigation
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

    // Returns your option price
    return tree[split + 1][0];
}




int main()
{
    // Declare Inputs
    double S = 100.0; // Stock Price
    double K = 105.0; // Strike Price
    double r = 0.05; // RiskFree Rate
    double t = 30.0/365.0; // Expiry
    double option_price = 1.60; // Option price input to compute IVol
    std::string opType = "call"; // Option type

    int nodes = 1000; // Size of tree
    double dt = t / (double) nodes;

    // Declare Tree
    int rows = 4*nodes + 2;
    int cols = nodes + 1;

    // Indexes middle of tree array which holds the initial stock price
    int split = rows / 2 - 1;

    // Used to calculate delta and gamma
    double ds = S*0.05;

    // Display the inputs
    std::cout << "Stock Price: " << S << std::endl;
    std::cout << "Strike Price: " << K << std::endl;
    std::cout << "RiskFree Rate: " << r << std::endl;
    std::cout << "Expiry: " << t << std::endl;
    std::cout << "Option Type: " << opType << std::endl;
    std::cout << "Nodes: " << nodes << std::endl;
    std::cout << std::endl;
    
    // Set the implied volatility calculation parameters
    double diff = 0, v0 = 0.1, v1 = 1.0, mkt = option_price, vega = 0;
    
    // Price and probability variables
    double pA, pB, pC;

    double U, D, D_UP, D_DOWN, D_M;
    
    // The change for vega
    double dh = 0.001;

    while(true){
        // Calculate UP/DOWN and Discount Probabilities for the regular option price
        U = u(v0, dt);
        D = d(v0, dt);

        D_UP = p_up(r, v0, dt);
        D_DOWN = p_dn(r, v0, dt);
        D_M = p_m(r, v0, dt);

        // Regular option price w/ v0 (volatility) as iterated
        pA = C(rows, cols, split, S, K, r, v0, dt, U, D, D_UP, D_M, D_DOWN, opType);

        // Calculate UP/DOWN and Probabilities for the upper vega formula
        U = u(v0+dh, dt);
        D = d(v0+dh, dt);

        D_UP = p_up(r, v0+dh, dt);
        D_DOWN = p_dn(r, v0+dh, dt);
        D_M = p_m(r, v0+dh, dt);

        // Upper vega price
        pB = C(rows, cols, split, S, K, r, v0+dh, dt, U, D, D_UP, D_M, D_DOWN, opType);

        // Calculate UP/DOWN and Probabilities for the lower vega formula
        U = u(v0-dh, dt);
        D = d(v0-dh, dt);

        D_UP = p_up(r, v0-dh, dt);
        D_DOWN = p_dn(r, v0-dh, dt);
        D_M = p_m(r, v0-dh, dt);

        // Lower vega price
        pC = C(rows, cols, split, S, K, r, v0-dh, dt, U, D, D_UP, D_M, D_DOWN, opType);

        // Calculate the differences and create a fixed point iterator formula
        diff = pA - mkt;
        
        // Compute vega
        vega = (pB - pC)/(2.0*dh);

        // Calculate each volatility iteration
        v1 = v0 - diff / vega;

        // Close loop if v0 & v1 are approximately the same
        if(abs(v1 - v0) <= 0.00001) {
            break;
        } else {
            v0 = v1;
        }

        
    }

    // Calculates the options price for the generated implied volatility measure
    double sim_price = C(rows, cols, split, S, K, r, v1, dt, U, D, D_UP, D_M, D_DOWN, opType);

    // Displays post results
    std::cout << "Market Price: " << mkt << std::endl;
    std::cout << "Option Price: " << sim_price << std::endl;
    std::cout << "Implied Vol: " << v1 << std::endl;

    return 0;
}

