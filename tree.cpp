#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <algorithm>

// Rounds two decimal places
double round_num(double x){
    return round(x*100) / 100;
}

// Prints options tree
void PrintM(double ** Z, int m, int n)
{
    for(int i = 0; i < m; ++i){
        for(int j = 0; j < n; ++j){
            if(Z[i][j] == 0){
                std::cout << "\t";
            } else {
                std::cout << Z[i][j] << "\t";
            }
        }
        std::cout << std::endl;
    }
}

// Calculates the simulated stock price at each node
double ** FeedForward(double ** Z, double S, int s, int m, int n, double U, double D)
{
    // Initial variables which loop through array
    int upx = 2;
    int dnx = 2;
    int cs = 0;

    while(cs <= n){
        // Setting the stock price at each iteration
        Z[s][cs] = S; 
        for(int i = cs + 1; i < n; ++i){
            // Values of up and down help create the treelike structure
            Z[s - (i - cs)*upx][i] = round_num(Z[s - (i - cs - 1)*upx][i - 1]*U);
            Z[s + (i - cs)*dnx][i] = round_num(Z[s + (i - cs -1)*dnx][i - 1]*D);
        }
        // Increment tree
        cs += 2;
    }

    return Z;
}

// Calculates the payoffs at the last column before discounting
double ** Payoff(double ** Z, double K, int m, int n, std::string opType){
    // Increment variable
    int b = 0;

    while(b <= m){
        // Stock - Strike if Call else Strike - Stock
        if(opType == "call"){
            Z[b + 1][n - 1] = round_num(std::max(Z[b][n - 1] - K, 0.0));
        } else {
            Z[b + 1][n - 1] = round_num(std::max(K - Z[b][n - 1], 0.0));
        }
        b += 4;
    }

    return Z;
}

// Discounts the simulated option payoffs and discounts values back to the current option price
double ** Discount(double ** Z, double P, double nP, double r, double dt, int m, int n, std::string opType){
    // Increment Variables
    int cx = 0;
    int ls = 0;
    int osx = 0;

    for(int i = 1; i < n; ++i){
        // Loops from last column to first
        cx = n - (i + 1);
        ls = 4 + osx;
        while(ls < m - osx){
            // Weighted discount function
            if(opType == "call"){
                Z[ls - 1][cx] = round_num(std::max(exp(-r*dt)*(P*Z[ls - 3][cx + 1] + nP*Z[ls + 1][cx + 1]) , 0.0));
            } else {
                Z[ls - 1][cx] = round_num(std::max(exp(-r*dt)*(nP*Z[ls - 3][cx + 1] + P*Z[ls + 1][cx + 1]) , 0.0));
            }
            ls += 4;
        }
        osx += 2;
    }
 
    return Z;
}

int main()
{
    // Initial Inputs

    double S = 100;              // Stock Price
    double K = 105;              // Strike Price
    double r = 0.02;             // Risk-Free Rate
    double v = 0.10;             // Volatility
    double T = 6.0/12.0;         // Maturity
    std::string opType = "call"; // Option Type


    int n = 10;                   // Number of Steps
    double dt = T/(double) n;    // Time divided by 'n'

    double U = exp(v*sqrt(dt));  // Upper Step Function
    double D = exp(-v*sqrt(dt)); // Lower Step Function

    double P = (exp(r*dt) - D) / (U - D); // Discount Function
    double nP = 1 - P;                    // Opposite Discount

    int Row = 4*n + 2; // Number of Rows in Array
    int Col = n + 1;   // Number of Columns in Array

    int split = Row / 2 - 1; // Fetches middle row

    // Build the double array which will hold the tree
    double ** X = new double*[Row];
    for(int i = 0; i < Row; ++i){
        X[i] = new double[Col];
    }

    // Place stock price in array
    X[split][0] = S;

    // Feed Forward Tree Calculations
    X = FeedForward(X, S, split, Row, Col, U, D);
    
    // Payoff Equations
    X = Payoff(X, K, Row, Col, opType);
    
    // Backwards Discounting
    X = Discount(X, P, nP, r, dt, Row, Col, opType);

    // Print the Option Price
    double opPrice = X[split + 1][0];
    std::cout << "Option Price: " << opPrice << std::endl;

    // Print the Options Tree
    PrintM(X, Row, Col);

    return 0;
}


