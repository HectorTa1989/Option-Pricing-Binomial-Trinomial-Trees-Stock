# Binomial/Trinomial Option Pricing Algos (no dividends)
These two programs calculate an options price by using the Binomial Tree and Trinomial Tree algorithms. A series of inputs are taken and the tree is calculated in a feed forward way. The tree then calculates the spread between the Stock/Strike or Strike/Stock depending on the option price. The payoffs are discounted backwards through the tree (time based) and equate to the options price. I got most of the equations from the paper https://warwick.ac.uk/fac/sci/maths/people/staff/oleg_zaboronski/fm/trinomial_tree_2009.pdf.

## Compiling
For the binomial tree you should compile it as ``` g++ binomial.cpp -std=c++14```. For the trinomial tree you should complie it as ``` g++ trinomial.cpp -std=c++14```. Once you have these compiled simply run ./a.out.

## Running Binomial Tree
![alt](https://github.com/marscolony2040/Option-Pricing-Binomial-Trinomial-Trees/blob/main/images/binomial.png)

## Running Trinomial Tree
![alt](https://github.com/marscolony2040/Option-Pricing-Binomial-Trinomial-Trees/blob/main/images/trinomial.png)
