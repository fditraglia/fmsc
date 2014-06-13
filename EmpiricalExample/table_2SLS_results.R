#Frank DiTraglia
#Last Updated: June 13th, 2014

#This script constructs a table of 2SLS Results for all instrument sets considered in my empirical example. It should not be run on its own: it is called by run_empirical_example.R

#There is a very slight discrepancy (a difference of 0.01 in a few places) between these results and those in the original version of my paper. This comes from the way I've chosen to handle missing observations: I now exclude Vietnam from the dataset completely, since we only observe the baseline instruments for this country and I want to hold everything constant except the instruments in this exercise. In constrast CG use Vietnam when it is available, as did I in the original version of the empirical example. Again, the difference in results is miniscule.

