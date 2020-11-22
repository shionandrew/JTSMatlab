#include <cmath>
#include <cstddef>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <sstream>
#include <stdio.h>
#include <sys/stat.h>
#include <vector>
#include <unistd.h>

using namespace std;
//string TIMESTEP1 = "4e-5";
//string TIMESTEP2 = "5e-4";
//string file1 = ".00004";
//string file2 = "0005";
string TIMESTEP1 = "Romanova Ldot";
string TIMESTEP2 = "Romanova Mass Accretion";
string file1 = "0";
string file2 = "MassAcc";

/*
 * Produces two overlaid histograms with log scale on x axis. Uses files created by log histogram
*/
void comparisonLogHistograms(string flargeRM, string flargeRStar, string fsmallRM, string fsmallRStar) {
    FILE* gp3=popen("gnuplot -persistent","w");

    fprintf(gp3, "%s \n", "set terminal postscript eps enhanced color font 'Helvetica,10'");
    fprintf(gp3, "%s%s%s%s%s \n", "set output \'comparison_",TIMESTEP1.c_str(),"_",TIMESTEP2.c_str(),".eps'");
    fprintf(gp3, "%s\n", "set logscale x");
    fprintf(gp3, "%s%s %s \n", "set multiplot layout 2,1 title \"","Period Distribution","\"");
    fprintf(gp3, "%s \n", "set ylabel \"Number of stars\"");
    fprintf(gp3, "%s \n", "unset xlabel");
    fprintf(gp3, "%s%s %s \n", "set label 1\"","m < 0.25 solar mass","\" at graph 0.03,0.9");
    fprintf(gp3, "%s%s%s%s%s %s%s%s%s%s \n", "plot \'",fsmallRM.c_str(),"\' using 1:2:3:($1*0.4) with boxerrorbars title \'",TIMESTEP1.c_str(),"\',",
                                          "\'",fsmallRStar.c_str(),"\' using 1:2:($1*0.4) with boxes title \'",TIMESTEP2.c_str(),"\'");
    fprintf(gp3, "%s \n", "set xlabel");
    fprintf(gp3, "%s \n", "set xlabel \"Period (days)\"");
    fprintf(gp3, "%s%s %s \n", "set label 1\"","m > 0.25 solar mass","\" at graph 0.03,0.9");
    fprintf(gp3, "%s%s%s%s%s %s%s%s%s%s \n", "plot \'",flargeRM.c_str(),"\' using 1:2:3:($1*0.4) with boxerrorbars title \'",TIMESTEP1.c_str(),"\',",
                                          "\'",flargeRStar.c_str(),"\' using 1:2:($1*0.4) with boxes title \'",TIMESTEP2.c_str(),"\'");
    fprintf(gp3, "%s \n", "unset multiplot");
}

void periodComparison(string fRM, string fRStar) {
    FILE* gp3=popen("gnuplot -persistent","w");

    fprintf(gp3, "%s \n", "set terminal postscript eps enhanced color font 'Helvetica,10'");
    fprintf(gp3, "%s \n", "set output 'period_comparison.eps'");
    fprintf(gp3, "%s \n", "set xlabel 'Time (Myr)'");
    fprintf(gp3, "%s \n", "set ylabel 'Period (days)'");

    fprintf(gp3, "%s'%s'%s '%s' \n", "plot ",fRM.c_str(),",",fRStar.c_str());
}

int main() {
    comparisonLogHistograms("Period" + file1 + "Large.temp","Period" + file2 + "Large.temp","Period" + file1 + "Small.temp","Period" + file2 + "Small.temp");
    //periodComparison("Period005.temp","Period001.temp");
}