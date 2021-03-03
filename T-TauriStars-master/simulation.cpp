/**
 * \file simulation.cpp
 *
 * \authors Cynthia Yan ('18'), Shion Andrew ('21), Mia Taylor ('22)
 *
 * \brief Provides the main() function for simulating a ttauri star cluster or one star; defines other helper functions not associated with TTauriStar class.
 */

#include "ttauristar.hpp"
#include <stdio.h>
using namespace std;
#include <string>
#include <algorithm>

bool const RAND_DIST = true; // False for (testing) NO random distribution of Bfield and Mdot factor, True (normal mode) for random distribution
int const TotalSimulations = 1;
static const vector<double> TIMESTEPS {0.05, 0.05}; // .005 the timestep we want
static const vector<bool> PropellerTests = {true, true};


/**
 * \brief holds period and phase (1,2,3 or 4) of each simulated star
 */
struct starData {
	double period;
	int phase;
};

/**
 * \brief computes number of stars in each period bucket
 */
void computeStats(vector<starData> periods, int bucket, int simulationNumber)
{
	// calculate average
	int phase1=0;
	int phase2=0;
	int phase3=0;
	int phase4=0;

  for(size_t bucket=0; bucket<periods.size() ; ++bucket){
		int phase = periods[bucket].phase;
		if(phase <= 1){
			++phase1;
		}
		else if(phase <= 2){
			++phase2;
		}
		else if(phase <= 3){
			++phase3;
		}
		else if(phase <= 4){
			++phase4;
		}
	}

	// save star population distribution per bucket in csv
	ofstream myfile;
	myfile.open ("starStages" + to_string(simulationNumber) + ".csv", ios_base::app);
	myfile << bucket << "," << phase1 << "," << phase2 << "," << phase3 << "," << phase4 << "\n";
	myfile.close();

	//cout << "BUCKET: " << bucket << endl;
	//cout << "SpinUp: " << phase1 << endl;
	//cout << "Propellar: " << phase2 << endl;
	//cout << "Locked: " << phase3 << endl;
	//cout << "Unlocked: " << phase4 << endl;
	//cout << "" << endl;

}


/**
 * \brief sorts periods into buckets and populates a vector with those periods for each bucket
 *
 * \returns vectors of vectors of star data)
 */
vector<vector<starData>> populateBuckets(vector<starData> periods, int numBuckets, int bucketSize, int simulationNumber)
{
	// initialize vector of vectors; each bucket contains vector of periods
	vector<vector<starData>> periodsData;
	for(int i = 0; i < numBuckets; ++i){
		vector<starData> bucket;
		periodsData.push_back(bucket);
	}
	// populate buckets
	for(size_t i = 0; i < periods.size();++i){
		// determine which bucket period fits into
		int index = (int)(periods[i].period/bucketSize);
		// ignoring stars with periods greater than 14 days
		if(index < numBuckets) {
		periodsData[index].push_back(periods[i]);
	}
}

for(size_t bucket = 0; bucket < periodsData.size(); ++bucket){
		if(!periodsData[bucket].empty()) {
			computeStats(periodsData[bucket], bucket, simulationNumber);
		}
}

	for(size_t bucket = 0; bucket < periodsData.size(); ++bucket){
		// save number of stars in each bucket to csv
		ofstream myfile;
		myfile.open ("bucketCount" +  to_string(simulationNumber) + ".csv", ios_base::app);
		myfile << bucket << "," << periodsData[bucket].size()  << "\n";
		myfile.close();
	}
	return periodsData;
}


/**
 * \brief returns stars in periods vector that have specified phase
 *
 * \param starData struct periods to be sorted, phaseNumber specifying desired phase, logDistribution = 0 if periods are to be plotted logarithmically, 1 if not logarithmic
 * \returns subset of starData struct periods that have specified phase
 */
 vector<starData> sortPhase(vector<starData> periods, int phaseNumber)
 {
	vector<starData> sortedPeriods;
	for(size_t star = 0; star<periods.size(); ++star) {
		if(periods[star].phase == phaseNumber){
			sortedPeriods.push_back(periods[star]);
		}
	 }
	return sortedPeriods;
 }



/**
 * \brief plot period histogram,
 *
 * \param starData struct periods1 for stars with mass < 0.25, periods2 for stars with mass > 0.25
 * \returns
 */
void plothistogram(vector<starData> periods1, vector<starData>periods2, int simulationNumber)
{

    vector<double> periodbins1; // stores number of low mass stars in each period bin
    vector<double> periodbins2; // stores number of high mass stars in each period bin
    vector<double> bins; // vector storing lower bound of period bin at corresponding index
    
    // iterate from minimum period to maximum period
    for(double i = 0; i < 20; i += 1) {
        periodbins1.push_back(0);
        periodbins2.push_back(0);
        bins.push_back(i);
    }
    // sorting periods into bins
    for(size_t k=0;k<periods1.size();k++) {
        size_t binIndex = floor(periods1[k].period); // find bin star falls into (rounding down) 
        ++periodbins1[binIndex]; // increase number of stars in this bin
    }
    // repeat for high mass stars
    for(size_t k=0;k<periods2.size();k++) {
        size_t binIndex =  floor(periods2[k].period);
        ++periodbins2[binIndex];
    }

    // printing to files
    FILE * temp1 = fopen("periods1.temp", "w");
    FILE * temp2 = fopen("periods2.temp", "w");
	FILE* gp3=popen("gnuplot -persistent","w");

    for(size_t k=0;k<periodbins1.size();k++) {
        double numPer1 = periodbins1[k]; // number of stars in bin
        // save period at bottom of bin, N stars in bin, and error (sqrt[N])
        fprintf(temp1,"%f %f %f \n",bins[k]+.5,numPer1,pow(numPer1,.5)); //.5 needed because of gnuplot centering...
    }
    for(size_t k=0;k<periodbins2.size();k++) {
        double numPer2 = periodbins2[k];
        fprintf(temp2,"%f %f %f \n",bins[k]+.5,numPer2,pow(numPer2,.5));
    }

	// write periods to file for gnuplot
    //for(size_t k=0;k<periods1.size();k++) {
    //    fprintf(temp1,"%f %x \n",periods1[k].period, pow(periods1[k].period,.5));
    //}
    //for(size_t k=0;k<periods2.size();k++) {
    //    fprintf(temp2,"%f %x \n",periods2[k].period,  pow(periods2[k].period,.5));
    //}
	//cout << "Number of Stars: " << periods1.size()+periods2.size() << endl;

    fprintf(gp3, "%s \n", "set terminal postscript eps enhanced color font 'Helvetica,10'");
    fprintf(gp3, "%s \n", "set style histogram errorbars linewidth 1");
    fprintf(gp3, "%s%s%s \n", "set output 'distribution", to_string(PropellerTests[simulationNumber]).c_str(), ".eps'");
    fprintf(gp3, "%s\n", "binwidth=1");
    fprintf(gp3, "%s\n", "set boxwidth binwidth");
    fprintf(gp3, "%s\n", "bin(x,width)=width*floor(x/width) + binwidth/2.0");
    fprintf(gp3, "%s%s %s \n", "set multiplot layout 2,1 title \"","Period Distribution","\"");
    fprintf(gp3, "%s \n", "set ylabel \"Number of stars\"");
    fprintf(gp3, "%s \n", "unset xlabel");
    fprintf(gp3, "%s \n", "set xrange [0:14]");
    fprintf(gp3, "%s%s %s \n", "set label 1\"","m < 0.25 solar mass","\" at graph 0.8,0.9");
    fprintf(gp3, "%s \n", "plot 'periods1.temp' using 1:2:3 with boxerrorbars");
    fprintf(gp3, "%s \n", "set xlabel \"Period (days)\"");
    fprintf(gp3, "%s%s %s \n", "set label 1\"","m > 0.25 solar mass","\" at graph 0.8,0.9");
    fprintf(gp3, "%s \n", "plot 'periods2.temp' using 1:2:3 with boxerrorbars" );
    fprintf(gp3, "%s \n", "unset multiplot");

    fclose(temp1);
    fclose(temp2);
    pclose(gp3);

}



/**
 * \brief plot period histogram with phases of stars color coded
 *
 * \param starData struct periods1 for stars with mass < 0.25, periods2 for stars with mass > 0.25, simulationNumber (if performing multiple simulations)
 * \returns
 */
void plothistogramWithPhase(vector<starData> periods1, vector<starData>periods2, int simulationNumber)
{
    //stringstream folderNameStream;
    //folderNameStream << fixed << setprecision(4);
    //folderNameStream << TIMESTEPS[simulationNumber]; // OUTPUTID*simulationNumber;
    //string folderName = folderNameStream.str();

	// make new folder and change directories into new folder
    //mkdir(folderName.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
    //chdir(folderName.c_str());

	  vector<vector<starData>> periods1Data = populateBuckets(periods1, 14 ,1, simulationNumber);
		vector<vector<starData>> periods2Data = populateBuckets(periods2, 14 ,1, simulationNumber);

		vector<starData> periods1_4 = sortPhase(periods1, 4);
		vector<starData> periods1_3 = sortPhase(periods1, 3);
		vector<starData> periods1_2 = sortPhase(periods1, 2);
		vector<starData> periods1_1 = sortPhase(periods1, 1);

		vector<starData> periods2_4 = sortPhase(periods2, 4);
		vector<starData> periods2_3 = sortPhase(periods2, 3);
		vector<starData> periods2_2 = sortPhase(periods2, 2);
		vector<starData> periods2_1 = sortPhase(periods2, 1);


    FILE * temp1 = fopen("periods1.temp", "w");
		FILE * temp1_4 = fopen("periods1_4.temp", "w");
		FILE * temp1_3 = fopen("periods1_3.temp", "w");
		FILE * temp1_2 = fopen("periods1_2.temp", "w");
		FILE * temp1_1 = fopen("periods1_1.temp", "w");

    FILE * temp2 = fopen("periods2.temp", "w");
		FILE * temp2_4 = fopen("periods2_4.temp", "w");
		FILE * temp2_3 = fopen("periods2_3.temp", "w");
		FILE * temp2_2 = fopen("periods2_2.temp", "w");
		FILE * temp2_1 = fopen("periods2_1.temp", "w");

		FILE* gp3=popen("gnuplot -persistent","w");

		// write periods to file for gnuplot
    for(size_t k=0;k<periods1.size();k++) {
        fprintf(temp1,"%f %x \n",periods1[k].period, periods1[k].phase);
    }

		for(size_t k=0;k<periods1_4.size();k++) {
				fprintf(temp1_4,"%f %x \n",periods1_4[k].period, periods1_4[k].phase);
		}

		for(size_t k=0;k<periods1_3.size();k++) {
				fprintf(temp1_3,"%f %x \n",periods1_3[k].period, periods1_3[k].phase);
		}

		for(size_t k=0;k<periods1_2.size();k++) {
				fprintf(temp1_2,"%f %x \n",periods1_2[k].period, periods1_2[k].phase);
		}

		for(size_t k=0;k<periods1_1.size();k++) {
				fprintf(temp1_1,"%f %x \n",periods1_1[k].period, periods1_1[k].phase);
		}

    for(size_t k=0;k<periods2.size();k++) {
        fprintf(temp2,"%f %x \n",periods2[k].period,  periods2[k].phase);
    }

		for(size_t k=0;k<periods2_4.size();k++) {
        fprintf(temp2_4,"%f %x \n",periods2_4[k].period,  periods2_4[k].phase);
    }

		for(size_t k=0;k<periods2_3.size();k++) {
        fprintf(temp2_3,"%f %x \n",periods2_3[k].period,  periods2_3[k].phase);
    }
		for(size_t k=0;k<periods2_2.size();k++) {
				fprintf(temp2_2,"%f %x \n",periods2_2[k].period,  periods2_2[k].phase);
		}
		for(size_t k=0;k<periods2_1.size();k++) {
				fprintf(temp2_1,"%f %x \n",periods2_1[k].period,  periods2_1[k].phase);
		}

		//cout << "Number of Stars: " << periods1.size()+periods2.size() << endl;

    fprintf(gp3, "%s \n", "set terminal postscript eps enhanced color font 'Helvetica,10'");
    fprintf(gp3, "%s \n", "set output 'distributionPhase.eps'");
    fprintf(gp3, "%s\n", "binwidth=1");
    fprintf(gp3, "%s\n", "set boxwidth binwidth");
    fprintf(gp3, "%s\n", "bin(x,width)=width*floor(x/width) + binwidth/2.0");
    fprintf(gp3, "%s%s %s \n", "set multiplot layout 2,1 title \"","Period Distribution","\"");
    fprintf(gp3, "%s \n", "set ylabel \"Number of stars\"");
    fprintf(gp3, "%s \n", "unset xlabel");
    fprintf(gp3, "%s \n", "set xrange [0:14]");
    fprintf(gp3, "%s%s %s \n", "set label 1\"","m < 0.25 solar mass","\" at graph 0.6,0.9");
    fprintf(gp3, "%s \n", "plot 'periods1.temp' using (bin($1,binwidth)):(1.0) smooth freq with boxes title 'Spin up', 'periods1_4.temp' using (bin($1,binwidth)):(1.0) smooth freq with boxes title 'Unlocked'");
    fprintf(gp3, "%s \n", "set xlabel \"Period (days)\"");
    fprintf(gp3, "%s%s %s \n", "set label 1\"","m > 0.25 solar mass","\" at graph 0.6,0.9");
    fprintf(gp3, "%s \n", "plot 'periods2.temp' using (bin($1,binwidth)):(1.0) smooth freq with boxes title 'Spin up', 'periods2_4.temp' using (bin($1,binwidth)):(1.0) smooth freq with boxes title 'Propeller', 'periods2_3.temp' using (bin($1,binwidth)):(1.0) smooth freq with boxes title 'Locked', 'periods2_2.temp' using (bin($1,binwidth)):(1.0) smooth freq with boxes title 'Unlocked'");
    fprintf(gp3, "%s \n", "unset multiplot");

    fclose(temp1);
    fclose(temp1_4);
    fclose(temp1_3);
    fclose(temp1_2);
    fclose(temp1_1);
    fclose(temp2);
    fclose(temp2_4);
    fclose(temp2_3);
    fclose(temp2_2);
    fclose(temp2_1);
    pclose(gp3);

}

/**
 * \brief plot period histogram with log scale on the x axis,
 *
 * \param periods1 for stars with mass < 0.25 and periods2 for stars with mass > 0.25
 * \returns
 * Added 4/5/20 by Mia Taylor
 * gnuplot commands largely based on this:
 * https://stackoverflow.com/questions/24207850/gnuplot-histogram-x-logscale
 */
void plotlogscalehistogram(vector<double> periods1, vector<double> periods2, int simulationNumber)
{
    string str1 = "Period" + to_string(PropellerTests[simulationNumber]) + "Small.temp";
    const char* F1 = str1.c_str();
    string str2 = "Period" + to_string(PropellerTests[simulationNumber]) + "Large.temp";
    const char* F2 = str2.c_str();


    //FILE * temp1 = fopen("periods1loghist.temp", "w");
    //FILE * temp2 = fopen("periods2loghist.temp", "w");
    FILE * temp1 = fopen(F1, "w");
    FILE * temp2 = fopen(F2, "w");
    FILE* gp3=popen("gnuplot -persistent","w");

    double binsPerDecade = 10;
    double intervalWidth = pow(10,1/binsPerDecade);

    // determining the range of the periods
    double minPeriod = min(*min_element(periods1.begin(), periods1.end()),
                          *min_element(periods2.begin(), periods2.end()));
    
    double maxPeriod = max(*max_element(periods1.begin(), periods1.end()),
                          *max_element(periods2.begin(), periods2.end()));

    // min and max OOM
    int minOOM = floor(log10(minPeriod));
    int maxOOM = floor(log10(maxPeriod))+1; // round down and then add 1
    // with bins of width 
    vector<double> periodbins1;
    vector<double> periodbins2;
    vector<double> bins;
    
    // iterate from minimum period to maximum period
    for(double i = pow(10,minOOM); i < pow(10,maxOOM); i *= intervalWidth) {
        periodbins1.push_back(0);
        periodbins2.push_back(0);
        bins.push_back(i);
    }

    // sorting periods into bins
    for(size_t k=0;k<periods1.size();k++) {
        size_t binIndex = floor((log10(periods1[k]) - minOOM)*binsPerDecade);
        ++periodbins1[binIndex];
    }

    for(size_t k=0;k<periods2.size();k++) {
        size_t binIndex = floor((log10(periods2[k]) - minOOM)*binsPerDecade);
        ++periodbins2[binIndex];
    }


    // printing to files
    for(size_t k=0;k<periodbins1.size();k++) {
        double numPer1 = periodbins1[k]; // number of stars in bin
        double numPer2 = periodbins2[k];
        // save period at bottom of bin, N stars in bin, and error (sqrt[N])
        fprintf(temp1,"%f %f %f \n",bins[k],numPer1,pow(numPer1,.5));
        fprintf(temp2,"%f %f %f \n",bins[k],numPer2,pow(numPer2,.5));

    }
    cout << "distributionlogscale" << to_string(PropellerTests[simulationNumber]).c_str() << endl;
    fprintf(gp3, "%s \n", "set terminal postscript eps enhanced color font 'Helvetica,10'");
    fprintf(gp3, "%s%s%s \n", "set output 'distributionlogscale", to_string(PropellerTests[simulationNumber]).c_str(), ".eps'");

    fprintf(gp3, "%s\n", "set logscale x");
    // fprintf(gp3, "%s\n", "set boxwidth 1.0");

    fprintf(gp3, "%s%s %s \n", "set multiplot layout 2,1 title \"","Period Distribution","\"");
    fprintf(gp3, "%s \n", "set style histogram errorbars linewidth 1");
    fprintf(gp3, "%s \n", "set ylabel \"Number of stars\"");
    fprintf(gp3, "%s \n", "unset xlabel");
    fprintf(gp3, "%s%s %s \n", "set label 1\"","m < 0.25 solar mass","\" at graph 0.8,0.9");
    fprintf(gp3, "%s%s%s \n", "plot '",F1,"' using 1:2:3:($1*0.2) with boxerrorbars"); //adjust boxwidth on plot here ($1*0.2)
    fprintf(gp3, "%s \n", "set xlabel");
    fprintf(gp3, "%s \n", "set xlabel \"Period (days)\"");
    fprintf(gp3, "%s%s %s \n", "set label 1\"","m > 0.25 solar mass","\" at graph 0.8,0.9");
    fprintf(gp3, "%s%s%s \n", "plot '",F2,"' using 1:2:3:($1*0.2) with boxerrorbars"); //adjust boxwidth on plot here
    fprintf(gp3, "%s \n", "unset multiplot");
    cout << "PLOTTING GRAPH" << endl;
    fclose(temp1);
    fclose(temp2);
    pclose(gp3);

}


/**
 * \brief plot period histogram with log scale on the x axis, color codes locked vs unlocked stars
 *
 * \param periods1 for stars with mass < 0.25 and periods2 for stars with mass > 0.25
 * \returns
 * Added 4/5/20 by Mia Taylor, edited by Shion Andrew
 * gnuplot commands largely based on this:
 * https://stackoverflow.com/questions/24207850/gnuplot-histogram-x-logscale
 */
 void plotlogscalehistogramWithPhase(vector<starData> periods1, vector<starData>periods2)
{

	// sort stars by phase
	vector<starData> periods1_4 = sortPhase(periods1, 4);
	vector<starData> periods1_3 = sortPhase(periods1, 3);
	vector<starData> periods1_2 = sortPhase(periods1, 2);
	vector<starData> periods1_1 = sortPhase(periods1, 1);
	vector<starData> periods2_4 = sortPhase(periods2, 4);
	vector<starData> periods2_3 = sortPhase(periods2, 3);
	vector<starData> periods2_2 = sortPhase(periods2, 2);
	vector<starData> periods2_1 = sortPhase(periods2, 1);

	// Create vector for locked stars (M < 0.25)
	vector<double> period1_lockedVector;
	for(size_t k = 0; k< periods1_1.size(); ++k){
		period1_lockedVector.push_back(periods1_1[k].period);
	}
	for(size_t k = 0; k< periods1_2.size(); ++k){
		period1_lockedVector.push_back(periods1_2[k].period);
	}
	for(size_t k = 0; k< periods1_3.size(); ++k){
		period1_lockedVector.push_back(periods1_3[k].period);
	}

	// Create vector for locked stars (M > 0.25)
	vector<double> period2_lockedVector;
	for(size_t k = 0; k< periods2_1.size(); ++k){
		period2_lockedVector.push_back(periods2_1[k].period);
	}
	for(size_t k = 0; k< periods2_2.size(); ++k){
		period2_lockedVector.push_back(periods2_2[k].period);
	}
	for(size_t k = 0; k< periods2_3.size(); ++k){
		period2_lockedVector.push_back(periods2_3[k].period);
	}

	// group unlocked stars (M < 0.25)
	vector<double> period1_unlockedVector;
	for(size_t k = 0; k< periods1_4.size(); ++k){
		period1_unlockedVector.push_back(periods1_4[k].period);
	}

	// group unlocked stars (M > 0.25)
	vector<double> period2_unlockedVector;
	for(size_t k = 0; k< periods2_4.size(); ++k){
		period2_unlockedVector.push_back(periods2_4[k].period);
	}

  // in case either vector is empty, populate it with a dummy star so minPeriod and maxPeriod doesn't return a segmentation error
	if(period1_lockedVector.size()==0){
		period1_lockedVector.push_back(1);
	}
  if(period2_lockedVector.size()==0){
		period2_lockedVector.push_back(1);
	}
	if(period1_unlockedVector.size() == 0){
		period1_unlockedVector.push_back(1);
	}
	if(period2_unlockedVector.size() == 0){
		period2_unlockedVector.push_back(1);
	}

	FILE * temp1_locked = fopen("lockedPeriods1.temp", "w");
	FILE * temp1_unlocked = fopen("unlockedPeriods1.temp", "w");

	FILE * temp2_locked = fopen("lockedPeriods2.temp", "w");
	FILE * temp2_unlocked = fopen("unlockedPeriods2.temp", "w");

	FILE* gp3=popen("gnuplot -persistent","w");

  double binsPerDecade = 5;
  double intervalWidth = pow(10,1/binsPerDecade);
  // determining the range of the periods

  double minPeriod = min(*min_element(period1_lockedVector.begin(), period1_lockedVector.end()),
                          min( min(*min_element(period1_unlockedVector.begin(), period1_unlockedVector.end()),
												*min_element(period2_unlockedVector.begin(), period2_unlockedVector.end())),
											*min_element(period2_unlockedVector.begin(), period2_unlockedVector.end())));

  double maxPeriod = max(*max_element(period1_lockedVector.begin(), period1_lockedVector.end()),
                        max( max(*max_element(period1_unlockedVector.begin(), period1_unlockedVector.end()),
											*min_element(period2_unlockedVector.begin(), period2_unlockedVector.end())),
										*min_element(period2_unlockedVector.begin(), period2_unlockedVector.end())));

  // min and max OOM
  int minOOM = floor(log10(minPeriod));
  int maxOOM = floor(log10(maxPeriod))+1; // round down and then add 1

  // with bins of width
  vector<double> periodbins1_locked;
  vector<double> periodbins1_unlocked;
  vector<double> periodbins2_locked;
  vector<double> periodbins2_unlocked;
  vector<double> bins;

  for(double i = pow(10,minOOM); i < pow(10,maxOOM); i *= intervalWidth) {
      periodbins2_locked.push_back(0);
      periodbins2_unlocked.push_back(0);
			periodbins1_locked.push_back(0);
			periodbins1_unlocked.push_back(0);
      bins.push_back(i);
  }
    // calculate bin that each star belongs to and increment stars in bin (periodsbins)
	for(size_t k=0;k<period1_lockedVector.size();k++) {
			size_t binIndex = floor((log10(period1_lockedVector[k]) - minOOM)*binsPerDecade);
            
			++periodbins1_locked[binIndex]; 
	}

	for(size_t k=0;k<period1_unlockedVector.size();k++) {
      size_t binIndex = floor((log10(period1_unlockedVector[k]) - minOOM)*binsPerDecade);

      ++periodbins1_unlocked[binIndex];
  }

  for(size_t k=0;k<period2_lockedVector.size();k++) {
      size_t binIndex = floor((log10(period2_lockedVector[k]) - minOOM)*binsPerDecade);
      ++periodbins2_locked[binIndex];

  }

  for(size_t k=0;k<period2_unlockedVector.size();k++) {
      size_t binIndex = floor((log10(period2_unlockedVector[k]) - minOOM)*binsPerDecade);
      ++periodbins2_unlocked[binIndex];
  }

	// printing to files
  for(size_t k=0;k<periodbins1_locked.size();k++) {
      fprintf(temp1_locked,"%f %f \n",bins[k],periodbins1_locked[k]);
  }
	// printing to files
	for(size_t k=0;k<periodbins1_unlocked.size();k++) {
			fprintf(temp1_unlocked,"%f %f \n",bins[k],periodbins1_unlocked[k]);
	}

  // printing to files
  for(size_t k=0;k<periodbins2_locked.size();k++) {
      fprintf(temp2_locked,"%f %f \n",bins[k],periodbins2_locked[k]);
  }
	// printing to files
	for(size_t k=0;k<periodbins2_unlocked.size();k++) {
			fprintf(temp2_unlocked,"%f %f \n",bins[k],periodbins2_unlocked[k]);
	}


  fprintf(gp3, "%s \n", "set terminal postscript eps enhanced color font 'Helvetica,10'");
  fprintf(gp3, "%s \n", "set output 'distributionlogscale.eps'");
  fprintf(gp3, "%s \n", "set terminal postscript eps enhanced color font 'Helvetica,10'");
  fprintf(gp3, "%s \n", "set output 'distributionlogscale.eps'");
  fprintf(gp3, "%s\n", "set logscale x");
  fprintf(gp3, "%s%s %s \n", "set multiplot layout 2,1 title \"","Period Distribution","\"");
  fprintf(gp3, "%s \n", "set ylabel \"Number of stars\"");
  fprintf(gp3, "%s \n", "unset xlabel");
  fprintf(gp3, "%s%s %s \n", "set label 1\"","m < 0.25 solar mass","\" at graph 0.85,0.8");
  fprintf(gp3, "%s \n", "plot 'lockedPeriods1.temp' using 1:2:($1*0.4) with boxes title 'locked', 'unlockedPeriods1.temp' using 1:2:($1*0.4) with boxes title 'unlocked'");
  fprintf(gp3, "%s \n", "set xlabel");
  fprintf(gp3, "%s \n", "set xlabel \"Period (days)\"");
  fprintf(gp3, "%s%s %s \n", "set label 1\"","m > 0.25 solar mass","\" at graph 0.85,0.8");
  fprintf(gp3, "%s \n", "plot 'lockedPeriods2.temp' using 1:2:($1*0.4) with boxes title 'locked', 'unlockedPeriods2.temp' using 1:2:($1*0.4) with boxes title 'unlocked'");
  fprintf(gp3, "%s \n", "unset multiplot");

  fclose(temp1_locked);
  fclose(temp1_unlocked);
  fclose(temp2_locked);
  fclose(temp2_unlocked);
  pclose(gp3);
  chdir(".."); 

}


/**
 * \brief plot propeller strength histogram,
 *
 * \param starData struct periods1 for stars with mass < 0.25, periods2 for stars with mass > 0.25, simulationNumber (if performing multiple simulations)
 * \returns
 */
void plotPropellerStrengthhistogram(vector<double> PropellerStrengths1, vector<double> PropellerStrengths2, int simulationNumber)
{


    FILE * temp1 = fopen("propellerStrengths1.temp", "w");
    FILE * temp2 = fopen("propellerStrengths2.temp", "w");


	FILE* gp3=popen("gnuplot -persistent","w");

	// write periods to file for gnuplot
    for(size_t k=0;k<PropellerStrengths1.size();k++) {
        fprintf(temp1,"%f \n",PropellerStrengths1[k]);
    }
    for(size_t k=0;k<PropellerStrengths2.size();k++) {
        fprintf(temp2,"%f \n",PropellerStrengths2[k]);
    }

    fprintf(gp3, "%s \n", "set terminal postscript eps enhanced color font 'Helvetica,10'");
    fprintf(gp3, "%s%s%s \n", "set output 'PropellerStrengthDistribution ", to_string(simulationNumber).c_str(), ".eps'");
    fprintf(gp3, "%s\n", "binwidth=1");
    fprintf(gp3, "%s\n", "set boxwidth binwidth");
    fprintf(gp3, "%s\n", "bin(x,width)=width*floor(x/width) + binwidth/2.0");
    fprintf(gp3, "%s%s %s \n", "set multiplot layout 2,1 title \"","Propeller Strength Distribution","\"");
    fprintf(gp3, "%s \n", "set ylabel \"Number of stars\"");
    fprintf(gp3, "%s \n", "unset xlabel");
    fprintf(gp3, "%s \n", "set xrange [0:50]");
    fprintf(gp3, "%s%s %s \n", "set label 1\"","m < 0.25 solar mass","\" at graph 0.8,0.9");
    fprintf(gp3, "%s \n", "plot 'propellerStrengths1.temp' using (bin($1,binwidth)):(1.0) smooth freq with boxes notitle");
    fprintf(gp3, "%s \n", "set xlabel \"PropellerStrength (v_s_t_a_r / v_d_i_s_k)\"");
    fprintf(gp3, "%s%s %s \n", "set label 1\"","m > 0.25 solar mass","\" at graph 0.8,0.9");
    fprintf(gp3, "%s \n", "plot 'propellerStrengths2.temp' using (bin($1,binwidth)):(1.0) smooth freq with boxes notitle");
    fprintf(gp3, "%s \n", "unset multiplot");

    fclose(temp1);
    fclose(temp2);
    pclose(gp3);
}

/**
 * \brief read cmk data file, Kunitomo data file
 *
 * \param fname      string representing the cmk data filename
 * \returns          a vector of vectors of doubles
 */
vector<vector<double>> readKunitomo()
{
	ifstream inputFile("KunitomoData.txt");
	if (!inputFile.good()) {
        throw invalid_argument( "Couldn't open cmk file for reading" );
    }

    // throw the first line
    string line;
    getline(inputFile, line);

    // allocate space on the heap to store the cmk data as a vector of vectors of doubles
    vector<vector<double>> cmktable;

    vector<double> cmkmasses;
    vector<double> cmkages;
    vector<double> cmkradii;
    vector<double> cmkxi;

    double mass; // Mass in M_sun
    double radius; 
    double age; // age in Myrs
    double xi; 

    // read one line at a time
    while (getline(inputFile, line)){
    	stringstream lineStream(line);
        // read the four doubles separated by spaces
        // and put them in the corresponding vector
        lineStream >> mass;
        lineStream >> age;
        lineStream >> radius;
        lineStream >> xi;
        cmkmasses.push_back(mass);
        cmkages.push_back(age);
        cmkradii.push_back(radius);
        cmkxi.push_back(xi);
        //cout << mass << " " << age << " " << radius << endl;
    }
    // store mass, age, radius vectors in cmktable
    cmktable.push_back(cmkmasses);
    cmktable.push_back(cmkages);
    cmktable.push_back(cmkradii);
    cmktable.push_back(cmkxi);
    // close file
    inputFile.close();
    return cmktable;
}


/**
 * \brief read cmk data file
 *
 * \param fname      string representing the cmk data filename
 * \returns          a vector of vectors of doubles
 */
vector<vector<double>> readcmk(string fname)
{
	ifstream inputFile(fname);

	if (!inputFile.good()) {
        throw invalid_argument( "Couldn't open cmk file for reading" );
    }

    // throw the first line
    string line;
    getline(inputFile, line);

    // allocate space on the heap to store the cmk data as a vector of vectors of doubles
    vector<vector<double>> cmktable;

    vector<double> cmkmasses;
    vector<double> cmkages;
    vector<double> cmkradii;
    vector<double> cmkxi;

    const double PI = 3.141592653589793;
    const double SIGMAB = 5.67e-5;   // stefan-boltzmann constatn in erg/K^4 cm^2 s
    const double SOLARRADIUS = 7e10; // Solar radius in cm
    const double SOLARLUMINOSITY = 3.86e33;  // Solar luminosity in ergs per second
    // Solar tempertaure in Kelvins
    const double SOLARTEMPERATURE = pow(SOLARLUMINOSITY/(4*PI*pow(SOLARRADIUS,2)*SIGMAB),0.25);


    double mass; // Mass in M_sun
    double logage; 
    double age; // age in Myrs
    double logL; // log(Luminosity)
    double logT; // log(Temperature)
    double radius; // R_sun

    // read one line at a time
    while (getline(inputFile, line)){
    	stringstream lineStream(line);
        // read the four doubles separated by spaces
        // and put them in the corresponding vector
        lineStream >> mass;
        lineStream >> logage;
        lineStream >> logL;
        lineStream >> logT;
        // calculate age
        age = pow(10.0,logage - 6.0);
        // calculate radius
        radius = pow(10,logL/2) * pow(pow(10,logT)/SOLARTEMPERATURE,-2);
        cmkmasses.push_back(mass);
        cmkages.push_back(age);
        cmkradii.push_back(radius);
        cmkxi.push_back(0.5);
        //cout << mass << " " << age << " " << radius << endl;
    }
    // store mass, age, radius vectors in cmktable
    cmktable.push_back(cmkmasses);
    cmktable.push_back(cmkages);
    cmktable.push_back(cmkradii);
    cmktable.push_back(cmkxi);
    // close file
    inputFile.close();
    return cmktable;
}



/**
 * \brief read cluster mass and age file
 *
 * \param fname      string representing the cluster filename
 * \returns          a vector of vectors of doubles
 */
vector<vector<double>> readcluster(string fname)
{
    ifstream inputFile(fname);

    if (!inputFile.good()) {
        throw invalid_argument( "Couldn't open cluster file for reading" );
    }

    // throw the first line
    string line;
    getline(inputFile, line);

    // allocate space on the heap to store the cmk data as a vector of vectors of doubles
    vector<vector<double>> startable;

    vector<double> logmasses;
    vector<double> logages;

    double mass;
    double age;

    // read one line at a time
    while (getline(inputFile, line)){
        stringstream lineStream(line);
        // read the four doubles separated by spaces
        // and put them in the corresponding vector
        lineStream >> mass;
        lineStream >> age;
        // push back log10 of masses and ages
        logmasses.push_back(log10(mass));
        if (log10(age) < -3) {
            cout << "mass " << mass << "age " << age << endl;
        }
        logages.push_back(log10(age));
    }
    // store mass, age, radius vectors in cmktable
    startable.push_back(logmasses);
    startable.push_back(logages);

    // close file
    inputFile.close();
    return startable;
}

/**
 * \brief generate the distribution corresponding to a cluster
 * \param fname      vectors of logmasses and logages
 * \param n          size of the simulated star table
 * \returns          a distribution
 */
vector<vector<double>> generatedistribution(vector<vector<double>> startable, double n)
{
    // create mass bins
    vector<double> logmassbins{-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1};
    // number of mass bins
    size_t nummassbin = logmassbins.size() - 1;
    // initialize weights with zeros
    vector<double> logmassweights(nummassbin,0.0);

    // create age bins
    vector<double> logagebins{-1.2,-1,-0.8,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0};
    // number of age bins
    size_t numagebin = logagebins.size() - 1;
    // initilize age weights with zeros
    vector<vector<double>> logageweights;
    for (size_t i = 0; i < nummassbin; ++i) {
        logageweights.push_back(vector<double>(numagebin,0.0));
    }
    // iterate through startable to get the weight
    // get logmass and logage from the table
    vector<double> logmasses = startable[0];
    vector<double> logages = startable[1];
    // the number of stars in the table
    size_t numstar = logmasses.size();
    // iterate over the stars in the cluster
    for (size_t i = 0; i < numstar; ++i) {
        // find the interval that the logmass belongs to
        for (size_t j = 0; j < nummassbin; ++j) {
            if (logmassbins[j] <= logmasses[i] && logmasses[i] < logmassbins[j+1]) {
                // increase the mass count in the appropriate interval
                logmassweights[j] += 1.0;
                // look at its age
                for (size_t k = 0; k < numagebin; ++k) {
                    if (logagebins[k] <= logages[i] && logages[i] < logagebins[k+1]) {
                        // increase the age weight
                        logageweights[j][k] += 1.0;
                    }
                }
            }
        }
    }

    // create distributions
    piecewise_linear_distribution<> logmassdist(logmassbins.begin(),logmassbins.end(),logmassweights.begin());
    vector<piecewise_linear_distribution<>> logagedists;
    for (size_t i = 0; i < nummassbin; ++i) {
        logagedists.push_back(piecewise_linear_distribution<>(logagebins.begin(),
        logagebins.end(),logageweights[i].begin()));
    }
    
    // vectors to store simulated mass
    vector<double> simlogmasses;
    vector<double> simlogages;

    // random number generator
    random_device rd;
    mt19937 gen(rd());


    // choose n stars (spcified size of simulated star table) from the distribution
    for (size_t i = 0; i < n; ++i){
        double logmass;
        double logage;
        for (size_t j = 0; j < nummassbin; ++j) {
            if (logmassbins[j] <= logmass && logmass < logmassbins[j+1]) {
                      // draw from distribution of relevent mass bin
                      logage = logagedists[j](gen);
                      logmass = logmassdist(gen);
            } 
        }
        if(RAND_DIST){
            simlogmasses.push_back(logmass);
            simlogages.push_back(logage);
        }
        else{
            simlogmasses.push_back(logmasses[i]);
            simlogages.push_back(logages[i]);
        }
    }
    //cout << simlogmasses[0] << endl;
    //cout << simlogages[0] << endl;

    // pack the two vectors as a table
    vector<vector<double>> simstartable;
    simstartable.push_back(simlogmasses);
    simstartable.push_back(simlogages);
    return simstartable;
}

/**
 * \brief plot startable, a 2-D vector of star values
 *
 * \param startable  vector of <logmass,logage>
 * \returns
 */
void plotstartable(vector<vector<double>> startable, size_t cluster, bool simulated)
{
    FILE * temp = fopen("star.temp", "w");
    FILE* gp=popen("gnuplot -persistent","w");
    for(size_t k=0;k<startable[0].size();k++) {
        fprintf(temp,"%f %f \n",startable[0][k],startable[1][k]);
    }
    fprintf(gp, "%s \n", "set terminal postscript eps enhanced color font 'Helvetica,10'");
    if (cluster == 1) {
        if (simulated == false) {
            fprintf(gp, "%s \n", "set output 'ONC.eps'");
            fprintf(gp, "%s%s %s \n", "set title \"","ONC","\"");
        } else {
            fprintf(gp, "%s \n", "set output 'simulatedONC.eps'");
            fprintf(gp, "%s%s %s \n", "set title \"","simulatedONC","\"");
        }
        fprintf(gp, "%s \n", "set yrange [-5:3]");
    } else {
        if (simulated == false) {
            fprintf(gp, "%s \n", "set output 'NGC2264.eps'");
            fprintf(gp, "%s%s %s \n", "set title \"","NGC2264","\"");
        } else {
            fprintf(gp, "%s \n", "set output 'simulatedNGC2264.eps'");
            fprintf(gp, "%s%s %s \n", "set title \"","simulatedNGC2264","\"");
        }
    }

    fprintf(gp, "%s \n", "set xlabel \"LogMass\"");
    fprintf(gp, "%s \n", "set ylabel \"LogAge\"");
    fprintf(gp, "%s \n", "plot 'star.temp'");

    fclose(temp);
}

/**
 * \brief plot observed ONC period distribution
 *
 * \returns
 */

void plotdistribution()
{
    // herbst (2002) data
    ifstream infile("herbst.txt");
    string line;
    // skip the first line
    getline(infile, line);
    double period, mass, disgard;
    FILE * temp1 = fopen("ONC_periods1.temp", "w");
	FILE * temp2 = fopen("ONC_periods2.temp", "w");
	FILE * gp3=popen("gnuplot -persistent","w");

	// write periods to file for gnuplot
    while (infile >> disgard >> period >> disgard >> disgard >> mass >> disgard) {
        //cout << period << endl;
        if (period > 0.001) {
            if (mass < 0.25) {
                fprintf(temp1,"%f \n",period);
            } else {
                fprintf(temp2,"%f \n",period);
            }
        }
    }
    // plot
    fprintf(gp3, "%s \n", "set terminal postscript eps enhanced color font 'Helvetica,10'");
    fprintf(gp3, "%s \n", "set output 'ONC_observed_distribution.eps'");
    fprintf(gp3, "%s\n", "binwidth=1");
    fprintf(gp3, "%s\n", "set boxwidth binwidth");
    fprintf(gp3, "%s\n", "bin(x,width)=width*floor(x/width) + binwidth/2.0");
    fprintf(gp3, "%s%s %s \n", "set multiplot layout 2,1 title \"","Period Distribution","\"");
    fprintf(gp3, "%s \n", "set ylabel \"Number of stars\""); 
    fprintf(gp3, "%s \n", "unset xlabel");
    fprintf(gp3, "%s \n", "set xrange [0:14]");
    fprintf(gp3, "%s%s %s \n", "set label 1\"","m < 0.25 solar mass","\" at graph 0.8,0.9"); 
    fprintf(gp3, "%s \n", "plot 'ONC_periods1.temp' using (bin($1,binwidth)):(1.0) smooth freq with boxes notitle");
    fprintf(gp3, "%s \n", "set xlabel");
    fprintf(gp3, "%s \n", "set xlabel \"Period (days)\"");
    fprintf(gp3, "%s%s %s \n", "set label 1\"","m > 0.25 solar mass","\" at graph 0.8,0.9");
    fprintf(gp3, "%s \n", "plot 'ONC_periods2.temp' using (bin($1,binwidth)):(1.0) smooth freq with boxes notitle");
    fprintf(gp3, "%s \n", "unset multiplot");

    fclose(temp1);
    fclose(temp2);
    pclose(gp3);

}

/**
 * \brief simulation of n stars
 *
 * \param cluster 1=ONC, 2=NGC and simulationNumber (to keep track if performing multiple iteration)
 * \returns
 */
void simulation(vector<vector<double>> simstartable, int simulationNumber)
{
    double timestep = TIMESTEPS[simulationNumber]; 
    //cout<< "TIMESTEP: " << timestep<<endl;
    
    // compute cmk table
    vector<vector<double>> cmktable = readcmk("cmkdata.txt"); //cmkdata.txt");
    //vector<vector<double>> cmktable;
    if(PropellerTests[simulationNumber]==true) {
        //cmktable = readKunitomo();
        cmktable = readcmk("cmkdata.txt");
        cout << "Kunitomo" << endl;
        }
    // random number generator
    default_random_engine genm;
    default_random_engine genb; 
    
    mt19937 genm_fixed(1); //Standard mersenne_twister_engine with fixed seed
    mt19937 genb_fixed(1); //Standard mersenne_twister_engine with fixed seed
    
	// normal distribution of log(massfactor)
    normal_distribution<double> logmdotfactordist(0.0,0.32);

	// normal distribution of bfieldstrength
    normal_distribution<double> bfielddist(1.67,0.1);

	// uniform distribution of bfieldstrength
	// uniform_real_distribution<double> bfielddist(1.2,2.0);

    clock_t start;
    double duration;
    start = clock();
    FILE * datafile;
    datafile = fopen("simulationONC.txt", "w");
    // write first line
    fprintf(datafile,"%s \n","mass(solarmass) age(Myr)        period(days)    log(mdotfactor) bfield(kG)");

    // keeping track of stars spinning faster than Keplerian speed
    FILE * datafile2 = fopen("fasterThanKeplerian.txt", "w");
    fprintf(datafile2,"%s \n","mass(solarmass) age(Myr)        period(days)    log(mdotfactor) bfield(kG)        Keplerian period");
    size_t numFasterThanKeplerian = 0;

	vector<starData> periods1_phase, periods2_phase;
	vector<double> periods1, periods2;
    vector<double> propellerStrengths1, propellerStrengths2;

	// iterate through distribution of stars
    for (size_t starNumber = 0; starNumber < simstartable[0].size(); ++starNumber) {
		// note: repeat variable for testing only
        for(size_t repeat = 0; repeat < 1; ++repeat){
		    double mass = pow(10,simstartable[0][starNumber]);
		    double age = pow(10,simstartable[1][starNumber]);
            double bfieldstrength = bfielddist(genb);
            double logmdotfactor = logmdotfactordist(genm);
            if (RAND_DIST == false){
                bfieldstrength = 1.67;
                logmdotfactor = 0.0;
            }
            TTauriStar star = TTauriStar(cmktable, mass, age, pow(10,logmdotfactor), bfieldstrength, timestep, PropellerTests[simulationNumber]);                        
            vector<double> dataVector = star.update();

            if (!dataVector.empty()) {
                double propellerStrength = star.propeller_strength;
                
                starData currentStar;
                currentStar.period = dataVector.at(0);
                currentStar.phase = dataVector.at(1);

                double period = currentStar.period;
                //if(period > 30){
                    //cout << "mass: " << mass << "age: " << age << "period: " << period << endl;
                //}
                // write to file
                fprintf(datafile,"%f        %f        %f        %f        %f \n",mass,age,period, logmdotfactor,bfieldstrength);
                    
                // Check if final period is greater than Keplerian period
                vector<double> kepler = star.getvector(12);
                double keplerianPeriod = kepler.back();
                
                if (keplerianPeriod > period) {
                    ++numFasterThanKeplerian;
                    fprintf(datafile2,"%f        %f        %f        %f        %f         %f \n",mass,age,period,logmdotfactor,bfieldstrength,keplerianPeriod);
                }
                if (mass < 0.25) {
                    periods1_phase.push_back(currentStar);
                    periods1.push_back(period);
                    propellerStrengths1.push_back(propellerStrength);
                    
                } else {
                    periods2_phase.push_back(currentStar);
                    periods2.push_back(period);
                    propellerStrengths2.push_back(propellerStrength);

                }
            }
            
		}
	}

    fclose(datafile);
    fclose(datafile2);
    duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
    cout<<"the program takes "<< duration << " s" << endl;
    cout << "there are " << numFasterThanKeplerian << " stars with periods shorter than break-up period" << endl;
    // plot
	plothistogram(periods1_phase, periods2_phase, simulationNumber);
    plotlogscalehistogram(periods1, periods2, simulationNumber);
    //plotPropellerStrengthhistogram(propellerStrengths1, propellerStrengths2, simulationNumber);
    //plotlogscalehistogramWithPhase(periods1_phase, periods1_phase);
}




/*The main program
 */
int main()
{
    // choose cluster simulation of single-star simulation
    size_t typeofsimulation;
    cout << "Type 1 for single-star simulation, type 2 for cluster simulation" << endl;
    cin >> typeofsimulation;
    cout << "The value you entered is " << typeofsimulation << endl;

    if (typeofsimulation == 1) {
        // mass
        double mass;
        cout << "mass: ";
        cin >> mass;
        cout << "The mass is " << mass << endl;

        // age
        double age; //Myrs
        cout << "age: ";
        cin >> age;
        cout << "The age is " << age << endl;

        // compute cmk table
        //vector<vector<double>> cmktable = readcmk("cmkdata.txt"); //cmkdata.txt");
        
        // if performing multiple simulations, iterate
        for(int simulationNumber = 0; simulationNumber < TotalSimulations; ++simulationNumber) {
            vector<vector<double>> cmktable;
            if(PropellerTests[simulationNumber]==true) {
                //cmktable = readKunitomo();
                cmktable = readcmk("cmkdata.txt");
            }
            else {
                cmktable = readcmk("cmkdata.txt"); //cmkdata.txt");
            }
            //for(int i = 0; i < 1; ++i){
                // Saving plots to folder
                stringstream folderNameStream0;
                folderNameStream0 << fixed << setprecision(2);
                //folderNameStream0 << "Star: Mass_" << mass << "_" << TIMESTEPS[simulationNumber];
                folderNameStream0 << "Star: Mass_" << mass << " Age_ " << age << "TF_" << PropellerTests[simulationNumber];
                string folderName0 = folderNameStream0.str();
                // make new folder and change directories into new folder
                mkdir(folderName0.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
                chdir(folderName0.c_str());
                // create a star with magnetic field of 1.67kG
                double timestep = .1*TIMESTEPS[simulationNumber]; //simulationNumber*OUTPUTID;
                TTauriStar star = TTauriStar(cmktable, mass, age, 1, 1.67, timestep, PropellerTests[simulationNumber]);
                star.update();
                // Plot all of the star's parameters as a function of time
                int TIME_INDEX = 1;
                int N_Y_INDEX = 13; // MT: changed from 10 to 12 so it plots perioddot and keplerian period
                for(int y_index = 2; y_index <= N_Y_INDEX; y_index += 1) {
                    star.plot(TIME_INDEX, y_index); 
                }
                // star.plot(2, 5);
                chdir(".."); 
            //}
        }


    } else if (typeofsimulation == 2) {
	        // Simuate a cluster: 1 = ONC, 2 = NGC 2264
	        size_t cluster;
	        cout << "Which cluster do you want to simulate?" << endl;
	        cout << "Type 1 for ONC, type 2 for NGC2264" << endl;
	        cin >> cluster;
	        cout << "The value you entered is " << cluster << endl;

            
	        // sample from a given cluster

	        vector<vector<double>> startable;
	        if (cluster == 1) {
	            startable = readcluster("hillenbrand.txt");
	        } else if (cluster == 2) {
	            startable = readcluster("dahm.txt");
	        } else {
	            throw invalid_argument( "Invalid input" );
	        }

          // plotstartable(startable, cluster, false);
          // simulate a startable of size n
          vector<vector<double>> simstartable = generatedistribution(startable,800);
          // plotstartable(simstartable, cluster, true);
          for(int simulationNumber = 0; simulationNumber < TotalSimulations; ++simulationNumber) {
                cout << simulationNumber << endl;
                simulation(simstartable, simulationNumber);
            }


    } else {
        throw invalid_argument( "Invalid input" );
    }

    // Plot the observed period distribution for ONC
    // plotdistribution();
    return 0;
}
