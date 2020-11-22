/**
 * \file ttauristar.cpp
 *
 * \authors Cynthia Yan ('18'), Shion Andrew
 *
 * \brief Implements the TTauriStar class.
 */
#include "ttauristar.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <vector>
#include <algorithm>
#include <random>

using namespace std;

double const TIMESTEP = 0.001;       /// timestep (Myrs) fraction of age
double const ALPHA = 0.01;           /// viscocity parameter
double const BETA = 1.35;            /// R_M/R_A
double const GAMMA = 1.0;            /// R_C/R_M
double const SHAPEFACTOR = 0.17;     /// f in the rotational inertia I=fMR^2
double const PROPEFF = 0.3;          /// propeller coefficient
double const CRITICALDENSITYCOEFF = 1.2e22;  /// Adjustable parameter "A" (T^1/2*cm^-1/2) used to determine critical density
double const PROPSTARTTIME = 0.05;   /// B-field turn-on time (Mrs); simulation starts at this time
//double const TURNONTIME = 0.049;      /// should be the same as above
//double const PROPTIMESPREAD = 0.002; /// introduce randomness for PROPSTARTTIME
double const DELTAM = 0.001;            /// deltam to determine when to stop mass iterations
double const ACCPOWER = -2.1;		// new parameter to investigate accretion efficiency
double const DEFAULTAGE = 1; 		/// default age of 1 million (for when user inputs age less than or equal to zero)
double const MAXITERATIONS = 50; // maximum number of iterations
double const GRAMS_TO_SOLAR_MASS = 5.02785*pow(10,-34); // conversion factor from grams to solar mass
double const SECONDS_TO_MYR = 3.17098*pow(10,-14); // conversion factor from seconds to Myr
double const CM_TO_SOLAR_RADIUS = 1.437*pow(10,-11); // conversion factor from cm to solar radius
double const G = 2937.5; // G in units of (solar radius^3)/(day^2 solarmass)

TTauriStar::TTauriStar(vector<vector<double>> cmktable,
	double mass, double age, double massdotfactor, double bfieldstrength, double timestep, bool Romanova)
	// initilizing private member variables
    :cmktable_(cmktable), mass_(mass), mass0_(mass), mass2_(0),
     age_(age), massdotfactor_(massdotfactor), bfieldstrength_(bfieldstrength),Romanova_(Romanova)
	{
	// set age to default age if input age is invalid
	if (age_ <= 0) {
			age_ = DEFAULTAGE;
	}
	// initialize validity of mass
	if (mass > 3) {
		valid_ = false;
	} else {
		valid_ = true;
	}
	timestep_ = timestep;

	//random_device rd1;  //Will be used to obtain a seed for the random number engine
	//mt19937 genp(rd1()); //Standard mersenne_twister_engine seeded with rd(
	//uniform_real_distribution<double> proptimeDist(0.048,0.052);
	//propstarttime_ = proptimeDist(genp);

	propstarttime_ = PROPSTARTTIME;

	// initially set propeller endtime to be the same as starttime (no phase 2)
	propendtime_ = propstarttime_;

  // assuming accretion is on at start
	acceff_ = 0;

	// save initial values of age and acceff in the corresponding vectors
	// go backwards in time
	while (age_ > 0) {

		// populate ages vector
		ages_.push_back(age_);

		// populate accretion efficiency vector
		acceffs_.push_back(acceff_);

		// calculate new age
		age_ -= timestep_;
	}

	// reverse the two vectors to be in forward time order
	reverse(ages_.begin(),ages_.end());
	reverse(acceffs_.begin(),acceffs_.end());
}


double TTauriStar::calculatemassdot()
{
	// calculate current value of Mdot in unts M_sun/yr
	// massdotfactor mimics scatter
		return 7.0e-8*massdotfactor_*pow(age_,ACCPOWER)*pow(mass_,2.4);
}

double TTauriStar::calculateradius()
{

	size_t index1 = 0; // natural number data type
	size_t index3 = 0;

	// find masslower and massupper such that masslower <= mass < massupper
	double masslower = cmktable_[0][0]; // [column][row], column = 0 is mass, 1 is age, 2 is radius
	double massupper = cmktable_[0][0];
	while (cmktable_[0][index3] <= mass_) {
	    if (cmktable_[0][index3] > masslower) {
        	index1 = index3;
        	masslower = cmktable_[0][index3];
            }
	    ++index3;
	}
	// At this point index1 is the first row # of the masslower block.
	// and index3 is the first row # of the massupper block.
	massupper = cmktable_[0][index3];

	// find agelower1 and ageupper1 such that agelower1 <= age < ageupper1
	double agelower1 = cmktable_[1][index1];
	double ageupper1 = cmktable_[1][index1];
	size_t index2 = index1;

	while (cmktable_[1][index2] <= age_) {
		++index2;
	}
	index1 = index2 - 1;
	agelower1 = cmktable_[1][index1];
	ageupper1 = cmktable_[1][index2];

	// find agelower2 and ageupper2 such that agelower2 <= age < ageupper2
	double agelower2 = cmktable_[1][index3];
	double ageupper2 = cmktable_[1][index3];
	size_t index4 = index3;
	while (cmktable_[1][index4] <= age_) {
		++index4;
	}
	index3 = index4 - 1;
	agelower2 = cmktable_[1][index3];
	ageupper2 = cmktable_[1][index4];

	// four coefficients
	double coeff1 = (massupper-mass_)*(ageupper1-age_)/(massupper-masslower)/(ageupper1-agelower1);
	double coeff2 = (massupper-mass_)*(age_-agelower1)/(massupper-masslower)/(ageupper1-agelower1);
	double coeff3 = (mass_-masslower)*(ageupper2-age_)/(massupper-masslower)/(ageupper2-agelower2);
	double coeff4 = (mass_-masslower)*(age_-agelower2)/(massupper-masslower)/(ageupper2-agelower2);

	// return the interpolated radius
	return coeff1*cmktable_[2][index1]+coeff2*cmktable_[2][index2]+coeff3*cmktable_[2][index3]+coeff4*cmktable_[2][index4];
}

double TTauriStar::calculatebfield()
{
	// turn on a constant dipole magnetic field at some input time

	if (age_ >= propstarttime_) {
		return bfieldstrength_;
	} else {
		return 0;
	}
}

double TTauriStar::calculaterm()
{
	//return .0396574*BETA*pow(massdot_,-2./7.)*pow(mass_,-1./7.)*pow(bfield_,4./7.)*pow(radius_,12./7.);
	return 14.4*BETA*pow(massdot_/(1.e-8),-2./7.)*pow(mass_/0.5,-1./7.)*pow(bfield_,4./7.)*pow(radius_/2.,12./7.);
	// radius at which ram pressure equals magnetic pressure*BETA (magnetispheric radius)
}



double TTauriStar::calculatediskdensity()
{
	// density of accretion disk at rm
	return 5.2*pow(ALPHA,-4.0/5.0)*pow(6.34*1e9*massdot_,0.7)*pow(mass_,0.25)*pow(rm_*6.957,-0.75)*pow(1.0-pow(radius_/rm_,0.5),7.0/10.0);

}

double TTauriStar::calculatecriticaldensity() {
	// calculates and returns current value of disk density in units of g/cm^2

	// calculate temperature of star at magentospheric radius in units of K
	double temprm = 1400*pow(.01/ALPHA, 0.2)*pow(massdot_/1e-8,.3)*pow(mass_/.5,.25)*pow(40/rm_,.75);
	// convert units of R_m into cm (1 Rdot = 6.957e+8m)
	return CRITICALDENSITYCOEFF/(pow(temprm,.5)*pow(rm_*6.957*1e10,3.0/2.0));
}

void TTauriStar::calculatemasses()
{
	// clear vectors involved
	masses_.clear();
	// initial values
	masses_.push_back(mass0_);
	// go backwards in time
	for (size_t i = ages_.size() - 1; i >= 1; --i) {
		// retrieve age and acceff (accreting or not acccreting)
		age_ = ages_[i];
		acceff_ = acceffs_[i];
		// calculate mass accretion rate
		massdot_ = calculatemassdot();
		// calculate new mass
		mass_ -= 1.0e6*massdot_*acceff_*timestep_;
		// push_back the mass into the masses vector
		// will reversed after all the push_backs are done for efficiency
		masses_.push_back(mass_);
	}
	// reverse the two vectors to be in forward time order
	reverse(masses_.begin(),masses_.end());
}

void TTauriStar::calculateperiods()
{
	// MT -- TEMPORARY -- open file:
	FILE * temp = fopen("period_change_RStar.temp","w");
	vector<double> periodChange;
	
	// clear vectors involved
	periods_Romanova_.clear();
	periods_.clear();
	massdots_.clear();
	radii_.clear();
	rms_.clear();
	diskdensities_.clear();
	acceffs_.clear();
	phase_.clear();
	rmPeriods_.clear();
	perioddots_.clear();
	keplerianPeriods_.clear();

	// initialze mass2_
	mass2_ = masses_[0];
	// initialize radius_ !!!!!
	mass_ = masses_[0];

	age_ = ages_[0];

	radius_ = calculateradius();



	massdot_ = calculatemassdot();
	// initial phase
	size_t phase = 1;
	// go forward in time
	for (size_t i = 0; i < ages_.size(); ++i) {
		// retrieve mass and age stored in vectors
		mass_ = masses_[i];
		// cannot handle mass > 3 solar mass
		if (mass_ > 3) {
			valid_ = false;
			break;
		}
		age_ = ages_[i];

		// update data members
		massdot_ = calculatemassdot();
		double radius = radius_;	
		radius_ = calculateradius();
		bfield_ = calculatebfield();
		rm_ = calculaterm();
		diskdensity_ = calculatediskdensity();
		radiusdot_ = (radius_-radius)/timestep_;

		// calculate the Keplerian period at rm
		periodrm_ = 0.1159*pow(rm_,3./2.)*pow(mass_,-1./2.);

		double criticaldensity = calculatecriticaldensity();


		// Calculate the periods
		// Phase 1: spin at break-up period.  SHOULD ALLOW FOR MORE THAT ONE POINT!!!!
		if (age_ < propstarttime_ && phase <= 1) {
	    	period_ = 0.1159*pow(radius_,3./2.)*pow(mass_,-1./2.);
	    	// assume accretion at start
	    	acceff_ = 1;
			period_Romanova_ = period_;
		// Phase 2: spin down due to propeller effect
		} else if (period_ < periodrm_ && phase <= 2 && diskdensity_ > criticaldensity) {
			// doesn't accrete
			//if (acceff_ = 1){
			//	cout << (rm_/radius_)/60 << endl;
			//	cout << periodrm_/period_ << endl;
			//	double mu = (rm_/radius_)/60;
			//	double omega = periodrm_/period_;
			//	cout << 3.5*pow(10,34)*pow(mu,-2)*pow(bfield_,2)*pow(radius_/2,3)*GRAMS_TO_SOLAR_MASS*pow(CM_TO_SOLAR_RADIUS,2)*86400/SECONDS_TO_MYR << endl;
			//}
			acceff_ = 0;
			//double PdotTerm1 = 2*radius_*radiusdot_*mass_*period_/(mass_*pow(radius_,2));
			double mu = (rm_/radius_)/60;
			double omega = periodrm_/period_;
			//acceff_ = 1-(0.57*pow(omega,.3));
			double PdotTerm1 = ((massdot_*acceff_*pow(radius_,2)+2*radius_*radiusdot_*mass_)*period_)/(mass_*pow(radius_,2));
			double R0 = 1.4*pow(10,11)*CM_TO_SOLAR_RADIUS;
			double romanovaL = pow(mu,-2)*pow(bfield_,2)*pow(R0/CM_TO_SOLAR_RADIUS,3)*pow(radius_/R0,6)*GRAMS_TO_SOLAR_MASS*pow(CM_TO_SOLAR_RADIUS,2)*86400/SECONDS_TO_MYR;
			//double romanovaL = 7.61*pow(10,35)*pow(mu,-2)*GRAMS_TO_SOLAR_MASS*pow(CM_TO_SOLAR_RADIUS,2)*86400/SECONDS_TO_MYR; // units of solar mass * (solar radius)^2/(days*Myr)
			//double romanovaL = 3.5*pow(10,34)*pow(mu,-2)*pow(bfield_,2)*pow(radius_/2,3)*GRAMS_TO_SOLAR_MASS*pow(CM_TO_SOLAR_RADIUS,2)*86400/SECONDS_TO_MYR;
			double romanovaTorque = 0.51*pow(omega,.99)*romanovaL;
			if (Romanova_) {
				period_ += timestep_*PdotTerm1 + timestep_*romanovaTorque*pow(period_,2)/(2*3.14*SHAPEFACTOR*mass_*pow(radius_,2));
			}
			else{
				period_ += timestep_*PdotTerm1 + timestep_*PROPEFF*0.972*pow(BETA,-3.)*pow(period_,2.)*pow(mass_,-4./7.)*pow(bfield_,2./7.)*pow(radius_,-8./7.)*pow(massdot_/1.e-8,6./7.);
			}
			//cout << "Romanava Torque:" << romanavaTorque*pow(period_,2)/(2*3.14*SHAPEFACTOR*mass_*pow(radius_,2)) << endl;
			//period_Romanava_ += timestep_*romanavaTorque*pow(period_Romanava_,2)/(2*3.14*SHAPEFACTOR*mass_*pow(radius_,2));
			//double ourTorque = 2*3.14*1.44385*pow(BETA,-3.)*pow(mass_,3./7.)*pow(bfield_,2./7.)*pow(radius_,6./7.)*pow(massdot_/1.e-8,6./7.);
			//if (romanavaTorque < ourTorque) {
			//	cout << romanavaTorque << endl;
			//cout << "Our Torque:" << ourTorque << endl;
			// keep track of the propeller endtime
			
			propendtime_ = age_;			
			phase = 2;
			
		// Phase 3: disk-locked
		} else if (diskdensity_ > criticaldensity && phase <= 3) {
			//if (phase == 2) {
			//	cout << age_ << endl;
			//}
			//cout << criticaldensity << endl;
			// condition for turning on propellar effect is w(Rc) = Req = gamma Rcm.
			period_ = 8.*pow(GAMMA*BETA/0.9288,3./2.)*pow(massdot_/1.0e-8,-3./7.)*pow(mass_/0.5,-5./7.)*pow(radius_/2.,18./7.)*pow(bfield_,6./7.);
			phase = 3;
			acceff_ = 1;
			period_Romanova_ = period_;

		// Phase 4: unlocked
		} else {
			acceff_ = 0;

			// version where accretion happens at magnetospheric radius
			// double periodChange = (TIMESTEP*acceff_*period_*massdot_/mass_
			// 	-50.74*TIMESTEP*acceff_*pow(period_,2)*massdot_*pow(mass_,-0.5)*pow(rm_,0.5)*pow(radius_,-2.))/TIMESTEP;

			// period_ += period_*2*(radius_-radius)/radius_
			//         +TIMESTEP*acceff_*period_*massdot_/mass_
			// 	-50.74*TIMESTEP*acceff_*pow(period_,2)*massdot_*pow(mass_,-0.5)*pow(rm_,0.5)*pow(radius_,-2.);



			// version where we assume accretion happens at the surface of the star
			//    double periodChange = (TIMESTEP*acceff_*period_*massdot_/mass_
			//    	-50.74*TIMESTEP*acceff_*pow(period_,2)*massdot_*pow(mass_,-0.5)*pow(radius_,-1.5))/TIMESTEP;

			   period_ += period_*2*(radius_-radius)/radius_  
			   		+timestep_*acceff_*period_*massdot_/mass_
			 	-50.74*timestep_*acceff_*pow(period_,2)*massdot_*pow(mass_,-0.5)*pow(radius_,-1.5)/(2*3.1415*SHAPEFACTOR);
				double breakupPeriod = 0.1159*pow(radius_,3./2.)*pow(mass_,-1./2.);
				if(period_ < breakupPeriod){
					period_ = breakupPeriod;
				}

				period_Romanova_ = period_;
			// fprintf(temp, "%f %f \n", age_, periodChange);

			// version where we neglect final term
			// period_ += period_*2*(radius_-radius)/radius_
			// 		   +TIMESTEP*acceff_*period_*massdot_/mass_;

		   	// phase = 4;
			
			//period_ += period_*2*(radius_-radius)/radius_
			//        +TIMESTEP*1e6*acceff_*period_*massdot_/mass_
			//	-50.74*TIMESTEP*1e6*acceff_*pow(period_,2)*massdot_*pow(mass_,-0.5)*pow(rm_,0.5)*pow(radius_,-2.);
			 phase = 4;
		}

		// calculate mass moving forward
		if (i < ages_.size() - 1) {
			mass2_ += 1.0e6*massdot_*acceff_*timestep_;
		}
		
		// finding and storing perioddot
		if (periods_.empty()) {
			perioddots_.push_back(0); // placeholder
		} else {
			perioddots_.push_back((period_ - periods_.back())/timestep_);
		}
		
		// store the period
		periods_.push_back(period_);
		periods_Romanova_.push_back(period_Romanova_);
		propellerStrength_ = periodrm_/period_;
		propellerStrengths_.push_back(propellerStrength_);
		massdots_.push_back(massdot_);
		radii_.push_back(radius_);
		rms_.push_back(rm_);
		diskdensities_.push_back(diskdensity_);
		acceffs_.push_back(acceff_);
		phase_.push_back(phase);
		rmPeriods_.push_back(periodrm_);

		keplerianPeriods_.push_back(0.1159*pow(radius_,3./2.)*pow(mass_,-1./2.));

	}
	fclose(temp);
}

vector<double> TTauriStar::update()
{
	vector<double> dataVector; // holds data to return 
	if (valid_) {
		// keep track of the number of iterations
		int numIterations = 0;
		while (abs((mass2_-mass0_)/mass0_) > DELTAM && numIterations < MAXITERATIONS) {
			calculatebfield();
			calculatediskdensity();
			calculatemassdot();
			calculateradius();
			calculatemasses();
		  calculateperiods();
			++numIterations;
		}

		if(!isnan(period_)) {
		dataVector.push_back(period_);
		dataVector.push_back(1.0*phase_.back());
		return dataVector;
		}
	}
	else {
		cout << "Age: " << age_ << endl;
		cout << "Star is not valid" << endl;
		dataVector.push_back(0);
		dataVector.push_back(0);
	}
	return dataVector;
}

vector<double> TTauriStar::getvector(int n)
{
	vector<double> output;
	if (n == 1) {
		output = ages_;
	}
	if (n == 2) {
		output = masses_;
	}
	if (n == 3) {
		output = periods_;
	}
	if (n == 4) {
		output = massdots_;
	}
	if (n == 5) {
		output = radii_;
	}
	if (n == 6) {
		output = rms_;
	}
	if (n == 7) {
		output = diskdensities_;
	}
	if (n == 8) {
		output = acceffs_;
	}
	if (n == 9) {
		output = phase_;
	}

	if (n == 10) {
		output = rmPeriods_;
	}

	if (n == 11) {
		output = perioddots_;
	}

	if (n == 12) {
		output = keplerianPeriods_;
	}
	if (n == 13) {
		output = propellerStrengths_;
	}
	if (n==14) {
		output = periods_Romanova_;
	}

	return output;
}

string TTauriStar::getname(int n)
{
	string output;
	if (n == 1) {
		output = "Age";
	}
	if (n == 2) {
		output = "Mass";
	}
	if (n == 3) {
		output = "Period";
	}
	if (n == 4) {
		output = "AccretionRate";
	}
	if (n == 5) {
		output = "Radius";
	}
	if (n == 6) {
		output = "R_M";
	}
	if (n == 7) {
		output = "DiskDensity";
	}
	if (n == 8) {
		output = "Acceff";
	}
	if (n== 9) {
		output = "Phase";
	}
	if (n==10) {
		output = "R_mPeriods";
	}
	if (n==11) {
		output = "Period_Dot";
	}
	if (n==12) {
		output = "Kepler_Period";
	}
	if (n==13) {
		output = "Propeller Strength";
	}
	if (n==14) {
		output = "Period_PropAccretion_";
	}
	return output;
}

string TTauriStar::getunit(int n)
{
	string output;
	if (n == 1) {
		output = "Myr";
	}
	if (n == 2) {
		output = "solar mass";
	}
	if (n == 3) {
		output = "day";
	}
	if (n == 4) {
		output = "solar mass/year";
	}
	if (n == 5) {
		output = "solar radius";
	}
	if (n == 6) {
		output = "solar radius";
	}
	if (n == 7) {
		output = "g/cm^2";
	}
	if (n == 8) {
		output = "";
	}
	if (n == 9) {
		output = "1 = spin up, 2 = propeller effect, 3 = disk locked, 4 = disk unlocked";
	}

	if (n==10) {
		output = "days";
	}

	if (n==11) {
		output = "days/Myr";
	}

	if (n==12) {
		output = "days";
	}
	if (n==13) {
		output = "v_s_t_a_r / v_d_i_s_k";
	}
	if (n==14) {
		output = "days";
	}
	return output;
}
/**
 * \brief plot parameters of a single star
 *
 * \inputs: m - x axis
 * 			n - y axis
 */
void TTauriStar::plot(int m, int n)

{
  vector<double> xVector = getvector(m); //gets row and columns
  vector<double> yVector = getvector(n);
	string title = getname(n); //creating title of file

	ofstream myfile; //creating file for data
	myfile.open(title);

	//stringstream plotNameStream;
	//plotNameStream << fixed << setprecision(2);
	//plotNameStream << title;
	//string plotName = plotNameStream.str();

	// construct log-log plot for accretion rate
	if(getname(n) == "AccretionRate1") {
		//generate a log-log plot
		vector<double> yVectorLog = yVector;
		//vector<double> xVectorLog = xVector;

		for(size_t k=0; k<yVector.size(); k++) {
			yVectorLog[k] = log10(yVector[k]);
			myfile<<to_string(xVector[k])<< " " << to_string(yVectorLog[k]) << endl; //save table data to myfile
		}
		//plotNameStream << " SLOPE_" << to_string(slope(xVector, yVectorLog));
		//plotName = plotNameStream.str();
		myfile.close();

		FILE* gp=popen("gnuplot -persistent","w"); //popen opens a plotting window

		fprintf(gp, "%s \n", "set terminal postscript eps enhanced color font 'Helvetica,20'");
		fprintf(gp, "%s'%s.%s' \n", "set output", title.c_str(), "eps");
		// To create a screen output remove previous two lines, but then this image cannot be saved.
		fprintf(gp, "%s%s %s %s%s\n", "set title \"",getname(n).data(),"vs",getname(m).data(),"\"");

		fprintf(gp, "%s%s %s%s%s\n", "set xlabel \"",getname(m).data(),"(",getunit(m).data(),")\"");
		fprintf(gp, "%s%s %s%s%s\n", "set ylabel \"Ln of ",getname(n).data(),"(",getunit(n).data(),")\"");
		fprintf(gp, "plot '%s'\n", title.c_str());

	}
	
	// construct overlaid plot for romanavaPeriods
	else if(n == 14) {
		for(size_t k=0;k<xVector.size();k++) {
			//myfile<< to_string(yVector[k]) << endl; //save table data to myfile
			myfile<< to_string(xVector[k])<<" " << to_string(yVector[k]) << endl; //save table data to myfile
		}
		myfile.close();

		string title2 = getname(3); //overlay regular period plot
		// FILE * temp = fopen("data.temp", "w"); // * is for a pointer; "w" means write, could also be "r" for read.
		FILE* gp=popen("gnuplot -persistent","w"); //popen opens a plotting window

		// for(size_t k=0;k<xVector.size();k++) {
		//         fprintf(temp,"%f %f \n",xVector[k],yVector[k]);
		// }


		fprintf(gp, "%s \n", "set terminal postscript eps enhanced color font 'Helvetica,20'");
		fprintf(gp, "%s'%s.%s' \n", "set output", title.c_str(), "eps");
		// To create a screen output remove previous two lines, but then this image cannot be saved.
		fprintf(gp, "%s%s %s %s%s\n", "set title \"",getname(n).data(),"vs",getname(m).data(),"\"");

		fprintf(gp, "%s%s %s%s%s\n", "set xlabel \"",getname(m).data(),"(",getunit(m).data(),")\"");
		fprintf(gp, "%s%s %s%s%s\n", "set ylabel \"",getname(n).data(),"(",getunit(n).data(),")\"");
		fprintf(gp, "plot '%s', '%s'\n", title.c_str(), title2.c_str());
	}
	
	
	
	
	else {
		for(size_t k=0;k<xVector.size();k++) {
			//myfile<< to_string(yVector[k]) << endl; //save table data to myfile
			myfile<< to_string(xVector[k])<<" " << to_string(yVector[k]) << endl; //save table data to myfile
		}
		myfile.close();

		// FILE * temp = fopen("data.temp", "w"); // * is for a pointer; "w" means write, could also be "r" for read.
		FILE* gp=popen("gnuplot -persistent","w"); //popen opens a plotting window

		// for(size_t k=0;k<xVector.size();k++) {
		//         fprintf(temp,"%f %f \n",xVector[k],yVector[k]);
		// }


		fprintf(gp, "%s \n", "set terminal postscript eps enhanced color font 'Helvetica,20'");
		//fprintf(gp, "%s \n", "set term png");
		fprintf(gp, "%s'%s.%s' \n", "set output", title.c_str(), "eps");
		// To create a screen output remove previous two lines, but then this image cannot be saved.
		fprintf(gp, "%s%s %s %s%s\n", "set title \"",getname(n).data(),"vs",getname(m).data(),"\"");

		fprintf(gp, "%s%s %s%s%s\n", "set xlabel \"",getname(m).data(),"(",getunit(m).data(),")\"");
		fprintf(gp, "%s%s %s%s%s\n", "set ylabel \"",getname(n).data(),"(",getunit(n).data(),")\"");
		fprintf(gp, "plot '%s'\n", title.c_str());
	}

}

double TTauriStar::slope( vector<double>& x, vector<double>& y){
    // if(x.size() != y.size()){
    //     throw exception("...");
    // }

    double n = x.size();
    double avgX = accumulate(x.begin(), x.end(), 0.0) / n;
    double avgY = accumulate(y.begin(), y.end(), 0.0) / n;

    double numerator = 0.0;
    double denominator = 0.0;

    for(int i=0; i<n; ++i){
        numerator += (x[i] - avgX) * (y[i] - avgY);
        denominator += (x[i] - avgX) * (x[i] - avgX);
    }

    if(denominator == 0){
        return 0;
    }

    return numerator / denominator;
}
