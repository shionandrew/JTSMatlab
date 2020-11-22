/**
 * \file ttauristar.hpp
 *
 * \authors Cynthia Yan (18'), Shion Andrew
 *
 * \brief Declares the TTauriStar class.
 */

#ifndef TTAURISTAR_HPP_INCLUDED
#define TTAURISTAR_HPP_INCLUDED

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
#include <tuple>
#include <string>
using namespace std;

/**
 * \class TTauriStar
 *
 * \brief Represents a T Tauri Star throughout its evolution.
 *
 * \details The first parameter is a pointer matrix to the mass-age-radius table,
 *          mass and age are the final values for a particular star,
 *          massdotfactor is the randomization Mdot parameter,
 *          bfieldstrength is the randomization magnetic field strength.
 */
class TTauriStar {
public:

	/**
	 * \brief Parameterized constructor for T Tauri Star. Initializes input parameters, as well as ages_ and
	 *  acceffs_ based on assumption that there is no accretion until propstarttime (no propellar phase)
	 *
	 * \param cmktable: a
	 * \param mass: final mass of star
	 * \param age: final age of star
	 * \param massdotfactor:
	 * \param bfieldstrength:
	 */
	TTauriStar(vector<vector<double>> cmktable, double mass, double age,
		double massdotfactor, double bfieldstrength, double timestep, bool Romanava);

	/**
	 * \brief update runs the simulation for each star until convergence.
	 *
	 * \returns vector with period and phase
	 */

	vector<double> update();

	/**
	 * \brief getvector fetches various evolutionary vectors for each star,
	 *
	 * \param integer code of the vector you want
	 *
	 * \returns the vector you want
	 */
	vector<double> getvector(int n);

	/**
	 * \brief fetch names of these vectors for use in plotting
	 *
	 * \param integer code of the vector you want
	 *
	 * \returns a string telling the name
	 */
	string getname(int n);

	/**
	 * \brief fetch units of these vectors for use in plotting
	 *
	 * \param integer code of the vector you want
	 *
	 * \returns a string telling the unit
	 */
	string getunit(int n);

	/**
	 * \brief makes plots through gnuplot vector m vs. vector n is plotted.
	 *
	 * \param integer codes of the vectors you want
	 *
	 * \returns
	 */
	void plot(int m, int n);


	/**
	 * \brief finds slope of two vectors, assuming linear relationships
	 *
	 * \param integer codes of the vectors you want
	 *
	 * \returns slope
	 */
	double slope( vector<double>& x, vector<double>& y);

	// public data members
	double propellerStrength_;							/// ratio of star velocity to keplerian velocity of disk
	vector<double> propellerStrengths_; /// ratio of star velocity to keplerian velocity of disk


private:
	double calculatecriticaldensity();
	double calculatemassdot();
	double calculateradius();
	double calculatebfield();
	double calculaterm();
	double calculatediskdensity();
	void calculatemasses();
	void calculateperiods();


	vector<vector<double>> cmktable_; // TTauriStar data members (Matrix)

	// data members at a given step in time
	double mass_;											/// (solar masses)
	double mass0_;            			  /// (solar mass) mass that we are trying to converge to, equal to the original input mass
	double mass2_;             				/// (solar mass) mass calculated going forward to test for convergence against mass0
	double age_;											/// Myears
	double massdotfactor_;      			/// equal to the original input parameter
	double bfieldstrength_;						/// (KGauss)
	double massdot_;            			/// M_sun/yr
	double period_;										/// days
	double period_Romanova_; 						/// TEMPORARY! Stores period calculated by results of Romanava et al with accretion during propeller effect 
	double radius_;										/// R_sun
	double bfield_;										/// KGauss
	double rm_;												/// R_sun
	double periodrm_;									/// days
	double diskdensity_;        			/// surface density at R_M, g/cm^2
	double propstarttime_;						/// start time for the propellar effect
	double propendtime_;        			/// end time for the propeller effect
	double acceff_;            				/// fraction of deposition by accretion
	bool valid_;  /// if the mass is > 3, there are not values in the table, so we drop stars with those masses.
	//This should be changed in the future with better data tables.
	double timestep_; 		// timestep (in Myrs)
	bool Romanova_; 	// TEMPORARY boolean determining whether to allow for mass accretion during prop. stage
	double radiusdot_; //R_sun/Myr


	vector<double> ages_;
	vector<double> masses_;         /// mass of the protostar changed backward
	vector<double> massdots_;       /// mass accretion rate of the protostar
	vector<double> periods_;        /// period of the protostar
	vector<double> periods_Romanova_;  /// TEMPORARY period of the protostar, allowing for accretion during prop.
	vector<double> acceffs_;
	vector<double> radii_;
	vector<double> rms_;
	vector<double> diskdensities_;
	vector<double> phase_;					// vector storing "phases" the star. Solely for the purpose of understanding the stages of the star's evolution in our graph
	vector<double> rmPeriods_;			// stores keplerian period calculated at magnetospheric radius (days)
	vector<double> perioddots_;
	vector<double> keplerianPeriods_;
};

#endif
