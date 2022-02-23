#pragma once
#include"Lattice.h"
#include <sstream>
#include<cmath>
namespace monteCarlo{
	class Exciton{
	private:	
		bool decayed;
		int xPos, yPos, zPos;
		int xhops,yhops,zhops;
		int index;
		int excitonStream;//set which run (stream) this exciton is in
		int startIndex;
		double averageLifeTime;
		double stationaryTime;
		double timeOfRecombination;
		int typeOfRecombination; /*0-radiatve,1-non-radiative (EEA), 2-non-radiative (ext),-1-Did not decay*/
		std::vector<int> ExcitonRoute;
		std::vector<double> LatticeEnergies;
		std::vector<double> dwellTimes;
		double exponentialRand(double,std::default_random_engine&);
		
		
	public:
		int getXpos();
		int getYpos();
		int getZpos();
		int getIndex();
		int getStream();
		double getLengthMoved();
		double getLifeTime();
		double getTimeOfRecombination();
		int getTypeOfRecombination();
		int getStartIndex();
		int getExcitonRouteLength();
		bool getDecayed();
		void setDecayed(bool);
		void resetExciton(monteCarlo::Lattice&,std::default_random_engine&);
		void EEA(monteCarlo::Lattice &,double);
		void advanceExciton(monteCarlo::Lattice&,double,std::default_random_engine&);
		double calculateDwellTime(double,std::default_random_engine&);
		Exciton(monteCarlo::Lattice&,double,int,bool,std::default_random_engine&);
		void construct(monteCarlo::Lattice&,double,int,bool,std::default_random_engine&);
		void Extraction(monteCarlo::Lattice&,double);
		void moveExciton(monteCarlo::Lattice&,std::default_random_engine&);
		
		
		/*printing and saving methods*/
		void saveExcitonRoute(std::string,monteCarlo::Lattice&,int);
		
		
};}