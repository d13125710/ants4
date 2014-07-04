#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>


#include <array> //for std::array
#include "Coords.h"
#include "Tokenizer.h"
#include "TokenFinder.h"


using namespace std;


typedef vector<Coords> Coords_vector;
typedef Tokenizer<TokenFinder> Tk;

typedef std::vector<std::vector<double> > MatrixArrayType;


class FileReader
{
private: 
	std::ifstream m_inFile; 
	std::string m_fileName;
	string sInput;
	Coords_vector m_strVector;
   // MatrixArrayType *m_theMatrix;
   
	std::default_random_engine m_generator;

public:
	bool Read(int setStartCity=0);
	MatrixArrayType *m_theMatrix;
	MatrixArrayType *m_nearestneighborMatrix;
	double ShowSolution();
	double ReadMatrix(int x, int y) const
	 {
		 double d=(*m_theMatrix) [x][y];
		 return  d;
	 }
	int getsize() const
	{
	     return  (int)(*m_theMatrix).size();
	}




	FileReader(string &str);
	~FileReader(void);

	
	 void swap2(double *v, double *v2, int i, int j)
    /*
     * FUNCTION: auxiliary routine for sorting an integer array
     * INPUT: two arraya, two indices
     * OUTPUT: none
     * (SIDE)EFFECTS: elements at position i and j of the two arrays are swapped
     */
    {
		double tmp;

		tmp = v[i];
		v[i] = v[j];
		v[j] = tmp;
		tmp = v2[i];
		v2[i] = v2[j];
		v2[j] = tmp;
    }

    void sort2(double *v, double *v2, int left, int right)
    /*
     * FUNCTION: recursive routine (quicksort) for sorting one array; second
     * arrays does the same sequence of swaps
     * INPUT: two arrays, two indices
     * OUTPUT: none
     * (SIDE)EFFECTS: elements at position i and j of the two arrays are swapped
     */
    {
		int k, last;

		if (left >= right)
			return;
		swap2(v, v2, left, (left + right) / 2);
		last = left;
		for (k = left + 1; k <= right; k++)
			if (v[k] < v[left])
			swap2(v, v2, ++last, k);
		swap2(v, v2, left, last);
		sort2(v, v2, left, last);
		sort2(v, v2, last + 1, right);
    }

  int getGaussianDistributionRowIndex(std::uniform_int_distribution<int> &distribution);
  MatrixArrayType& calculateNearestNeigbhor(int NUMOFANTS);
	
  
private:
	//open
	bool Openfile();
	double calculateRoundDistance(double x1, double y1,double x2, double y2) ;
	double calculateEuclidianDistance(double x1, double y1,double x2, double y2); 
	void matrix_calculate(MatrixArrayType& theMatrix);
	
	

};

