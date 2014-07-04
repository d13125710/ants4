
#include "StdAfx.h"
#include "FileReader.h"


std::string::size_type sz;     // alias of size_t
 

FileReader::FileReader(string &str)
{
	m_fileName = str;
	m_inFile.open(str + ".tsp");

}
double FileReader::ShowSolution()
{
	std::ifstream m_solutionFile; 
	m_solutionFile.open(m_fileName + ".opt.tour");
	bool ReadAhead = false;
	std::vector<int>  strVector;
	double distance = 0;

	if(m_solutionFile.is_open())
	{
		m_solutionFile.setf(ios::skipws);

			while(!m_solutionFile.eof() && m_solutionFile.good())
			{
				//we need to read up to the newline
				getline(m_solutionFile , sInput , '\n');
						
				if(ReadAhead)
				{
					//split up the input into tokens of whit spaces
					Tk ti(sInput.c_str() , TokenFinder(" ")); // find white spaces
						
					while( ti != Tk()) // skip the first as format in "place" , distance , distance , dont need the first
					{
						//convert the string to a double
						int position = std::stoi(*ti,&sz);
						if(position == -1)
						{
							ReadAhead = false;
							break;
						}
						strVector.push_back(position - 1); //our arrays start at 0
					
				
						++ti;
					}


					//m_strVector.push_back(sInput);
				}
				//start of the COORDS in the TSP format , so read the coords from here....
				if(sInput == "TOUR_SECTION")
				{
					ReadAhead = true;
				}
		
			}

				m_solutionFile.close();	

			// strVector 1 14 13 12 7 6 15 5 11 9 10 16 3 2 4 8
	
			//int lastNode = strVector[0];
			//int next = strVector[1];
			
			int i, start, end;

			for(unsigned int i = 1; i < strVector.size() ; i++)
			{
				
				start = strVector[i-1];
				end = strVector[i];
				distance += ReadMatrix(start , end);
				
			}

	
	distance += ReadMatrix(end, strVector[0]);
	
	}
	printf("Shortest tour found: distance =%.2f \n" , distance);
				for (int j=0; j<strVector.size(); j++) 
				printf("%i|",strVector[j]+1);


	return distance;
}

bool FileReader::Read(int setStartCity)
{
	bool status = true;
	bool ReadAhead = false;

	if(m_inFile.is_open())
	{
		m_inFile.setf(ios::skipws);

		while(!m_inFile.eof() && m_inFile.good())
		{
			//we need to read up to the newline
			getline(m_inFile , sInput , '\n');
						
			if(ReadAhead)
			{
				//split up the input into tokens of whit spaces
				Tk ti(sInput.c_str() , TokenFinder(" ")); // find white spaces
				Coords cd;

				
				while( ti++ != Tk()) // skip the first as format in "place" , distance , distance , dont need the first
				{
					//convert the string to a double
					cd.x = std::stod (*ti,&sz);
					//cd.x = std::stoi (*ti,&sz);
					++ti; //get the next token
					cd.y = std::stod (*ti,&sz);
					//add the struct to the vector to to be stored
					m_strVector.push_back(cd);
				
					++ti;
				}


				//m_strVector.push_back(sInput);
			}
			//start of the COORDS in the TSP format , so read the coords from here....
			if(sInput == "NODE_COORD_SECTION" || sInput == "EDGE_WEIGHT_SECTION")

			{
				ReadAhead = true;
			}
				
		}
		
		int Vsize = m_strVector.size();


		

	   	static MatrixArrayType theMatrix(Vsize , Vsize);
	  
	    matrix_calculate(theMatrix);
	  
		

		m_theMatrix = &theMatrix;

	 	
		//calculateNearestNeigbhor();
		//std::uniform_int_distribution<int> distribution(0,Vsize -1);
	}
	cout << "init..";  
	return status;

}

int FileReader::getGaussianDistributionRowIndex(std::uniform_int_distribution<int> &distribution)
{
	//std::uniform_int_distribution<int> distribution1(0,51);
	int dist = distribution(m_generator);
    return dist;
}


double FileReader::calculateEuclidianDistance(double x1, double y1,double x2, double y2) 
{
      double x = x1 - x2;
   double y = y1 - y2;
	double dist=(x*x) + (y*y);           //calculating distance by euclidean formula
   return  sqrt(dist);    
}

double FileReader::calculateRoundDistance(double x1, double y1,double x2, double y2) 
{
   double x = x1 - x2;
   double y = y1 - y2;

 //  cout << x << " y " << y << ",";
   double dist=((x*x) + (y*y));           //calculating distance by euclidean formula
   double round =  sqrt(dist) + 0.5;  
   return  round;
}


MatrixArrayType& FileReader::calculateNearestNeigbhor(int NUMOFANTS)
{
	int n = m_strVector.size();
	int i = 0;
		
	double *pHelpArray = new double[m_strVector.size() ];
	double *pdistanceArray = new double[m_strVector.size() ];
	*(pdistanceArray) = std::numeric_limits<double>::max();


	static MatrixArrayType neigbhors(n , NUMOFANTS);

		for (int node = 0; node < n; node++) 
		{ 

			for (i = 0; i < n; i++) 
			{
				*(pdistanceArray+i)=(*m_theMatrix)[node][i];
				*(pHelpArray+i)=  i;
			}
			double max = std::numeric_limits<double>::max() - 1; 
			*(pdistanceArray+node) = 100000000;  // set to a large value .. 
			sort2(pdistanceArray, pHelpArray, 0, n - 1);
			for (i = 0; i < NUMOFANTS; i++) 
			{
				neigbhors[node][i] = *(pHelpArray+i);
			}
	}
	this->m_nearestneighborMatrix = &neigbhors;
	return *m_nearestneighborMatrix;
}
void FileReader::matrix_calculate(MatrixArrayType& theMatrix)
{

	int i, j;
	int N = (int) m_strVector.size();
	
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++)
		{
			//theMatrix[i][j] = 0;
			theMatrix[i][j] = (int) calculateRoundDistance(m_strVector[i].x, m_strVector[i].y, m_strVector[j].x, m_strVector[j].y);
		//	theMatrix[i][j] += (coords[i][Y_VAL] - coords[j][Y_VAL]) * 	(coords[i][Y_VAL] - coords[j][Y_VAL]);
			//theMatrix[i][j] = sqrt(ddistance[i][j]);	
		//	theMatrix[j][i] = theMatrix[i][j];
		}
	}
	return ; 


	
}








FileReader::~FileReader(void)
{
	m_inFile.close();
}

