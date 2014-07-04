// Ant4.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"


#include <vector>


#include <time.h>
#include <math.h>
#include <random>

#include "FileReader.h"
#include "Helpers.h"



#define START_NODE 0
double RHO=0.5; //decay rate
int N_GENERATIONS=50000;
int ALPHA= 1;
int BETA= 2;
int Q= 52;
bool EAS=false;
bool ACS=!EAS;
int num_ants =52;

int neg = 20;

int elitist_ants = 52;

#define randomize() (srand((unsigned int)time(0)))

void calculate_city_distances(void);
double eval_tour(std::vector<int> &tour);
void adjust_pheremone(std::vector<int> &tour , double tour_length);
void evaporate_pheremone(void);

bool reset_ant=false;

int num_cities;
typedef std::vector<std::vector<double> > MatrixArrayType;
MatrixArrayType pherom;
MatrixArrayType ddistance;

MatrixArrayType nn_list;
void ddrdoobs();
double node_branching(double l);


//std::vector<std::vector<int> > nn_list(52);

///ranked ant system
int ras_ranks; /*
			   * additional parameter for rank-based version
			   * of ant system
			   */

double trail_max; /* maximum pheromone trail in MMAS */
double trail_min; /* minimum pheromone trail in MMAS */

double  trail_0;

double lambda = 0.05;

//std::vector<int>  restart_best_ant(52);
//std::vector<int> best_so_far_ant(52);
std::vector<int> nntour(52);


template <class RandomAccessIterator, class URNG>
  void shuffle (RandomAccessIterator first, RandomAccessIterator last, URNG&& g)
{
  for (auto i=(last-first)-1; i>0; --i) {
    std::uniform_int_distribution<decltype(i)> d(0,i);
    swap (first[i], first[d(g)]);
  }
}

void ea_adjust_pheremone_weihhted(std::vector<int> &tour , double tour_length, int weight)
{
	
	double d_tau = weight / tour_length; 
	
	//double tour_length1 = eval_tour (tour); 
	for(int city = 1; city < tour.size(); city++)
	 {
	      int from = tour[city-1];
		  int to = tour[city];
 
		  // eq 14.2 / 14.3
		  pherom[from][to] += d_tau ;
		  pherom[to][from] += d_tau ; 
	
  }


}


double eval_tour (std::vector<int> &tour)
{
	int i, start, end;
	double d = 0;
	
	for (i=1; i<tour.size(); i++) {
		start = tour[i-1];
		end = tour[i];
		d += ddistance[start][end];	
	}

	d += ddistance[tour[num_cities-1]][tour[0]];
	
	return d;
}

void adjust_pheremone(std::vector<int> &tour , double tour_length)
{
	//double tour_length1 = eval_tour (tour); 
	for(int city = 1; city < tour.size(); city++)
	 {
	      int from = tour[city-1];
		  int to = tour[city];
 
		  // eq 14.2 / 14.3
		  pherom[from][to] += (Q/tour_length) ;
		  pherom[to][from] += (Q/tour_length) ; 
	
  }

}


void evaporate_pheremone (void)
{
	int from, to;
	
	for (from=0; from < num_cities; from++) {
		for (to=0; to<from; to++) {
			pherom[from][to] *= (1-RHO);
			pherom[to][from] *= (1-RHO);

		}
	}
}



void Init()
{
	printf("\nWelcome to ACO");
	printf("\nLast compiled on %s, %s", __DATE__, __TIME__);
	printf("\n");

    randomize();

	//string file("berlin52");
	//string file("ulysses16");
	string file("eil51");
	//string file("brg180");
	//string file("fl1400");
	//pcb1173
	//	pcb442
	//	pr2392
	//	rat783
	//	lin318
	//string file("lin318");
	FileReader fileReader(file);
	fileReader.Read();
	cout << "\nShorest optimized distance is = " << fileReader.ShowSolution();

	//get nn negbour list
	 num_cities = fileReader.getsize();
	 MatrixArrayType theMatrix1(num_cities , num_cities);
	 MatrixArrayType theMatrix2= *(fileReader.m_theMatrix);
	// fileReader.calculateNearestNeigbhor();
	// MatrixArrayType negibours= *(fileReader.m_nearestneighborMatrix);
	 if(neg >= num_cities )
		 neg = num_cities-1;

	 nn_list = fileReader.calculateNearestNeigbhor(neg);

	 pherom = theMatrix1 ;
	 ddistance = theMatrix2;
	 
    

}
void CreateInitialPhero()
{
	/* Set up small, random amounts of pheremone to trails */
	for (int i=0; i<num_cities; i++)	{
		for (int j=0; j<num_cities; j++)
			pherom[i][j] =  0.1 * rand() / (double)RAND_MAX;
	}
	
}

void ddrdoobs()
{
	int  k , epoch=0;
	double x, d ;

	Init();

	
	
	//std::default_random_engine m_generator;
	//std::uniform_int_distribution<int> distribution(0,num_cities);
	
		
	
	std::vector<double> t_prob(num_cities);
	
	
	std::vector<int> best_tour(num_cities);
	std::vector<bool> visited(num_cities);

	double min_tour = 10000000000;
		std::vector<double> strength(num_cities+1);
		CreateInitialPhero(); 

//			std::vector<double> strength(num_cities+1);
		
    
	

    /* Main loop starts here */
	for (int t1=0; t1<N_GENERATIONS; t1++)
	{
			
		
		std::random_device random_dev;
		std::mt19937       generator(random_dev());

		/* 1 ant for 5 city */
		for (int t=0; t<num_ants; t++)
		{	
			
				double bestdistance = 10000000000;
			//cout << "init" << t;
			/* Initialise tour */
				std::vector<int> tour(num_cities);
			for (int i=0; i<num_cities; i++) {
				visited[i] = false;
				tour[i] = -1;
			}
			
			/* Tour always starts at position 0 */
			
			std::vector<int> vec(num_cities);
				for (int i=0; i<num_cities; i++) {
					vec[i]=i;
				}
		   
			shuffle(vec.begin(), vec.end(), generator);
			std::vector<int>::iterator it1 = vec.begin();
			int start = *it1;// distribution(m_generator);
			tour[0] = start;
			visited[start] = true;

		
			
			/* j marks where we are in constructing the tour, 
			   k is the next node to visit */
			for (int j=1; j<num_cities; j++) {
		

			
			/* Get transition probabilities - cannot choose a 
			   node that has previously been used */

				for (std::vector<int>::iterator it = vec.begin() ; it != vec.end(); it++ )
				{
								int i = *it;
				
				//for (int i=0; i<num_cities; i++) {
				
					t_prob[i] = 0.0;
					if (visited[i] == false)	
					{
						int tx = tour[j-1];
						double x = pherom[i][tx];
						//try nearest neb..
						
						
						double y = ddistance[i][tx];
					    
						t_prob[i] = pow(x , ALPHA) / pow(y,BETA);
						
						
					}
				}
			
				/* Roulette wheel selection - equivalent to probability-
			   	   based selection of next city. The index of the next 
				   city on the tour is given in k */
				strength[0] = 0;			
				for (int z =0; z<num_cities; z++)
					strength[z+1] = t_prob[z] + strength[z];

				long _rnd = rand();

				if (_rnd == 0)
					_rnd++;

				x = strength[num_cities] * (double)_rnd / RAND_MAX;	
				
				k = 0;
				while (!((strength[k] <= x) && (x <= strength[k+1])))
					k++;
			

				tour[j] = k;				
				visited[k] = true;								
			}
					
			/* Save best tour so far */
			double d = eval_tour(tour);
		
	
			if (d < min_tour)
			{
				min_tour = d;
				for (int i=0; i<num_cities; i++)
					best_tour[i] = tour[i];


				printf("\n");
				printf("Shortest tour to date: distance =%.2f at Epoch %d \n" , eval_tour(best_tour) , t1);
				for (int j=0; j<num_cities; j++) 
				printf("%i|",best_tour[j]+1);
				printf("\n tour length: %d" , best_tour.size());
				printf("\n");
				
			}
			if(EAS)
			{
				if(d < bestdistance)
				{
					bestdistance = d;
					ea_adjust_pheremone_weihhted(tour , d , elitist_ants);

				}
			}
			if(ACS)
			{
				adjust_pheremone (tour , d);
			}
		
			
		}

		epoch = t1;

			
				//for (int j=0; j<num_cities; j++) 
				//printf("%i|",best_tour[j]+1);
				//printf("\n tour length: %d" , best_tour.size());
				//printf("\n");

	if(t1%1000 == 0)
	{
			printf("\n");
				printf("Shortest tour to date: distance =%.2f at Epoch %d \n" , eval_tour(best_tour) , t1);
				for (int j=0; j<num_cities; j++) 
				printf("%i|",best_tour[j]+1);
				printf("\n tour length: %d" , best_tour.size());
				printf("\n");

	}
		evaporate_pheremone ();	

		
		
		
	}
	
	printf("\n");
	printf("Shortest tour found: distance =%.2f at Epoch %d \n" , eval_tour(best_tour) , epoch);
	for (int j=0; j<num_cities; j++) 
		printf("%i|",best_tour[j]+1);
	printf("\n tour length: %d" , best_tour.size());
	printf("\n");
	//printf("Length of shortest tour found: %.2f\n", eval_tour(best_tour));
	///printf("This is %.2f%% longer than the best solution known\n",	100*((eval_tour(best_tour) / BEST_SOLUTION)-1));
}


void init_pheromone_trails(double initial_trail)
 {
	int i, j;

	// TRACE ( System.out.println(" init trails with %.15f\n",initial_trail); );

	/* Initialize pheromone trails */
		for (i = 0; i < num_cities; i++)
		{
			for (j = 0; j <= i; j++) 
			{
				pherom[i][j] = initial_trail;
				pherom[j][i] = initial_trail;
			}
		}
}



void empty_memory(std::vector<bool> &antsvisted)
{
	for (int i=0; i<num_cities; i++)
	{
		antsvisted[i] = false;
	}
}

void choose_closest_next(std::vector<bool> &antsvisted , std::vector<int> &nntour)
{
	int city, current_city, next_city, min_distance;
	next_city =num_cities;

	for(int phase = 1; phase<num_cities; phase++)
	{
		current_city = nntour[phase - 1];
		min_distance = 1000000000; /* Search shortest edge */
		
			for (city = 0; city < num_cities; city++)
			{
				if (antsvisted[city])
				; /* city already visited */
				else
				{
					if (ddistance[current_city][city] < min_distance)
					{
						next_city = city;
						min_distance = ddistance[current_city][city];
					}
				}
			}
	
	
		nntour[phase] = next_city;
		antsvisted[next_city] = true;
	}

}

double nn_tour(std::vector<bool> &antsvisted)
{
	std::vector<int> nntour(num_cities);
	int phase, help;
	empty_memory(antsvisted);
	phase = 0; /* counter of the construction steps */
	int rnd= (rand()%(antsvisted.size()-1))+1;
	nntour[0] =rnd;
	antsvisted[rnd] = true;
	choose_closest_next(antsvisted, nntour);
	return eval_tour(nntour);
}
void global_update_pheromone(std::vector<int> &AntTour)
{

	double tour_length = eval_tour (AntTour); 
	for(int city = 1; city < AntTour.size(); city++)
	 {
	      int from = AntTour[city-1];
		  int to = AntTour[city];
 
		  // eq 14.2 / 14.3
		  pherom[from][to] += (Q/tour_length) ;
		  pherom[to][from] += (Q/tour_length) ; 
	
  }

}
void check_pheromone_trail_limits()
    /*
     * FUNCTION: only for MMAS without local search:
     * keeps pheromone trails inside trail limits
     * INPUT: none
     * OUTPUT: none
     * (SIDE)EFFECTS: pheromones are forced to interval [trail_min,trail_max]
     */
    {
	int i, j;

	// TRACE ( System.out.println("mmas specific: check pheromone trail limits\n"); );

	for (i = 0; i < num_cities ; i++) {
	    for (j = 0; j < i; j++) {
		if (pherom[i][j] < trail_min) {
		    pherom[i][j] = trail_min;
		    pherom[j][i] = trail_min;
		} else if (pherom[i][j] > trail_max) {
		    pherom[i][j] = trail_max;
		    pherom[j][i] = trail_max;
		}
	    }
	}
    }


//	mmas_update(t1 , iteration_best_ant , restart_best_ant_path , best_to_date_tour);	
void mmas_update(int iteration , std::vector<int> &iteration_best_ant , std::vector<int> &restart_best_ant_path ,  std::vector<int> &best_to_date_tour )
{
	

	if (iteration % 20 == 0) //even
	{
	   global_update_pheromone(best_to_date_tour);  //best to date for even
	} 
	else //odd
	{
	    if(reset_ant)
			global_update_pheromone(restart_best_ant_path);
	    else
			global_update_pheromone(iteration_best_ant);
	}

}

//There are variants in the selection of the ants allowed to update pheromones: the best to date ant, 
//or the best for current iteration, or the best after latest reset ant, or the best to date ant for even iterations, 
//and the best for iteration for odd iterations.
void mmAnts(void)
{
	int  k , epoch=0;
	double x, d ;

	double  branching_factor;

	Init();
	std::vector<double> t_prob(num_cities);
	std::vector<bool> visited(num_cities);

	std::vector<double> strength(num_cities+1);


	//calculate min max values inital
	double tourlength = nn_tour(visited);
    trail_max = 1. / ((RHO) * tourlength);
	trail_min = trail_max / (2. * num_cities);
	init_pheromone_trails(trail_max);

	
	//the ant paths
	std::vector<int> best_to_date_tour(num_cities);
	std::vector<int> iteration_best_ant(num_cities);
	std::vector<int> restart_best_ant_path(num_cities);


	double best_ant_to_date = 10000000000;


    /* Main loop starts here */
	for (int t1=0; t1<N_GENERATIONS; t1++)
	{
			
		double restart_best_ant = 10000000000;

		std::random_device random_dev;
		std::mt19937       generator(random_dev());

		/* 1 ant for 5 city */
		for (int t=0; t<num_ants; t++)
		{	
			double iteration_best_ant_distance = 10000000000;
			
		
			//cout << "init" << t;
			/* Initialise tour */
			std::vector<int> tour(num_cities);
		

			for (int i=0; i<num_cities; i++) {
				visited[i] = false;
				tour[i] = -1;
			}
			
			/* Tour always starts at position 0 */
			
			std::vector<int> vec(num_cities);
				for (int i=0; i<num_cities; i++) {
					vec[i]=i;
				}
		   
		//	shuffle(vec.begin(), vec.end(), generator);
			std::vector<int>::iterator it1 = vec.begin();
			int start = *it1;// distribution(m_generator);
			tour[0] = start;
			visited[start] = true;

		
			
			/* j marks where we are in constructing the tour, 
			   k is the next node to visit */
			for (int j=1; j<num_cities; j++)
			{
				for (std::vector<int>::iterator it = vec.begin() ; it != vec.end(); it++ )
				{
						int i = *it;
						t_prob[i] = 0.0;
						if (visited[i] == false)	
						{
							int tx = tour[j-1];
							double x = pherom[i][tx];
							double y = ddistance[i][tx];
							t_prob[i] = pow(x , ALPHA) / pow(y,BETA);
						
						
						}
				}
			
					/* Roulette wheel selection - equivalent to probability-
			   		   based selection of next city. The index of the next 
					   city on the tour is given in k */
					strength[0] = 0;			
					for (int z =0; z<num_cities; z++)
						strength[z+1] = t_prob[z] + strength[z];

					long _rnd = rand();

					if (_rnd == 0)
						_rnd++;

					x = strength[num_cities] * (double)_rnd / RAND_MAX;	
				
					k = 0;
					while (!((strength[k] <= x) && (x <= strength[k+1])))
						k++;
			

					tour[j] = k;				
					visited[k] = true;								
				}
					
				/* Save best tour so far */
				double d = eval_tour(tour);
				//best ant so far for itteration
				if (d < iteration_best_ant_distance)
				{
					//copy the tour
					iteration_best_ant_distance = d;
					for (int i=0; i<num_cities; i++)
					{
						iteration_best_ant[i] = tour[i];
					}
			     	branching_factor = node_branching(lambda);
					double p_x = exp(log(0.05) / num_cities);
					trail_min = 1. * (1. - p_x) / (p_x * (double) ((num_ants + 1) / 2));
					trail_max = 1. / ((RHO) * eval_tour(iteration_best_ant));
					trail_0 = trail_max;
					trail_min = trail_max * trail_min;
				}
			if (d < restart_best_ant)
			{
				restart_best_ant = d;
				for (int i=0; i<num_cities; i++)
				{
					restart_best_ant_path[i] = tour[i];
				}
		
			}
			if (d < best_ant_to_date)
			{
				best_ant_to_date = d;
				for (int i=0; i<num_cities; i++)
					best_to_date_tour[i] = tour[i];
	
			}
		
	  	}
		
			//mmas_update( , min_tour);
		mmas_update(t1 , iteration_best_ant , restart_best_ant_path , best_to_date_tour);	
	
       //end of ants walk
	    evaporate_pheremone ();	
		check_pheromone_trail_limits();
			

		if ((t1 % 100) == 0) 
		{
			//create nn_list ..
			branching_factor = node_branching(lambda);
			double  branch_fac = 1.00001;

	
			  if ( (branching_factor < branch_fac)
					&& (t1  > 150)) 
			  {
				cout << "INIT TRAILS!!!\n";
				restart_best_ant = 10000000000;
				init_pheromone_trails(trail_max);
				reset_ant = true;
		 	  }

		}
		else
		{
			reset_ant = false;
		}
		printf("\n");
		printf("Shortest tour found: distance =%.2f at Epoch %d \n" , eval_tour(best_to_date_tour) , t1);
		printf("iteration_best_ant tour found: distance =%.2f at Epoch %d \n" , eval_tour(iteration_best_ant) , t1);
		printf("restart_best_ant tour found: distance =%.2f at Epoch %d" , eval_tour(restart_best_ant_path) , t1);
	

	//end itteration	
	}
	
	printf("\n");
	printf("Shortest tour found: distance =%.2f  \n" , eval_tour(best_to_date_tour) );
	for (int j=0; j<num_cities; j++) 
		printf("%i|",best_to_date_tour[j]+1);
	printf("\n tour length: %d" , best_to_date_tour.size());
	printf("\n");
	//printf("Length of shortest tour found: %.2f\n", eval_tour(best_tour));
	///printf("This is %.2f%% longer than the best solution known\n",	100*((eval_tour(best_tour) / BEST_SOLUTION)-1));

}
 
 double node_branching(double l)
 {
	int i, m;
	double min, max, cutoff;
	double avg;

	std::vector<double> num_branches(num_cities);

	
	for (m = 0; m < neg; m++) 
	{
	   
	    min = pherom[m][nn_list[m][1]];
	    max = pherom[m][nn_list[m+1][1]];


	    for (i = 1; i < neg; i++) 
		{
			if (pherom[m][nn_list[m][i]] > max)
				max = pherom[m][nn_list[m][i]];
			if (pherom[m][nn_list[m][i]] < min)
				min = pherom[m][nn_list[m][i]];
	    }
	    cutoff = min + l * (max - min);

	    for (i = 1; i < neg; i++) {
		if (pherom[m][nn_list[m][i]] > cutoff)
		    num_branches[m] += 1.;
	    }
	}
	avg = 0.;
	for (m = 0; m < num_cities; m++) {
	    avg += num_branches[m];
	}
	
	return (avg / (double) (num_cities * 2));
}




int _tmain(int argc, _TCHAR* argv[])
{
	
	//ddrdoobs();
	mmAnts();
	return 0;
}

