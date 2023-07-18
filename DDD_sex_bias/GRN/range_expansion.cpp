/*
	Copyright (C) 2023  Jhelam N. Deshpande
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
*/

//============================================================================


#include <cstdlib>		//standard C library
#include <iostream>		//standard input output function library
#include <fstream>		//file stream library
#include <numeric>
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>    //random number generator gsl library      
#include <gsl/gsl_randist.h> //gsl random distribution library
#include <algorithm>
#include <math.h>
using namespace std;

//----------------------------------------------------------------------------------------
//------------------------------------------------------------Gene-regulatory network size
const int INPUT_SIZE=2;	//size of input layer of the gene-regulatory network
const int REGULATORY_SIZE=4;	//number of regulatory genes in the GRN
const int GRN_TIME=20;	//Time for which the GRN is iterated for an individual

//________________________________________________________________________________________
//----------------------------------------------------------Class defining the individuals
class TInd{       
public:
	TInd();
	double input_weights[2][INPUT_SIZE][REGULATORY_SIZE];          //Weights input to regulatory layer
	double regulatory_weights[2][REGULATORY_SIZE][REGULATORY_SIZE];  //Weights regulatory layer within itself    
	double output_weights[2][REGULATORY_SIZE];      //Weights regulatory layer to output layer
	double regulatory_threshold[2][REGULATORY_SIZE];    //thresholds of regulatory layer
	double regulatory_slope[2][REGULATORY_SIZE];    //slope of regulatory layer
};

//----------------------------------------------------------Constructor for the class TInd
TInd::TInd(){     
	for (int i = 0; i < 2; i++)
	{
		for(int j=0;j<INPUT_SIZE;j++)
		{
			for(int k=0; k<REGULATORY_SIZE;k++)
			{
				input_weights[i][j][k]=0;
			}
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for(int j=0;j<REGULATORY_SIZE;j++)
		{
			for(int k=0; k<REGULATORY_SIZE;k++)
			{
				regulatory_weights[i][j][k]=0;
			}
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for(int j=0;j<REGULATORY_SIZE;j++)
		{
			output_weights[i][j]=0;
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for(int j=0;j<REGULATORY_SIZE;j++)
		{
			regulatory_threshold[i][j]=0;
		}
	}
	for (int i = 0; i < 2; i++)
	{
		for(int j=0;j<REGULATORY_SIZE;j++)
		{
			regulatory_slope[i][j]=0;
		}
	}

	//dispersal_probability=0;
}

//________________________________________________________________________________________
//--------------------------------------------------------------Class defining the patches
class TPatch      
{
public:
	TPatch();
	vector <TInd> females;    //females in the simulation
	vector <TInd> newfemales;   //vector to store new disperses or newborn females
	vector <TInd> males;        
	vector <TInd> newmales;
	double measured_dispersal;    //measures fraction of individuals dispersing
	double grn_mortality;
};

//--------------------------------------------------------Constructor for the class TPatch
TPatch::TPatch(){    
	females.clear();
	males.clear();
	newfemales.clear();
	newmales.clear();
	measured_dispersal=0;
	grn_mortality=0;
}

//________________________________________________________________________________________
//------------------------------------------------------------------------Model Parameters 

int No;         //Initial number of individuals per patch
double LAMBDA;      //Intrinsic growth rate or mean female fecundity---Beverton-Holt model
double ALPHA;     //Intraspecific competition coefficient---Beverton-Holt model
double DISPERSAL_PROB;  // Dispersal probability
int BURN_IN_TIME;       //Number of time steps before range expansions begin (duration of burn-in period)
int REPLICATES;         //Number of replicates
double DISP_MORT; //dispersal mortality
double EXTINCTION_PROB;		//probability of local patch extinction
double VARIANCE_MIN; 	//variance of mutation
double MUT_RATE; 	//mutation rate at the beginning of the simulation
double MUT_RATE_MIN; 	//mutation rate after 5000 time steps
double GENETIC_VARIATION; 	//standing genetic variation in network parameters---weights, thresholds.etc
double VARIANCE_MAX; 	//standard deviation of mutation effects in the beginning of the simulation
double MUT_SLOW;		//slow mutation rate, set to 0
double VAR_SLOW;		//slow standard deviation, set to 0

const int world_size_x = 500;//size of landscape in x direction
const int world_size_y = 5;//size of landscape in y direction
const int burn_in_x=10;//size of landscape before range expansion
TPatch world[world_size_x][world_size_y];//Creating patches

const gsl_rng *gBaseRand;

//________________________________________________________________________________________
//------------------------------------------------------Initialize Random Number Generator

void specify_rng(unsigned long randSeed)
{
	gBaseRand = gsl_rng_alloc(gsl_rng_rand);

	srand(randSeed);
	unsigned long r = rand();
	gsl_rng_set(gBaseRand, r);
}

//________________________________________________________________________________________
//-------------------------------------------------------------------------Simplifications

//-------------------------------------------------Simplify Random Drawing between 0 and 1

double ran()
{
	return gsl_rng_uniform(gBaseRand);
}

//---------------------------------------------------------------Simplify Gaussian Randoms

double gauss(double sd)
{
	return gsl_ran_gaussian(gBaseRand,sd);
}

//-----------------------------------------------------------------Simplify Poisson Random

int poisson(double sd)
{
	return gsl_ran_poisson(gBaseRand,sd);
}


const int RS = 100;                 // random seed


// random number generator
		//specify_rng(time(NULL));
		//specify_rng(RS);


//________________________________________________________________________________________
//--------------------------------------------------------Miscellaneous required functions

vector<double> quartiles(vector<double> trait)
{
  sort(trait.begin(), trait.end());//sort vector
  double q1;  //stores quartile
  double q2;  //stores quartile
  double q3;  //stores quartile
  int n=trait.size();  //stores size of vector
  q2=trait.at((n+1)/2)+double((n+1)%2)*(trait.at(((n+1)/2)+1)-trait.at((n+1)/2))/double(2);
  q1=trait.at((n+1)/4)+double((n+1)%4)*(trait.at(((n+1)/4)+1)-trait.at((n+1)/4))/double(4);
  q3=trait.at(3*(n+1)/4)+double(3*(n+1)%4)*(trait.at((3*(n+1)/4)+1)-trait.at(3*(n+1)/4))/double(4);
  vector<double> quartile;
  quartile.push_back(q1);
  quartile.push_back(q2);
  quartile.push_back(q3);
  return quartile;

}


//---------------------------------------------return a subset of individuals given patches
vector <TInd> subset_metapopulation(int x_min, int x_max)
{
  vector <TInd> all_individuals;//stores all individuals in patches between x_min to x_max 
  for(int x=x_min;x<x_max;x++)
  {
  	for(int y=0; y<world_size_y; y++)
  	{
  		for(int f=0;f<world[x][y].females.size();f++)
  		{
  			all_individuals.push_back(world[x][y].females.at(f));
  		}
  		for(int m=0;m<world[x][y].males.size();m++)
  		{
  			all_individuals.push_back(world[x][y].males.at(m));
  		}
  	}
  }
  int no_of_individuals=100;
  if(no_of_individuals>all_individuals.size())
  	no_of_individuals=all_individuals.size();
  vector<TInd> random_individuals;
  for(int ind=0;ind<no_of_individuals;ind++)//draw no_of_individuals individuals
  { 
    int position_random=floor(ran()*all_individuals.size());//individuals drawn at random
    random_individuals.push_back(all_individuals.at(position_random));
}
return random_individuals;
}


//----------------------------------------------------function returns mean of two numbers
double mean(double a, double b) 
{
	return ((a+b)/double(2));
}




//-------------------------------------------function returns alleles with mutation added
double mutate(double d,int t) 
{
	double new_trait;			//stores value of the allele post mutation
	double mut;						//mutation rate
	double var;						//standard deviation of mutation effect
	if(t<5000)						//for the first 5000 time steps
		mut=-(t*(MUT_RATE-MUT_RATE_MIN)/double(5000))+MUT_RATE; //mutation rate decreases linearly with time
	else mut=MUT_RATE_MIN;			//mutation rate is constant after the first 5000 time steps
	if(t<5000)						//same as mutation rate
		var=-(t*(VARIANCE_MAX-VARIANCE_MIN)/double(5000))+VARIANCE_MAX;
	else var=VARIANCE_MIN;
	if(ran()<mut)					       //mutation is introduced to the allele with probability mut
		new_trait= d+ gauss(var);  //mutation size is drawn from a normal distribution with mean 0 and standard deviation var
	else new_trait=d;
	if(ran()<MUT_SLOW)
		new_trait= new_trait+ gauss(VAR_SLOW);  
	return new_trait;							//return allele value with mutation added

}

//________________________________________________________________________________________
//--------------------------------------------------------Inputting and setting parameters
void set_parameters()  
{
	string para[15];	//stores model parameter values
	string line;
	ifstream myfile ("input.txt");	//read parameter input file
	int count=0;
	if (myfile.is_open())
	{
		while ( getline (myfile,line))
		{
			if(count%2==1)
			{
				para[count/2]=line;	//read every alternate line
			}
			count++;

		}
		myfile.close();
	}
	else cout << "Unable to open file";
	No = (int) std::atof(para[0].c_str());							//sets initial population size per patch
	LAMBDA=std::atof(para[1].c_str());									//sets intrinsic growth rate mean fecundity of the female---Beverton-Holt model
	ALPHA=std::atof(para[2].c_str());										//sets intra-specific competition coefficient of the Beverton-Holt model
	BURN_IN_TIME= (int) std::atof(para[3].c_str());			//sets the number of time steps before the beginning of range expansion
	REPLICATES=(int) std::atof(para[4].c_str());				//sets number of replicate simulations that are run
	DISPERSAL_PROB=std::atof(para[5].c_str());					//sets dispersal probability if it is not genetically encoded
	DISP_MORT=std::atof(para[6].c_str());								//sets dispersal costs
	EXTINCTION_PROB=std::atof(para[7].c_str());					//sets probability of random patch extinction per time step per patch
	VARIANCE_MIN=std::atof(para[8].c_str());						//sets mutation effects after 5000 time steps
	MUT_RATE=std::atof(para[9].c_str());								//sets maximum mutation probability before 5000 time steps
	MUT_RATE_MIN=std::atof(para[10].c_str());						//sets mutation probability after 5000 time steps
	GENETIC_VARIATION=std::atof(para[11].c_str());			//sets standing genetic variation
	VARIANCE_MAX=std::atof(para[12].c_str());						//sets maximum mutation effects before 5000 time steps
	MUT_SLOW=std::atof(para[13].c_str());								//slow mutation rate---set to 0
	VAR_SLOW=std::atof(para[14].c_str());								//effects of slow mutations---set to 0
}



//________________________________________________________________________________________
//----------------------------------------------------------------------Initialising model
void init_world()
{
	for(int x=0;x<world_size_x;x++)
	{
		for(int y=0;y<world_size_y;y++) //clear all individuals in a patch for all patches
		{

			world[x][y].females.clear(); 
			world[x][y].males.clear();
			world[x][y].newfemales.clear(); 
			world[x][y].newmales.clear();
		}
	}//initialise individuals only in the range core
	for(int x=world_size_x/2-burn_in_x/2;x<world_size_x/2+burn_in_x/2;x++)	//go through all patches in the range core
	{
		for(int y=0;y<world_size_y ;y++)
		{
		world[x][y].females.clear(); //clear all individuals in a patch
		world[x][y].males.clear();
		for(int n=0; n<No; n++)	//create No initial individuals per patch
		{
			TInd newind;      //creat a new individual
			//initialising the genotype of the individual
			 //initialising network parameters for dispersal GRN with standing genetic variation
			for (int i = 0; i < 2; i++) //initialising  network parameters 
			{
				for(int j=0;j<INPUT_SIZE;j++)
				{
					for(int k=0; k<REGULATORY_SIZE;k++)
					{
						newind.input_weights[i][j][k]=gauss(GENETIC_VARIATION);
					}
				}
			}
			for (int i = 0; i < 2; i++)
			{
				for(int j=0;j<REGULATORY_SIZE;j++)
				{
					for(int k=0; k<REGULATORY_SIZE;k++)
					{
						newind.regulatory_weights[i][j][k]=gauss(GENETIC_VARIATION);
					}
				}
			}
			for (int i = 0; i < 2; i++)
			{
				for(int j=0;j<REGULATORY_SIZE;j++)
				{
					newind.output_weights[i][j]=gauss(GENETIC_VARIATION);
				}
			}
			for (int i = 0; i < 2; i++)
			{
				for(int j=0;j<REGULATORY_SIZE;j++)
				{
					newind.regulatory_threshold[i][j]=gauss(GENETIC_VARIATION);
				}
			}
			for (int i = 0; i < 2; i++)
			{
				for(int j=0;j<REGULATORY_SIZE;j++)
				{
					newind.regulatory_slope[i][j]=gauss(GENETIC_VARIATION);
				}
			}
			if(ran()<0.5)	//individual has equal chance of being male or female
			world[x][y].females.push_back(newind); //adding new females
			else world[x][y].males.push_back(newind);//adding new males
		}

	}
}
}

//________________________________________________________________________________________
//------------------------------------------------------------------------Output genotypes 
void output_genotype(ofstream& op5, int x_min,int x_max, int r, int t)
{
	for(int x=x_min; x<x_max;x++)
	{
		for(int y=0; y<world_size_y;y++)
		{

			for(int f=0;f<world[x][y].females.size();f++)
			{
				//output a sample half of the individuals in a specified region of the landscape
				if(ran()<0.5){
					op5 <<r<<" "<<t<<" "<<x<<" "<<y<<" ";
					for(int i=0;i<REGULATORY_SIZE;i++)
					{

						for(int j=0; j<INPUT_SIZE;j++)
						{

							op5<<mean(world[x][y].females.at(f).input_weights[0][j][i],
								world[x][y].females.at(f).input_weights[1][j][i])<<" ";
						}

						for(int j=0; j<REGULATORY_SIZE;j++)
						{

							op5<<mean(world[x][y].females.at(f).regulatory_weights[0][i][j],
								world[x][y].females.at(f).regulatory_weights[1][i][j])<<" ";
						}


						op5<<mean(world[x][y].females.at(f).regulatory_threshold[0][i],
							world[x][y].females.at(f).regulatory_threshold[1][i])<<" ";
						op5<<mean(world[x][y].females.at(f).output_weights[0][i],
							world[x][y].females.at(f).output_weights[1][i])<<" ";
						op5<<mean(world[x][y].females.at(f).regulatory_slope[0][i],
							world[x][y].females.at(f).regulatory_slope[1][i])<<" ";
					}


					op5<<endl;
				}
			}

			for(int m=0;m<world[x][y].males.size();m++)
			{
				if(ran()<0.5){
					op5 <<r<<" "<<t<<" "<<x<<" "<<y<<" ";
					for(int i=0;i<REGULATORY_SIZE;i++)
					{

						for(int j=0; j<INPUT_SIZE;j++)
						{

							op5<<mean(world[x][y].males.at(m).input_weights[0][j][i],
								world[x][y].males.at(m).input_weights[1][j][i])<<" ";
						}

						for(int j=0; j<REGULATORY_SIZE;j++)
						{

							op5<<mean(world[x][y].males.at(m).regulatory_weights[0][i][j],
								world[x][y].males.at(m).regulatory_weights[1][i][j])<<" ";
						}


						op5<<mean(world[x][y].males.at(m).regulatory_threshold[0][i],
							world[x][y].males.at(m).regulatory_threshold[1][i])<<" ";
						op5<<mean(world[x][y].males.at(m).output_weights[0][i],
							world[x][y].males.at(m).output_weights[1][i])<<" ";
						op5<<mean(world[x][y].males.at(m).regulatory_slope[0][i],
							world[x][y].males.at(m).regulatory_slope[1][i])<<" ";
					}


					op5<<endl;
				}
			}
		}  

	}
}

//________________________________________________________________________________________
//--------------------Output population density,range front position.etc for a given patch

void output_metapopulation(ofstream& op, int x,int  y, int r, int t, int margin_x_left, int margin_x_right)
{

	op <<r<<" "<<t<<" "<<x<<" "<<y<<" "<<world[x][y].females.size()+world[x][y].males.size()<<" "<<
	double(world[x][y].females.size())/double(world[x][y].females.size()+world[x][y].males.size())<<" "<<
	world[x][y].measured_dispersal<<" "<<world[x][y].grn_mortality<<" ";
	op<<margin_x_left<<" "<<margin_x_right;
	op<<endl;


}



//________________________________________________________________________________________
//--------------------------------------Deciding coordinates of new patch while dispersing
vector<int> decide_patch(int x, int y,int t) 
{ 
	//nearest neighbor 8 dispersal
	vector<int> coordinates;
	int newx=0;
	int newy=0;
	int decider=floor(ran()*double(8));
	switch(decider){
		case 0: newx++; 
		break;
		case 1: newy++;
		break;
		case 2: newx++;
		newy++;
		break;
		case 3:newx--;
		break;
		case 4:newy--;
		break;
		case 5:newx--;
		newy--;
		break;
		case 6: newx++;
		newy--;
		break;
		case 7:newx--;
		newy++;
		break;
		default: cout<<"error in nn8"<<endl;
	}
	newx=newx+x;
	newy=newy+y;
	//before range expansion begins reflecting boundary condition in x direction
	if(t<BURN_IN_TIME)
	{
		if (newx<world_size_x/2-burn_in_x/2)	//if the new patch is beyond the range core in the left direction
			newx=world_size_x/2-burn_in_x/2+1;	//send the individual to one patch to the right
		if(newx>world_size_x/2+burn_in_x/2-1)
			newx=world_size_x/2+burn_in_x/2-2;
	}
	if(t>=BURN_IN_TIME)	 //when range expansions begin
	{
		if (newx<0) //no conditions, does not matter at the landscape end as the simulations stop
			newx=0;
		if(newx==world_size_x)
			newx=world_size_x-1;
	}
	//torus in y direction
	if(newy<0)	//if the individual is at the bottom
		newy=world_size_y - 1;	//send it to the top
	if(newy==world_size_y)
		newy=0;
	coordinates.push_back(newx);	//store the new coordinates of the individual
	coordinates.push_back(newy);
	return coordinates;
}
//________________________________________________________________________________________
//---------------------------------function returns the dispersal rate calculated from GRN
double  grn_output(double inp[INPUT_SIZE],
	double iw[2][INPUT_SIZE][REGULATORY_SIZE], 
	double rw[2][REGULATORY_SIZE][REGULATORY_SIZE],
	double ow[2][REGULATORY_SIZE], 
	double rt[2][REGULATORY_SIZE], double rs[2][REGULATORY_SIZE])
{
	double input_weights[INPUT_SIZE][REGULATORY_SIZE];						//stores mid allelic value of input matrix weights
	double regulatory_weights[REGULATORY_SIZE][REGULATORY_SIZE];	//stores mid allelic value of regulatory matrix weights
	double output_weights[REGULATORY_SIZE];												//stores mid allelic value of output matrix weights
	double regulatory_threshold[REGULATORY_SIZE];									//stores mid allelic value of regulatory thresholds
	double regulatory_slope[REGULATORY_SIZE];											//stores mid allelic value of slopes of regulatory genes

		double genes[REGULATORY_SIZE][GRN_TIME]; //stores values of genes at each iteration
		double output[GRN_TIME];//stores values of output
		output[0]=0;
		double epsilon=0.0001;//tolerance
		double var=0; //stores the variation in the equilibrium phenotype
		for (int i = 0; i < REGULATORY_SIZE; i++) // averaging between the two alleles 
		{
			for(int j=0; j<REGULATORY_SIZE; j++)
			{
				regulatory_weights[i][j]=mean(rw[0][i][j],rw[1][i][j]);
			}
			for(int j=0; j<INPUT_SIZE; j++)
			{
				input_weights[j][i]=mean(iw[0][j][i],iw[1][j][i]);
			}
			regulatory_threshold[i]=mean(rt[0][i],rt[1][i]);
			regulatory_slope[i]=mean(rs[0][i],rs[1][i]);
			output_weights[i]=mean(ow[0][i],ow[1][i]);
			genes[i][0]=2*ran()-1; //initialising the genes randomly with a number between -1 and 1
		}
		for(int time=1;time<GRN_TIME;time++)	//iterate through the gene-regulatory network
		{ 
			double output_store=0.0; //stores GRN output at each time step
			for(int i=0;i<REGULATORY_SIZE;i++)	//go through all the genes in the regulatory layer
			{
				double regulatory_output=0.0;	//stores the input to a gene
				for(int j=0;j<INPUT_SIZE;j++)
				{
					regulatory_output=regulatory_output+input_weights[j][i]*inp[j]; //sums up all the inputs from input layer
				}
				for(int j=0;j<REGULATORY_SIZE;j++)
				{
					regulatory_output=regulatory_output+regulatory_weights[j][i]*genes[j][time-1]; //sums up the input from other regulatory genes
				}
				genes[i][time]=(2.0/(1.0+exp(-(regulatory_slope[i])*(regulatory_output-regulatory_threshold[i]))))-1; //gene expression state for a given iteration
				output_store=output_store+(genes[i][time]*output_weights[i]);  //store the output
			}
			output[time]=output_store;
			if(time >10)
			{
				var=var+(output[time]-output[time-1])*(output[time]-output[time-1]);
			}
		}
		var=sqrt(var)/(GRN_TIME-10);
		if(var<epsilon)	//check for steady gene expression
			return output[GRN_TIME-1] ;
		else
		{ 
			return 42.0;   
		}  
	}

//________________________________________________________________________________________
//-------------------------------------------------------------Density regulation function

	double densReg(double a) 
	{
		return(1 /(1+(a)));
	}



//calculate perturbed phenotype for a given individual
	double perturbed_trait(
		double iw[2][INPUT_SIZE][REGULATORY_SIZE], 
		double rw[2][REGULATORY_SIZE][REGULATORY_SIZE],
		double ow[2][REGULATORY_SIZE], 
		double rt[2][REGULATORY_SIZE], double rs[2][REGULATORY_SIZE])
	{
	//randomly draw the position of the perturbation
		double perturbed=0;
		double inp[1];
		inp[0]=0.5;
		int pos_i=floor(ran()*2);  //choose an allele to perturb
		int pos_j=floor(ran()*(REGULATORY_SIZE+INPUT_SIZE+3));	//choose the position and what part of the GRN to perturb
		int pos_k=floor(ran()*REGULATORY_SIZE);	//choose the position in the matrix
      if(pos_j<REGULATORY_SIZE) //perturb regulatory weights if the number drawn is 0,..,REGULATORY_SIZE-1
      {
      	double rw_new[2][REGULATORY_SIZE][REGULATORY_SIZE];//copy regulatory weights into a new matrix
      	for(int i=0;i<2; i++)
      	{
      		for(int j=0;j<REGULATORY_SIZE;j++)
      		{
      			for(int k=0;k<REGULATORY_SIZE;k++)
      			{
      				rw_new[i][j][k]=rw[i][j][k];
      			}
      		}
      	}
      	rw_new[pos_i][pos_j][pos_k]=rw[pos_i][pos_j][pos_k]+gauss(VARIANCE_MIN); //perturb regulatory weights at the chosen position
      	perturbed= grn_output(inp,iw,rw_new,ow,rt,rs); //calculate perturbed phenotype
      }
      else if(pos_j<REGULATORY_SIZE+INPUT_SIZE)	//if the random number drawn is between REGULATORY_SIZE and REGULATORY_SIZE+INPUT_SIZE-1 then perturb the input weights
      {
      	pos_j=pos_j-REGULATORY_SIZE;	//get the position within the input matrix to perturb
      	double iw_new[2][INPUT_SIZE][REGULATORY_SIZE];//copy the input weights into a new matrix
      	for(int i=0;i<2; i++)
      	{
      		for(int j=0;j<INPUT_SIZE;j++)
      		{
      			for(int k=0;k<REGULATORY_SIZE;k++)
      			{
      				iw_new[i][j][k]=iw[i][j][k];
      			}
      		}
      	}
      	iw_new[pos_i][pos_j][pos_k]=iw[pos_i][pos_j][pos_k]+gauss(VARIANCE_MIN);	//add a perturbation to the chosen input weight
      	perturbed= grn_output(inp,iw_new ,rw,ow,rt,rs);
      }
      else if(pos_j<REGULATORY_SIZE+INPUT_SIZE+1)
      {
        double ow_new[2][REGULATORY_SIZE];//stores perturbed output weights  
        for(int i=0;i<2; i++)
        {
        	for(int j=0;j<REGULATORY_SIZE;j++)
        	{
        		ow_new[i][j]=ow[i][j];
        	}
        }
        ow_new[pos_i][pos_k]=ow[pos_i][pos_k]+gauss(VARIANCE_MIN);
        perturbed= grn_output(inp,iw ,rw,ow_new,rt,rs);
    }
    else if(pos_j<REGULATORY_SIZE+INPUT_SIZE+2)
    {
        double rs_new[2][REGULATORY_SIZE];//stores perturbed slope of the gene        
        for(int i=0;i<2; i++)
        {
        	for(int j=0;j<REGULATORY_SIZE;j++)
        	{
        		rs_new[i][j]=rs[i][j];
        	}
        }
        rs_new[pos_i][pos_k]=rs[pos_i][pos_k]+gauss(VARIANCE_MIN);
        perturbed= grn_output(inp,iw ,rw,ow,rt,rs_new);
    }
    else if(pos_j<REGULATORY_SIZE+INPUT_SIZE+3)
    {
        double rt_new[2][REGULATORY_SIZE];//stores perturbed  threshold of the gene
        for(int i=0;i<2; i++)
        {
        	for(int j=0;j<REGULATORY_SIZE;j++)
        	{
        		rt_new[i][j]=rt[i][j];
        	}
        }
        rt_new[pos_i][pos_k]=rt[pos_i][pos_k]+gauss(VARIANCE_MIN);
        perturbed= grn_output(inp,iw ,rw,ow,rt_new,rs);
    }
    return perturbed;
}



vector <double> sensitivity_to_mutation(int x_min,int x_max)
{
	vector<TInd> individuals;
	vector <double> sensitivity;
	sensitivity.push_back(0);
	sensitivity.push_back(0);
	individuals=subset_metapopulation(x_min,x_max);
	int no_perturbations=100;


	individuals=subset_metapopulation(x_min,x_max);
	int count_disp=0;
int max_density=15;
	for(int p=0;p<individuals.size();p++)
	{
		for(int s=0;s<2;s++){
			  for(int d=0;d<max_density;d++){//all possible values of normalized density

			  	double inp[2];
			  	inp[0]=double(d)/double(10)+0.05;
			  	inp[1]=double(s);
			  	double unperturbed_disp=grn_output(inp,individuals.at(p).input_weights ,individuals.at(p).regulatory_weights,
			  		individuals.at(p).output_weights,individuals.at(p).regulatory_threshold,individuals.at(p).regulatory_slope);
			  	for(int q=0;q<no_perturbations;q++){
			  		double perturbed_disp=perturbed_trait(individuals.at(p).input_weights ,individuals.at(p).regulatory_weights,
			  			individuals.at(p).output_weights,individuals.at(p).regulatory_threshold,individuals.at(p).regulatory_slope);

			  		if(perturbed_disp !=42 && unperturbed_disp !=42){
			  			if(unperturbed_disp>1)
			  				unperturbed_disp=1;
			  			if(unperturbed_disp<0)
			  				unperturbed_disp=0;
			  			if(perturbed_disp>1)
			  				perturbed_disp=1;
			  			if(perturbed_disp<0)
			  				perturbed_disp=0;
			  			sensitivity.at(1)=sensitivity.at(1)+(perturbed_disp-unperturbed_disp)*(perturbed_disp-unperturbed_disp);
			  			count_disp++;
			  		}

			  	}
			  }
			}
		}
			sensitivity.at(1)=sqrt(sensitivity.at(1)/double(count_disp));


			return sensitivity;
		}



vector< vector <double> > trait_calulator(int x_min,int x_max, double sex) //takes the boundaries for which 
{
  vector <TInd> all_individuals;//stores all individuals in patches between x_min to x_max 
  for(int x=x_min;x<x_max;x++)
  {
  	for(int y=0; y<world_size_y; y++)
  	{
  		for(int f=0;f<world[x][y].females.size();f++)
  		{
  			all_individuals.push_back(world[x][y].females.at(f));
  		}
  		for(int m=0;m<world[x][y].males.size();m++)
  		{
  			all_individuals.push_back(world[x][y].males.at(m));
  		}
  	}
  }
  int no_of_individuals=1000;
  int max_density=15;
  vector<TInd> random_individuals;
  for(int ind=0;ind<no_of_individuals;ind++)//draw no_of_individuals individuals
  { 
    int position_random=floor(ran()*all_individuals.size());//individuals drawn at random
    random_individuals.push_back(all_individuals.at(position_random));
}
all_individuals.clear();
vector<vector<double> > quartile;
  for(int d=0;d<max_density;d++){//all possible values of normalised density
    vector<double> disp_all; // stores all dispersal probabilities corresponding to that density
    for(int p=0;p<no_of_individuals;p++){
    	double inp[2];
    	inp[0]=double(d)/double(10)+0.05;
    	inp[1]=sex;
    	double disp=grn_output(inp,
    		random_individuals.at(p).input_weights,
    		random_individuals.at(p).regulatory_weights,random_individuals.at(p).output_weights,
            random_individuals.at(p).regulatory_threshold,random_individuals.at(p).regulatory_slope);//calculates dispersal
    	if(disp<0)
    		disp=0;
    	if(disp>1)
    		disp=1;
    	disp_all.push_back(disp);
    }
    quartile.push_back(quartiles(disp_all));
}
return quartile;
}





//________________________________________________________________________________________
//-------------------------------------------------------------------Life cycle procedures
//---------------------------------------------------------------------Dispersal procedure
void disperse(int t)
{
	for(int x=0; x<world_size_x;x++)    //clears newmales and newfemales vectors
	{
		for(int y=0; y<world_size_y;y++)
		{
			world[x][y].newfemales.clear();	
			world[x][y].newmales.clear();
		}
	}
	for(int x=0;x<world_size_x;x++)   //dispersal loop---go through all patches
	{
		for(int y=0;y<world_size_y;y++)
		{
			int count_dispersers=0;	//stores number of individuals leaving a patch
			//inputs of the GRN stored in inp
			double population_density=double(world[x][y].males.size()+world[x][y].females.size());
			double normalised_density=population_density*ALPHA/(LAMBDA-1);	//calculate local population density normalised by expected equilibrium density of the Beverton-Holt model
			int grn_mortality_count=0;	//count number dying because of unviable GRN
			for(int f=0;f<world[x][y].females.size();f++)	//go through all the females
			{
				double dispersal_probability;	//stores individual dispersal probability
				double inp[2];	//stores input to the GRN
				inp[0]=normalised_density;	//the GRN takes the local normalised population density as an input
				inp[1]=0.0;
				//dispersal probability calculation
				dispersal_probability=grn_output(inp,
					world[x][y].females.at(f).input_weights,
					world[x][y].females.at(f).regulatory_weights,world[x][y].females.at(f).output_weights,
					world[x][y].females.at(f).regulatory_threshold,world[x][y].females.at(f).regulatory_slope);	//calculate dispersal probability
				//dispersal
				if(int(dispersal_probability) != 42)	//if the dispersal GRN in viable
				{
					if(ran()< dispersal_probability)	//disperse with specified probability
					{
						std::vector<int> coor=decide_patch(x,y,t);	//stores the target patch 
						//dispersal mortalilty
						if(ran()>DISP_MORT)	//if the individual survives dispersal
							world[coor.at(0)][coor.at(1)].newfemales.push_back(world[x][y].females.at(f));	//append it to the target patch newfemales vector
						world[x][y].females.erase(world[x][y].females.begin()+f);	//remove the individual from natal patch
						f--;	//the next individual in the vector is now in the same position as the removed individual
						count_dispersers++;	//count number of individuals dispersing
					}
				}
				else	//remove the individual if it has a non-viable GRN
				{
					world[x][y].females.erase(world[x][y].females.begin()+f);
					f--;
					grn_mortality_count++;
				}
				
			}
			//same as above for males
			for(int m=0;m<world[x][y].males.size();m++)
			{
			double dispersal_probability;	//stores individual dispersal probability
				double inp[2];	//stores input to the GRN
				inp[0]=normalised_density;	//the GRN takes the local normalised population density as an input
				inp[1]=1.0;
				dispersal_probability=grn_output(inp,
					world[x][y].males.at(m).input_weights,
					world[x][y].males.at(m).regulatory_weights,world[x][y].males.at(m).output_weights,
					world[x][y].males.at(m).regulatory_threshold,world[x][y].males.at(m).regulatory_slope);
				if(int(dispersal_probability)!= 42)
				{
					if(ran()< dispersal_probability)
					{
						std::vector<int> coor=decide_patch(x,y,t);
						if(ran()>DISP_MORT)
							world[coor.at(0)][coor.at(1)].newmales.push_back(world[x][y].males.at(m));
						world[x][y].males.erase(world[x][y].males.begin()+m);
						m--;
						count_dispersers++;
					}
				}
				else
				{
					world[x][y].males.erase(world[x][y].males.begin()+m);
					m--;
					grn_mortality_count++;
				}
				
			}
			
			world[x][y].measured_dispersal=double(count_dispersers)/population_density;
			world[x][y].grn_mortality =double(grn_mortality_count)/population_density;


		}
		
	}
	for(int x=0;x<world_size_x;x++)	//go through all the patches
	{
		for(int y=0;y<world_size_y;y++)
		{
			if(world[x][y].newfemales.size()>0)	 //if there are dispersers that are arriving in the patch
			{
				for(int f=0;f<world[x][y].newfemales.size();f++)
				{
					world[x][y].females.push_back(world[x][y].newfemales.at(f));	//add these dispersers to the new patch
				}
			}
			if(world[x][y].newmales.size()>0) //same for males
			{
				for(int m=0;m<world[x][y].newmales.size();m++)
				{
					world[x][y].males.push_back(world[x][y].newmales.at(m));
				}
			}
			world[x][y].newfemales.clear();	//clear the newfemales and newmales vectors
			world[x][y].newmales.clear();

		}
	}

}


//---------------------------------------------------------------------reproduction procedure
void reproduce(int t) //reproduction loop
{ 
	for(int x=0; x<world_size_x;x++)	//go through all the patches
	{
		for(int y=0; y<world_size_y;y++)
		{
			world[x][y].newfemales.clear();	//clear the newfemales vector just in case
			world[x][y].newmales.clear();
			if(world[x][y].males.size() !=0 && world[x][y].females.size() != 0)	//if there are individuals present in the patch
			{
				int Nf=world[x][y].females.size();	//calculate number of females
				int Nm=world[x][y].males.size();	//calculate number of males
				double alpha_net=0;	//stores the amount of density regulation
				for(int f=0; f<Nf;f++)      //calculating net alpha for males and females
				{
					alpha_net=alpha_net+ALPHA;
				}
				for(int m=0; m<Nm;m++)
				{
					alpha_net=alpha_net+ALPHA;
				}
				for(int f=0;f<Nf;f++)	//go through all the females in a patch
				{
					int mate_position=floor(ran()*world[x][y].males.size());	//randomly choose a male to mate with
						double mean_lambda= LAMBDA;// setting mean as trait value
						double mean_offspring = 2*mean_lambda * densReg(alpha_net);//calculate mean number of offspring from the Beverton-Holt model
						int no_of_babies= poisson(mean_offspring); //number of offspring Poisson distributed
						for(int b=0;b<no_of_babies;b++)         //setting parental characteristics
						{
							TInd newind;	//create the new offspring
							//inheritance of parental traits
							for(int i=0;i<REGULATORY_SIZE;i++)
							{
								for(int j=0;j<INPUT_SIZE;j++)
								{
									if(ran()<0.5)
										newind.input_weights[0][j][i]=mutate(world[x][y].females.at(f).input_weights[0][j][i],t);
									else newind.input_weights[0][j][i]=mutate(world[x][y].females.at(f).input_weights[1][j][i],t);
									if(ran()<0.5)
										newind.input_weights[1][j][i]=mutate(world[x][y].males.at(mate_position).input_weights[0][j][i],t);
									else newind.input_weights[1][j][i]=mutate(world[x][y].males.at(mate_position).input_weights[1][j][i],t); 
								}
								for(int j=0;j<REGULATORY_SIZE;j++)
								{
									if(ran()<0.5)
										newind.regulatory_weights[0][i][j]=mutate(world[x][y].females.at(f).regulatory_weights[0][i][j],t);
									else newind.regulatory_weights[0][i][j]=mutate(world[x][y].females.at(f).regulatory_weights[1][i][j],t);
									if(ran()<0.5)
										newind.regulatory_weights[1][i][j]=mutate(world[x][y].males.at(mate_position).regulatory_weights[0][i][j],t);
									else newind.regulatory_weights[1][i][j]=mutate(world[x][y].males.at(mate_position).regulatory_weights[1][i][j],t);
								}
								if(ran()<0.5)
									newind.output_weights[0][i]=mutate(world[x][y].females.at(f).output_weights[0][i],t);
								else newind.output_weights[0][i]=mutate(world[x][y].females.at(f).output_weights[1][i],t);
								if(ran()<0.5)
									newind.output_weights[1][i]=mutate(world[x][y].males.at(mate_position).output_weights[0][i],t);
								else newind.output_weights[1][i]=mutate(world[x][y].males.at(mate_position).output_weights[1][i],t); 
								if(ran()<0.5)
									newind.regulatory_threshold[0][i]=mutate(world[x][y].females.at(f).regulatory_threshold[0][i],t);
								else newind.regulatory_threshold[0][i]=mutate(world[x][y].females.at(f).regulatory_threshold[1][i],t);
								if(ran()<0.5)
									newind.regulatory_threshold[1][i]=mutate(world[x][y].males.at(mate_position).regulatory_threshold[0][i],t);
								else newind.regulatory_threshold[1][i]=mutate(world[x][y].males.at(mate_position).regulatory_threshold[1][i],t); 
								if(ran()<0.5)
									newind.regulatory_slope[0][i]=mutate(world[x][y].females.at(f).regulatory_slope[0][i],t);
								else newind.regulatory_slope[0][i]=mutate(world[x][y].females.at(f).regulatory_slope[1][i],t);
								if(ran()<0.5)
									newind.regulatory_slope[1][i]=mutate(world[x][y].males.at(mate_position).regulatory_slope[0][i],t);
								else newind.regulatory_slope[1][i]=mutate(world[x][y].males.at(mate_position).regulatory_slope[1][i],t); 
								
							}

							if(ran()<0.5)	//offspring has equal chance of being male or female
								world[x][y].newfemales.push_back(newind);	//add newborn females
							else world[x][y].newmales.push_back(newind);	//add newborn males

						}
					}

				}
				else 
				{
					world[x][y].females.clear();	//clear the patches 
					world[x][y].males.clear();
				}

			}
		}
		
	}





void death()	//death procedure
{
	for(int x=0;x<world_size_x;x++)	//go through all the patches
	{
		for(int y=0;y<world_size_y;y++)
		{
			world[x][y].females.clear();	//clear the previous generation males and females
			world[x][y].males.clear();
			world[x][y].females=world[x][y].newfemales;	//put the offspring in the males and females vector
			world[x][y].males=world[x][y].newmales;
			world[x][y].newfemales.clear();	//clear the newfemales vector
			world[x][y].newmales.clear();
		}
	}
}

void patch_extinction()	//random patch extinction
{
	for(int x=0;x<world_size_x;x++)	//go through all the patches
	{
		for(int y=0; y<world_size_y;y++)
		{
			if(ran()<EXTINCTION_PROB)	//the patch is cleared with a probability EXTINCTION_PROB
			{
				world[x][y].females.clear();
				world[x][y].males.clear();
			}
		}
	}
}


//------------------------------------------------------------------------------Main function
int main()
	{
	//output the population size, sex ratio, measured dispersal at each patch
		ofstream op;
		op.open("output.txt");
		ofstream op1;
		op1.open("genotype_properties.txt");	//output sensitivity to mutation
		ofstream op2;
		op2.open("phenotype_properties_core.txt");	//output phenotype in the range core
		ofstream op3;
		op3.open("phenotype_properties_front_left.txt");
		ofstream op4;
		op4.open("phenotype_properties_front_right.txt");
		ofstream op5;
		op5.open("output1.txt");	//output genotypes at the end of burn-in and at the end of range expansion
		specify_rng(RS);
		set_parameters();	//read and set the model parameters
		//headers to the output files
		op5 <<"rep"<<" "<<"t"<<" "<<"x"<<" "<<"y ";
		for(int i=0;i<REGULATORY_SIZE;i++)
		{
			for(int j=0;j<INPUT_SIZE;j++)
			{
				op5<<"iw_"<<j+1<<i+1<<" ";
			}
			for(int j=0;j<REGULATORY_SIZE;j++)
			{
				op5<<"rw_"<<i+1<<j+1<<" ";
			}
			op5<<"rt_"<<i+1<<" ";
			op5<<"ow_"<<i+1<<" ";
			op5<<"rs_"<<i+1<<" ";
		}
		op5<<endl;
		op <<"rep"<<" "<<"t"<<" "<<"x"<<" "<<"y"<<" "<<"N"<<" "<<"sex_ratio"<<" "<<"disp_rate grn_mortality ";
		/*for(int i=0;i<REGULATORY_SIZE;i++)
		{
			for(int j=0;j<INPUT_SIZE;j++)
			{
				op<<"iw_"<<j+1<<i+1<<" ";
			}
			for(int j=0;j<REGULATORY_SIZE;j++)
			{
				op<<"rw_"<<i+1<<j+1<<" ";
			}
			op<<"rt_"<<i+1<<" ";
			op<<"ow_"<<i+1<<" ";
			op<<"rs_"<<i+1<<" ";
		}*/
		op<<"margin_x_left margin_x_right";
		op<<endl;
		op1<<"rep t sensitivity_core sensitivity_front_left sensitivity_front_right"<<endl;
		op2<<"rep t ";
		op3<<"rep t ";
		op4<<"rep t ";

		int max_density=15;
		for(int d=0;d<max_density;d++)
		{
			op2<<"d_f"<<d<<"_1 "<<"d_f"<<d<<"_2 "<<"d_f"<<d<<"_3 ";
		}
		
		for(int d=0;d<max_density;d++)
		{
			op3<<"d_f"<<d<<"_1 "<<"d_f"<<d<<"_2 "<<"d_f"<<d<<"_3 ";
		}

		for(int d=0;d<max_density;d++)
		{
			op4<<"d_f"<<d<<"_1 "<<"d_f"<<d<<"_2 "<<"d_f"<<d<<"_3 ";
		}
		for(int d=0;d<max_density;d++)
		{
			op2<<"d_m"<<d<<"_1 "<<"d_m"<<d<<"_2 "<<"d_m"<<d<<"_3 ";
		}
		op2<<endl;
		for(int d=0;d<max_density;d++)
		{
			op3<<"d_m"<<d<<"_1 "<<"d_m"<<d<<"_2 "<<"d_m"<<d<<"_3 ";
		}
		op3<<endl;
		for(int d=0;d<max_density;d++)
		{
			op4<<"d_m"<<d<<"_1 "<<"d_m"<<d<<"_2 "<<"d_m"<<d<<"_3 ";
		}


		op4<<endl;
	for(int r=0; r<REPLICATES; r++)     //replicates
	{
		init_world();	//initialise the landscape
		int t=0;	//set time count to 0
		int margin_x_left=world_size_x/2-burn_in_x/2;	//set initial range front position
		int margin_x_right=world_size_x/2+burn_in_x/2-1;
		do{	//loop throug time
			//output
			for(int x=0; x<world_size_x;x++)
			{
				for(int y=0; y<world_size_y ;y++)
				{

				 if(x>=world_size_x/2-burn_in_x/2 && x<=world_size_x/2+burn_in_x/2-1 && t>BURN_IN_TIME-50 )//output the central 10x5 patches
				 	output_metapopulation(op,x,y,r,t,margin_x_left,margin_x_right);
				 if(x<margin_x_left && world[x][y].females.size()+world[x][y].males.size()>0)
				 	margin_x_left=x;	//calculate position of range front every time step
				 if(x>margin_x_right && world[x][y].females.size()+world[x][y].males.size()>0)
				 	margin_x_right=x;
				} 
			}
			if(t>BURN_IN_TIME)	//output range front
			{
				for(int x=margin_x_left; x<margin_x_left+5;x++)
				{
					for(int y=0; y<world_size_y ;y++)
					{

						output_metapopulation(op,x,y,r,t,margin_x_left,margin_x_right);

					}
				} 

				for(int x=margin_x_right-4; x<margin_x_right;x++)
				{
					for(int y=0; y<world_size_y ;y++)
					{

						output_metapopulation(op,x,y,r,t,margin_x_left,margin_x_right);

					}
				}
			}

			if(t<=BURN_IN_TIME)	//output sensitivity to mutation and dispersal trait quartiles
			{
				if(t%2500==0)	//every 2500 time steps
				{
					vector<double> rob;	//stores the sensitivity to mutation
					rob=sensitivity_to_mutation(world_size_x/2-burn_in_x/2,world_size_x/2+burn_in_x/2);
					op1<<r<<" "<<t<<" ";
					op1<<rob.at(1)<<" NA NA ";
					//op1<<heterozygosity(world_size_x/2-burn_in_x/2,world_size_x/2+burn_in_x/2)<<" NA NA"<<endl;
					op2<<r<<" "<<t<<" ";
					vector<vector<double> > quartile=trait_calulator(world_size_x/2-burn_in_x/2, world_size_x/2+burn_in_x/2,0.0);	//storefemale reaction norm median and quartiles
					int max_density=15;	//output for 15 population densities
					for(int d=0;d<max_density;d++)
					{
						op2<<quartile.at(d).at(0)<<" "<<quartile.at(d).at(1)<<" "<<quartile.at(d).at(2)<<" ";
					}
					quartile=trait_calulator(world_size_x/2-burn_in_x/2, world_size_x/2+burn_in_x/2,1.0);	//store male reaction norms
					for(int d=0;d<max_density;d++)
					{
						op2<<quartile.at(d).at(0)<<" "<<quartile.at(d).at(1)<<" "<<quartile.at(d).at(2)<<" ";
					}
					op2<<endl;


				}
			}
			if(t>BURN_IN_TIME && t%10==0)	//output sensitivity to mutation and trait quartiles every 10 time steps
			{
				vector<double> rob;	//stores sensitivity to mutation
				rob=sensitivity_to_mutation(world_size_x/2-burn_in_x/2,world_size_x/2+burn_in_x/2);
				op1<<r<<" "<<t<<" ";
				op1<<rob.at(1)<<" ";
				rob=sensitivity_to_mutation(margin_x_left,margin_x_left+5);

				op1<<rob.at(1)<<" ";
				rob=sensitivity_to_mutation(margin_x_right-4,margin_x_right+1);
				op1<<rob.at(1)<<" ";
				op2<<r<<" "<<t<<" ";
				op3<<r<<" "<<t<<" ";
				op4<<r<<" "<<t<<" ";
				vector<vector<double> > quartile=trait_calulator(world_size_x/2-burn_in_x/2, world_size_x/2+burn_in_x/2,0.0);	//store female reaction norm for different population densities
				int max_density=15;	//number of population densities tested
				for(int d=0;d<max_density;d++)	//go through all the population densities
				{
					op2<<quartile.at(d).at(0)<<" "<<quartile.at(d).at(1)<<" "<<quartile.at(d).at(2)<<" ";	//output female reaction norm
				}

				quartile=trait_calulator(margin_x_left,margin_x_left+5,0.0);
				for(int d=0;d<max_density;d++)
				{
					op3<<quartile.at(d).at(0)<<" "<<quartile.at(d).at(1)<<" "<<quartile.at(d).at(2)<<" ";
				}

				quartile=trait_calulator(margin_x_right-4,margin_x_right+1,0.0);
				for(int d=0;d<max_density;d++)
				{
					op4<<quartile.at(d).at(0)<<" "<<quartile.at(d).at(1)<<" "<<quartile.at(d).at(2)<<" ";
				}
				quartile=trait_calulator(world_size_x/2-burn_in_x/2, world_size_x/2+burn_in_x/2,1.0);//output male reaction norm
				for(int d=0;d<max_density;d++)
				{
					op2<<quartile.at(d).at(0)<<" "<<quartile.at(d).at(1)<<" "<<quartile.at(d).at(2)<<" ";
				}
				op2<<endl;
				quartile=trait_calulator(margin_x_left,margin_x_left+5,1.0);
				for(int d=0;d<max_density;d++)
				{
					op3<<quartile.at(d).at(0)<<" "<<quartile.at(d).at(1)<<" "<<quartile.at(d).at(2)<<" ";
				}
				op3<<endl;
				quartile=trait_calulator(margin_x_right-4,margin_x_right+1,1.0);
				for(int d=0;d<max_density;d++)
				{
					op4<<quartile.at(d).at(0)<<" "<<quartile.at(d).at(1)<<" "<<quartile.at(d).at(2)<<" ";
				}
				op4<<endl;
			}

			if(t==BURN_IN_TIME-1)
				output_genotype(op5,world_size_x/2-burn_in_x/2,world_size_x/2+burn_in_x/2,r,t);

//output the genotypes at the range front at the end of range expansion
			if(margin_x_right==world_size_x-1 && margin_x_left==0){
				output_genotype(op5,margin_x_right-4,margin_x_right+1,r,t);
				output_genotype(op5,margin_x_left,margin_x_left+5,r,t);
			}

			//life cycle
			disperse(t);
			reproduce(t);
			death();
			patch_extinction();
			t++;
			if(t>BURN_IN_TIME+10000)
				break;

		}
		while(margin_x_right!=world_size_x-1 || margin_x_left!=0);
	}
	op.close();
	op1.close();
	op2.close();
	op3.close();
	op4.close();
	op5.close();
	return 0;
}







