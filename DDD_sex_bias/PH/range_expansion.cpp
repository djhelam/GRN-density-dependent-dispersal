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


#include <cstdlib>    //standard C library
#include <iostream>   //standard input output function library
#include <fstream>    //file stream library
#include <numeric>
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>    //random number generator gsl library      
#include <gsl/gsl_randist.h> //gsl random distribution library
#include <algorithm>
#include <math.h>
using namespace std;

const int NO_OF_LOCI=1; //number of loci

//________________________________________________________________________________________
//----------------------------------------------------------Class defining the individuals
class TInd{       
public:
  TInd();
  double lambda[2]; //growth rate
  double alpha[2];  //competition coefficient
  double C_thresh_female[2][NO_OF_LOCI]; //threshold from Poethke and Hovestdtat (2002) if female
  double C_thresh_male[2][NO_OF_LOCI]; //threshold from Poethke and Hovestdtat (2002) if male
};

//----------------------------------------------------------Constructor for the class TInd
TInd::TInd(){     
  lambda[0]=0;      
  lambda[1]=0;        
  alpha[0]=0;
  alpha[1]=0;
  for(int i=0;i<NO_OF_LOCI;i++)
  {
    C_thresh_female[0][i]=0;
    C_thresh_female[1][i]=0;
    C_thresh_male[0][i]=0;
    C_thresh_male[1][i]=0;
  }

}

//________________________________________________________________________________________
//--------------------------------------------------------------Class defining the patches
class TPatch      
{
public:
  TPatch();
  vector <TInd> females;    //females in the simulation
  vector <TInd> newfemales;   //vector to store new disperses of new born females
  vector <TInd> males;        //males in the simulation
  vector <TInd> newmales; //vector stores newborn males or dispersers
  double measured_dispersal;    //fraction of individuals dispersing
};

//--------------------------------------------------------Constructor for the class TPatch
TPatch::TPatch(){    
  females.clear();
  males.clear();
  newfemales.clear();
  newmales.clear();
  measured_dispersal=0;
}

//________________________________________________________________________________________
//------------------------------------------------------------------------Model Parameters 

int No;         //Initial number of individuals
double LAMBDA;      //Intrinsic growth rate---Beverton-Holt model
double ALPHA;     //Intra-specific competition coefficient
double DISPERSAL_PROB;  // Dispersal probability if it does not evolve
int BURN_IN_TIME;       //Number of time steps before beginning range expansion
int REPLICATES;         //Number of replicates
double DISP_MORT; //dispersal mortality
double EXTINCTION_PROB;//probability of local patch extinction
double VARIANCE; //standard deviation of mutation effects
double MUT_RATE; //mutation rate for dispersal trait

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

//----------------------------------------------------function returns 1,2, and 3 quartile
vector<double> quartiles(vector<double> trait)
{
  sort(trait.begin(), trait.end());//sort vector
  double q1;  //stores quartile
  double q2;  //stores quartile
  double q3;  //stores quartile
  int n=trait.size();  //stores size of vector
  q2=trait.at((n+1)/2)+double((n+1)%2)*(trait.at(((n+1)/2)+1)-trait.at((n+1)/2))/double(2); //median
  q1=trait.at((n+1)/4)+double((n+1)%4)*(trait.at(((n+1)/4)+1)-trait.at((n+1)/4))/double(4); //25%
  q3=trait.at(3*(n+1)/4)+double(3*(n+1)%4)*(trait.at((3*(n+1)/4)+1)-trait.at(3*(n+1)/4))/double(4);//75%
  vector<double> quartile;
  quartile.push_back(q1);
  quartile.push_back(q2);
  quartile.push_back(q3);
  return quartile;

}


//----------------------------------------------------function returns mean of two numbers
double mean(double a, double b) 
{
  return ((a+b)/double(2));
}


//------------------------------------------function returns mutated dispersal probability
double mutate(double d) 
{
  double new_trait;
  if(ran()<MUT_RATE)
    new_trait= d+ gauss(VARIANCE);  
  else new_trait=d;
  return new_trait;

}


//________________________________________________________________________________________
//--------------------------------------------------------Inputting and setting parameters
void set_parameters()  
{
  string para[10];  //stores  model parameters
  string line;
  ifstream myfile ("input.txt");  //open parameter input file
  int count=0;
  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      if(count%2==1)
      {
        para[count/2]=line; //store every only parameter values and not their names
      }
      count++;

    }
    myfile.close();
  }
  else cout << "Unable to open file";
  No = (int) std::atof(para[0].c_str());  //store initial population density per patch
  LAMBDA=std::atof(para[1].c_str());      //store intrinsic growth rate---Beverton-Holt model
  ALPHA=std::atof(para[2].c_str());       //store intra-specific competition coefficient
  BURN_IN_TIME= (int) std::atof(para[3].c_str()); //store time before range expansions begin
  REPLICATES=(int) std::atof(para[4].c_str());  //store number of replicates
  DISPERSAL_PROB=std::atof(para[5].c_str());  //store dispersal probability 
  DISP_MORT=std::atof(para[6].c_str()); //store dispersal mortality
  EXTINCTION_PROB=std::atof(para[7].c_str()); //store probability of random patch extinction
  VARIANCE=std::atof(para[8].c_str());  //store standard deviation of mutation effects
  MUT_RATE=std::atof(para[9].c_str());  //store mutation rate
}


//________________________________________________________________________________________
//----------------------------------------------------------------------Initialising model
void init_world()
{
  for(int x=0;x<world_size_x;x++) //go through all patches
  {
    for(int y=0;y<world_size_y ;y++)
    {
    //clear all individuals in a patch
      world[x][y].females.clear();  //clear all individuals in the patch
      world[x][y].males.clear();
      world[x][y].newfemales.clear(); 
      world[x][y].newmales.clear();
    }
  }
 //seed population only in the central 10 patches in x direction- 245 to 254
  for(int x=world_size_x/2-burn_in_x/2;x<world_size_x/2+burn_in_x/2;x++)
  {
    for(int y=0;y<world_size_y ;y++)
    {
    //Initialise No individuals in each patch
      for(int n=0; n<No; n++)
      {
          TInd newind;      //stores new individuals
          newind.lambda[0]=LAMBDA;  
          newind.lambda[1]=LAMBDA;
          newind.alpha[0]=ALPHA;
          newind.alpha[1]=ALPHA;
          for(int i=0;i<NO_OF_LOCI;i++)
          {
            newind.C_thresh_female[0][i]=ran();  //initialise genotype with standing genetic variation
            newind.C_thresh_female[1][i]=ran();
            newind.C_thresh_male[0][i]=ran();  //initialise genotype with standing genetic variation
            newind.C_thresh_male[1][i]=ran();
          }
    
          if(ran()<0.5) //individual has equal chance of being male or female
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
    for(int x=x_min; x<x_max;x++) //output half the genotypes in a specified region
    {
      for(int y=0; y<world_size_y ;y++)
      {
       for(int f=0;f<world[x][y].females.size();f++)
       {
        if(ran()<0.5){
          op5 <<r<<" "<<t<<" "<<x<<" "<<y<<" " ;
          for(int i=0;i<NO_OF_LOCI;i++)
          {
            op5<<mean(world[x][y].females.at(f).C_thresh_female[0][i],world[x][y].females.at(f).C_thresh_female[1][i])<<" "<<mean(world[x][y].females.at(f).C_thresh_male[0][i],world[x][y].females.at(f).C_thresh_male[1][i])<<" ";
          }
          op5<<endl;
        }
      }
      for(int m=0;m<world[x][y].males.size();m++)
      {
        if(ran()<0.5){
          op5 <<r<<" "<<t<<" "<<x<<" "<<y<<" " ;
          for(int i=0;i<NO_OF_LOCI;i++)
          {
            op5<<mean(world[x][y].males.at(m).C_thresh_female[0][i],world[x][y].males.at(m).C_thresh_female[1][i])<<" "<<mean(world[x][y].males.at(m).C_thresh_male[0][i],world[x][y].males.at(m).C_thresh_male[1][i])<<" ";
          }
          op5<<endl;
        }
      }

    }
  }

}

//________________________________________________________________________________________
//----------------------------------------------------------Output population density.etc

void output_metapopulation(ofstream& op, int x,int  y, int r, int t, int margin_x_left, int margin_x_right)
{
  op <<r<<" "<<t<<" "<<x<<" "<<y<<" "<<world[x][y].females.size()+world[x][y].males.size()<<" "<<
  double(world[x][y].females.size())/double(world[x][y].females.size()+world[x][y].males.size())<<" "<<
  world[x][y].measured_dispersal<<" "<<margin_x_left<<" "<<margin_x_right<<" ";
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
    if (newx<world_size_x/2-burn_in_x/2)  //if the new patch is beyond the range core in the left direction
      newx=world_size_x/2-burn_in_x/2+1;  //send the individual to one patch to the right
    if(newx>world_size_x/2+burn_in_x/2-1)
      newx=world_size_x/2+burn_in_x/2-2;
  }
  if(t>=BURN_IN_TIME)  //when range expansions begin
  {
    if (newx<0) //no conditions, does not matter at the landscape end as the simulations stop
      newx=0;
    if(newx==world_size_x)
      newx=world_size_x-1;
  }
  //torus in y direction
  if(newy<0)  //if the individual is at the bottom
    newy=world_size_y - 1;  //send it to the top
  if(newy==world_size_y)
    newy=0;
  coordinates.push_back(newx);  //store the new coordinates of the individual
  coordinates.push_back(newy);
  return coordinates;
}

//________________________________________________________________________________________
//----------------------------function returns the dispersal rate calculated from genotype

double  dispersal_probability_calc(double c[2][NO_OF_LOCI], double density)
{
  double c_mean=0;  //store mean threshold population density
  double disp;  //store dispersal probability
    for(int i=0;i<NO_OF_LOCI;i++) //go through the loci
    {
      c_mean=c_mean+mean(c[0][i],c[1][i]);  
    }
    c_mean=c_mean/double(NO_OF_LOCI); //threshold calculated from genotype
    if(density<c_mean)  //if population density is below the threshold
      disp=0; //no dispersal
    else 
      disp=1.0-(c_mean/density);  //otherwise disperse with specified probability
    return disp;  //return dispersal probability
}





//________________________________________________________________________________________
//-------------------------------------------------------------Density regulation function

  double densReg(double a) 
  {
    return(1 /(1+(a)));
  }




vector< vector <double> > trait_calulator(int x_min,int x_max,int sex) //takes the boundaries for which 
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
  vector<TInd> random_individuals;
  for(int ind=0;ind<no_of_individuals;ind++)//draw no_of_individuals individuals
  { 
    int position_random=floor(ran()*all_individuals.size());//individuals drawn at random
    random_individuals.push_back(all_individuals.at(position_random));
  }
  all_individuals.clear();
  vector<vector<double> > quartile;
    vector<double> disp_all; // stores all dispersal probabilities corresponding to that density
    for(int p=0;p<no_of_individuals;p++){
      double disp=0;
      for(int i=0; i<NO_OF_LOCI;i++)
      {
        if(sex==0)
          disp=disp+mean(random_individuals.at(p).C_thresh_female[0][i],random_individuals.at(p).C_thresh_female[1][i]);
        if(sex==1)
           disp=disp+mean(random_individuals.at(p).C_thresh_male[0][i],random_individuals.at(p).C_thresh_male[1][i]);
      }
      disp=disp/NO_OF_LOCI;
      disp_all.push_back(disp);
    }
    quartile.push_back(quartiles(disp_all));
    return quartile;
  }



//________________________________________________________________________________________
//-------------------------------------------------------------------Life cycle procedures


//---------------------------------------------------------------------Dispersal procedure
  void disperse(int t)
  {
  for(int x=0; x<world_size_x;x++)    //clears newmales and newfemales vectors
  {
    for(int y=0; y<world_size_y ;y++)
    {
      world[x][y].newfemales.clear();
      world[x][y].newmales.clear();
    }
  }
  for(int x=0;x<world_size_x;x++)   //dispersal loop
  {
    for(int y=0;y<world_size_y ;y++)  //go through all the patches
    {
      int count_dispersers=0; //track number of individuals dispersing
      //inputs of the GRN stored in inp
      double population_density=double(world[x][y].males.size()+world[x][y].females.size());  //stores population density
      
      for(int f=0;f<world[x][y].females.size();f++)
      {
        //dispersal probability calculation according to function in Poethke and Hovesdtat (2002)
        double dispersal_probability;
        dispersal_probability=dispersal_probability_calc(world[x][y].females.at(f).C_thresh_female,population_density*ALPHA/(LAMBDA-1)); 
        //dispersal
        if(ran()< dispersal_probability)  //if the individual emigrates
        {
          std::vector<int> coor=decide_patch(x,y,t);  //choose target patch
            //dispersal mortalilty
          if(ran()>DISP_MORT) //if the individual survives dispersal
            world[coor.at(0)][coor.at(1)].newfemales.push_back(world[x][y].females.at(f));  //send it to target patch
          world[x][y].females.erase(world[x][y].females.begin()+f); //remove the individual from natal patch
          f--;  //the next individual is now at its position in the vector
          count_dispersers++; //count number of individuals dispersing
        }
        
      }
      for(int m=0;m<world[x][y].males.size();m++) //same as above for maes
      {
        double dispersal_probability;
        dispersal_probability=dispersal_probability_calc(world[x][y].males.at(m).C_thresh_male,population_density*ALPHA/(LAMBDA-1));
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

      
      world[x][y].measured_dispersal=double(count_dispersers)/population_density;

    }
    
  }
for(int x=0;x<world_size_x;x++) //go through all the patches
  {
    for(int y=0;y<world_size_y;y++)
    {
      if(world[x][y].newfemales.size()>0)  //if there are dispersers that are arriving in the patch
      {
        for(int f=0;f<world[x][y].newfemales.size();f++)
        {
          world[x][y].females.push_back(world[x][y].newfemales.at(f));  //add these dispersers to the new patch
        }
      }
      if(world[x][y].newmales.size()>0) //same for males
      {
        for(int m=0;m<world[x][y].newmales.size();m++)
        {
          world[x][y].males.push_back(world[x][y].newmales.at(m));
        }
      }
      world[x][y].newfemales.clear(); //clear the newfemales and newmales vectors
      world[x][y].newmales.clear();

    }
  }

}

//---------------------------------------------------------------------reproduction procedure
void reproduce() //reproduction loop
{ 
  for(int x=0; x<world_size_x;x++)  //go through all the patches
  {
    for(int y=0; y<world_size_y ;y++)
    {
      world[x][y].newfemales.clear(); //clear newfemales and newmales vector
      world[x][y].newmales.clear();
      if(world[x][y].males.size() !=0 && world[x][y].females.size() != 0) //if there are individuals 
      {
        int Nf=world[x][y].females.size();  //calculate female population size
        int Nm=world[x][y].males.size();  //calculate male population size
        double alpha_net=0; //stores net effect of population density 
        for(int f=0; f<Nf;f++)      //calculating net alpha for males and females
        {
          alpha_net=alpha_net+mean(world[x][y].females.at(f).alpha[0],world[x][y].females.at(f).alpha[1]);
        }
        for(int m=0; m<Nm;m++)
        {
          alpha_net=alpha_net+mean(world[x][y].males.at(m).alpha[0],world[x][y].males.at(m).alpha[1]);
        }
        for(int f=0;f<Nf;f++) //go through all the females
        {
          int mate_position=floor(ran()*world[x][y].males.size());  //choose a male to mate with
            double mean_lambda= mean(world[x][y].females.at(f).lambda[0],world[x][y].females.at(f).lambda[1]);// setting mean as trait value
            double mean_offspring = 2*mean_lambda * densReg(alpha_net);//beverton holt model
            int no_of_babies= poisson(mean_offspring); //number of offspring poisson distributed
            for(int b=0;b<no_of_babies;b++)         //setting parental characteristics
            {
              TInd newind;
              if(ran()<0.5)
                newind.lambda[0]=world[x][y].females.at(f).lambda[0];
              else newind.lambda[0]=world[x][y].females.at(f).lambda[1];
              if(ran()<0.5)
                newind.lambda[1]=world[x][y].males.at(mate_position).lambda[0];
              else newind.lambda[1]=world[x][y].males.at(mate_position).lambda[1];
              if(ran()<0.5)
                newind.alpha[0]=world[x][y].females.at(f).alpha[0];
              else newind.alpha[0]=world[x][y].females.at(f).alpha[1];
              if(ran()<0.5)
                newind.alpha[1]=world[x][y].males.at(mate_position).alpha[0];
              else newind.alpha[1]=world[x][y].males.at(mate_position).alpha[1];
              for(int i=0;i<NO_OF_LOCI;i++) //inheritance of threshold
              {
                if(ran()<0.5)
                  newind.C_thresh_female[0][i]=mutate(world[x][y].females.at(f).C_thresh_female[0][i]);
                else newind.C_thresh_female[0][i]=mutate(world[x][y].females.at(f).C_thresh_female[1][i]);
                if(ran()<0.5)
                  newind.C_thresh_female[1][i]=mutate(world[x][y].males.at(mate_position).C_thresh_female[0][i]);
                else newind.C_thresh_female[1][i]=mutate(world[x][y].males.at(mate_position).C_thresh_female[1][i]);
              }
              for(int i=0;i<NO_OF_LOCI;i++) //inheritance of threshold
              {
                if(ran()<0.5)
                  newind.C_thresh_male[0][i]=mutate(world[x][y].females.at(f).C_thresh_male[0][i]);
                else newind.C_thresh_male[0][i]=mutate(world[x][y].females.at(f).C_thresh_male[1][i]);
                if(ran()<0.5)
                  newind.C_thresh_male[1][i]=mutate(world[x][y].males.at(mate_position).C_thresh_male[0][i]);
                else newind.C_thresh_male[1][i]=mutate(world[x][y].males.at(mate_position).C_thresh_male[1][i]);
              }
              if(ran()<0.5)
                world[x][y].newfemales.push_back(newind);
              else world[x][y].newmales.push_back(newind);
            }
          }
        }
        else 
        {
          world[x][y].females.clear();
          world[x][y].males.clear();
        }

      }

    }

  }

  void death()
  {
    for(int x=0;x<world_size_x;x++)
    {
      for(int y=0;y<world_size_y ;y++)
      {
        world[x][y].females.clear();
        world[x][y].males.clear();
        world[x][y].females=world[x][y].newfemales;
        world[x][y].males=world[x][y].newmales;
        world[x][y].newfemales.clear();
        world[x][y].newmales.clear();
      }
    }
  }

  void patch_extinction()
  {
    for(int x=0;x<world_size_x;x++)
    {
      for(int y=0; y<world_size_y ;y++)
      {
        if(ran()<EXTINCTION_PROB)
        {
          world[x][y].females.clear();
          world[x][y].males.clear();
        }
      }
    }
  }


 int main()
  {
  //output the population size, sex ratio, measured dispersal at each patch
    ofstream op;
    op.open("output.txt");
    //ofstream op1;
    //op1.open("genotype_properties.txt");
    ofstream op2;
    op2.open("phenotype_properties_core.txt");
    ofstream op3;
    op3.open("phenotype_properties_front_left.txt");
    ofstream op4;
    op4.open("phenotype_properties_front_right.txt");
    ofstream op5;
    op5.open("output1.txt");

    //seed 
    specify_rng(RS);

    //read parameters from input file
    set_parameters();

    //column names of output files
    op <<"rep"<<" "<<"t"<<" "<<"x"<<" "<<"y"<<" "<<"N"<<" "<<"sex_ratio"<<" "<<"disp_rate margin_x_left margin_x_right ";

    op<<endl;

    op5 <<"rep"<<" "<<"t"<<" "<<"x"<<" "<<"y"<<" ";
    for(int i=0;i<NO_OF_LOCI;i++)
    {
      op5<<"d_f"<<i<<" "<<"d_m"<<i<<" ";
    }
    op5<<endl;

    //op1<<"rep t hetero_core hetero_front_left hetero_front_right"<<endl;

    op2<<"rep t ";
    op3<<"rep t ";
    op4<<"rep t ";

    op2<<"d_f"<<0<<"_1 "<<"d_f"<<0<<"_2 "<<"d_f"<<0<<"_3 "<<"d_m"<<0<<"_1 "<<"d_m"<<0<<"_2 "<<"d_m"<<0<<"_3 ";

    op2<<endl;

    op3<<"d_f"<<0<<"_1 "<<"d_f"<<0<<"_2 "<<"d_f"<<0<<"_3 "<<"d_m"<<0<<"_1 "<<"d_m"<<0<<"_2 "<<"d_m"<<0<<"_3 ";

    op3<<endl;
    op4<<"d_f"<<0<<"_1 "<<"d_f"<<0<<"_2 "<<"d_f"<<0<<"_3 "<<"d_m"<<0<<"_1 "<<"d_m"<<0<<"_2 "<<"d_m"<<0<<"_3 ";
    op4<<endl;
    
  for(int r=0; r<REPLICATES; r++)     //replicates
  {
    init_world(); //initialise landscape
    int t=0;
    int margin_x_left=world_size_x/2-burn_in_x/2;
    int margin_x_right=world_size_x/2+burn_in_x/2-1;;
    do{
      //output
      for(int x=0; x<world_size_x;x++)
      {
        for(int y=0; y<world_size_y ;y++)
        {

         if(x>=world_size_x/2-burn_in_x/2 && x<=world_size_x/2+burn_in_x/2-1 && t>BURN_IN_TIME-50 )//output the central 10x5 patches
          output_metapopulation(op,x,y,r,t,margin_x_left,margin_x_right);
        if(x<margin_x_left && world[x][y].females.size()+world[x][y].males.size()>0)
          margin_x_left=x;
        if(x>margin_x_right && world[x][y].females.size()+world[x][y].males.size()>0)
          margin_x_right=x;
      } 
    }

    if(t>BURN_IN_TIME)
    {
          //output
      for(int x=0; x<world_size_x;x++)
      {
        for(int y=0; y<world_size_y ;y++)
        {

         if(x>=world_size_x/2-burn_in_x/2 && x<=world_size_x/2+burn_in_x/2-1 && t>BURN_IN_TIME-50 )//output the central 10x5 patches
          output_metapopulation(op,x,y,r,t,margin_x_left,margin_x_right);
        if(x<margin_x_left && world[x][y].females.size()+world[x][y].males.size()>0)
          margin_x_left=x;
        if(x>margin_x_right && world[x][y].females.size()+world[x][y].males.size()>0)
          margin_x_right=x;
      } 
    }
  }
  if(t>BURN_IN_TIME)
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

  if(t<=BURN_IN_TIME)
  {
    if(t%2500==0)
    {
      //op1<<r<<" "<<t<<" ";
      //op1<<heterozygosity(world_size_x/2-burn_in_x/2,world_size_x/2+burn_in_x/2)<<" NA NA"<<endl;
      vector<vector<double> > quartile=trait_calulator(world_size_x/2-burn_in_x/2, world_size_x/2+burn_in_x/2,0);
      op2<<r<<" "<<t<<" ";
      op2<<quartile.at(0).at(0)<<" "<<quartile.at(0).at(1)<<" "<<quartile.at(0).at(2)<<" ";
      quartile=trait_calulator(world_size_x/2-burn_in_x/2, world_size_x/2+burn_in_x/2,1);
      op2<<quartile.at(0).at(0)<<" "<<quartile.at(0).at(1)<<" "<<quartile.at(0).at(2)<<" ";
      op2<<endl;
    }
  }

  if(t>BURN_IN_TIME && t%10==0)
  {
    //op1<<r<<" "<<t<<" ";
    //op1<<heterozygosity(world_size_x/2-burn_in_x/2,world_size_x/2+burn_in_x/2)<<
    //" "<<heterozygosity(margin_x_left,margin_x_left+5)<<" "<<heterozygosity(margin_x_right-4,margin_x_right+1)<<endl;
    op2<<r<<" "<<t<<" ";
    op3<<r<<" "<<t<<" ";
    op4<<r<<" "<<t<<" ";
    vector<vector<double> > quartile=trait_calulator(world_size_x/2-burn_in_x/2, world_size_x/2+burn_in_x/2,0);
    op2<<quartile.at(0).at(0)<<" "<<quartile.at(0).at(1)<<" "<<quartile.at(0).at(2)<<" ";
    quartile=trait_calulator(world_size_x/2-burn_in_x/2, world_size_x/2+burn_in_x/2,1);
    op2<<quartile.at(0).at(0)<<" "<<quartile.at(0).at(1)<<" "<<quartile.at(0).at(2)<<" ";
    op2<<endl;
    
    quartile=trait_calulator(margin_x_left,margin_x_left+5,0);
    op3<<quartile.at(0).at(0)<<" "<<quartile.at(0).at(1)<<" "<<quartile.at(0).at(2)<<" ";
    quartile=trait_calulator(margin_x_left,margin_x_left+5,1);
    op3<<quartile.at(0).at(0)<<" "<<quartile.at(0).at(1)<<" "<<quartile.at(0).at(2)<<" ";
    op3<<endl;
    
    quartile=trait_calulator(margin_x_right-4,margin_x_right+1,0);
    op4<<quartile.at(0).at(0)<<" "<<quartile.at(0).at(1)<<" "<<quartile.at(0).at(2)<<" ";
    quartile=trait_calulator(margin_x_right-4,margin_x_right+1,1);
    op4<<quartile.at(0).at(0)<<" "<<quartile.at(0).at(1)<<" "<<quartile.at(0).at(2)<<" ";
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
  reproduce();
  death();
  patch_extinction();

  t++;
  if(t>BURN_IN_TIME+10000)
    break;

}
while(margin_x_right!=world_size_x-1 || margin_x_left!=0);

}
//close output files
op.close();

//op1.close();

op2.close();

op3.close();

op4.close();


op5.close();
return 0;

}







