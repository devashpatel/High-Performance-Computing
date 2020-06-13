//Candidate Number: 099500 - Solving the Advection equation, when placed on the University of Exeter HPC System
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>

int main(){

    double u[1000][1000];
    double updatedArray[1000][1000];

    int x,y,i; 
    float dUdt; 

    const float deltaX = 0.001;
    const float deltaY = 0.001; 
    const float sigma = 0.03;
    const float x0= 0.1;
    const float y0 = 0.1; 
    const float velocity = 0.01;   
    const float timestep = 0.05;
    const bool doIO = true;

    float denominator = 2*pow(sigma,2); 


    //Inital Conditions
    #pragma omp parallel for default(none) shared(u,denominator) private(x,y) schedule(dynamic)
    for (x=0; x<1000; x++) 
    {    
        for (y=0; y<1000; y++) 
        {
        //Inital Conditions Calculation
        u[x][y] = exp(- (pow((float)x*deltaX - x0, 2)
                    + pow((float)y*deltaY - y0, 2))
                    /denominator);
        }  
    }
        

    //Rate of change, 1500 time step
    for (i=0; i<1500; i++) 
    { 
    #pragma omp parallel for default(none) shared(u,updatedArray) private(dUdt,x,y) schedule(dynamic)
        for (x=0; x<1000; x++)
        {
            for (y=0; y<1000; y++)
            {
            //Finite Difference Calculation
            dUdt = ((-velocity * (u[x][y]-u[x-1][y]))/deltaX)
                    -((velocity  * (u[x][y]-u[x][y-1]))/deltaY); 
            updatedArray[x][y]=u[x][y]+dUdt*timestep; 
            }
        }

    #pragma omp parallel for default(none) shared(u,updatedArray) private(x,y) schedule(dynamic)
        for (x=0; x<1000; x++)
        {
            for (y=0; y<1000; y++) 
            {
             u[x][y]= updatedArray[x][y]; 
            }
        }    
    }


    //Writing Values to file    
    if (doIO)
    {
    FILE *outfile; 
    outfile=fopen("u.dat", "w");
        

    //Write U values to file
    for (x=0; x<1000; x++)
    {
        for (y=0; y<1000; y++)
        {
        float i = (float)x/(float)1000; 
        float j = (float)y/(float)1000;
        fprintf(outfile,"%f %f %f\n",i, j,u[x][y]);
        }
    }  
    fclose(outfile);
    }
    return 0; 
}