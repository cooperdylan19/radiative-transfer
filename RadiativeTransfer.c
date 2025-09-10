// two slashes used to comment on code
/// three slashes used to comment on question

// the world explodes error message means file hasnt opened succesfully

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdint.h> 
#include <stdbool.h>
#include <stdlib.h> // relevant libraries for all question parts

/// functions used throughout the code

int64_t random_number(int64_t *seed, int64_t a, int64_t c, int64_t m) // code from chapter 2.1 - PRNG
{
    *seed = (a * (*seed) + c) % m;
    return *seed; // modified, seed is now a pointer
}

double P(double tau) // P(tau) defined in question 1a
{ 
    return exp(-tau);
}

double inverseCDF(double y) // inverse CDF of P (pen and paper calculation)
{
    return -log(1 - y);
}

/// functions defined for Q1a begin here

void rawData(int N, int64_t *seed, int64_t a, int64_t c, int64_t m) // to generate a plot of data before apply sampling methods
{
    clock_t start, end;
    double CPUtime;

    start = clock(); //time how long it takes to run through functions

    FILE *fp;
    fp = fopen("rawData.csv", "w"); // store in a file

    if (!fp) 
    {
        printf("the world explodes\n"); // ensure file is actually open
        return;
    }

    for (int i = 0; i < N - 1; i++) // loop through N points, randomly generate x,y values in [0,1]
    {
        double y = ((double)random_number(seed, a, c, m) / (double)m);
        double x = ((double)random_number(seed, a, c, m) / (double)m);

        fprintf(fp, "%f,%f\n", x, y); // print values into x, y
    }

    printf("rawData.csv file has been created with %d data points. Function runtime: %f seconds\n", N, CPUtime); 

    fclose(fp); // close file

    end = clock(); // End timer
    CPUtime = ((double) (end - start)) / CLOCKS_PER_SEC; 
}

void rejectionMethod(int N, int64_t *seed, int64_t a, int64_t c, int64_t m) // take same arguments as rawData for convenience in main()
{
    clock_t start, end;
    double CPUtime;

    start = clock(); //time how long it takes to run through functions

    FILE *fp;
    fp = fopen("rejectionData.csv", "w");

    if (!fp) 
    {
        printf("the world explodes\n");
        return;
    }

    for (int i = 0; i < N - 1; i++) 
    {
        double x, y;
        do 
        {
            x = ((double)random_number(seed, a, c, m) / (double)m); 
            y = ((double)random_number(seed, a, c, m) / (double)m); 
        } while (y > P(x)); // if coefficient A is not 1, multiply y by 1 / A (as in sin^2 class example)
        
        fprintf(fp, "%f,%f\n", x, y); // as before...
    }

    printf("rejectionData.csv file has been created with %d data points. Function runtime: %f seconds\n", N, CPUtime); 

    fclose(fp); // close file

    end = clock(); // End timer
    CPUtime = ((double) (end - start)) / CLOCKS_PER_SEC; 
}

void CDFMethod(int N, int64_t *seed, int64_t a, int64_t c, int64_t m) // same arguments again
{
    clock_t start, end;
    double CPUtime;

    start = clock(); //time how long it takes to run through functions

    FILE *fp;
    fp = fopen("CDFData.csv", "w");

    if (!fp) 
    {
        printf("the world explodes\n");
        return;
    }

    for (int i = 0; i < N - 1; i++) // same loop
    {
        double y = ((double)random_number(seed, a, c, m) / (double)m);
        double x = inverseCDF(y); 

        fprintf(fp, "%f\n", x); // only care for x values (histogram)
    }

    printf("CDFData.csv file has been created with %d data points. Function runtime: %f seconds\n", N, CPUtime); 

    fclose(fp); // close file

    end = clock(); // End timer
    CPUtime = ((double) (end - start)) / CLOCKS_PER_SEC; 
}

/// functions defined for Q1a end here

/// functions defined for Q1b begin here

typedef struct 
{
    double r;
    double theta;
    double phi; // SPCs

    double z; // atmosphere depth [zmin = 0, zmax = 1]
} variablesVector;

variablesVector initializeVariables(int64_t *seed, int64_t a, int64_t c, int64_t m) // takes same arguments as random_number
{
    variablesVector photon;
    double tau = 10.0; // optical depth (in question)

    photon.r = inverseCDF((double)random_number(seed, a, c, m) / m) / tau; // inverseCDF / optical depth
    photon.theta = (M_PI * (double)random_number(seed, a, c, m) / m) - (M_PI / 2); // arccos scales theta to [-pi/2, +pi/2] as expected for SPC
    photon.phi = 2 * acos(2 * ((double)random_number(seed, a, c, m) / m) - 1); // [0,+2pi]

    photon.z = photon.r * cos(photon.theta); // z = r * cos(theta)

    return photon; // need for simulatePhotons
}

FILE** photonTraceFiles(int NtrackedPhotons, FILE *mainFilePointer) 
{
    FILE **fps = malloc(sizeof(FILE*) * NtrackedPhotons);
    if (!fps) 
    {
        printf("memory allocation error\n"); // allocate some memory for files
        return NULL;
    }

    for (int i = 0; i < NtrackedPhotons; i++) 
    {
        char filename[256];
        sprintf(filename, "photon%d.csv", i + 1); // naming files
        fps[i] = fopen(filename, "w");
        if (!fps[i]) 
        {
            printf("the world explodes. File: %s\n", filename);
            for (int j = 0; j < i; j++) 
            {
                fclose(fps[j]);
            }
            fclose(mainFilePointer); // remember to close for error
            free(fps);
            return NULL;
        }
    }
    return fps;
}

void simulatePhotons(int *binCount, int Nphotons, int Nbins, double zMax, double albedo, double binWidth, int64_t *seed, int64_t a, int64_t c, int64_t m) 
{
    int successfullyBinnedPhotons = 0;
    variablesVector photon;
    int NtrackedPhotons = 5; // no. photons I want to track (1 million expectedly made my computer comparable to a jet engine)

    FILE *fp = fopen("finalCoordinates.csv", "w"); // plot final positions to see scattering effects
    if (!fp) 
    {
        printf("the world explodes\n");
        return;
    }

    FILE **fps = photonTraceFiles(NtrackedPhotons, fp); // open trace files
    if (!fps) 
    {
        return; // errors dealt with in pTF function
    }

    int photonIndex = 0;
    while (successfullyBinnedPhotons < Nphotons) 
    {
        photon = (variablesVector){0.0, 0.0, 0.0, 0.0};
        if (photonIndex < NtrackedPhotons) 
        {
            fprintf(fps[photonIndex], "%f,%f,%f\n", 0.0, 0.0, 0.0); // print inital position into files
        }

        photon = initializeVariables(seed, a, c, m);

        while (photon.z <= zMax && photon.z >= 0) // while in the mediuum
        {
            if (albedo >= random_number(seed, a, c, m) / m) // absorbtion check
            {
                variablesVector step = initializeVariables(seed, a, c, m);
                photon.r += step.r; 
                photon.theta = step.theta;
                photon.phi = step.phi;
                photon.z += step.z; // moving the photon randomly

                if (photonIndex < NtrackedPhotons) // tracing first 5 photons
                {
                    double x = photon.r * sin(photon.theta) * cos(photon.phi);
                    double y = photon.r * sin(photon.theta) * sin(photon.phi);
                    double z = photon.z;
                    fprintf(fps[photonIndex], "%f,%f,%f\n", x, y, z);
                }
            } else 
            {
                photon.z = -1; // absorb photon (place outside medium)
                break;
            }
        }

        if ((photon.z > zMax || photon.z < 0) && photon.z != -1) // also want examples of photons exiting wrong way (z=0)
        {
            double x = photon.r * sin(photon.theta) * cos(photon.phi);
            double y = photon.r * sin(photon.theta) * sin(photon.phi);
            double z = photon.z;
            fprintf(fp, "%f,%f,%f\n", x, y, z); // convert to CC for file as before
            
            int i_bin = (int)(photon.theta / binWidth);
            if (i_bin >= 0 && i_bin < Nbins) 
            {
                binCount[i_bin]++;
                successfullyBinnedPhotons++; // binning photons according to theta
            }
        }
        
        photonIndex++;
    }

    fclose(fp);
    for (int i = 0; i < NtrackedPhotons; i++) 
    {
        fclose(fps[i]); // close files...
    }
    free(fps); // ... and free memory
}

void initializeBins(int *arrayPoint, int Nbins) 
{
    for (int i = 0; i < Nbins; i++) 
    {
        arrayPoint[i] = 0; // convenience function, sets all array elements to 0
    }
}

void printResults(int *arrayPoint, int Nbins, double binWidth, int Nphotons) 
{
    printf("Bin range        | Midpoint μ     | No. photons in bin | Normalized count\n");
    printf("------------------------------------------------------------------------\n"); // making a table
    
    FILE *fp = fopen("photonBinning.csv", "w"); 
    if (!fp) 
    {
        printf("the world explodes\n"); 
        return;
    }
    
    int totalBinnedPhotons = 0; // count the number of photons binned (ensure all 1 million binned)
    
    for (int i = 0; i < Nbins; i++) 
    {
        double lowerBoundFactor = (double)i * binWidth / M_PI; // express the bin ranges in terms of pi, easier to check correct calculation
        double upperBoundFactor = (double)(i + 1) * binWidth / M_PI;

        double midpoint = ((i + 0.5) * binWidth); // midpoints expressed normally as not sure if pi character corresponds to an actual number in Maple

        double normalizedCount = (double)arrayPoint[i] / Nphotons; // normalize
        totalBinnedPhotons += arrayPoint[i]; // sum photons in bins

        // print results to console
        printf("[%1.3fπ, %1.3fπ] | %-14.3f | %-18d | %-16f\n", lowerBoundFactor, upperBoundFactor, midpoint, arrayPoint[i], normalizedCount);

        // write the same data into the CSV
        fprintf(fp, "[%1.3fπ, %1.3fπ],%f,%d,%f\n", lowerBoundFactor, upperBoundFactor, midpoint, arrayPoint[i], normalizedCount);
    }

    printf("\nTotal number of photons binned: %d\n", totalBinnedPhotons); // print total number of photons binned to console
    
    fclose(fp);
}

void runFlowchart(int Nphotons, int Nbins, double zMax, double albedo, double binWidth, int64_t seed, int64_t a, int64_t c, int64_t m) 
{
    clock_t start, end;
    double CPUtime;

    start = clock();

    int arrayPoint[Nbins]; // convenience function, means only need to run one function in main()

    initializeBins(arrayPoint, Nbins);
    simulatePhotons(arrayPoint, Nphotons, Nbins, zMax, albedo, binWidth, &seed, a, c, m);
    printResults(arrayPoint, Nbins, binWidth, Nphotons);

    end = clock(); // End timer
    CPUtime = ((double) (end - start)) / CLOCKS_PER_SEC; 

    printf("\n"); // gap from table
    printf("Function runtime: %f seconds\n", CPUtime); 
}

/// functions defined for Q1b end here

/// functions defined for Q1c start here

double RayleighPDF(double x) 
{
    return (3 / 4) * (1 + cos(x) * cos(x)) * (1 / (4 * M_PI)); // Rayleigh scattering equation for rejection method sampling
}

double SPCarray[3] = {0, 0, 0}; // same logic as previous question, using arrays now for matrix multiplication later [r, theta, phi]

void initializeVariablesMatrix(int64_t *seed, int64_t a, int64_t c, int64_t m, double tau) // reworked iV from earlier, want tau variable now
{
    double Pmax = 3 / (8 * M_PI); // maximum value since MAX{cos^2} = 1 (normalization factor)

    SPCarray[0] = inverseCDF((double)random_number(seed, a, c, m) / m) / tau; // r

    double x = 2 * ((double)random_number(seed, a, c, m) / m) - 1;
    double theta = 2 * acos(x);
    double y = Pmax * (random_number(seed, a, c, m) / m); // using rejection method to sample theta values [1]

    int iterationCount = 0;

    while (y > RayleighPDF(x)) // && iterationCount < 10) 
    {
        double y = Pmax * (random_number(seed, a, c, m) / m);
        double x = 2 * ((double)random_number(seed, a, c, m) / m) - 1;
        double theta = 2 * acos(x);
        // iterationCount++;
    }
    //if (iterationCount >= 10) // TEST
    // printf("Max iterations reached, breaking loop\n");
    

    SPCarray[1] = theta; // theta
    SPCarray[2] = ((double)random_number(seed, a, c, m) / m) * 2 * M_PI; // phi



    //printf("Initial r: %f, theta: %f, phi: %f\n", SPCarray[0], SPCarray[1], SPCarray[2]); TEST CORRECT RANGES
}

double E[4][4]; // initialize 4 x 4 rotation matrix

void scatteringMatrix(double theta) // defined according to ref [1] {bottom of page}
{
    double M = cos(theta);
    double a = 3.0 / 4.0;

    E[0][0] = a * (M * M + 1); // E just represents an element of the scattering matrix (P used in [1], but already used here in line 19)
    E[0][1] = a * (M * M - 1);
    E[0][2] = 0;
    E[0][3] = 0;

    E[1][0] = a * (M * M - 1);
    E[1][1] = a * (M * M + 1);
    E[1][2] = 0;
    E[1][3] = 0;

    E[2][0] = 0;
    E[2][1] = 0;
    E[2][2] = a * 2 * M;
    E[2][3] = 0;

    E[3][0] = 0;
    E[3][1] = 0;
    E[3][2] = 0;
    E[3][3] = a * 2 * M;
}

void RotateVector(double r, double theta, double phi) 
{
    double xNew = r * sin(theta) * cos(phi); // SPC -> CC (required to multiply with rotation matrix)
    double yNew = r * sin(theta) * sin(phi);
    double zNew = r * cos(theta);

    scatteringMatrix(theta); // determine rotation matrix for CURRENT theta (axis is therefore variable)

    double inputVector[4][1] = {{xNew}, {yNew}, {zNew}, {1}}; // set up initial vector for multiplication ({1} required dimensionally)

    double outputVector[4][1] = {0}; // initialize output coordinates

    for (int i = 0; i < 4; i++) 
    {
        for (int k = 0; k < 4; k++) 
        {
            outputVector[i][0] += E[i][k] * inputVector[k][0]; // '1 x 4' * '4 x 4' = '1 x 4'
        }
    }

    for (int i = 0; i < 3; i++) // dont care about 4th element, want x', y', z'
    {
        SPCarray[i] = outputVector[i][0]; // update SPC array with new rotated coordinates
    }
}

void simulatePhotons2(int Nphotons, int Nbins, double zMax, double albedo, double binWidth, int64_t *seed, int64_t a, int64_t c, int64_t m, double tau) 
{
    int successfullyBinnedPhotons = 0; // modified simulatePhotons to work with rotated matrix coordinates
    int photonsRemoved = 0;

    char fileName[256]; 
    sprintf(fileName, "finalCoordinates2_tau%.2f.csv", tau); // want file name to include tau value for clarity (will call for unique tau)

    FILE *fp = fopen(fileName, "w");
    if (!fp) {
        printf("the world explodes\n");
        return;
    }

    for (int photonIndex = 0; photonIndex < Nphotons; photonIndex++) 
    {
        double x = 0.0, y = 0.0, z = 0.0; // photons begin at the origin

        initializeVariablesMatrix(seed, a, c, m, tau); // generate initial movements...

        SPCarray[1] = 0; // ... but set theta = 0 since we want to eject photons straight upwards (along z)

        x += SPCarray[0] * sin(SPCarray[1]) * cos(SPCarray[2]);
        y += SPCarray[0] * sin(SPCarray[1]) * sin(SPCarray[2]);
        z += SPCarray[0] * cos(SPCarray[1]);

        while (z <= zMax && z >= 0) // same scattering logic as in simulatePhotons...
        {
            if (albedo >= (double)random_number(seed, a, c, m) / m) 
            {
                initializeVariablesMatrix(seed, a, c, m, tau);
                RotateVector(SPCarray[0], SPCarray[1], SPCarray[2]); // here the coordinates are rotated according to sM

                x += SPCarray[0] * sin(SPCarray[1]) * cos(SPCarray[2]);
                y += SPCarray[0] * sin(SPCarray[1]) * sin(SPCarray[2]); // update coordinates in CC
                z += SPCarray[0] * cos(SPCarray[1]);
            } else 
            {
                z = -1; // absorb
                photonsRemoved++;
                break;
            }
        }

        if (z >= zMax || z <= 0) 
        {
            successfullyBinnedPhotons++;
            fprintf(fp, "%f, %f, %f\n", x, y, z); // track final coordinates as before
        }
    }
    printf("%s file has been created with %d data points for an optical depth of %f\n", fileName, successfullyBinnedPhotons, tau);
    fclose(fp);
}

/// functions defined for Q1c end here



int main() // want to keep main() as simple as possible
{
    /// Q1a
    printf("QUESTION 1A\n");

    int N = 10000; //10,000 
    int64_t seed = time(NULL); // time elapsed since the unich epoch (always different = random)
    int64_t a = 16807;
    int64_t c = 0;
    int64_t m = 147483647; // defined in chapter 2 (minimal standard generator)

    rawData(N, &seed, a, c, m);
    rejectionMethod(N, &seed, a, c, m);
    CDFMethod(N, &seed, a, c, m);

    /// Q1b
    printf("\n");
    printf("QUESTION 1B\n");

    int Nphotons = 1000000; // 1 million photons (1,000 for maple plots due to runtime)
    int Nbins = 10;
    double zMax = 1.0;
    double albedo = 1.0;
    double binWidth = M_PI / (2 * Nbins); // all defined in question

    runFlowchart(Nphotons, Nbins, zMax, albedo, binWidth, seed, a, c, m);

    /// Q1c
    printf("\n");
    printf("QUESTION 1C\n");

    double tauBlue = 10;
    double tauOther = 0.1;

    // initializeVariablesMatrix(&seed, a, c, m);
    simulatePhotons2(Nphotons, Nbins, zMax, albedo, binWidth, &seed, a, c, m, tauBlue);
    simulatePhotons2(Nphotons, Nbins, zMax, albedo, binWidth, &seed, a, c, m, tauOther);

    return 0;
}

// REFERENCES
// [1] https://adsabs.harvard.edu/pdf/2011BASI...39..101W - rotation/ scattering matrix for Rayleigh scattering
