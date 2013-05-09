#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

//  Monster hack so we can print comd iteration number.
static int comd_count = 0;
static int loop_count = 0;

// If the line below is commented, use a linear stress-strain relationship
// If uncommented, calculate the stress from CoMD
#define ZERO_TEMP_COMD "./CoMD -x 6 -y 6 -z 6"

#define GNUPLOT "/usr/local/bin/gnuplot -persist"

int msleep(unsigned long milisec) {
    struct timespec req={0};
    time_t sec=(int)(milisec/1000);
    milisec=milisec-(sec*1000);
    req.tv_sec=sec;
    req.tv_nsec=milisec*1000000L;
    while(nanosleep(&req,&req)==-1)
        continue;
    return 1;
}

#define L 50
const double dx   = 1.0;
const double dt   = 0.5;

const double c    = 0.5; // wave speed for linear stress-strain relationship
const double rho0 = 0.5; // mass / undeformed MD volume

double t    = 0.0; // time (evolving)


// Notation: A = A_11 -- strain (deformation gradient) along the x axis
//           p        -- momentum density
//           e        -- energy density

double A0[L]; // initial condition for A

double A[L];
double p[L];
double e[L];

double temp1[L], temp2[L], temp3[L];
double temp4[L], temp5[L], temp6[L];


// The elastic energy corresponding to the deformation gradient A at zero temperature
// (for the linear stress-strain relationship)
double zeroTempEnergyDensity(double A) {
    return rho0*c*c*(A-1)*(A-1)/2;
}

double stressFn(double A, double e) {
  #ifdef ZERO_TEMP_COMD


//  Instead of opening a pipe, we write the 
  // open pipe to CoMD
  FILE *fPipe = popen(ZERO_TEMP_COMD,"r");
  if (fPipe == NULL) {
    printf("Error launching CoMD < %s >!\n", ZERO_TEMP_COMD);
    exit(0);
  }

  double stressErlang;
  fscanf(fPipe, "%lf", &stressErlang);
  fprintf(stdout, "stressErlang{%3d}[%2d]: %lf\n", loop_count, comd_count, stressErlang);
  /*

  // throw away all stdout from CoMD
  int bufferSize = 256;
  char buffer[bufferSize];
  while (fgets(buffer, sizeof(char)*bufferSize, fPipe)) {}
*/

  // close the pipe
  int exitStatus = pclose(fPipe);
  if (exitStatus != 0) {
    printf("CoMD failed with status %d\n", exitStatus);
    exit(0);
  }
/*
  // read the first element from the output stress file
  const char *stressFilename = "stressTensor.txt";
  FILE *fStress = fopen (stressFilename, "r");
  if (fStress == NULL) {
    printf("Error reading stress file < %s >!\n", stressFilename);
    exit(0);
  }
  double stressXX = 0;
  int nelems = fscanf(fStress, "%lf", &stressXX);
  if (nelems != 1) {
    printf("Error reading first element of stress file!\n");
  }
  fclose(fStress);
*/
  return stressErlang;

  #else
  
  return rho0*c*c*(A-1); // For this toy problem, just the derivative of the zero temp energy wrt A
  
  #endif
}

int mod(int x, int y) {
    int rem = x % y;
    return rem < 0 ? rem + y : rem;
}

// Exact solution for A(t), assuming linear stress/strain relation
double exactA(int i) {
    double del_i = (int) ((t * c) / dx);
    int i1 = mod(i+del_i, L);
    int i2 = mod(i-del_i, L);
    return 0.5*(A0[i1] + A0[i2]);
}


// inputs: conserved fields A, p, e
// outputs: fluxes f_A, f_p, f_e
void fluxes(double *A, double *p, double *e, double *f_A, double *f_p, double *f_e) {
    int i;
    comd_count = 0;
    for (i = 0; i < L; i++) {
        double stress = stressFn(A[i], e[i]);
        double v = p[i] / rho0;
        
        f_A[i] = -v;
        f_p[i] = -stress;
        f_e[i] = -stress*v;
        comd_count++;
    }
}

void initializedConservedFields() {
    int i;
    for (i = 0; i < L; i++) {
        A0[i] = A[i] = (i < L/2) ? 1.01 : 1.0; // small initial step in deformation gradient
        p[i] = 0;
        e[i] = zeroTempEnergyDensity(A[i]);
    }
}

double netEnergy() {
    double acc = 0.0;
    int i;
    for (i = 0; i < L; i++) {
        acc += e[i];
    }
    return acc;
}

double minmod(double x, double y) {
    return (x > 0 == y > 0) ? 0.5*(x+y) : 0;
}

void modulatedDerivative(double *w, double *ret) {
    int i;
    for (i = 0; i < L; i++) {
        int im = (i-1+L)%L;
        int ip = (i+1)%L;
        ret[i] = minmod(w[ip]-w[i], w[i]-w[im]) / dx;
    }
}

void halfStep2ndOrder() {    
    int i,j;
    double *ws[3] = {A, p, e};
    double *fs[3] = {temp1, temp2, temp3};
    double *wps[3] = {temp4, temp5, temp6};
    
    fluxes(ws[0], ws[1], ws[2], fs[0], fs[1], fs[2]);
    
    // calculate w^(n+1/2)
    for (j = 0; j < 3; j++) {
        double *w  = ws[j];
        double *wp = wps[j];
        double *f  = fs[j];
        modulatedDerivative(f, wp);
        for (i = 0; i < L; i++)
            wp[i] = w[i] - (0.5*dt/(2*dx)) * wp[i];
    }
    
    // calculate f^(n+1/2)
    fluxes(wps[0], wps[1], wps[2], fs[0], fs[1], fs[2]);
        
    // calculate w^(n+1)
    for (j = 0; j < 3; j++) {
        double *w  = ws[j];
        double *f = fs[j];
        
        double *dw = temp4;
        modulatedDerivative(w, dw);
        
        double w0 = w[0];
        for (i = 0; i < L; i++) {
            int ip = (i+1)%L;
            double wn = (i == L-1) ? w0 : w[i+1];
            w[i] = 0.5*(w[i] + wn) + (dx/8) * (dw[i] - dw[ip]) - (0.5*dt/dx) * (f[ip] - f[i]);
        }
    }
}


void fullStep() {
    int i,j;
    // two half steps on a staggered grid
    halfStep2ndOrder();
    halfStep2ndOrder();
    
    // rotate back to original index <-> coordinate map
    double *ws[3] = {A, p, e};
    for (j = 0; j < 3; j++) {
        double *w = ws[j];
        double w_last = w[L-1];
        for (i = L-1; i >= 0; i--)
            w[i] = w[i-1];
        w[0] = w_last;
    }
    
    t += dt;
}

void plotFields(FILE *gp) {
    int i,j;
    fprintf(gp, "set yrange [0.99:1.02]\n");
    fprintf(gp, "plot '-' u 1:2 ti 'A' w lp, '-' u 1:3 ti 'A exact' w lp\n\n");
    for (j = 0; j < 2; j++) { // gnuplot forces us to repeat ourselves
        for (i = 0; i < L; i++) {
            fprintf(gp, "%f %f %f\n", i*dx, A[i], exactA(i));
        }
        fprintf(gp, "e\n");
    }
    fflush(gp);
}

int main(int argc, char **argv) {
    int i,j;
    // open pipe to gnuplot instance
    FILE *gp;
    gp = popen(GNUPLOT,"w");
    if (gp==NULL) {
      printf("Error opening pipe to GNU plot < %s > t! \n", GNUPLOT);
        exit(0);
    }
    
    initializedConservedFields();
    
    loop_count = 0;
    for (i = 0; i < 100; i++) {
        for (j = 0; j < 1; j++) {
            fullStep();
        }
        plotFields(gp);
        //msleep(100);
        loop_count++;
    }
    
    pclose(gp);
    
    return 0;
}

