#include "scarf.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>

#define ROUNDS 8
#define INPUT_SP (1<<N)
// nb of 5 bit words
#define HALF_INPUT_SP (1<<(N/2))
// for computing the ddt
#define NB_SAMPLE_DDT (1<<15)
// for computing the ddt on 5bit differentials
#define NB_SAMPLE_5bit_DDT (1<<20)
// for testing a fixed differential
#define NB_SAMPLE_TEST (1<<20)
#define tenb_MASK 0b1111111111
#define sixtyb_MASK 0xFFFFFFFFFFFFFFF


const float prob = 1/((1<<N)-1.0);

double log2(double arg);

float * DDT(int rounds);

float test_differential(int alpha, int beta, int rounds);

void sample_ddt(int mode);

uint64_t getRandom(void);

float getmax(float * ddt, int * alpha, int * beta, int size);

void write_ddt_to_files(FILE * fp, float* ddt, int size);

float * fivebitDDT(int rounds);

double log2(double arg){
    return log(arg)/log(2);
}

//get best differential in the ddt, puts the differential in alpha_m and beta_m
float getmax(float * ddt, int * alpha_m, int * beta_m, int size){
    float max = 0.0; 
    for (int alpha = 1; alpha < size; alpha ++){
        for (int beta = 1; beta < size; beta++){
            if (max < ddt[alpha+(size*beta)]){
                max = ddt[alpha+(size*beta)];
                *alpha_m = alpha;
                *beta_m = beta;
            }
        }
    }
    return max;
}

// returns a pseudo random 60bits value (4 upper bits set to 0)
uint64_t getRandom(void){
    uint64_t x, y;
    x = arc4random();
    y = arc4random();
    uint64_t ret = ((x<<32)+y);
    return ret&sixtyb_MASK;
}

float * DDT(int rounds){
    float * ddt = malloc(sizeof(float)*INPUT_SP*INPUT_SP);
   
    for (int c = 0; c < NB_SAMPLE_DDT; c++){
        uint64_t k3, k2, k1, k0, t1, t2;
        k3 = getRandom(), k2 = getRandom(), k1 = getRandom(), k0 = getRandom(), 
        t1 = getRandom(), t2 = getRandom(); 

        for (int alpha = 1; alpha < INPUT_SP; alpha++){
            for(int x=0; x < INPUT_SP; x++){
                int xx = (x^alpha)&tenb_MASK;
                int y = etilda(x,k3,k2,k1,k0,t1,t2,rounds);
                int yy = etilda(xx,k3,k2,k1,k0,t1,t2,rounds);
                int beta = (y ^ yy)&tenb_MASK;
                ddt[alpha + (INPUT_SP*beta)]++;
            }
        
        }
    }
    for (int alpha = 1; alpha < INPUT_SP; alpha++){
        for (int beta=1; beta < INPUT_SP; beta++){
            // averaging the differential probabilites
            ddt[alpha+(INPUT_SP*beta)] = ddt[alpha+(INPUT_SP*beta)]/(INPUT_SP*NB_SAMPLE_DDT);
        }
    }
    return ddt;
}

void write_ddt_to_file(FILE * fp, float* ddt, int size){
    for (int alpha = 0; alpha < size; alpha++){
        for (int beta = 0; beta < size; beta++){
            fprintf(fp, "%f ", ddt[alpha + (size*beta)]);
        }
        fprintf(fp,"\n");
    }
}

float * fivebitDDT(int rounds){
    float * ddt = malloc(sizeof(float)*HALF_INPUT_SP*HALF_INPUT_SP);
    for (int c = 0; c < NB_SAMPLE_5bit_DDT; c++){
        uint64_t k3, k2, k1, k0, t1, t2;
        k3 = getRandom(), k2 = getRandom(), k1 = getRandom(), k0 = getRandom(), 
        t1 = getRandom(), t2 = getRandom();

        for (int alpha = 1; alpha < HALF_INPUT_SP; alpha++){
            for(int x=0; x < INPUT_SP; x++) {
                int xx = (x^alpha)&tenb_MASK;
                int y = etilda(x,k3,k2,k1,k0,t1,t2,rounds);
                int yy = etilda(xx,k3,k2,k1,k0,t1,t2,rounds);
                int beta = (y ^ yy)&tenb_MASK;
                if (beta == alpha) {
                    ddt[alpha+(HALF_INPUT_SP*beta)]++;
                }
            }
        }
    }
    for (int alpha = 1; alpha < HALF_INPUT_SP; alpha++){
        for (int beta=1; beta < HALF_INPUT_SP; beta++){
            // averaging the differential probabilites
            ddt[alpha+(HALF_INPUT_SP*beta)] = ddt[alpha+(HALF_INPUT_SP*beta)]/(INPUT_SP*NB_SAMPLE_5bit_DDT);
        }
    }
    return ddt;
}

float test_differential(int alpha, int beta, int rounds){
    int counter = 0;
    for (int c = 0; c < NB_SAMPLE_TEST; c++){
        uint64_t k3, k2, k1, k0, t1, t2;
        k3 = getRandom(), k2 = getRandom(), k1 = getRandom(), k0 = getRandom(), 
        t1 = getRandom(), t2 = getRandom();

        for(int x=0; x < INPUT_SP; x++){
            int xx = (x^alpha)&tenb_MASK;
            int y = etilda(x,k3,k2,k1,k0,t1,t2,rounds);
            int yy = etilda(xx,k3,k2,k1,k0,t1,t2,rounds);
            int bbeta = (y ^ yy)&tenb_MASK;
            if (bbeta == beta) {
                counter++;
            } 
        }
    }
    return ((float)counter)/(NB_SAMPLE_TEST*INPUT_SP);
}

// if mode = 0 then we do full ddt
// else if mode = 1 then we do 5bit ddt
void sample_ddt(int mode){
    FILE * fp1, * fp2;
    switch (mode){
        case 0:
            fp1 = fopen("best_log_results.txt", "w+");

            fp2 = fopen("avg_ddt_result.txt", "w+");

            #pragma omp parallel for ordered
            for (int r = 2; r < 9; r++){
                int alpha_m = 0, beta_m = 0;
                float * ddt = DDT(r);
                float max = getmax(ddt, &alpha_m, &beta_m, INPUT_SP);
                float bias = max - prob;
                double logbias = -log2(bias);
                #pragma omp ordered
                {
                fprintf(fp2, "rounds : %d\n", r);
                printf("nb rounds : %d, logbias = %f, alpha = %b, beta = %b \n", r, logbias, alpha_m, beta_m);
                fprintf(fp1, "nb rounds : %d, logbias = %f, alpha = %b, beta = %b \n", r, logbias, alpha_m, beta_m);
                write_ddt_to_file(fp2,ddt, INPUT_SP);
                }
                
                free(ddt);
            }
            break;
        case 1: 
            fp1 = fopen("best_log_results_5bits.txt", "w+");

            fp2 = fopen("avg_ddt_result_5bits.txt", "w+");

            #pragma omp parallel for ordered
            for (int r = 1; r < 9; r++){
                int alpha_m = 0, beta_m = 0;
                float * ddt = fivebitDDT(r);
                float max = getmax(ddt, &alpha_m, &beta_m, HALF_INPUT_SP);
                float bias = max - prob;
                double logbias = -log2(bias);
                #pragma omp ordered
                {
                fprintf(fp2, "rounds : %d\n", r);
                printf("nb rounds : %d, logbias = %f, alpha = %b, beta = %b \n", r, logbias, alpha_m, beta_m);
                fprintf(fp1, "nb rounds : %d, logbias = %f, alpha = %b, beta = %b \n", r, logbias, alpha_m, beta_m);
                write_ddt_to_file(fp2,ddt, HALF_INPUT_SP);
                }

                free(ddt); 
            }
            break;
    }
    
}

int main(){
    sample_ddt(1);
    return 0;
}