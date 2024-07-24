#include "scarf.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>
#include <time.h>

#define ROUNDS 8
#define INPUT_SP (1<<N)
// nb of 5 bit words
#define HALF_INPUT_SP (1<<(N/2))
// for computing the ddt
#define NB_SAMPLE_DDT (1<<15)
// for computing the ddt on 5bit differentials
#define NB_SAMPLE_5bit_DDT (1<<20)
// for testing a fixed differential
#define NB_SAMPLE_TEST (1<<15)
// for computing the DDT for the round function 
#define NB_SAMPLE_ROUND_DDT (1<<15)
#define tenb_MASK 0b1111111111
#define sixtyb_MASK 0xFFFFFFFFFFFFFFF


const double prob = 1/((1<<N)-1.0);

double log2(double arg);

double * DDT(int rounds);

double * G_DDT(void);

double test_differential(int alpha, int beta, int rounds,int mode);

void sample_ddt(int mode);

uint64_t getRandom_(void);

uint64_t getRandom(void);

double getmax(double * ddt, int * alpha, int * beta, int size);

void write_ddt_to_files(FILE * fp, double* ddt, int size);

double * fivebitDDT(int rounds);

double log2(double arg){
    return log(arg)/log(2);
}

double * round_DDT(int mode);

//get best differential in the ddt, puts the differential in alpha_m and beta_m
double getmax(double * ddt, int * alpha_m, int * beta_m, int size){
    double max = 0.0; 
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

uint64_t getRandom(void){
    uint64_t x, y;
    x = rand();
    y = rand();
    uint64_t ret = ((x<<30)+y);
    return ret&sixtyb_MASK;
}


// returns a pseudo random 60bits value (4 upper bits set to 0) using arc4random
// uint64_t getRandom_(void){
//     uint64_t x, y;
//     x = arc4random();
//     y = arc4random();
//     uint64_t ret = ((x<<32)+y);
//     return ret&sixtyb_MASK;
// }

double * G_DDT(void){
    double * ddt = malloc(sizeof(double)*HALF_INPUT_SP);
    for (int c = 0; c < NB_SAMPLE_ROUND_DDT; c++){
                uint64_t k = getRandom(); k = extract_bits(k,0,30);
                for (int alpha = 1; alpha < HALF_INPUT_SP; alpha++){
                    for(int x=0; x < HALF_INPUT_SP; x++){
                        int xx = (x^alpha)&tenb_MASK;
                        int y = Gfunc(x,k);
                        int yy = Gfunc(xx,k);
                        if (y == yy){ddt[alpha]++;}
                    }    
                }
            }
            for (int alpha = 1; alpha < HALF_INPUT_SP; alpha++){
                ddt[alpha] = ddt[alpha]/(HALF_INPUT_SP*NB_SAMPLE_ROUND_DDT);
            }
    return ddt;
}



// mode = 0 : 10bit
// mode = 1 : 5bit
double * round_DDT(int mode){
    double * ddt = malloc(sizeof(double)*INPUT_SP*INPUT_SP);

    switch (mode){
        case 0 :
            for (int c = 0; c < NB_SAMPLE_ROUND_DDT; c++){
                uint64_t k = getRandom(); k = extract_bits(k,0,30);
                for (int alpha = 1; alpha < INPUT_SP; alpha++){
                    for(int x=0; x < INPUT_SP; x++){
                        int xx = (x^alpha)&tenb_MASK;
                        int y = R1(x,k);
                        int yy = R1(xx,k);
                        int beta = (y ^ yy)&tenb_MASK;
                        ddt[alpha + (INPUT_SP*beta)]++;
                    }    
                }
            }
            for (int alpha = 1; alpha < INPUT_SP; alpha++){
                for (int beta=1; beta < INPUT_SP; beta++){
                // averaging the differential probabilites
                ddt[alpha+(INPUT_SP*beta)] = ddt[alpha+(INPUT_SP*beta)]/(INPUT_SP*NB_SAMPLE_ROUND_DDT);
                }
            }
            break;

        case 1 : 
            {
            for (int c = 0; c < NB_SAMPLE_ROUND_DDT; c++){
                uint64_t k = getRandom(); uint64_t k1 = extract_bits(k,0,30); uint64_t k2 = extract_bits(k,30,60);
                for (int alpha = 1; alpha < HALF_INPUT_SP; alpha++){
                    for(int x=0; x < INPUT_SP; x++){
                        int xx = (x^alpha)&tenb_MASK;
                        int y = R1(x,k1);
                        int yy = R1(xx,k1);
                        y = R1(y,k2);
                        yy = R1(yy,k2);
                        int beta = ((y ^ yy)&tenb_MASK);
                        if (beta == alpha){
                            ddt[alpha + (HALF_INPUT_SP*beta)]++;
                        }
                        // printf("%f %d %x %x\n", ddt[alpha+(HALF_INPUT_SP*beta)], c, alpha, beta);
                    }    
                }
            }
            for (int alpha = 1; alpha < HALF_INPUT_SP; alpha++){
                for (int beta=1; beta < HALF_INPUT_SP; beta++){
                // averaging the differential probabilites
                ddt[alpha+(HALF_INPUT_SP*beta)] = ddt[alpha+(HALF_INPUT_SP*beta)]/(INPUT_SP*NB_SAMPLE_ROUND_DDT);
                }
            }
            break;
            }
    }
    return ddt;
}

double * DDT(int rounds){
    double * ddt = malloc(sizeof(double)*INPUT_SP*INPUT_SP);
   
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

void write_ddt_to_file(FILE * fp, double* ddt, int size){
    for (int alpha = 0; alpha < size; alpha++){
        for (int beta = 0; beta < size; beta++){
            fprintf(fp, "%f ", ddt[alpha + (size*beta)]);
        }
        fprintf(fp,"\n");
    }
}

double * fivebitDDT(int rounds){
    double * ddt = malloc(sizeof(double)*HALF_INPUT_SP*HALF_INPUT_SP);
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
// mode 0 test for full rounds
// mode 1 test for 2r
double test_differential(int alpha, int beta, int rounds, int mode){
    int counter = 0;
    for (int c = 0; c < NB_SAMPLE_TEST; c++){
        uint64_t k3, k2, k1, k0, t1, t2;
        k3 = getRandom(), k2 = getRandom(), k1 = getRandom(), k0 = getRandom(), 
        t1 = getRandom(), t2 = getRandom();

        switch (mode) {
            case 0 :
                for(int x=0; x < INPUT_SP; x++){
                    int xx = (x^alpha)&tenb_MASK;
                    int y = etilda(x,k3,k2,k1,k0,t1,t2,rounds);
                    int yy = etilda(xx,k3,k2,k1,k0,t1,t2,rounds);
                    int bbeta = (y ^ yy)&tenb_MASK;
                    if (bbeta == beta) {
                        counter++;
                    } 
                }
                break;
            case 1 :
                for(int x=0; x < INPUT_SP; x++){
                    int xx = (x^alpha)&tenb_MASK;
                    int y = R1(x,k1);
                    int yy = R1(xx,k1);
                    y = R1(y,k2);
                    yy = R1(yy,k2);
                    int bbeta = (y ^ yy)&tenb_MASK;
                    if (bbeta == beta) {
                        counter++;
                    } 
                }
        }
    }
    return ((double)counter)/(NB_SAMPLE_TEST*INPUT_SP);
}

// if mode = 0 then we do full ddt
// mode = 1 then we do 5bit ddt
// mode = 2 one round full ddt
// mode = 3 one round 5bit ddt
void sample_ddt(int mode){
    FILE * fp1, * fp2;
    switch (mode){
        case 0:
            fp1 = fopen("best_log_results.txt", "w+");

            fp2 = fopen("avg_ddt_result.txt", "w+");

            #pragma omp parallel for ordered
            for (int r = 2; r < 9; r++){
                int alpha_m = 0, beta_m = 0;
                double * ddt = DDT(r);
                double max = getmax(ddt, &alpha_m, &beta_m, INPUT_SP);
                double bias = max - prob;
                double logbias = -log2(bias);
                #pragma omp ordered
                {
                fprintf(fp2, "rounds : %d\n", r);
                printf("nb rounds : %d, logbias = %f, alpha = %x, beta = %x \n", r, logbias, alpha_m, beta_m);
                fprintf(fp1, "nb rounds : %d, logbias = %f, alpha = %x, beta = %x \n", r, logbias, alpha_m, beta_m);
                write_ddt_to_file(fp2,ddt, INPUT_SP);
                }
                
                free(ddt);
            }
            break;
        case 1: 
            fp1 = fopen("best_log_results_5bits.txt", "w+");

            fp2 = fopen("avg_ddt_result_5bits.txt", "w+");

            #pragma omp parallel for ordered
            for (int r = 2; r < 9; r++){
                int alpha_m = 0, beta_m = 0;
                double * ddt = fivebitDDT(r);
                double max = getmax(ddt, &alpha_m, &beta_m, HALF_INPUT_SP);
                double bias = max - prob;
                double logbias = -log2(bias);
                #pragma omp ordered
                {
                fprintf(fp2, "rounds : %d\n", r);
                printf("nb rounds : %d, logbias = %f, alpha = %x, beta = %x \n", r, logbias, alpha_m, beta_m);
                fprintf(fp1, "nb rounds : %d, logbias = %f, alpha = %x, beta = %x \n", r, logbias, alpha_m, beta_m);
                write_ddt_to_file(fp2,ddt, HALF_INPUT_SP);
                }
                free(ddt); 
            }
            break;
        case 2: 
            {
            fp1 = fopen("best_log_results_round_full.txt", "w+");

            fp2 = fopen("avg_ddt_result_round_full.txt", "w+");

            int alpha_m = 0, beta_m = 0;
            double * ddt = round_DDT(0);
            double max = getmax(ddt, &alpha_m, &beta_m, INPUT_SP);
            double bias = max - prob;
            double logbias = -log2(bias);
            printf("logbias = %f, alpha = %x, beta = %x \n", logbias, alpha_m, beta_m);
            fprintf(fp1, "logbias = %f, alpha = %x, beta = %x \n", logbias, alpha_m, beta_m);
            write_ddt_to_file(fp2,ddt, INPUT_SP);
            break;
            }
        case 3:
            {
            fp1 = fopen("best_log_results_round_5bits.txt", "w+");

            fp2 = fopen("avg_ddt_result_round_5bits.txt", "w+");

            int alpha_m = 0, beta_m = 0;
            double * ddt = round_DDT(1);
            double max = getmax(ddt, &alpha_m, &beta_m, HALF_INPUT_SP);
            double bias = max - prob;
            double logbias = -log2(bias);
            printf("logbias = %f, alpha = %x, beta = %x \n", logbias, alpha_m, beta_m);
            fprintf(fp1, "logbias = %f, alpha = %x, beta = %x \n", logbias, alpha_m, beta_m);
            write_ddt_to_file(fp2,ddt, HALF_INPUT_SP);
            }
    }
    
}

int main(){
    srand ( time(NULL) );
    // for (int i = 0; i < 20; i++){
    //     sample_ddt(3);
    // }
    double * ddt = G_DDT();
    for (int alpha = 1; alpha < HALF_INPUT_SP; alpha++){
        printf("%f\n", -log2(ddt[alpha]));
    }

    uint64_t k = getRandom(); k = extract_bits(k,0,30);
    int tab[HALF_INPUT_SP];
    for(int x=0; x < HALF_INPUT_SP; x++){
        tab[x] = Gfunc(x,k);
    }
    int count = 0;
    for (int val=0; val < HALF_INPUT_SP; val++){
        for (int x=0; x < HALF_INPUT_SP; x++){
            if (tab[x]==val){
                count++;
                break;
            }
        }
    }
    printf("%d\n",count);
    for (int x = 0; x < HALF_INPUT_SP; x++){
        printf("%d %d\n", x, tab[x]);
    }

    // double pr = test_differential(7,7,8,1);
    // pr = pr - prob;
    // double logbias = -log2(pr);
    // printf("bias : %f, logbias : %f\n",pr, logbias);
    return 0;
}