#include "LinearComplexity.h"
#include "cephes.h"
#include "LC_Random_Test.h"
#include "my_struct.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define length 1000000
#define samplesize 1000



unsigned char buffer[length/8];
const double alpha = 0.01;
const double alphau = 0.0001;
const unsigned int mask[15] = {0x0001, 0x0002, 0x0004, 0x0008,
                               0x0010, 0x0020, 0x0040, 0x0080, 
                               0x0100, 0x0200, 0x0400, 0x0800, 
                               0x1000, 0x2000, 0x4000};
int main(int argc, char ** argv) {
	
    Initialfenweidian();
	if (argc != 3) {
		printf("Usage: prog  position  seedFile\n");
		return -1;
	}
	int pos = atoi(argv[1]);
    if( pos < 0 || pos > 15) {
        printf("the input position is invalid. legal interval: [1, 15].");
        printf("Usage: prog  position  seedFile\n");
		return -1;
    }
	FILE *seed_fp = fopen(argv[2], "r");
	if (seed_fp == NULL) {
		printf("file %s is not exist\n");
		printf("Usage: prog  position  seedFile\n");
		return -1;
	}
	
    char outputFileName[100];
    sprintf(outputFileName, "C_Rand_pos%d_result.txt", pos);
    FILE *fp;
	fp = fopen(outputFileName, "w");
	if(fp == NULL) {
        printf("error: fp is NULL!\n");
        return -1;
    }
    fprintf(fp, 
"N=%d ,  n=%d \n\
NIST, HSY, DM1, DM2, HSY-CS7, DM1-CS7, DM2-CS7\n\
Applying erfc function\n", samplesize, length);
	
    int RejectNum[7] = {0};
	int pass[7];
	int p_value[7][10];
	double t[7];
	double2 d;
    char TestName[7][13] = {"NIST LC", "HSY", "DM1", "DM2", "HSY-CS7", "DM1-CS7", "DM2-CS7"}; 
    int numOfTimes = 100;
	int seed;
	for (int wai = 0; wai < numOfTimes; wai++) {
		/*Initialization*/
		fscanf(seed_fp, "%d,", &seed);
		srand(seed); 
		memset(pass, 0, sizeof(pass)); 
		memset(p_value, 0, sizeof(p_value));//assign all the values of the 2-dimension array p_value 0 
		int j ;
		for (int i = 0; i < samplesize; i++) {
            j = -1; 
            memset(buffer, 0, sizeof(buffer));
			for(int i1 = 0; i1 < length/8; i1++) {
                for(int j1 = 0; j1 < 8; j1++) {
                    buffer[i1] <<= 1;
                    if((unsigned int)rand() & mask[pos - 1]) {
                        buffer[i1] ^= 1;
                    }
                }
            }
			
			j ++;
			t[j] = NIST_Test(buffer, length);
			if (t[j] >= alpha) pass[j]++;
			if      (t[j] < 0.1)p_value[j][0]++;
			else if (t[j] < 0.2)p_value[j][1]++;
			else if (t[j] < 0.3)p_value[j][2]++;
			else if (t[j] < 0.4)p_value[j][3]++;
			else if (t[j] < 0.5)p_value[j][4]++;
			else if (t[j] < 0.6)p_value[j][5]++;
			else if (t[j] < 0.7)p_value[j][6]++;
			else if (t[j] < 0.8)p_value[j][7]++;
			else if (t[j] < 0.9)p_value[j][8]++;
			else                p_value[j][9]++;

			
			j ++;
			t[j] = Modified_HSY_Test(buffer, length, 6);
			if (t[j] >= alpha) pass[j]++;
			if      (t[j] < 0.1)p_value[j][0]++;
			else if (t[j] < 0.2)p_value[j][1]++;
			else if (t[j] < 0.3)p_value[j][2]++;
			else if (t[j] < 0.4)p_value[j][3]++;
			else if (t[j] < 0.5)p_value[j][4]++;
			else if (t[j] < 0.6)p_value[j][5]++;
			else if (t[j] < 0.7)p_value[j][6]++;
			else if (t[j] < 0.8)p_value[j][7]++;
			else if (t[j] < 0.9)p_value[j][8]++;
			else                p_value[j][9]++;

			j ++;
			d = Our_Test(buffer, length, 6);
			t[j] = d.d1;
			if (t[j] >= alpha) pass[j]++;
			if      (t[j] < 0.1)p_value[j][0]++;
			else if (t[j] < 0.2)p_value[j][1]++;
			else if (t[j] < 0.3)p_value[j][2]++;
			else if (t[j] < 0.4)p_value[j][3]++;
			else if (t[j] < 0.5)p_value[j][4]++;
			else if (t[j] < 0.6)p_value[j][5]++;
			else if (t[j] < 0.7)p_value[j][6]++;
			else if (t[j] < 0.8)p_value[j][7]++;
			else if (t[j] < 0.9)p_value[j][8]++;
			else                p_value[j][9]++;

			j ++;
			t[j] = d.d2;
			if (t[j] >= alpha) pass[j]++;
			if      (t[j] < 0.1)p_value[j][0]++;
			else if (t[j] < 0.2)p_value[j][1]++;
			else if (t[j] < 0.3)p_value[j][2]++;
			else if (t[j] < 0.4)p_value[j][3]++;
			else if (t[j] < 0.5)p_value[j][4]++;
			else if (t[j] < 0.6)p_value[j][5]++;
			else if (t[j] < 0.7)p_value[j][6]++;
			else if (t[j] < 0.8)p_value[j][7]++;
			else if (t[j] < 0.9)p_value[j][8]++;
			else                p_value[j][9]++;

			j ++;
			t[j] = HSYChiSquare7(buffer, length);
			if (t[j] >= alpha) pass[j]++;
			if      (t[j] < 0.1)p_value[j][0]++;
			else if (t[j] < 0.2)p_value[j][1]++;
			else if (t[j] < 0.3)p_value[j][2]++;
			else if (t[j] < 0.4)p_value[j][3]++;
			else if (t[j] < 0.5)p_value[j][4]++;
			else if (t[j] < 0.6)p_value[j][5]++;
			else if (t[j] < 0.7)p_value[j][6]++;
			else if (t[j] < 0.8)p_value[j][7]++;
			else if (t[j] < 0.9)p_value[j][8]++;
			else                p_value[j][9]++;
			
            j ++;
			t[j] = DM1ChiSquare7(buffer, length);
			if (t[j] >= alpha) pass[j]++;
			if      (t[j] < 0.1)p_value[j][0]++;
			else if (t[j] < 0.2)p_value[j][1]++;
			else if (t[j] < 0.3)p_value[j][2]++;
			else if (t[j] < 0.4)p_value[j][3]++;
			else if (t[j] < 0.5)p_value[j][4]++;
			else if (t[j] < 0.6)p_value[j][5]++;
			else if (t[j] < 0.7)p_value[j][6]++;
			else if (t[j] < 0.8)p_value[j][7]++;
			else if (t[j] < 0.9)p_value[j][8]++;
			else                p_value[j][9]++;
			
            j ++;
			t[j] = DM2ChiSquare7(buffer, length);
			if (t[j] >= alpha) pass[j]++;
			if      (t[j] < 0.1)p_value[j][0]++;
			else if (t[j] < 0.2)p_value[j][1]++;
			else if (t[j] < 0.3)p_value[j][2]++;
			else if (t[j] < 0.4)p_value[j][3]++;
			else if (t[j] < 0.5)p_value[j][4]++;
			else if (t[j] < 0.6)p_value[j][5]++;
			else if (t[j] < 0.7)p_value[j][6]++;
			else if (t[j] < 0.8)p_value[j][7]++;
			else if (t[j] < 0.9)p_value[j][8]++;
			else                p_value[j][9]++;
		}
		fprintf(fp, "Round  %d  has been completed!\n" , wai + 1);
		double KF[7];
        for (int i = 0; i < j + 1; i++) {
			KF[i] = 0;
			for (int j = 0; j < 10; j++) {
				fprintf(fp, "%4d  ", p_value[i][j]);
				KF[i] += (p_value[i][j] - samplesize / 10)*(p_value[i][j] - samplesize / 10);
			}
			KF[i] /= (samplesize / 10);
			double temp = cephes_igamc(9.0 / 2, KF[i] / 2);
			if (temp < alphau) RejectNum[i]++;  // temp is the p-value of p-values
			fprintf(fp, "%.6f  %d/%d\n", temp, pass[i], samplesize); 		
		}
		for (int i = 0; i < j + 1; i++)
			fprintf(fp, "%d  ", RejectNum[i]);
		fprintf(fp, "\n");
		fprintf(fp, "------------------------------------------------------------------------\n");
	}
	fclose(fp);
	fclose(seed_fp);
	return 0;
}
