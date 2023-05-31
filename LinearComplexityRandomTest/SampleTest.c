#include "LinearComplexity.h"
#include "cephes.h"
#include "LC_Random_Test.h"
#include "my_struct.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<math.h>
#define length 1000000
#define samplesize 1000
unsigned int next1 = 1;
void my_srand(unsigned short seed){
	next1 = (unsigned int)seed;
}
char  my_rand(){
	next1 = next1 * 214013L + 2531011L;
	if(next1&0x00010000)
		    return 1;
	else    return 0;
}
unsigned char buffer[length/8];
const double alpha = 0.01;
const double alphau = 0.0001;
const int FW = 12;
int main(int argc, char** argv) {
	FILE *fp, *datafp;
	char writename[100];
	strcpy(writename, argv[1]);
	fp = fopen(strcat(writename, ".txt"), "w");
	int RejectNum[14]={0};
	/*generate file names*/
	for (int wai = 0; wai < 1; wai++) {
		datafp = fopen(argv[1], "rb");
		/*Initialization*/ 
		int pass[7]={0};
		int p_value[7][10]={0};
		double t[7];
		double2 d;
		for (int i = 0; i < samplesize; i++) {
			
			fread(buffer, sizeof(unsigned char), length/8, datafp);
			
			
			int j ;
			j = 0;
			/*the linear complexity random test included in NIST test suite*/
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
			else if (t[j] < 1)  p_value[j][9]++;

			
			j = 1;
			/*HSY test, which has been skipped steps H4 and H5 to generate a P-value validly*/
			t[j] = Modified_HSY_Test(buffer, length, 1);
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
			else if (t[j] < 1)  p_value[j][9]++;


			j = 2;
			/*Our DM-1 test*/
			d = Our_Test(buffer, length, 6);
			t[j]=d.d1;
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
			else if (t[j] < 1)  p_value[j][9]++;

			j = 3;
			/*Our DM-2 test*/
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
			else if (t[j] < 1)  p_value[j][9]++;

		}
		fclose(datafp);
		fprintf(fp, "Round  %d  has been completed!\n" , wai + 1);
		double KF[7];
		for (int i = 0; i < 4; i++) {
			KF[i] = 0;
			for (int j = 0; j < 10; j++) {
				KF[i] += (p_value[i][j] - samplesize / 10)*(p_value[i][j] - samplesize / 10);
			}
			KF[i] /= (samplesize / 10);
			/*the P-value of the 1000 generated P-values*/
			double temp = cephes_igamc(9.0 / 2, KF[i] / 2);
			if (temp < alphau) RejectNum[i]++;
			fprintf(fp, "%.6e  %d\n", temp, pass[i]);	
		}
		for (int i = 0; i < 4; i++)
			fprintf(fp, "%d  ", RejectNum[i]);
		fprintf(fp, "\n");
		fprintf(fp, "------------------------------------------------------------------------\n");
	}
	fclose(fp);
	return 0;
}
