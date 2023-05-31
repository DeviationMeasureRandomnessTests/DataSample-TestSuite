#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<math.h>
#include<string.h>
#define length 1000000  //10^6 bits for a Sequence
#define samplesize 1000
/*calculate the linear complexity in  bits*/
char B_[501];
char C_[501];
char T_[501];
int N_, m_, L_;
char d_;


int main() {
	FILE *fp1,*fp2,*fp3;
	char basename1[40] = "MT32T64PRGGAPA.data";
	fp1 = fopen(basename1, "rb");  //Read Only
	char basename2[40] = "MT32PRGValid.data";
	fp2 = fopen(basename2, "rb");  //Read Only
	char basename3[40] = "Mix66-MT&PLCPGAPA.data";
	fp3 = fopen(basename3, "wb+");  //Write
	
	int i;
	int Len = 66000000; //Number of bits to Handle  10^8 bits
	
	int FullLen=1000000000; //Total bits to Handle  10^9 bits<2^30
	 
	unsigned char buffer1, buffer2;
	
	for (i = 0; i < Len / 8; i++) 
	{
		fread(&buffer1, sizeof(unsigned char), 1, fp1);
		fread(&buffer2, sizeof(unsigned char), 1, fp2);
		//buffer1 = 0;	
		fwrite(&buffer1, sizeof(unsigned char), 1, fp3);
	}
	
	for (i = Len / 8; i < FullLen/8; i++) 
	{
		fread(&buffer2, sizeof(unsigned char), 1, fp2);
		fwrite(&buffer2, sizeof(unsigned char), 1, fp3);
	}
	
		
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	return 0;
}
