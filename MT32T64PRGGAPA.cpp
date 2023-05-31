#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <bitset>
#include <stdint.h>

using namespace std;

enum
{
    W = 32,
    N = 624,
    M = 397,
    R = 31,
    A = 0x9908B0DF,

    F = 1812433253,

    U = 11,
    

    S = 7,
    B = 0x9D2C5680,//2636928640

    T = 15,
    C = 0xEFC60000,//4022730752

    L = 18,

    MASK_LOWER = 0x7fffffff,
    MASK_UPPER = 0x80000000
};

bool isInit;
static int index = N + 1;
static uint32_t MT[N];  //624 * 32 - 31 = 19937

void srand(const uint32_t seed)
{
    uint32_t  i, t;
    isInit = 1;
    MT[0] = seed;
    //Initialization
    for (i = 1; i < N; i++)
    {
        t = F * (MT[i - 1] ^ (MT[i - 1] >> (W - 2))) + i;
        MT[i] = t & 0xffffffff;   
    }
    index = N;
}

void generate()
{
    uint32_t  x, xA;
    for (int i = 0; i < N; i++)
    {
        // 2^31 = 0x80000000
        // 2^31-1 = 0x7fffffff
        x = (MT[i] & MASK_UPPER) + (MT[(i + 1) % N] & MASK_LOWER);
        xA = x >> 1;
        if (x & 0x1)
        {
            xA ^= A;
        }
        MT[i] = MT[(i + M) % N] ^ xA;
    }
    index = 0;
}

uint32_t Rand()
{
    if (!isInit)
        srand(5489);
    if (index >= N)
    {
        generate();
    }
    uint32_t y;
    y = MT[index];
    y = y ^ (y >> U);                
    y = y ^ ((y << S) & B);   
    y = y ^ ((y << T) & C); 
    y = y ^ (y >> L);                
    index = index + 1;
    return y;
}
int main()
{
    srand(20230526);  //Set a seed
    int cnt = 0;
    uint32_t res;
   
   
    FILE *fp1;
	char basename1[20] = "MT32T64PRGGAPA.data";
	fp1 = fopen(basename1, "wb+");  
	
	int  SeqLength,Blocklength;
	SeqLength=1000000;//10^6
	Blocklength=32;
	int i,j,k,l;
	
	int S1[32]; // A Seed block	
	int S2[32]; // A Seed block	
	int S1E[64]; // Expanded block
	
	unsigned char Seq[8]; //8 Bytes corresponding to S[64]
	
 for(l=0;l<15625000;l++)  //10^9/64 
 {	
  //S1---------------------------------------Confirmed Correctness for the first Blocklength=500 bits!
  
        res = Rand();
		//printf("%d:%u\n",i, res);  //Test Code
		for(k=0;k<32;k++)
		{
			S1[k]= (res & (1<<31) )>>31;  // A boolean value
			res=res<<1;
		}
	
	if(l%125==0) //Every 8000 bits=16*500bits
	{
	
						      
    // SeedBlock 32---> Block 64    Perfect LCPs
    S1E[0]=1; 
	for(i=1;i<Blocklength;i++)
	{
		S1E[2*i-1]=S1[i-1];
		S1E[2*i]=(S1E[2*i-1]+S1E[i-1])%2;
	}
	    S1E[2*Blocklength-1]=S1[Blocklength-1];
	}
	
	else
	{
		
		
		res = Rand();
		//printf("%d:%u\n",i, res);  //Test Code
		for(k=0;k<32;k++)
		{
			S2[k]= (res & (1<<31) )>>31;  // A boolean value
			res=res<<1;
		}
		//S1+S2-->S1E
		for(i=0;i<32;i++)
		{
			S1E[i]=S1[i];
			S1E[i+32]=S2[i];		
		}
	}
	    

  
 //Convert to Seq in bytes	
	for(j=0;j<8;j++)
	{
		Seq[j]=(S1E[8*j]<<7)+(S1E[8*j+1]<<6)+(S1E[8*j+2]<<5)+ (S1E[8*j+3]<<4)+ (S1E[8*j+4]<<3) + (S1E[8*j+5]<<2) + (S1E[8*j+6]<<1) +S1E[8*j+7];
	}
	
	 fwrite(Seq, sizeof(unsigned char), 8, fp1);  
}

   
   /* for (int i = 0; i < 1000; i++)//Test code for comparison
    {
        res = Rand();
        printf("%d:%u\n",i, res);
    }*/ 
    /*
0:1791095845
1:4282876139
2:3093770124
3:4005303368
4:491263
622:3906854836
623:2006116153
624:1104314680
625:939235918
997:2752667134
998:978824302
999:548926898*/
    
    	fclose(fp1); 
	    return 0;
	    
}
