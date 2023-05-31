#ifndef _LinearComplexity_H_
#define _LinearComplexity_H_
#include "my_struct.h"
int  NumOneBits(unsigned i);
double Jueduizhi(double a);


unsigned char* BoolArrayToByteArray(char s[], int start, int n);

char * ByteArrayToBoolArray(unsigned char s[], int start, int n);


int bitL(char *s, int n);
int intL(char *s, int n);
int byteL(unsigned char* S, int start, int n);

double bitA(char *s, int n);
double intA(char *s, int n);
double byteA(unsigned char* S, int start, int n);
void byteLA(unsigned char *S, int start, int n, int *wL, double *wA);

double bitB(char *s, int n);
double intB(char *s, int n);
double byteB(unsigned char* S, int start, int n);

double2 bitAB(char *s, int n);
double2 intAB(char *s, int n);
double2 byteAB(unsigned char* S, int start, int n);

int byteJumps(unsigned char*S, int start, int n);
int byteOddHopSum(unsigned char* S, int start, int n);
int byteEvenHopSum(unsigned char* S, int start, int n);

int bitJumps(char*s, int n);
int bitOddHopSum(char*s, int n);
int bitEvenHopSum(char*s, int n);
#endif
