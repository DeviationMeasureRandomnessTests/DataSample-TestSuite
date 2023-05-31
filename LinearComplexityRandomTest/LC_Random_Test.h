#ifndef LC_RANDOM_TEST
#define LC_RANDOM_TEST
#include "my_struct.h"
void Initialfenweidian();


double NISToriginal_Test(unsigned char *s, int n);
void LinearBlock_Test_num(unsigned char*s, int n, long long* N, long long* N1) ;
void Our_Test_num(unsigned char *s, long long n, int fw, long long *N, long long *N1, long long *N2);
void Modified_HSY_Test_num(unsigned char*s, int n, int fw, long long *N1_, long long*N2_); 
void NIST_Test_num(unsigned char *s, int n, long long *N, double **nu); 


double NIST_Test(unsigned char *s, int n);
double HSY_Test(unsigned char *s, int n, int fw);
double Modified_HSY_Test(unsigned char *s, int n, int fw);
double2 Our_Test(unsigned char *s, int n, int fw);
double DM1_Test(unsigned char *s, int n, int fw);
double DM2_Test(unsigned char *s, int n, int fw);
double LinearBlock_Test(unsigned char *s, int n);
double Jump_Test(unsigned char *s, int n, int fw);
double Ohs_Test(unsigned char *s, int n, int fw);




double DM1ChiSquare5(unsigned char *s, int n);
double DM1ChiSquare7(unsigned char *s, int n);
double DM2ChiSquare7(unsigned char *s, int n);
double HSYChiSquare7(unsigned char *s, int n);
double OhsChiSquare7(unsigned char *s, int n);
double EhsChiSquare7(unsigned char *s, int n);
double JumpChiSquare7(unsigned char *s, int n);

double JumpNmbers_Test(unsigned char *s, int n);
#endif
