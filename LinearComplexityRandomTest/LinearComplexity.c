#include "LinearComplexity.h"
#include <stdlib.h>
unsigned char BytePos[8] = {
	0x80, 0x40, 0x20, 0x10,
	0x08, 0x04, 0x02, 0x01
};
int NumOneBits(unsigned n)
{
	n = n - ((n >> 1) & 0x55555555);
	n = (n & 0x33333333) + ((n >> 2) & 0x33333333);
	return (((n + (n >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

double Jueduizhi(double a) {
	if (a < 0)return -a;
	else return a;
}

/*calculate the linear complexity in  bits*/
int bitL(char *s, int n) {
	char *C, *B, *T;
	B = (char*)calloc(n + 1, sizeof(char));
	C = (char*)calloc(n + 1, sizeof(char));
	T = (char*)calloc(n + 1, sizeof(char));
	int N, m, L;
	char d;
    int i;
	m = -1, L = 0;
	B[0] = C[0] = T[0] = 1;
	for (i = 1; i <= n; i++)
		B[i] = C[i] = T[i] = 0;
	char *rs = (char*)calloc(n, sizeof(char));
	for (i = 0; i < n; i++)
		rs[i] = s[n - 1 - i];
	for (N = 0; N < n; N++){
		d = s[N];
		for (i = 1; i <= L; i++)d ^= C[i] & rs[n - 1 - (N - i)];
		if (d) {
			if (N + 1 - L > L) {
				for (i = 1; i <= L; i++)T[i] = C[i];
				for (i = N - m; i <= N + 1 - L; i++)C[i] ^= B[i - (N - m)];
				for (i = 1; i <= L; i++)B[i] = T[i];
				L = N + 1 - L;
				m = N;
			}
			else {
				for (i = N - m; i <= N + 1 - L; i++)C[i] ^= B[i - (N - m)];
			}
		}
	}
	free (B);
	free (C);
	free (T);
	free (rs);

	return L;
}

/*calculate the one order derivation sum in bits*/
double bitA(char *s, int n) {
	double result = 0.0;
	char *C, *B, *T;
	B = (char*)calloc(n + 1, sizeof(char));
	C = (char*)calloc(n + 1, sizeof(char));
	T = (char*)calloc(n + 1, sizeof(char));
	int N, m, L;
	char d;

	m = -1, L = 0;
	B[0] = C[0] = T[0] = 1;
	int i;
	for (i = 1; i <= n; i++)
		B[i] = C[i] = T[i] = 0;
	char *rs = (char*)calloc(n, sizeof(char));
	for (i = 0; i < n; i++)
		rs[i] = s[n - 1 - i];
	for (N = 0; N < n; N++){
		d = s[N];
		for (i = 1; i <= L; i++)d ^= C[i] & rs[n - 1 - (N - i)];
		if (d) {
			if (N + 1 - L > L) {
				for (i = 1; i <= L; i++)T[i] = C[i];
				for (i = N - m; i <= N + 1 - L; i++)C[i] ^= B[i - (N - m)];
				for (i = 1; i <= L; i++)B[i] = T[i];
				L = N + 1 - L;
				m = N;
			}
			else {
				for (i = N - m; i <= N + 1 - L; i++)C[i] ^= B[i - (N - m)];
			}
		}
		result += Jueduizhi(L - (N + 1) / 2.0);
	}
	free (B);
	free (C);
	free (T);
	free (rs);
	return result;
}

/*calculate the second order derivation sum in bits*/
double bitB(char *s, int n) {
	double result = 0.0;
	char *C, *B, *T;
	B = (char*)calloc(n + 1, sizeof(char));
	C = (char*)calloc(n + 1, sizeof(char));
	T = (char*)calloc(n + 1, sizeof(char));
	int N, m, L;
	char d;

	m = -1, L = 0;
	B[0] = C[0] = T[0] = 1;
	int i;
	for (i = 1; i <= n; i++)
		B[i] = C[i] = T[i] = 0;
	char *rs = (char*)calloc(n, sizeof(char));
	for (i = 0; i < n; i++)
		rs[i] = s[n - 1 - i];
	for (N = 0; N < n; N++){
		d = s[N];
		for (i = 1; i <= L; i++)d ^= C[i] & rs[n - 1 - (N - i)];
		if (d) {
			if (N + 1 - L > L) {
				for (i = 1; i <= L; i++)T[i] = C[i];
				for (i = N - m; i <= N + 1 - L; i++)C[i] ^= B[i - (N - m)];
				for (i = 1; i <= L; i++)B[i] = T[i];
				L = N + 1 - L;
				m = N;
			}
			else {
				for (i = N - m; i <= N + 1 - L; i++)C[i] ^= B[i - (N - m)];
			}
		}
		result += (L - (N + 1) / 2.0) * (L - (N + 1) / 2.0);
	}
	free (B);
	free (C);
	free (T);
	free (rs);
	return result;
}

/*calculate the one and two order derivation sum in bits in the same time*/
double2 bitAB(char *s, int n) {
	double2 result;
	result.d1=result.d2= 0.0;
	char *C, *B, *T;
	B = (char*)calloc(n + 1, sizeof(char));
	C = (char*)calloc(n + 1, sizeof(char));
	T = (char*)calloc(n + 1, sizeof(char));
	int N, m, L;
	char d;
    int i;
	m = -1, L = 0;
	B[0] = C[0] = T[0] = 1;
	for (i = 1; i <= n; i++)
		B[i] = C[i] = T[i] = 0;
	char *rs = (char*)calloc(n, sizeof(char));
	for (i = 0; i < n; i++)
		rs[i] = s[n - 1 - i];
	for (N = 0; N < n; N++){
		d = s[N];
		for (i = 1; i <= L; i++)d ^= C[i] & rs[n - 1 - (N - i)];
		if (d) {
			if (N + 1 - L > L) {
				for (i = 1; i <= L; i++)T[i] = C[i];
				for (i = N - m; i <= N + 1 - L; i++)C[i] ^= B[i - (N - m)];
				for (i = 1; i <= L; i++)B[i] = T[i];
				L = N + 1 - L;
				m = N;
			}
			else {
				for (i = N - m; i <= N + 1 - L; i++)C[i] ^= B[i - (N - m)];
			}
		}
		result.d1 += Jueduizhi(L - (N + 1) / 2.0);
		result.d2 += (L - (N + 1) / 2.0)*(L - (N + 1) / 2.0);
	}
	free (B);
	free (C);
	free (T);
	free (rs);
	return result;
}


/*calculate the linear complexity in unsigned int*/
int intL(char *s, int n) {
	unsigned *B, *C, *T;
	B = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	C = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	T = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	B[0] = C[0] = T[0] = (1 << 31);
	unsigned *RS = (unsigned *)calloc(n, sizeof(unsigned));

	RS[0] = s[n - 1];
	int i, j;
	for (j = 1; j < 32; j++){
		RS[0] <<= 1;
		if (n - 1 - j >= 0 && s[n - 1 - j])RS[0] ^= 1;
	}
	for (i = 1; i < n; i++){
		RS[i] = RS[i - 1] << 1;
		if (n - 1 - i - 31 >= 0 && s[n - 1 - i - 31])RS[i] ^= 1;
	}


	int N = 0, L = 0, m = -1;
	unsigned d;
	while (N < n) {
		d = 0;

		for ( i = 0; i <= L / 32; i++) 
            d ^= C[i] & RS[n - 1 - N + i * 32];
		
        if ((NumOneBits(d) & 1) == 1) {
			int ShiftB = (N - m) % 32;
			int StartC = (N - m) / 32;
			int StopC = (N + 1 - L) / 32;
			if (N + 1 - L > L) {
				for ( i = 0; i <= L / 32; i++)
                    T[i] = C[i];
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)
                        C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC+1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
				for ( i = 0; i <= L / 32; i++)
                    B[i] = T[i];
				L = N + 1 - L;
				m = N;
			}
			else {
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)
                        C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC + 1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
			}
		}
		N++;
	}
	free (B);
	free (C);
	free (T);
	free (RS);
	return L;
}
/*calculate the one order derivation sum in unsigned int*/
double intA(char *s, int n) {
	double result = 0.0;
	unsigned *B, *C, *T;
	B = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	C = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	T = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	B[0] = C[0] = T[0] = (1 << 31);
	unsigned *RS = (unsigned *)calloc(n, sizeof(unsigned));
    int i, j;
	RS[0] = s[n - 1];
	for (j = 1; j < 32; j++){
		RS[0] <<= 1;
		if (n - 1 - j >= 0 && s[n - 1 - j])RS[0] ^= 1;
	}
	for (i = 1; i < n; i++){
		RS[i] = RS[i - 1] << 1;
		if (n - 1 - i - 31 >= 0 && s[n - 1 - i - 31])RS[i] ^= 1;
	}


	int N = 0, L = 0, m = -1;
	unsigned d;
	while (N < n) {
		d = 0;
		for ( i = 0; i <= L / 32; i++)d ^= C[i] & RS[n - 1 - N + i * 32];
		if ((NumOneBits(d) & 1) == 1) {
			int ShiftB = (N - m) % 32;
			int StartC = (N - m) / 32;
			int StopC = (N + 1 - L) / 32;
			if (N + 1 - L > L) {
				for ( i = 0; i <= L / 32; i++)T[i] = C[i];
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC+1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
				for ( i = 0; i <= L / 32; i++)B[i] = T[i];
				L = N + 1 - L;
				m = N;
			}
			else {
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC + 1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
			}
		}
		N++;
		result += Jueduizhi(L-N/2.0);
	}
	free (B);
	free (C);
	free (T);
	free (RS);

	return result;
}
/*calculate the two order derivation sum in unsigned int*/
double intB(char *s, int n) {
	double result = 0.0;
	unsigned *B, *C, *T;
	B = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	C = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	T = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	B[0] = C[0] = T[0] = (1 << 31);
	unsigned *RS = (unsigned *)calloc(n, sizeof(unsigned));
    int i, j;
	RS[0] = s[n - 1];
	for (j = 1; j < 32; j++){
		RS[0] <<= 1;
		if (n - 1 - j >= 0 && s[n - 1 - j])RS[0] ^= 1;
	}
	for (i = 1; i < n; i++){
		RS[i] = RS[i - 1] << 1;
		if (n - 1 - i - 31 >= 0 && s[n - 1 - i - 31])RS[i] ^= 1;
	}


	int N = 0, L = 0, m = -1;
	unsigned d;
	while (N < n) {
		d = 0;
		for ( i = 0; i <= L / 32; i++)d ^= C[i] & RS[n - 1 - N + i * 32];
		if ((NumOneBits(d) & 1) == 1) {
			int ShiftB = (N - m) % 32;
			int StartC = (N - m) / 32;
			int StopC = (N + 1 - L) / 32;
			if (N + 1 - L > L) {
				for ( i = 0; i <= L / 32; i++)T[i] = C[i];
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC+1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
				for ( i = 0; i <= L / 32; i++)B[i] = T[i];
				L = N + 1 - L;
				m = N;
			}
			else {
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC + 1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
			}
		}
		N++;
		result += (L-N/2.0)*(L-N/2.0);
	}
	free (B);
	free (C);
	free (T);
	free (RS);
	return result;
}
/*calculate the one and two order derivation sum in unsigned int in the same time*/
double2 intAB(char *s, int n) {
	double2 result = {0.0, 0.0};
	unsigned *B, *C, *T;
	B = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	C = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	T = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	B[0] = C[0] = T[0] = (1 << 31);
	unsigned *RS = (unsigned *)calloc(n, sizeof(unsigned));
    int i, j;
	RS[0] = s[n - 1];
	for (j = 1; j < 32; j++){
		RS[0] <<= 1;
		if (n - 1 - j >= 0 && s[n - 1 - j])RS[0] ^= 1;
	}
	for (i = 1; i < n; i++){
		RS[i] = RS[i - 1] << 1;
		if (n - 1 - i - 31 >= 0 && s[n - 1 - i - 31])RS[i] ^= 1;
	}


	int N = 0, L = 0, m = -1;
	unsigned d;
	while (N < n) {
		d = 0;
		for ( i = 0; i <= L / 32; i++)d ^= C[i] & RS[n - 1 - N + i * 32];
		if ((NumOneBits(d) & 1) == 1) {
			int ShiftB = (N - m) % 32;
			int StartC = (N - m) / 32;
			int StopC = (N + 1 - L) / 32;
			if (N + 1 - L > L) {
				for ( i = 0; i <= L / 32; i++)T[i] = C[i];
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC+1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
				for ( i = 0; i <= L / 32; i++)B[i] = T[i];
				L = N + 1 - L;
				m = N;
			}
			else {
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC + 1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
			}
		}
		N++;
		result.d1 += Jueduizhi(L-N/2.0);
		result.d2 += (L-N/2.0)*(L-N/2.0);

	}
	free (B);
	free (C);
	free (T);
	free (RS);
	return result;
}


/*calculate linear complexity with the sequence stored in bytes*/
int byteL(unsigned char* S, int start, int n){
	unsigned  *B, *C, *T;
	B = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	C = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	T = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	B[0] = C[0] = T[0] = (1 << 31);
	unsigned *RS = (unsigned *)calloc(n, sizeof(unsigned));

	char *s = (char *)calloc(n, sizeof(char));
	int i, j;
	for (i = 0; i < n; i++){
		if (S[(start + i) >> 3] & BytePos[(start + i) & 7])s[i] = 1;
		else s[i] = 0;
	}

	RS[0] = s[n - 1];
	for (j = 1; j < 32; j++){
		RS[0] <<= 1;
		if (n - 1 - j >= 0 && s[n - 1 - j])RS[0] ^= 1;
	}
	for (i = 1; i < n; i++){
		RS[i] = RS[i - 1] << 1;
		if (n - 1 - i - 31 >= 0 && s[n - 1 - i - 31])RS[i] ^= 1;
	}


	int N = 0, L = 0, m = -1;
	unsigned d;
	while (N < n) {
		d = 0;
		for ( i = 0; i <= L / 32; i++)d ^= C[i] & RS[n - 1 - N + i * 32];
		if ((NumOneBits(d) & 1) == 1) {
			int ShiftB = (N - m) % 32;
			int StartC = (N - m) / 32;
			int StopC = (N + 1 - L) / 32;
			if (N + 1 - L > L) {
				for ( i = 0; i <= L / 32; i++)T[i] = C[i];
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC+1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
				for ( i = 0; i <= L / 32; i++)B[i] = T[i];
				L = N + 1 - L;
				m = N;
			}
			else {
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC + 1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
			}
		}
		N++;
	}
	free (B);
	free (C);
	free (T);
	free (RS);
	free (s);

	return L;

}

/*calculate one order deviation sum with the sequence stored in bytes*/
double byteA(unsigned char* S, int start, int n){
	double result=0.0;
	unsigned  *B, *C, *T;
	B = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	C = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	T = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	B[0] = C[0] = T[0] = (1 << 31);
	unsigned *RS = (unsigned *)calloc(n, sizeof(unsigned));

    int i, j;
	char *s = (char *)calloc(n, sizeof(char));
	for (i = 0; i < n; i++){
		if (S[(start + i) >> 3] & BytePos[(start + i) & 7])s[i] = 1;
		else s[i] = 0;
	}

	RS[0] = s[n - 1];
	for (j = 1; j < 32; j++){
		RS[0] <<= 1;
		if (n - 1 - j >= 0 && s[n - 1 - j])RS[0] ^= 1;
	}
	for (i = 1; i < n; i++){
		RS[i] = RS[i - 1] << 1;
		if (n - 1 - i - 31 >= 0 && s[n - 1 - i - 31])RS[i] ^= 1;
	}


	int N = 0, L = 0, m = -1;
	unsigned d;
	while (N < n) {
		d = 0;
		for ( i = 0; i <= L / 32; i++)d ^= C[i] & RS[n - 1 - N + i * 32];
		if ((NumOneBits(d) & 1) == 1) {
			int ShiftB = (N - m) % 32;
			int StartC = (N - m) / 32;
			int StopC = (N + 1 - L) / 32;
			if (N + 1 - L > L) {
				for ( i = 0; i <= L / 32; i++)T[i] = C[i];
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC+1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
				for ( i = 0; i <= L / 32; i++)B[i] = T[i];
				L = N + 1 - L;
				m = N;
			}
			else {
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC + 1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
			}
		}
		N++;
		result+=Jueduizhi(L-N/2.0);
	}
	free (B);
	free (C);
	free (T);
	free (RS);
	free (s);
	return result;

}

/*calculate the one order deviation sum and the linear complexity for the binary sequence stored in bytes,
the output will be written in wL and wA */
void byteLA(unsigned char* S, int start, int n, int *wL, double *wA) {
    double result = 0.0;
    unsigned  *B, *C, *T;
    B = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
    C = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
    T = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
    B[0] = C[0] = T[0] = (1 << 31);
    unsigned *RS = (unsigned *)calloc(n, sizeof(unsigned));

    int i, j;
    char *s = (char *)calloc(n, sizeof(char));
    for (i = 0; i < n; i++) {
        if (S[(start + i) >> 3] & BytePos[(start + i) & 7])s[i] = 1;
        else s[i] = 0;
    }

    RS[0] = s[n - 1];
    for (j = 1; j < 32; j++) {
        RS[0] <<= 1;
        if (n - 1 - j >= 0 && s[n - 1 - j])RS[0] ^= 1;
    }
    for (i = 1; i < n; i++) {
        RS[i] = RS[i - 1] << 1;
        if (n - 1 - i - 31 >= 0 && s[n - 1 - i - 31])RS[i] ^= 1;
    }


    int N = 0, L = 0, m = -1;
    unsigned d;
    while (N < n) {
        d = 0;
        for (i = 0; i <= L / 32; i++)d ^= C[i] & RS[n - 1 - N + i * 32];
        if ((NumOneBits(d) & 1) == 1) {
            int ShiftB = (N - m) % 32;
            int StartC = (N - m) / 32;
            int StopC = (N + 1 - L) / 32;
            if (N + 1 - L > L) {
                for (i = 0; i <= L / 32; i++)T[i] = C[i];
                if (ShiftB == 0) {
                    for (i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
                }
                else {
                    C[StartC] ^= B[0] >> ShiftB;
                    for (i = StartC + 1; i <= StopC; i++) {
                        C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
                        C[i] ^= B[i - StartC] >> ShiftB;
                    }
                }
                for (i = 0; i <= L / 32; i++)B[i] = T[i];
                L = N + 1 - L;
                m = N;
            }
            else {
                if (ShiftB == 0) {
                    for (i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
                }
                else {
                    C[StartC] ^= B[0] >> ShiftB;
                    for (i = StartC + 1; i <= StopC; i++) {
                        C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
                        C[i] ^= B[i - StartC] >> ShiftB;
                    }
                }
            }
        }
        N++;
        result += Jueduizhi(L - N / 2.0);
    }
    free(B);
    free(C);
    free(T);
    free(RS);
    free(s);
    *wA = result;
    *wL = L;
}


/*calculate two order deviation sum with the sequence stored in bytes*/
double byteB(unsigned char* S, int start, int n){
	double result=0.0;
	unsigned  *B, *C, *T;
	B = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	C = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	T = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	B[0] = C[0] = T[0] = (1 << 31);
	unsigned *RS = (unsigned *)calloc(n, sizeof(unsigned));

	int i, j;
	char *s = (char *)calloc(n, sizeof(char));
	for (i = 0; i < n; i++){
		if (S[(start + i) >> 3] & BytePos[(start + i) & 7])s[i] = 1;
		else s[i] = 0;
	}

	RS[0] = s[n - 1];
	for (j = 1; j < 32; j++){
		RS[0] <<= 1;
		if (n - 1 - j >= 0 && s[n - 1 - j])RS[0] ^= 1;
	}
	for (i = 1; i < n; i++){
		RS[i] = RS[i - 1] << 1;
		if (n - 1 - i - 31 >= 0 && s[n - 1 - i - 31])RS[i] ^= 1;
	}


	int N = 0, L = 0, m = -1;
	unsigned d;
	while (N < n) {
		d = 0;
		for ( i = 0; i <= L / 32; i++)d ^= C[i] & RS[n - 1 - N + i * 32];
		if ((NumOneBits(d) & 1) == 1) {
			int ShiftB = (N - m) % 32;
			int StartC = (N - m) / 32;
			int StopC = (N + 1 - L) / 32;
			if (N + 1 - L > L) {
				for ( i = 0; i <= L / 32; i++)T[i] = C[i];
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC+1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
				for ( i = 0; i <= L / 32; i++)B[i] = T[i];
				L = N + 1 - L;
				m = N;
			}
			else {
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC + 1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
			}
		}
		N++;
		result+=(L-N/2.0)*(L-N/2.0);
	}
	free (B);
	free (C);
	free (T);
	free (RS);
	free (s);
	return result;

}

/*calculate one and two order derivation sum with the sequence stored in bytes in the same time*/
double2 byteAB(unsigned char* S, int start, int n){
	double2 result={0.0, 0.0};
	unsigned  *B, *C, *T;
	B = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	C = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	T = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	B[0] = C[0] = T[0] = (1 << 31);
	unsigned *RS = (unsigned *)calloc(n, sizeof(unsigned));

	int i, j;
	char *s = (char *)calloc(n, sizeof(char));
	for (i = 0; i < n; i++){
		if (S[(start + i) >> 3] & BytePos[(start + i) & 7])s[i] = 1;
		else s[i] = 0;
	}

	RS[0] = s[n - 1];
	for (j = 1; j < 32; j++){
		RS[0] <<= 1;
		if (n - 1 - j >= 0 && s[n - 1 - j])RS[0] ^= 1;
	}
	for (i = 1; i < n; i++){
		RS[i] = RS[i - 1] << 1;
		if (n - 1 - i - 31 >= 0 && s[n - 1 - i - 31])RS[i] ^= 1;
	}


	int N = 0, L = 0, m = -1;
	unsigned d;
	while (N < n) {
		d = 0;
		for ( i = 0; i <= L / 32; i++)d ^= C[i] & RS[n - 1 - N + i * 32];
		if ((NumOneBits(d) & 1) == 1) {
			int ShiftB = (N - m) % 32;
			int StartC = (N - m) / 32;
			int StopC = (N + 1 - L) / 32;
			if (N + 1 - L > L) {
				for ( i = 0; i <= L / 32; i++)T[i] = C[i];
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC+1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
				for ( i = 0; i <= L / 32; i++)B[i] = T[i];
				L = N + 1 - L;
				m = N;
			}
			else {
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC + 1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
			}
		}
		N++;
		result.d1+=Jueduizhi(L-N/2.0);
		result.d2+=(L-N/2.0)*(L-N/2.0);
	}
	free (B);
	free (C);
	free (T);
	free (RS);
	free (s);
	return result;

}


/*translate a bit string into a byte string from the position of 'start' for n bits*/
unsigned char* BoolArrayToByteArray(char s[], int start, int numOfBits){
	int size, i, j;
	size = (numOfBits - 1) / 8 + 1;
	unsigned char *t = (unsigned char*)calloc (size, sizeof(unsigned char));
	for (i = 0; i < size - 1; i++){
		for (j = 0; j < 8; j++){
			t[i] <<= 1;
			if (s[start + i * 8 + j])
				t[i] ^= 1;
		}
	}
	for (j = 0; j < 8; j++){
		t[i] <<= 1;
		if (8 * i + j < numOfBits && s[start + i * 8 + j])
			t[i] ^= 1;
	}
	return t;
}
/*translate a byte string into a bit string from the position of 'start' in the
format of bit for totally n bits*/
char * ByteArrayToBoolArray(unsigned char s[], int start, int numOfBits){
	char *t = (char *)calloc(numOfBits, sizeof(char));
	int i;
	for (i = 0; i < numOfBits; i++){
		t[i] = 0;
		if (s[(start + i) / 8] & (1 << (7 - (start + i) % 8)))
			t[i] = 1;
	}
	return t;
}

/*calculate the jump times with the sequence  stored in bytes*/
int byteJumps(unsigned char* S, int start, int n){
	int result = 0;
	unsigned  *B, *C, *T;
	B = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	C = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	T = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	B[0] = C[0] = T[0] = (1 << 31);
	unsigned *RS = (unsigned *)calloc(n, sizeof(unsigned));

	int i, j;
	char *s = (char *)calloc(n, sizeof(char));
	for (i = 0; i < n; i++){
		if (S[(start + i) >> 3] & BytePos[(start + i) & 7])s[i] = 1;
		else s[i] = 0;
	}

	RS[0] = s[n - 1];
	for (j = 1; j < 32; j++){
		RS[0] <<= 1;
		if (n - 1 - j >= 0 && s[n - 1 - j])RS[0] ^= 1;
	}
	for (i = 1; i < n; i++){
		RS[i] = RS[i - 1] << 1;
		if (n - 1 - i - 31 >= 0 && s[n - 1 - i - 31])RS[i] ^= 1;
	}


	int N = 0, L = 0, m = -1;
	unsigned d;
	while (N < n) {
		d = 0;
		for ( i = 0; i <= L / 32; i++)d ^= C[i] & RS[n - 1 - N + i * 32];
		if ((NumOneBits(d) & 1) == 1) {
			int ShiftB = (N - m) % 32;
			int StartC = (N - m) / 32;
			int StopC = (N + 1 - L) / 32;
			if (N + 1 - L > L) {
				for ( i = 0; i <= L / 32; i++)T[i] = C[i];
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC+1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
				for ( i = 0; i <= L / 32; i++)B[i] = T[i];
				L = N + 1 - L;
				m = N;
				result++;
			}
			else {
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC + 1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
			}
		}
		N++;
	}
	free (B);
	free (C);
	free (T);
	free (RS);
	free (s);
	return result;

}
/* calculate the odd hop sum of a binary sequence with the format of 8 bits one byte */
int byteOddHopSum(unsigned char* S, int start, int n){
	int result = 0;
    int odd = 1;
	unsigned  *B, *C, *T;
	B = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	C = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	T = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	B[0] = C[0] = T[0] = (1 << 31);
	unsigned *RS = (unsigned *)calloc(n, sizeof(unsigned));

	int i, j;
	char *s = (char *)calloc(n, sizeof(char));
	for (i = 0; i < n; i++){
		if (S[(start + i) >> 3] & BytePos[(start + i) & 7])s[i] = 1;
		else s[i] = 0;
	}

	RS[0] = s[n - 1];
	for (j = 1; j < 32; j++){
		RS[0] <<= 1;
		if (n - 1 - j >= 0 && s[n - 1 - j])RS[0] ^= 1;
	}
	for (i = 1; i < n; i++){
		RS[i] = RS[i - 1] << 1;
		if (n - 1 - i - 31 >= 0 && s[n - 1 - i - 31])RS[i] ^= 1;
	}


	int N = 0, L = 0, m = -1;
	unsigned d;
	while (N < n) {
		d = 0;
		for ( i = 0; i <= L / 32; i++)d ^= C[i] & RS[n - 1 - N + i * 32];
		if ((NumOneBits(d) & 1) == 1) {
			int ShiftB = (N - m) % 32;
			int StartC = (N - m) / 32;
			int StopC = (N + 1 - L) / 32;
			if (N + 1 - L > L) {
				for ( i = 0; i <= L / 32; i++)T[i] = C[i];
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC+1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
				for ( i = 0; i <= L / 32; i++)B[i] = T[i];
				if(odd==1){
                    result += (N + 1 - 2*L);
                    odd = 0;
                }
                else
                    odd = 1;
                L = N + 1 - L;
				m = N;
			}
			else {
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC + 1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
			}
		}
		N++;
	}
	free (B);
	free (C);
	free (T);
	free (RS);
	free (s);
	return result;
}

/* calculate the even hop sum of a binary sequence with the format of 8 bits one byte */
int byteEvenHopSum(unsigned char* S, int start, int n){
	int result = 0;
    int even = 0;
	unsigned  *B, *C, *T;
	B = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	C = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	T = (unsigned *)calloc(n / 32 + 1, sizeof(unsigned));
	B[0] = C[0] = T[0] = (1 << 31);
	unsigned *RS = (unsigned *)calloc(n, sizeof(unsigned));
	char *s = (char *)calloc(n, sizeof(char));
	int i, j;
	for (i = 0; i < n; i++){
		if (S[(start + i) >> 3] & BytePos[(start + i) & 7])s[i] = 1;
		else s[i] = 0;
	}

	RS[0] = s[n - 1];
	for (j = 1; j < 32; j++){
		RS[0] <<= 1;
		if (n - 1 - j >= 0 && s[n - 1 - j])RS[0] ^= 1;
	}
	for (i = 1; i < n; i++){
		RS[i] = RS[i - 1] << 1;
		if (n - 1 - i - 31 >= 0 && s[n - 1 - i - 31])RS[i] ^= 1;
	}


	int N = 0, L = 0, m = -1;
	unsigned d;
	while (N < n) {
		d = 0;
		for ( i = 0; i <= L / 32; i++)d ^= C[i] & RS[n - 1 - N + i * 32];
		if ((NumOneBits(d) & 1) == 1) {
			int ShiftB = (N - m) % 32;
			int StartC = (N - m) / 32;
			int StopC = (N + 1 - L) / 32;
			if (N + 1 - L > L) {
				for ( i = 0; i <= L / 32; i++)T[i] = C[i];
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC+1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
				for ( i = 0; i <= L / 32; i++)B[i] = T[i];
				if(even==1){
                    result += (N + 1 - 2*L);
                    even = 0;
                }
                else
                    even = 1;
                L = N + 1 - L;
				m = N;
			}
			else {
				if (ShiftB == 0) {
					for ( i = StartC; i <= StopC; i++)C[i] ^= B[i - StartC];
				}
				else {
					C[StartC] ^= B[0] >> ShiftB;
					for ( i = StartC + 1; i <= StopC; i++) {
						C[i] ^= B[i - StartC - 1] << (32 - ShiftB);
						C[i] ^= B[i - StartC] >> ShiftB;
					}
				}
			}
		}
		N++;
	}
	free (B);
	free (C);
	free (T);
	free (RS);
	free (s);
	return result;
}

/*calculate the jump times with the sequence stored in bits*/
int bitJumps(char *s, int n) {
	int Jump = 0;
	char *C, *B, *T;
	B = (char*)calloc(n + 1, sizeof(char));
	C = (char*)calloc(n + 1, sizeof(char));
	T = (char*)calloc(n + 1, sizeof(char));
	int N, m, L;
	char d;

	m = -1, L = 0;
	B[0] = C[0] = T[0] = 1;
	char *rs = (char*)malloc(n*sizeof(char));
	int i;
	for (i = 0; i < n; i++)
		rs[i] = s[n - 1 - i];
	for (N = 0; N < n; N++){
		d = s[N];
		for (i = 1; i <= L; i++)d ^= C[i] & rs[n - 1 - (N - i)];
		if (d) {
			if (N + 1 - L > L) {
				for (i = 1; i <= L; i++)T[i] = C[i];
				for (i = N - m; i <= N + 1 - L; i++)C[i] ^= B[i - (N - m)];
				for (i = 1; i <= L; i++)B[i] = T[i];
				L = N + 1 - L;
				m = N;
				Jump++;
			}
			else {
				for (i = N - m; i <= N + 1 - L; i++)C[i] ^= B[i - (N - m)];
			}
		}
	}
	free (B);
	free (C);
	free (T);
	free (rs);

	return Jump;
}

/*calculate the odd hop sum with the sequence stored in bits*/
int bitOddHopSum(char *s, int n) {
	int Ohs = 0;
	int odd = 1;
	char *C, *B, *T;
	B = (char*)calloc(n + 1, sizeof(char));
	C = (char*)calloc(n + 1, sizeof(char));
	T = (char*)calloc(n + 1, sizeof(char));
	int N, m, L;
	char d;

	m = -1, L = 0;
	B[0] = C[0] = T[0] = 1;
	char *rs = (char*)malloc(n*sizeof(char));

	int i;
	for (i = 0; i < n; i++)
		rs[i] = s[n - 1 - i];
	for (N = 0; N < n; N++){
		d = s[N];
		for (i = 1; i <= L; i++)d ^= C[i] & rs[n - 1 - (N - i)];
		if (d) {
			if (N + 1 - L > L) {
				for (i = 1; i <= L; i++)T[i] = C[i];
				for (i = N - m; i <= N + 1 - L; i++)C[i] ^= B[i - (N - m)];
				for (i = 1; i <= L; i++)B[i] = T[i];
				if(odd==1){
                    Ohs += (N+1-L-L);
				}
                odd = 1-odd;
				L = N + 1 - L;
				m = N;

			}
			else {
				for (i = N - m; i <= N + 1 - L; i++)C[i] ^= B[i - (N - m)];
			}
		}
	}
	free (B);
	free (C);
	free (T);
	free (rs);

	return Ohs;
}

/*calculate the even hop sum with the sequence stored in bits*/
int bitEvenHopSum(char *s, int n) {
	int Ehs = 0;
	int even = 0;
	char *C, *B, *T;
	B = (char*)calloc(n + 1, sizeof(char));
	C = (char*)calloc(n + 1, sizeof(char));
	T = (char*)calloc(n + 1, sizeof(char));
	int N, m, L;
	char d;

	m = -1, L = 0;
	B[0] = C[0] = T[0] = 1;
	char *rs = (char*)malloc(n*sizeof(char));

	int i;
	for (i = 0; i < n; i++)
		rs[i] = s[n - 1 - i];
	for (N = 0; N < n; N++){
		d = s[N];
		for (i = 1; i <= L; i++)d ^= C[i] & rs[n - 1 - (N - i)];
		if (d) {
			if (N + 1 - L > L) {
				for (i = 1; i <= L; i++)T[i] = C[i];
				for (i = N - m; i <= N + 1 - L; i++)C[i] ^= B[i - (N - m)];
				for (i = 1; i <= L; i++)B[i] = T[i];
				if(even==1){
                    Ehs += (N+1-L-L);
				}
                even = 1 - even;
				L = N + 1 - L;
				m = N;

			}
			else {
				for (i = N - m; i <= N + 1 - L; i++)C[i] ^= B[i - (N - m)];
			}
		}
	}
	free (B);
	free (C);
	free (T);
	free (rs);

	return Ehs;
}

