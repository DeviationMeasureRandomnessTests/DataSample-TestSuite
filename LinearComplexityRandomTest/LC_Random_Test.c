#include "LC_Random_Test.h"
#include "LinearComplexity.h"
#include "my_struct.h"
#include "cephes.h"
#include <math.h>
#include <stdio.h>
static const int M = 500;
//FILE *fp1, *fp2, *fp3, *fp4;
/*-------------------the percentiles of DM-1-----------------------*/
double T_Hamo1[13] = { 480, 441, 423, 403, 390, 379,369, 360, 351, 341, 327, 317, 299 };//have limit in the end
double fenwei_Hamo1[13] = { 0.00993521, 0.0496111, 0.0985057, 0.196668, 0.292549, 0.393766,
						0.498258, 0.597295, 0.694731, 0.792997, 0.899212, 0.94834, 0.989752 };
double T_DM1[13] = { 482, 443, 425, 404, 391, 380, 371, 362, 352, 342, 328, 318, 300 };//no limit in the end
double fenwei_DM1[13] = { 0.00985309, 0.0489237, 0.0969601, 0.199704, 0.296058, 0.397407,
						0.490929, 0.589356, 0.697183, 0.794679, 0.899932, 0.948615, 0.989736 };
/*-------------------the percentiles of DM-2-----------------------*/
double T_Hamo2[13] = { 1044.5, 825.5, 738.5, 653.5, 601.5, 562.5, 529.5, 500.5, 471.5, 440.5, 402.5, 375.5, 331.5 };//have limit in the end
double fenwei_Hamo2[13] = { 0.00999305, 0.0496213, 0.0998859, 0.199246, 0.299574, 0.399401, 0.499772, 0.596826, 0.696678
, 0.797886, 0.899286, 0.948936, 0.989596 };
double T_DM2[13] = { 1051.5, 830.5, 743.5, 657.5, 605.5, 566.5, 533.5, 503.5, 473.5, 443.5, 405.5, 377.5, 333.5 };//no limit in the end
double fenwei_DM2[13] = { 0.00997701, 0.0497821, 0.099758, 0.199665, 0.299349, 0.398346, 0.497813, 0.597438, 0.699957, 0.797106, 0.89811, 0.949325, 0.989571 };
/*-------------------jump times---------------------*/
int  J_T[13] = { 143, 138, 135, 131, 129, 127, 125, 123, 121, 118, 115, 112, 107 };
double j_fenwei[13] = { 0.01057644, 0.04752037, 0.09877715, 0.21730577, 0.29879308, 0.39185401,
             0.49156801, 0.59181476, 0.68637532, 0.80670430, 0.89367072, 0.94814558, 0.98818247 };
/*--------------------odd hop sum------------------*/
double T_ohs[13] = {145, 140, 137, 133, 131, 129, 127, 125, 122, 120, 116, 114, 108};
double fenwei_ohs[13]={0.00613880, 0.03044762, 0.06724289, 0.16027063, 0.22912656, 0.31203323, 0.40577431, 0.50530516,
    0.65212319, 0.73990461, 0.87441753, 0.91925262, 0.98458395 };


/*-----------------------------------------------------------------------*/
double pi_1[5]={0.205321,0.194332,0.192944,0.201227,0.206175};

/*----------------------------------------------------------------------------------------------------
--------------------------------------卡方测试所分区间------------------------------------------------
-----------------------------------------------------------------------------------------------------*/

double kafang_fwd_Japan[6];
double pi_Japan[7];

double kafang_fwd_DM1[6];
double pi_DM1[7];

double kafang_fwd_DM2[6];
double pi_DM2[7];

double kafang_fwd_Ohs[6];
double pi_Ohs[7];

double kafang_fwd_Ehs[6];
double pi_Ehs[7];

double kafang_fwd_Jump[6];
double pi_Jump[7];
/*----------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------*/
void Initialfenweidian() {
    FILE *fp = fopen("fwd_Japan_13.data", "rb");
    fread(T_Hamo1, sizeof(double), 13, fp);
    fread(fenwei_Hamo1, sizeof(double), 13, fp);

    fp = freopen("fwd_DM1_13.data", "rb", fp);
    fread(T_DM1, sizeof(double), 13, fp);
    fread(fenwei_DM1, sizeof(double), 13, fp);

    fp = freopen("fwd_DM2_13.data", "rb", fp);
    fread(T_DM2, sizeof(double), 13, fp);
    fread(fenwei_DM2, sizeof(double), 13, fp);

    fp = freopen("fwd_Ohs_13.data", "rb", fp);
    fread(T_ohs, sizeof(double), 13, fp);
    fread(fenwei_ohs, sizeof(double), 13, fp);

    fp = freopen("japanPercentile7.data", "rb", fp);
    fread(kafang_fwd_Japan, sizeof(double), 6, fp);
    fread(pi_Japan, sizeof(double), 7, fp);

    fp = freopen("dm1Percentile7.data", "rb", fp);
    fread(kafang_fwd_DM1, sizeof(double), 6, fp);
    fread(pi_DM1, sizeof(double), 7, fp);

    fp = freopen("dm2Percentile7.data", "rb", fp);
    fread(kafang_fwd_DM2, sizeof(double), 6, fp);
    fread(pi_DM2, sizeof(double), 7, fp);

    fp = freopen("ohsPercentile7.data", "rb", fp);
    fread(kafang_fwd_Ohs, sizeof(double), 6, fp);
    fread(pi_Ohs, sizeof(double), 7, fp);

    fp = freopen("ehsPercentile7.data", "rb", fp);
    fread(kafang_fwd_Ehs, sizeof(double), 6, fp);
    fread(pi_Ehs, sizeof(double), 7, fp);

    fp = freopen("jumpPercentile7.data", "rb", fp);
    fread(kafang_fwd_Jump, sizeof(double), 6, fp);
    fread(pi_Jump, sizeof(double), 7, fp);

    fclose(fp);
}




void LinearBlock_Test_num(unsigned char*s, int n, long long* N, long long* N1) {
	long long block_num = n / M;
	long long temp_N1 = 0;
	#pragma omp parallel
	{
		#pragma omp for reduction(+:temp_N1)
		for (long long i = 0; i < block_num ; i++) {
			if (byteL(s, M*i, M) == M / 2)
				temp_N1++;
		}
	}
	*N += block_num;
	*N1 += temp_N1;
}



double NISToriginal_Test(unsigned char *s, int n)
{
	int       i, L, N,  parity,  sign;
    int	       K = 6;
	double    p_value, nu[7], chi2,  mean,  T_;
	double    pi[7] = { 0.010417, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833 };

	for(i=0;  i<K+1;  i++)
        nu[i] = 0.00;
	N = n / M;

	for (i = 0; i < N; i++) {

		L = byteL(s , M*i, M);
        if((parity = (M+1)%2) == 0)
            sign = -1;
        else
            sign = 1;
		mean = M / 2.0 + (9.0 + sign) / 36.0 -1.0/pow(2,  M)* (M / 3.0 + 2.0 / 9.0);
		if((parity = M%2) == 0)
            sign = 1;
        else
            sign = -1;
		T_ =sign* (L - mean) + 2.0 / 9.0;

		if (T_ <= -2.5)
			nu[0]++;
		else if (T_ > -2.5 && T_ <= -1.5)
			nu[1]++;
		else if (T_ > -1.5 && T_ <= -0.5)
			nu[2]++;
		else if (T_ > -0.5 && T_ <= 0.5)
			nu[3]++;
		else if (T_ > 0.5 && T_ <= 1.5)
			nu[4]++;
		else if (T_ > 1.5 && T_ <= 2.5)
			nu[5]++;
		else
			nu[6]++;
	}

	chi2 = 0.00;
    for(i=0;  i<K+1;  i++)
        chi2 += pow(nu[i] - N*pi[i], 2)/(N*pi[i]);
	p_value = cephes_igamc(K / 2.0, chi2 / 2.0);
	//fwrite(&p_value, sizeof(double), 1, fp1);

	return p_value;
}


void NIST_Test_num(unsigned char *s, int n, long long *N, double **nu)
{
	int       block_Num;
    int	       K = 6;
	double    p_value, nu1, nu2, nu3, nu4, nu5, nu6, nu0, chi2;
	double    pi[7] = { 0.010417, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833 };

	nu0 = 0.00;
	nu1 = 0.00;
	nu2 = 0.00;
	nu3 = 0.00;
	nu4 = 0.00;
	nu5 = 0.00;
	nu6 = 0.00;
	block_Num = n / M;
	#pragma omp parallel
	{
		#pragma omp for reduction(+:nu0,nu1,nu2,nu3,nu4,nu5,nu6)
		for (int i = 0; i < block_Num ; i++) {

			int L = byteL(s , M*i, M);
			double mean = M / 2.0 + (9.0 -1) / 36.0 - (M / 3.0 + 2.0 / 9.0) / pow(2, M);

			double T_ = (L - mean) + 2.0 / 9.0;

			if (T_ <= -2.5)
				nu0+=1;
			else if (T_ > -2.5 && T_ <= -1.5)
				nu1+=1;
			else if (T_ > -1.5 && T_ <= -0.5)
				nu2+=1;
			else if (T_ > -0.5 && T_ <= 0.5)
				nu3+=1;
			else if (T_ > 0.5 && T_ <= 1.5)
				nu4+=1;
			else if (T_ > 1.5 && T_ <= 2.5)
				nu5+=1;
			else
				nu6+=1;
		}
	}
    *N += block_Num ;

    (*nu)[0] += nu0;
    (*nu)[1] += nu1;
    (*nu)[2] += nu2;
    (*nu)[3] += nu3;
    (*nu)[4] += nu4;
    (*nu)[5] += nu5;
    (*nu)[6] += nu6;
}
double HSY_Test(unsigned char *s, int n, int fw) {
	double p_value;
	int N = n / M, N1 = 0, N2 = 0;
	#pragma omp parallel
	{
		#pragma omp for reduction(+: N1, N2)
		for (int i = 0; i < N; i++) {
            int L;
            double A;
            byteLA(s, M*i, M, &L, &A);
            if (L != M / 2)continue;
			if (A > T_Hamo1[fw])N2 += 1;
			N1+=1;
		}
	}


	double p0 = N1*1.0/N;
	if(Jueduizhi(p0 - 0.5) > 1.5/sqrt(N)){

		return 0;
	}
	//cout << N << "," << N1 << "," << N1 << endl;
	//if (N1*Phi < 5)cout << "Too little samples!" << endl;
	double p1 = (N2*1.0) / N1;
	double z1 = (p1 - fenwei_Hamo1[fw]) / sqrt(fenwei_Hamo1[fw]*(1 - fenwei_Hamo1[fw]) / N1);
	//fwrite(&z1, sizeof(double), 1, fp2);
	p_value = cephes_erfc(Jueduizhi(z1) / sqrt(2));
	//p_value = cephes_normal(z1);

	return p_value;

}


void Modified_HSY_Test_num(unsigned char*s, int n, int fw, long long *N1_, long long*N2_) {
	int N = n / M, N1 = 0, N2 = 0;
	#pragma omp parallel
	{
		#pragma omp for reduction(+: N1, N2)
		for (int i = 0; i < N; i++) {
            int L;
            double A;
            byteLA(s, M*i, M, &L, &A);
            if (L != M / 2)continue;
            if (A > T_Hamo1[fw])N2 += 1;
            N1 += 1;
		}
	}

    *N1_ += N1;
    *N2_ += N2;

	//fwrite(&z1, sizeof(double), 1, fp2);

}


void Our_Test_num(unsigned char *s, long long n, int fw, long long *N, long long *N1, long long *N2) {
	long long block_num = n / M;
	long long N11 = 0, N12 = 0;

	#pragma omp parallel
	{
		#pragma omp for reduction(+:N11, N12)
		for (long long i = 0; i < block_num;  i++) {
			double2 D = byteAB(s , M*i, M);
			if (D.d1 > T_DM1[fw])  N11 += 1;
			if (D.d2 > T_DM2[fw])  N12 += 1;
		}
	}
	*N += block_num;
	*N1 += N11;
	*N2 += N12;
}

//按：此测试方法与其他大有不同，以其未尝分块也。然则序列之长，BM算法耗时亦久矣
double JumpNmbers_Test(unsigned char *s, int n) {
	double temp;
	double mu, sigma;
	double numJ;
	double z1;
	double p_value;
	int Parite;
	if (n >= 1024)
		temp = 0.0;
	else
		temp = pow(2.0, -(double)n);
	Parite = n & 1;
	mu = n / 4.0 + (4 + Parite) / 12.0 - temp / 3.0;
	sigma = n / 8.0 - (2 - Parite) / (9.0 - Parite) + n * temp / 6.0
		+ (6 + Parite) * temp / 18.0 - temp * temp / 9.0;
	sigma = sqrt(sigma);
	numJ = byteJumps(s, 0, n);
	z1 = (numJ - mu) / sigma;
	p_value = cephes_erfc(Jueduizhi(z1) / sqrt(2));
	//p_value = cephes_normal(z1);
    return p_value;
}

/*-------------------------------------------------------------------------------------------------------
-----------------------------------------只利用一个分位点计算p-value--------------------------------------
--------------------------------------------------------------------------------------------------------*/
double LinearBlock_Test(unsigned char*s, int n) {
    double p_value;
    int N = n / M, N1 = 0;
#pragma omp parallel
    {
#pragma omp for reduction(+:N1)
        for (int i = 0; i < N; i++) {
            if (byteL(s, M*i, M) == M / 2)
                N1++;
        }
    }
    double p0 = N1*1.0 / N;
    double z1 = 2 * (p0 - 0.5) * sqrt(1.0*N);
    p_value = cephes_erfc(Jueduizhi(z1) / sqrt(2));
    //p_value = cephes_normal(z1);
    return p_value;
}

double Modified_HSY_Test(unsigned char*s, int n, int fw) {
	double p_value;
	int N = n / M, N1 = 0, N2 = 0;
	#pragma omp parallel
	{
		#pragma omp for reduction(+: N1, N2)
		for (int i = 0; i < N; i++) {
            int L;
            double A;
            byteLA(s, M*i, M, &L, &A);
            if (L != M / 2)continue;
            if (A > T_Hamo1[fw])N2 += 1;
            N1 += 1;
		}
	}


	double p1 = (N2*1.0) / N1;
	double z1 = (p1 - fenwei_Hamo1[fw]) / sqrt(fenwei_Hamo1[fw]*(1 - fenwei_Hamo1[fw]) / N1);
	p_value = cephes_erfc(Jueduizhi(z1) / sqrt(2));
	//p_value = cephes_normal(z1);
	//fwrite(&z1, sizeof(double), 1, fp2);

	return p_value;
}

double2 Our_Test(unsigned char *s, int n, int fw) {
    double2 p_value;
    int N = n / M;
    int N11 = 0, N12 = 0;

#pragma omp parallel
    {
#pragma omp for reduction(+:N11,N12)
        for (int i = 0; i < N; i++) {
            double2 D = byteAB(s, M*i, M);
            if (D.d1 > T_DM1[fw])  N11 += 1;
            if (D.d2 > T_DM2[fw])  N12 += 1;
        }
    }


    //if (N1*Phi < 5)cout << "Too little samples!"<<endl;
    //cout << N << "," << N1 << endl;
    double p1 = N11*1.0 / N;
    double z1 = (p1 - fenwei_DM1[fw]) / sqrt(fenwei_DM1[fw] * (1 - fenwei_DM1[fw]) / N);
    //fwrite(&z1, sizeof(double), 1, fp3);
    p_value.d1 = cephes_erfc(Jueduizhi(z1) / sqrt(2));
    //p_value.d1 = cephes_normal(z1);
    double p2 = N12*1.0 / N;
    double z2 = (p2 - fenwei_DM2[fw]) / sqrt(fenwei_DM2[fw] * (1 - fenwei_DM2[fw]) / N);
    //fwrite(&z2, sizeof(double), 1, fp4);
    p_value.d2 = cephes_erfc(Jueduizhi(z2) / sqrt(2));
    //p_value.d2 = cephes_normal(z2);
    return p_value;
}

double Ohs_Test(unsigned char *s, int n, int fw) {
    double p_value;
    int N = n / M;
    int N1 = 0;


#pragma omp parallel
    {
#pragma omp for reduction(+:N1)
        for (int i = 0; i < N; i++) {
            if (byteOddHopSum(s, M*i, M) > T_ohs[fw])
                N1++;
        }
    }
    //if (N1*Phi < 5)cout << "Too little samples!"<<endl;
    //cout << N << "," << N1 << endl;
    double p1 = N1*1.0 / N;
    double z1 = (p1 - fenwei_ohs[fw]) / sqrt(fenwei_ohs[fw] * (1 - fenwei_ohs[fw]) / N);
    p_value = cephes_erfc(Jueduizhi(z1) / sqrt(2));
    //p_value = cephes_normal(z1);
    return p_value;
}

double Jump_Test(unsigned char *s, int n, int fw) {
    double p_value;
    int N = n / M;
    int N1 = 0;

#pragma omp parallel
    {
#pragma omp for reduction(+:N1)
        for (int i = 0; i < N; i++) {
            if (byteJumps(s, M*i, M) > J_T[fw])
                N1++;
        }
    }

    //if (N1*Phi < 5)cout << "Too little samples!"<<endl;
    //cout << N << "," << N1 << endl;
    double p1 = N1*1.0 / N;
    double z1 = (p1 - j_fenwei[fw]) / sqrt(j_fenwei[fw] * (1 - j_fenwei[fw]) / N);
    p_value = cephes_erfc(Jueduizhi(z1) / sqrt(2));
    //p_value = cephes_normal(z1);
    return p_value;
}

double DM1_Test(unsigned char *s, int n, int fw) {
    double   D;
    int N = n / M;
    int N1 = 0;
    int i;
#pragma omp parallel
    {
#pragma omp for reduction(+:N1)
        for (i = 0; i < N; i++) {
            D = byteA(s, M*i, M);
            if (D > T_DM1[fw])  N1++;
        }
    }
    //if (N1*Phi < 5)cout << "Too little samples!"<<endl;
    //cout << N << "," << N1 << endl;
    double p = N1*1.0 / N;
    double z = (p - fenwei_DM1[fw]) / sqrt(fenwei_DM1[fw] * (1 - fenwei_DM1[fw]) / N);
    //fwrite(&z1, sizeof(double), 1, fp3);
    double p_value = cephes_erfc(Jueduizhi(z) / sqrt(2));
    //p_value.d1 = cephes_normal(z1);
    return p_value;
}


double DM2_Test(unsigned char *s, int n, int fw) {
    double   D;
    int N = n / M;
    int N1 = 0;
    int i;
#pragma omp parallel
    {
#pragma omp for reduction(+:N1)
        for (i = 0; i < N; i++) {
            D = byteB(s, M*i, M);
            if (D > T_DM2[fw])  N1++;
        }
    }
    //if (N1*Phi < 5)cout << "Too little samples!"<<endl;
    //cout << N << "," << N1 << endl;
    double p = N1*1.0 / N;
    double z = (p - fenwei_DM2[fw]) / sqrt(fenwei_DM2[fw] * (1 - fenwei_DM2[fw]) / N);
    //fwrite(&z1, sizeof(double), 1, fp3);
    double p_value = cephes_erfc(Jueduizhi(z) / sqrt(2));
    //p_value.d1 = cephes_normal(z1);
    return p_value;
}

/*-------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------------------------------------
------------------------------------分多个区间利用卡方统计量得出p-value-----------------------------------
--------------------------------------------------------------------------------------------------------*/
double NIST_Test(unsigned char *s, int n)
{
    int       N;
    int	       K = 6;
    double    p_value, nu1, nu2, nu3, nu4, nu5, nu6, nu0, chi2;
    double    pi[7] = { 0.010417, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833 };

    nu0 = 0.00;
    nu1 = 0.00;
    nu2 = 0.00;
    nu3 = 0.00;
    nu4 = 0.00;
    nu5 = 0.00;
    nu6 = 0.00;
    N = n / M;
#pragma omp parallel
    {
#pragma omp for reduction(+:nu0,nu1,nu2,nu3,nu4,nu5,nu6)
        for (int i = 0; i < N; i++) {

            int L = byteL(s, M*i, M);
            double mean = M / 2.0 + (9.0 - 1) / 36.0 - (M / 3.0 + 2.0 / 9.0) / pow(2, M);

            double T_ = (L - mean) + 2.0 / 9.0;

            if (T_ <= -2.5)
                nu0 += 1;
            else if (T_ > -2.5 && T_ <= -1.5)
                nu1 += 1;
            else if (T_ > -1.5 && T_ <= -0.5)
                nu2 += 1;
            else if (T_ > -0.5 && T_ <= 0.5)
                nu3 += 1;
            else if (T_ > 0.5 && T_ <= 1.5)
                nu4 += 1;
            else if (T_ > 1.5 && T_ <= 2.5)
                nu5 += 1;
            else
                nu6 += 1;
        }
    }
    chi2 = 0.00;
    chi2 += pow(nu0 - N*pi[0], 2) / (N*pi[0]);
    chi2 += pow(nu1 - N*pi[1], 2) / (N*pi[1]);
    chi2 += pow(nu2 - N*pi[2], 2) / (N*pi[2]);
    chi2 += pow(nu3 - N*pi[3], 2) / (N*pi[3]);
    chi2 += pow(nu4 - N*pi[4], 2) / (N*pi[4]);
    chi2 += pow(nu5 - N*pi[5], 2) / (N*pi[5]);
    chi2 += pow(nu6 - N*pi[6], 2) / (N*pi[6]);
    p_value = cephes_igamc(K / 2.0, chi2 / 2.0);
    //fwrite(&p_value, sizeof(double), 1, fp1);

    return p_value;
}

double DM1ChiSquare5(unsigned char *s, int n) {
    int       N;
    int	      K = 4;
    double    p_value, nu1, nu2, nu3, nu4, nu0, chi2;

    nu0 = 0.00;
    nu1 = 0.00;
    nu2 = 0.00;
    nu3 = 0.00;
    nu4 = 0.00;

    N = n / M;
#pragma omp parallel
    {
#pragma omp for reduction(+:nu0,nu1,nu2,nu3,nu4)
        for (int i = 0; i < N; i++) {

            double T_ = byteA(s, M*i, M);

            if (T_ < 343)
                nu0 += 1;
            else if (T_ < 362)
                nu1 += 1;
            else if (T_ < 380)
                nu2 += 1;
            else if (T_ < 404)
                nu3 += 1;
            else
                nu4 += 1;

        }
    }
    chi2 = 0.00;
    chi2 += pow(nu0 - N*pi_1[0], 2) / (N*pi_1[0]);
    chi2 += pow(nu1 - N*pi_1[1], 2) / (N*pi_1[1]);
    chi2 += pow(nu2 - N*pi_1[2], 2) / (N*pi_1[2]);
    chi2 += pow(nu3 - N*pi_1[3], 2) / (N*pi_1[3]);
    chi2 += pow(nu4 - N*pi_1[4], 2) / (N*pi_1[4]);

    p_value = cephes_igamc(K / 2.0, chi2 / 2.0);
    //fwrite(&p_value, sizeof(double), 1, fp1);

    return p_value;
}

double DM1ChiSquare7(unsigned char *s, int n) {
    int       N;
    int	      K = 6;
    double    p_value, nu0, nu1, nu2, nu3, nu4, nu5, nu6, chi2;

    nu0 = 0.00;
    nu1 = 0.00;
    nu2 = 0.00;
    nu3 = 0.00;
    nu4 = 0.00;
    nu5 = 0.00;
    nu6 = 0.00;
    N = n / M;
#pragma omp parallel
    {
#pragma omp for reduction(+:nu0,nu1,nu2,nu3,nu4, nu5, nu6)
        for (int i = 0; i < N; i++) {

            double T_ = byteA(s, M*i, M);

            if (T_ < kafang_fwd_DM1[0])
                nu0 += 1;
            else if (T_ < kafang_fwd_DM1[1])
                nu1 += 1;
            else if (T_ < kafang_fwd_DM1[2])
                nu2 += 1;
            else if (T_ < kafang_fwd_DM1[3])
                nu3 += 1;
            else if (T_ < kafang_fwd_DM1[4])
                nu4 += 1;
            else if (T_ < kafang_fwd_DM1[5])
                nu5 += 1;
            else
                nu6 += 1;

        }
    }
    chi2 = 0.00;
    chi2 += pow(nu0 - N*pi_DM1[0], 2) / (N*pi_DM1[0]);
    chi2 += pow(nu1 - N*pi_DM1[1], 2) / (N*pi_DM1[1]);
    chi2 += pow(nu2 - N*pi_DM1[2], 2) / (N*pi_DM1[2]);
    chi2 += pow(nu3 - N*pi_DM1[3], 2) / (N*pi_DM1[3]);
    chi2 += pow(nu4 - N*pi_DM1[4], 2) / (N*pi_DM1[4]);
    chi2 += pow(nu5 - N*pi_DM1[5], 2) / (N*pi_DM1[5]);
    chi2 += pow(nu6 - N*pi_DM1[6], 2) / (N*pi_DM1[6]);



    p_value = cephes_igamc(K / 2.0, chi2 / 2.0);
    //fwrite(&p_value, sizeof(double), 1, fp1);

    return p_value;
}

double DM2ChiSquare7(unsigned char *s, int n) {
    int       N;
    int	      K = 6;
    double    p_value, nu0, nu1, nu2, nu3, nu4, nu5, nu6, chi2;

    nu0 = 0.00;
    nu1 = 0.00;
    nu2 = 0.00;
    nu3 = 0.00;
    nu4 = 0.00;
    nu5 = 0.00;
    nu6 = 0.00;
    N = n / M;
#pragma omp parallel
    {
#pragma omp for reduction(+:nu0,nu1,nu2,nu3,nu4, nu5, nu6)
        for (int i = 0; i < N; i++) {

            double T_ = byteB(s, M*i, M);

            if (T_ < kafang_fwd_DM2[0])
                nu0 += 1;
            else if (T_ < kafang_fwd_DM2[1])
                nu1 += 1;
            else if (T_ < kafang_fwd_DM2[2])
                nu2 += 1;
            else if (T_ < kafang_fwd_DM2[3])
                nu3 += 1;
            else if (T_ < kafang_fwd_DM2[4])
                nu4 += 1;
            else if (T_ < kafang_fwd_DM2[5])
                nu5 += 1;
            else
                nu6 += 1;

        }
    }
    chi2 = 0.00;
    chi2 += pow(nu0 - N*pi_DM2[0], 2) / (N*pi_DM2[0]);
    chi2 += pow(nu1 - N*pi_DM2[1], 2) / (N*pi_DM2[1]);
    chi2 += pow(nu2 - N*pi_DM2[2], 2) / (N*pi_DM2[2]);
    chi2 += pow(nu3 - N*pi_DM2[3], 2) / (N*pi_DM2[3]);
    chi2 += pow(nu4 - N*pi_DM2[4], 2) / (N*pi_DM2[4]);
    chi2 += pow(nu5 - N*pi_DM2[5], 2) / (N*pi_DM2[5]);
    chi2 += pow(nu6 - N*pi_DM2[6], 2) / (N*pi_DM2[6]);



    p_value = cephes_igamc(K / 2.0, chi2 / 2.0);
    //fwrite(&p_value, sizeof(double), 1, fp1);

    return p_value;
}


double HSYChiSquare7(unsigned char *s, int n) {
    int       N,N1=0;
    int	      K = 6;
    double    p_value, nu0, nu1, nu2, nu3, nu4, nu5, nu6, chi2;

    nu0 = 0.00;
    nu1 = 0.00;
    nu2 = 0.00;
    nu3 = 0.00;
    nu4 = 0.00;
    nu5 = 0.00;
    nu6 = 0.00;
    N = n / M;
    #pragma omp parallel
    {
    #pragma omp for reduction(+:N1,nu0,nu1,nu2,nu3,nu4, nu5, nu6)
        for (int i = 0; i < N; i++) {
            int L;
            double T_;
            byteLA(s, M*i, M, &L, &T_);
            if (L != M / 2)continue;
            N1 += 1;


            if (T_ < kafang_fwd_Japan[0])
                nu0 += 1;
            else if (T_ < kafang_fwd_Japan[1])
                nu1 += 1;
            else if (T_ < kafang_fwd_Japan[2])
                nu2 += 1;
            else if (T_ < kafang_fwd_Japan[3])
                nu3 += 1;
            else if (T_ < kafang_fwd_Japan[4])
                nu4 += 1;
            else if (T_ < kafang_fwd_Japan[5])
                nu5 += 1;
            else
                nu6 += 1;

        }
    }
    chi2 = 0.00;
    chi2 += pow(nu0 - N1*pi_Japan[0], 2) / (N1*pi_Japan[0]);
    chi2 += pow(nu1 - N1*pi_Japan[1], 2) / (N1*pi_Japan[1]);
    chi2 += pow(nu2 - N1*pi_Japan[2], 2) / (N1*pi_Japan[2]);
    chi2 += pow(nu3 - N1*pi_Japan[3], 2) / (N1*pi_Japan[3]);
    chi2 += pow(nu4 - N1*pi_Japan[4], 2) / (N1*pi_Japan[4]);
    chi2 += pow(nu5 - N1*pi_Japan[5], 2) / (N1*pi_Japan[5]);
    chi2 += pow(nu6 - N1*pi_Japan[6], 2) / (N1*pi_Japan[6]);



    p_value = cephes_igamc(K / 2.0, chi2 / 2.0);
    //fwrite(&p_value, sizeof(double), 1, fp1);

    return p_value;
}

double OhsChiSquare7(unsigned char *s, int n) {
    int       N;
    int	      K = 6;
    double    p_value, nu0, nu1, nu2, nu3, nu4, nu5, nu6, chi2;

    nu0 = 0.00;
    nu1 = 0.00;
    nu2 = 0.00;
    nu3 = 0.00;
    nu4 = 0.00;
    nu5 = 0.00;
    nu6 = 0.00;
    N = n / M;
#pragma omp parallel
    {
#pragma omp for reduction(+:nu0,nu1,nu2,nu3,nu4, nu5, nu6)
        for (int i = 0; i < N; i++) {

            double T_ = byteOddHopSum(s, M*i, M);

            if (T_ < kafang_fwd_Ohs[0])
                nu0 += 1;
            else if (T_ < kafang_fwd_Ohs[1])
                nu1 += 1;
            else if (T_ < kafang_fwd_Ohs[2])
                nu2 += 1;
            else if (T_ < kafang_fwd_Ohs[3])
                nu3 += 1;
            else if (T_ < kafang_fwd_Ohs[4])
                nu4 += 1;
            else if (T_ < kafang_fwd_Ohs[5])
                nu5 += 1;
            else
                nu6 += 1;

        }
    }
    chi2 = 0.00;
    chi2 += pow(nu0 - N*pi_Ohs[0], 2) / (N*pi_Ohs[0]);
    chi2 += pow(nu1 - N*pi_Ohs[1], 2) / (N*pi_Ohs[1]);
    chi2 += pow(nu2 - N*pi_Ohs[2], 2) / (N*pi_Ohs[2]);
    chi2 += pow(nu3 - N*pi_Ohs[3], 2) / (N*pi_Ohs[3]);
    chi2 += pow(nu4 - N*pi_Ohs[4], 2) / (N*pi_Ohs[4]);
    chi2 += pow(nu5 - N*pi_Ohs[5], 2) / (N*pi_Ohs[5]);
    chi2 += pow(nu6 - N*pi_Ohs[6], 2) / (N*pi_Ohs[6]);



    p_value = cephes_igamc(K / 2.0, chi2 / 2.0);
    //fwrite(&p_value, sizeof(double), 1, fp1);

    return p_value;
}

double EhsChiSquare7(unsigned char *s, int n) {
    int       N;
    int	      K = 6;
    double    p_value, nu0, nu1, nu2, nu3, nu4, nu5, nu6, chi2;

    nu0 = 0.00;
    nu1 = 0.00;
    nu2 = 0.00;
    nu3 = 0.00;
    nu4 = 0.00;
    nu5 = 0.00;
    nu6 = 0.00;
    N = n / M;
#pragma omp parallel
    {
#pragma omp for reduction(+:nu0,nu1,nu2,nu3,nu4, nu5, nu6)
        for (int i = 0; i < N; i++) {

            double T_ = byteEvenHopSum(s, M*i, M);

            if (T_ < kafang_fwd_Ehs[0])
                nu0 += 1;
            else if (T_ < kafang_fwd_Ehs[1])
                nu1 += 1;
            else if (T_ < kafang_fwd_Ehs[2])
                nu2 += 1;
            else if (T_ < kafang_fwd_Ehs[3])
                nu3 += 1;
            else if (T_ < kafang_fwd_Ehs[4])
                nu4 += 1;
            else if (T_ < kafang_fwd_Ehs[5])
                nu5 += 1;
            else
                nu6 += 1;

        }
    }
    chi2 = 0.00;
    chi2 += pow(nu0 - N*pi_Ehs[0], 2) / (N*pi_Ehs[0]);
    chi2 += pow(nu1 - N*pi_Ehs[1], 2) / (N*pi_Ehs[1]);
    chi2 += pow(nu2 - N*pi_Ehs[2], 2) / (N*pi_Ehs[2]);
    chi2 += pow(nu3 - N*pi_Ehs[3], 2) / (N*pi_Ehs[3]);
    chi2 += pow(nu4 - N*pi_Ehs[4], 2) / (N*pi_Ehs[4]);
    chi2 += pow(nu5 - N*pi_Ehs[5], 2) / (N*pi_Ehs[5]);
    chi2 += pow(nu6 - N*pi_Ehs[6], 2) / (N*pi_Ehs[6]);



    p_value = cephes_igamc(K / 2.0, chi2 / 2.0);
    //fwrite(&p_value, sizeof(double), 1, fp1);

    return p_value;
}


double JumpChiSquare7(unsigned char *s, int n) {
    int       N;
    int	      K = 6;
    double    p_value, nu0, nu1, nu2, nu3, nu4, nu5, nu6, chi2;

    nu0 = 0.00;
    nu1 = 0.00;
    nu2 = 0.00;
    nu3 = 0.00;
    nu4 = 0.00;
    nu5 = 0.00;
    nu6 = 0.00;
    N = n / M;
#pragma omp parallel
    {
#pragma omp for reduction(+:nu0,nu1,nu2,nu3,nu4, nu5, nu6)
        for (int i = 0; i < N; i++) {

            double T_ = byteJumps(s, M*i, M);

            if (T_ < kafang_fwd_Jump[0])
                nu0 += 1;
            else if (T_ < kafang_fwd_Jump[1])
                nu1 += 1;
            else if (T_ < kafang_fwd_Jump[2])
                nu2 += 1;
            else if (T_ < kafang_fwd_Jump[3])
                nu3 += 1;
            else if (T_ < kafang_fwd_Jump[4])
                nu4 += 1;
            else if (T_ < kafang_fwd_Jump[5])
                nu5 += 1;
            else
                nu6 += 1;

        }
    }
    chi2 = 0.00;
    chi2 += pow(nu0 - N*pi_Jump[0], 2) / (N*pi_Jump[0]);
    chi2 += pow(nu1 - N*pi_Jump[1], 2) / (N*pi_Jump[1]);
    chi2 += pow(nu2 - N*pi_Jump[2], 2) / (N*pi_Jump[2]);
    chi2 += pow(nu3 - N*pi_Jump[3], 2) / (N*pi_Jump[3]);
    chi2 += pow(nu4 - N*pi_Jump[4], 2) / (N*pi_Jump[4]);
    chi2 += pow(nu5 - N*pi_Jump[5], 2) / (N*pi_Jump[5]);
    chi2 += pow(nu6 - N*pi_Jump[6], 2) / (N*pi_Jump[6]);



    p_value = cephes_igamc(K / 2.0, chi2 / 2.0);
    //fwrite(&p_value, sizeof(double), 1, fp1);

    return p_value;
}

/*-------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------*/
