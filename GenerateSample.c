#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#define length 1000000
#define samplesize 1000
/*calculate the linear complexity in bits*/
char B_[501];
char C_[501];
char T_[501];
int N_, m_, L_;
char d_;

const int Len = length*samplesize;
char s[length*samplesize];
void OneRound(char *q) {
	d_ = q[N_];
	for (int i = 1; i <= L_; i++)d_ ^= C_[i] & q[N_ - i];
	if (d_) {
		if (N_ + 1 - L_ > L_) {
			for (int i = 1; i <= L_; i++)T_[i] = C_[i];
			for (int i = N_ - m_; i <= N_ + 1 - L_; i++)C_[i] ^= B_[i - (N_ - m_)];
			for (int i = 1; i <= L_; i++)B_[i] = T_[i];
			L_ = N_ + 1 - L_;
			m_ = N_;
		}
		else {
			for (int i = N_ - m_; i <= N_ + 1 - L_; i++)C_[i] ^= B_[i - (N_ - m_)];
		}
	}
	N_++;
}
/*pull down the point that is above the line y=0.5x*/
void AdjustOneRound(char *q) {
	d_ = q[N_];
	for (int i = 1; i <= L_; i++)d_ ^= C_[i] & q[N_ - i];
	if (d_) {
		/*for (int i = 1; i <= L_; i++)T_[i] = C_[i];
		for (int i = N_ - m_; i <= N_ + 1 - L_; i++)C_[i] ^= B_[i - (N_ - m_)];
		for (int i = 1; i <= L_; i++)B_[i] = T_[i];
		L_ = N_ + 1 - L_;
		m_ = N_;*/
		q[N_] ^= 1;
	}
	N_++;
}
/*push up the point that is below the line y=0.5x*/
void ReverseAdjustOneRound(char*q) {
	d_ = q[N_];
	for (int i = 1; i <= L_; i++)d_ ^= C_[i] & q[N_ - i];
	if (d_) {
		if (N_ + 1 - L_ > L_) {
			for (int i = 1; i <= L_; i++)T_[i] = C_[i];
			for (int i = N_ - m_; i <= N_ + 1 - L_; i++)C_[i] ^= B_[i - (N_ - m_)];
			for (int i = 1; i <= L_; i++)B_[i] = T_[i];
			L_ = N_ + 1 - L_;
			m_ = N_;
		}
		else {
			for (int i = N_ - m_; i <= N_ + 1 - L_; i++)C_[i] ^= B_[i - (N_ - m_)];
		}
	}
	else {
		q[N_] ^= 1;
		for (int i = 1; i <= L_; i++)T_[i] = C_[i];
		for (int i = N_ - m_; i <= N_ + 1 - L_; i++)C_[i] ^= B_[i - (N_ - m_)];
		for (int i = 1; i <= L_; i++)B_[i] = T_[i];
		L_ = N_ + 1 - L_;
		m_ = N_;
	}
	N_++;
}
/*
 * [n1] is the number of PCT T1s, [n2] is the number of PCT T2s, 
 * mm is the horizontal length of the big PCT,  [gap] is the tail 
 * length.
 *
 */
void MakeSequence(int n1, int n2, int mm, int gap) {
	FILE *fillfp = fopen("rdseed.data", "rb");
	unsigned char readbuffer[63], temp;
	int latestpoint;
	char q[500];
	int wai = 0;
	int frontlen = 500 - n1 * 2 - n2 * 4 - 2 * mm - gap;


	fread(readbuffer, sizeof(char), 63, fillfp);
	for (int i = 0; i < 62; i++) {
		temp = readbuffer[i];
		for (int j = 0; j<8; j++) {
			if (temp & 0x80)
				q[8 * i + j] = 1;
			else q[8 * i + j] = 0;
			temp <<= 1;
		}
	}
	temp = readbuffer[62];
	for (int j = 0; j < 4; j++) {
		if (temp & 0x80)
			q[8 * 62 + j] = 1;
		else q[8 * 62 + j] = 0;
		temp <<= 1;
	}

	/*initialize the Berlekamp-Massey algorithm*/
	m_ = -1, L_ = 0;
	memset(B_, 0, 501);
	memset(C_, 0, 501);
	memset(T_, 0, 501);
	B_[0] = C_[0] = T_[0] = 1;
	for (N_ = 0; N_ < frontlen; N_++) {
		d_ = q[N_];
		for (int i = 1; i <= L_; i++)d_ ^= C_[i] & q[N_ - i];
		if (d_) {
			if (N_ + 1 - L_ > L_) {
				for (int i = 1; i <= L_; i++)T_[i] = C_[i];
				for (int i = N_ - m_; i <= N_ + 1 - L_; i++)C_[i] ^= B_[i - (N_ - m_)];
				for (int i = 1; i <= L_; i++)B_[i] = T_[i];
				L_ = N_ + 1 - L_;
				m_ = N_;
			}
			else {
				for (int i = N_ - m_; i <= N_ + 1 - L_; i++)C_[i] ^= B_[i - (N_ - m_)];
			}
		}
		if (L_ == (N_ + 1) / 2)
			latestpoint = N_ + 1;
	}

	N_ = frontlen / 2 + latestpoint / 2 - 1;
	q[N_] ^= 1;

	int rounds = frontlen - N_;
	for (int i = 0; i < rounds; i++)OneRound(q);

	for (int i = 0; i < n2; i++) {
		AdjustOneRound(q);
		ReverseAdjustOneRound(q);
		OneRound(q);
		OneRound(q);
	}
	for (int i = 0; i < n1; i++) {
		ReverseAdjustOneRound(q);
		OneRound(q);
	}


	for (int i = 0; i < mm - 1; i++) {
		d_ = 0;
		for (int i = 1; i <= L_; i++)d_ ^= C_[i] & q[N_ - i];
		q[N_] = d_;
		N_++;
	}
	ReverseAdjustOneRound(q);


	for (int i = 0; i < 500; i++)s[wai * 500 + i] = q[i];

	for (wai = 1; wai< Len / 500; wai++) {
		fread(readbuffer, sizeof(char), 63, fillfp);
		for (int i = 0; i < 62; i++) {
			temp = readbuffer[i];
			for (int j = 0; j<8; j++) {
				if (temp & 0x80)
					q[8 * i + j] = 1;
				else q[8 * i + j] = 0;
				temp <<= 1;
			}
		}
		temp = readbuffer[62];
		for (int j = 0; j < 4; j++) {
			if (temp & 0x80)
				q[8 * 62 + j] = 1;
			else q[8 * 62 + j] = 0;
			temp <<= 1;
		}
		
        m_ = -1, L_ = 0;
		memset(B_, 0, 501);
		memset(C_, 0, 501);
		memset(T_, 0, 501);
		B_[0] = C_[0] = T_[0] = 1;
		for (N_ = 0; N_ < frontlen; N_++) {
			d_ = q[N_];
			for (int i = 1; i <= L_; i++)d_ ^= C_[i] & q[N_ - i];
			if (d_) {
				if (N_ + 1 - L_ > L_) {
					for (int i = 1; i <= L_; i++)T_[i] = C_[i];
					for (int i = N_ - m_; i <= N_ + 1 - L_; i++)C_[i] ^= B_[i - (N_ - m_)];
					for (int i = 1; i <= L_; i++)B_[i] = T_[i];
					L_ = N_ + 1 - L_;
					m_ = N_;
				}
				else {
					for (int i = N_ - m_; i <= N_ + 1 - L_; i++)C_[i] ^= B_[i - (N_ - m_)];
				}
			}
			if (L_ == (N_ + 1) / 2)
				latestpoint = N_ + 1;
		}

		N_ = frontlen / 2 + latestpoint / 2 - 1;
		q[N_] ^= 1;

		rounds = frontlen - N_;
		for (int i = 0; i < rounds; i++)OneRound(q);

		for (int i = 0; i < n2; i++) {
			AdjustOneRound(q);
			ReverseAdjustOneRound(q);
			OneRound(q);
			OneRound(q);
		}
		for (int i = 0; i < n1; i++) {
			ReverseAdjustOneRound(q);
			OneRound(q);
		}

		for (int i = 0; i < mm - 1; i++) {
			d_ = 0;
			for (int i = 1; i <= L_; i++)d_ ^= C_[i] & q[N_ - i];
			q[N_] = d_;
			N_++;
		}
		ReverseAdjustOneRound(q);

		for (int i = 0; i < 500; i++)s[wai * 500 + i] = q[i];
	}
	fclose(fillfp);
}



int main() {
	FILE *fp;
	char basename[20] = "biggap10_4_1";
	int adj = 4;
	fp = fopen(basename, "wb");
	MakeSequence(111, 46, 22, 10);
	unsigned char buffer2;
	for (int i = 0; i < Len / 8; i++) {
		buffer2 = 0;
		for (int j = 0; j<8; j++) {
			buffer2 <<= 1;
			if (s[i * 8 + j])
				buffer2 ^= 1;
		}
		fwrite(&buffer2, sizeof(unsigned char), 1, fp);
	}
	fclose(fp);
	return 0;
}
