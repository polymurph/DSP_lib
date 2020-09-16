#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>

using namespace std;

void shuffle(vector<complex<float>>& x);

bool IsPowerOfTwo(ulong x);

unsigned int reverseBits(unsigned int num);

unsigned int bitrev(unsigned int n,
	       	    const unsigned int& B);

void swapper(vector<int>& v)
{
	iter_swap(v.begin(), v.begin() + 3);
}

void dftmerge(vector<complex<float>>& x);

void fft(const unsigned int& N,
	 vector<complex<float>>& x);

void ifft(vector<complex<float>>& X);

void zero_padd(vector<complex<float>>& v,
	       const unsigned int& n);


void convolve(vector<complex<float>>& x,
      const vector<complex<float>>& h);

template<typename T>
void vect_add(vector<T>& a, const vector<T>& b);

template<typename T>
void vect_sub(vector<T>& a, const vector<T>& b);



int main(void)
{
	unsigned int i = 0;
	int indx = 0;
	complex<float> li(0,1);
	/*
	vector<complex<float>> v = {0 + 0 * li,1 + 0 * li,2 + 0 * li,
				    3 + 0 * li,4 + 0 * li,5 + 0 * li,
				    6 + 0 * li,7 + 0 * li,8 + 0 * li,
				    9 + 0 * li,10 + 0 * li,11 + 0 * li,
				    12 + 0 * li,13 + 0 * li,14 + 0 * li,
				    15 + 0 * li};
	*/

	vector<complex<float>> v;
	vector<complex<float>> h;
	/*
	v.push_back(complex<float>(indx++,0));
	v.push_back(complex<float>(indx++,0));
	v.push_back(complex<float>(indx++,0));
	v.push_back(complex<float>(indx++,0));
	v.push_back(complex<float>(indx++,0));
	v.push_back(complex<float>(indx++,0));
	v.push_back(complex<float>(indx++,0));
	v.push_back(complex<float>(indx++,0));
	*/
	
	v.push_back(complex<float>(1,0));
	v.push_back(complex<float>(0,0));
	v.push_back(complex<float>(-1,0));
	v.push_back(complex<float>(0,0));
	v.push_back(complex<float>(1,0));
	v.push_back(complex<float>(0,0));
	v.push_back(complex<float>(-1,0));
	v.push_back(complex<float>(0,0));
	
	h.push_back(complex<float>(1,0));
	h.push_back(complex<float>(0,0));
	h.push_back(complex<float>(0,0));
	h.push_back(complex<float>(-1,0));
	h.push_back(complex<float>(0,0));
	h.push_back(complex<float>(0,0));
	h.push_back(complex<float>(1,0));
	h.push_back(complex<float>(0,0));



	cout << "size = " << v.size() << endl;


	for(int i=0; i < v.size(); i++)
		std::cout << v.at(i) << ' ';
	cout << endl;
	
	//iter_swap(v.begin(), v.begin() + 2);	
	
	//swapper(v);

	//shuffle(v);
	
	fft(v.size(),v);

	
	for(i=0; i < v.size(); i++)
		std::cout << v.at(i) << ' ';
	cout << endl;

	for(i=0; i < v.size(); i++)
		std::cout << abs(v.at(i)) / v.size() << ' ';
	cout << endl;

	ifft(v);


	for(i=0; i < v.size(); i++)
		std::cout << v.at(i) << ' ';
	cout << endl;



	vect_add<complex<float>>(v,h);

	for(i=0; i < v.size(); i++)
		std::cout << v.at(i) << ' ';
	cout << endl;

	return 0;
}


void shuffle(vector<complex<float>>& x)
{
	int N = x.size();
	// check if vector size is multiple of 2. If not zero padd one zero;
	if(!IsPowerOfTwo(N)) return;
	
	unsigned int n = 0;
	unsigned int B = 1;
	unsigned int r = 0;

	while((N >> B) > 0) B++;
	
	B--;

	for(n = 0; n < N; n++) {
		r = bitrev(n, B);
		if(r < n) continue;
		iter_swap(x.begin() + n, x.begin() + r);	
	}
}

bool IsPowerOfTwo(ulong x)
{
    return (x != 0) && ((x & (x - 1)) == 0);
}


unsigned int reverseBits(unsigned int num)
{
    unsigned int  NO_OF_BITS = sizeof(num) * 8;
    unsigned int reverse_num = 0;
    int i;
    for (i = 0; i < NO_OF_BITS; i++)
    {
        if((num & (1 << i)))
           reverse_num |= 1 << ((NO_OF_BITS - 1) - i);
   }
    return reverse_num;
}

unsigned int bitrev(unsigned int n,
	       	    const unsigned int& B)
{
	int m;
        unsigned int r;

	for(r = 0, m = B - 1; m >= 0; m--) {
		if((n >> m) == 1) {
			r += 1 << (B - 1 - m);
			n -= 1 << (m);	
		}
	}
	return r;
}

void dftmerge(vector<complex<float>>& x)
{
	unsigned int k;
	unsigned int ii;
	unsigned int p;
	unsigned int q;
	unsigned int M = 2;
	
	complex<float> A;
	complex<float> B;
	complex<float> V;
	complex<float> W;

	complex<float> li(0,1);
	complex<float> temp(0,0);

	unsigned int N = x.size();
	
	while(M <= N) {
		// calculating twiddle factor
		temp.real(0);
		temp.imag(-2 * M_PI / M);
		W = exp(temp);
		V.real(1);
		V.imag(0);

		for(k = 0;k < M / 2; k++) {
			for(ii = 0; ii < N; ii += M) {
				p = k + ii;
				q = p + M / 2;
				A = x.at(p);
				B = x.at(q) * V;
				// butterfly operations
				x.at(p) = A + B;
				x.at(q) = A - B;	
			}
			V = V * W;
		}
		M = 2 * M;
	}
}

void fft(const unsigned int& N,
	 vector<complex<float>>& x)
{
	if(x.size() < N) {
		zero_padd(x, N - x.size());
	}

	if(!IsPowerOfTwo(N)) return;

	shuffle(x);
	dftmerge(x);
}


void ifft(vector<complex<float>>& X)
{
	unsigned int k;
	unsigned int N = X.size();

	complex<float> N_c(N,0);

	// conjugate input
	for(k = 0; k < N; k++) {
		X.at(k) = conj(X.at(k));
	}

	fft(N,X);

	for(k = 0; k < N; k++) {
		X.at(k) = conj(X.at(k)) / N_c;
	}
}

void zero_padd(vector<complex<float>>& v, 
	       const unsigned int& n)
{
	unsigned int i = 0;

	complex<float> c(0,0);

	for(i = n; i > 0; i--) {
		v.push_back(c);
	}
}

void convolve(vector<complex<float>>& x,
      	      const vector<complex<float>>& h)
{


}

template<typename T>
void vect_add(vector<T>& a, const vector<T>& b)
{
	if(a.size() != b.size()) return;
	
	for(unsigned int i = 0; i < a.size(); i++) {
		a.at(i) = a.at(i) + b.at(i);
	}
}

template<typename T>
void vect_sub(vector<T>& a, const vector<T>& b)
{
	if(a.size() != b.size()) return;
	
	for(unsigned int i = 0; i < a.size(); i++) {
		a.at(i) = a.at(i) - b.at(i);
	}
}



