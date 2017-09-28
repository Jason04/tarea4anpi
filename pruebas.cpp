#include <cmath>

/*struct LUdcmp
{
	Int n;
	MatDoub lu;
	VecInt indx;
	Doub d;
	LUdcmp(MatDoub_I &a);
	void solve(VecDoub_I &b, VecDoub_O &x);
	void solve(MatDoub_I &b, MatDoub_O &x);
	void inverse(MatDoub_O &ainv);
	Doub det();
	void mprove(VecDoub_I &b, VecDoub_IO &x);
	MatDoub_I &aref;
};*/

void lu(anpi::Matrix<T>& A, anpi::Matrix<T>& LU){
	const double TINY=1.0e-40;
	int n = A.rows();
	int imax;
	double d = 1.0;
	double big,temp;
	
	std::vector <double> vv(n); 
	std::vector <float> indx(n); //vector b

	for (int i=0;i<n;i++) {
		big=0.0;
		for (int j=0;j<n;j++){

			if ((temp=std::abs(LU[i][j])) > big){
				big = temp;
			}
		}
		if (big == 0.0){ 
			throw("Singular matrix in LUdcmp");
		}

		vv[i]=1.0/big;
	}

	for (int k=0; k<n; k++) { 

		big=0.0;
		imax=k;
		
		for (int i = k; i<n; i++) {

			temp = vv[i]*std::abs(LU[i][k]);
			if (temp > big) {
				big=temp;
				imax=i;
			}
		}
		if (k != imax) {
			for (int j=0;j<n;j++) {

				temp=LU[imax][j];
				LU[imax][j]=LU[k][j];
				LU[k][j]=temp;
			}
			d = -d;
			vv[imax]=vv[k];
		}
		indx[k]=imax;

		if (LU[k][k] == 0.0){ 
			LU[k][k]=TINY;
		}
		for (int i=k+1;i<n;i++) {

			temp=LU[i][k] /= LU[k][k];
			for (int j=k+1;j<n;j++)
				LU[i][j] -= temp*LU[k][j];
		}
	}
}



void LUdcmp::solve(std::vector<T>& b, std::vector<T>& x){
	int i,ii=0,ip,j;
	double sum;
	if (b.size() != n || x.size() != n){
		throw("LUdcmp::solve bad sizes");
	}
	for (i=0;i<n;i++){
	 x[i] = b[i];
	}
	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=x[ip];
		x[ip]=x[i];
		if (ii != 0){
			for (j=ii-1;j<i;j++){
			 sum -= lu[i][j]*x[j];
			}
		}
		else if (sum != 0.0){
			ii=i+1;
		}
		x[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=x[i];
		for (j=i+1;j<n;j++){
			sum -= lu[i][j]*x[j];
		}
		x[i]=sum/lu[i][i];
	}
}

void solve(anpi::Matrix<T>& b, anpi::Matrix<T>& x){
	int n = b.rows();
	int i,j,m=b.cols();
	if (b.rows() != n || x.rows() != n || b.cols() != x.cols()){
		throw("LUdcmp::solve bad sizes");
	}
	std::vector <double> xx(n); 
	for (j=0;j<m;j++) {
		for (i=0;i<n;i++){
		 xx[i] = b[i][j];
		}
		solve(xx,xx);
		for (i=0;i<n;i++){ 
			x[i][j] = xx[i];
		}
	}
}
/*
void LUdcmp::inverse(MatDoub_O &ainv)
{
	Int i,j;
	ainv.resize(n,n);
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) ainv[i][j] = 0.;
		ainv[i][i] = 1.;
	}
	solve(ainv,ainv);
}
Doub LUdcmp::det()
{
	Doub dd = d;
	for (Int i=0;i<n;i++) dd *= lu[i][i];
	return dd;
}
void LUdcmp::mprove(VecDoub_I &b, VecDoub_IO &x)
{
	Int i,j;
	VecDoub r(n);
	for (i=0;i<n;i++) {
		Ldoub sdp = -b[i];
		for (j=0;j<n;j++)
			sdp += (Ldoub)aref[i][j] * (Ldoub)x[j];
		r[i]=sdp;
	}
	solve(r,r);
	for (i=0;i<n;i++) x[i] -= r[i];
}
*/
