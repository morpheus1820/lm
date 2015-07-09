#include <iostream>
#include <vector>
#include <math.h>
#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"

using namespace cv;
using namespace std;

Mat t0,R0,w;

Mat rodriguez(Mat R)
{
	Mat w(1,3,CV_64F);
	double theta = acos( (trace(R)[0]-1)/2.0 );
	// vconcat( R.at<double>(2,1)-R.at<double>(1,2),   
	// 		 R.at<double>(0,2)-R.at<double>(2,0),
	// 		 R.at<double>(1,0)-R.at<double>(0,1),
	// 		 w);

	w.at<double>(0,0)=R.at<double>(2,1)-R.at<double>(1,2);
	w.at<double>(0,1)=R.at<double>(0,2)-R.at<double>(2,0);
	w.at<double>(0,2)=R.at<double>(1,0)-R.at<double>(0,1);

	w*=theta/(2*sin(theta));

	return w;
}

void lm(int n, Mat params, int iters=200, double lambda=0.01)
{
	bool updateJ=true;
	int Ndata=2*n; // corner points???
	int Nparams=11;

	Mat K(4,4,CV_64F);
	double wx=w.at<double>(0,0);
	double wy=w.at<double>(0,1);
	double wz=w.at<double>(0,2);
	
	double e;
	Mat H,J,d;

	// main cycle
	for (int it=1;it<iters;it++)
	{
		// create intr. param matrix
		if(updateJ)
		{
			K.at<double>(0,0)=params.at<double>(0,0);
			K.at<double>(1,1)=params.at<double>(0,1);
			K.at<double>(2,2)=params.at<double>(0,2);
			K.at<double>(3,3)=1.0;
			K.at<double>(0,3)=params.at<double>(0,3);
			K.at<double>(1,3)=params.at<double>(0,4);
			K.at<double>(2,3)=params.at<double>(0,5);

			cout << "K: " << K << endl;

			float theta2=wx*wx + wy*wy + wz*wz;
			double theta=sqrt(theta2);
			Mat omega(3,3,CV_64F);
			omega.at<double>(0,1)=-wz;
			omega.at<double>(0,2)=wy;
			omega.at<double>(1,0)=wz;
			omega.at<double>(1,2)=-wx;
			omega.at<double>(2,0)=-wy;
			omega.at<double>(2,1)=wx;

			Mat R=Mat::eye(3,3,CV_64F) + (sin(theta)/theta)*omega+((1-cos(theta)/(theta*theta)))*(omega*omega);
			Mat t(1,3,CV_64F);
			t=t0;
				
			J=Mat::zeros(Ndata,Nparams,CV_64F);
			d=Mat::zeros(Ndata,1,CV_64F);

			// evaluate jacobian with current params
			// for(int i=0;i<xe.cols;i++)
			// {


			// }
							
			//compute hessian
			H=J.t()*J;
			if(it==0)
			{
				e=d.dot(d);
				cout << "e: "<< e << endl;
			}				
		}

		// apply damping factor
		Mat H_lm=H+(lambda*Mat::eye(Nparams,Nparams,CV_64F));

		// update params
		Mat dp= -H_lm.inv()*(J.t()*d);	
		for(int j=0;j<6;j++)
			params.at<double>(0,j)=dp.at<double>(0,j);


		// evaluate total distance (error) at new params


	}


}

int main( int argc, char **argv )
{


	t0=Mat(1,3,CV_64F);
	R0=Mat(3,3,CV_64F);

	// initial guess 
	t0.at<double>(0,0)=0; t0.at<double>(0,1)=0; t0.at<double>(0,2)=0;
	R0.at<double>(0,0)=0.1; R0.at<double>(1,1)=0.2; R0.at<double>(2,2)=0.3;

	w=rodriguez(R0);
	Mat params(1,6,CV_64F);
	params.at<double>(0,0)=1;
	params.at<double>(0,1)=2;
	params.at<double>(0,2)=3;
	params.at<double>(0,3)=4;
	params.at<double>(0,4)=5;
	params.at<double>(0,5)=6;


	lm(100,params);

	cout << t0 << endl << R0 << endl;

	return 0;

}
