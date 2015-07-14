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

void lm(int n, Mat K, Mat params, int iters=200, double lambda=0.01)
{
	bool updateJ=true;
	int Ndata=2*n; // corner points???
	int Nparams=6;

	//Mat K(4,4,CV_64F);  // settare K...
    Mat omega(3,3,CV_64F);


	double wx=params.at<double>(0,0); //w.at<double>(0,0);
	double wy=params.at<double>(0,1); //w.at<double>(0,1);
	double wz=params.at<double>(0,2); //w.at<double>(0,2);
    double tx=params.at<double>(0,3);
    double ty=params.at<double>(0,4);
    double tz=params.at<double>(0,5);

    double wx_lm,wy_lm,wz_lm,tx_lm,ty_lm,tz_lm;

	double e;
	Mat H,J,d;

	// main cycle
	for (int it=1;it<iters;it++)
	{
		// create intr. param matrix
		if(updateJ)
		{

			float theta2=wx*wx + wy*wy + wz*wz;
			double theta=sqrt(theta2);

            // update w (rot matrix rodriguez)
			omega.at<double>(0,1)=-wz;
			omega.at<double>(0,2)=wy;
			omega.at<double>(1,0)=wz;
			omega.at<double>(1,2)=-wx;
			omega.at<double>(2,0)=-wy;
			omega.at<double>(2,1)=wx;

            // update R
			Mat R=Mat::eye(3,3,CV_64F) + (sin(theta)/theta)*omega+((1-cos(theta)/(theta*theta)))*(omega*omega);
			Mat t(1,3,CV_64F);
            // update t
			t.at<double>(0,0)=tx;
            t.at<double>(0,1)=ty;
            t.at<double>(0,2)=tz;

            cout << "K: " << K << " w: " << w << " t: " <<t << endl;


			// evaluate jacobian with current params w and t
            J=Mat::zeros(Ndata,Nparams,CV_64F);
            d=Mat::zeros(Ndata,1,CV_64F);

            //check depth over points
			for(int i=0;i<n;i++)
			{
                Mat x(1,3,CV_64F); // i-th point
                Mat Rt=R;
                // Rt.at<double>(0,2)=t.at<double>(0,0); // perchè??
                // Rt.at<double>(1,2)=t.at<double>(0,1); // perchè??
                // Rt.at<double>(2,2)=t.at<double>(0,2); // perchè??

                // //x.at<double>(0,0)=   point x coord
                // //x.at<double>(0,1)=   point y coord
                // //x.at<double>(0,2)= 1.0;

                // Mat uvs=K*Rt*x;
                // double up=uvs.at<double>(0,0)/uvs.at<double>(0,2);
                // double vp=uvs.at<double>(0,1)/uvs.at<double>(0,2);
                // compute distance
                //...

			}

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
		// for(int j=0;j<6;j++)
		// 	params.at<double>(0,j)=dp.at<double>(0,j);
        wx_lm=wx+dp.at<double>(0,0);
        wy_lm=wx+dp.at<double>(0,1);
        wz_lm=wx+dp.at<double>(0,2);
        tx_lm=tx+dp.at<double>(0,3);
        ty_lm=ty+dp.at<double>(0,4);
        tz_lm=tz+dp.at<double>(0,5);


		// evaluate total distance (error) at new params
        omega.at<double>(0,1)=-wz_lm;
        omega.at<double>(0,2)=wy_lm;
        omega.at<double>(1,0)=wz_lm;
        omega.at<double>(1,2)=-wx_lm;
        omega.at<double>(2,0)=-wy_lm;
        omega.at<double>(2,1)=wx_lm;
        double theta2=wx_lm*wx_lm + wy_lm*wy_lm + wz_lm*wz_lm;
        double theta=sqrt(theta2);

        Mat R=Mat::eye(3,3,CV_64F) + (sin(theta)/theta)*omega+((1-cos(theta)/(theta*theta)))*(omega*omega);
        Mat t(1,3,CV_64F);
        t.at<double>(0,0)=tx_lm;
        t.at<double>(0,1)=ty_lm;
        t.at<double>(0,2)=tz_lm;
        Mat d_lm=Mat::zeros(Ndata,1,CV_64F);

        //check depth over points
        // for(int i=0;i<...)
        // {


        // }

        double e_lm=d_lm.dot(d_lm);
        if(e_lm<e)
        {
            lambda/=10.0;
            wx=wx_lm;
            wy=wy_lm;
            wz=wz_lm;
            tx=tx_lm;
            ty=ty_lm;
            tz=tz_lm;
            e=e_lm;
            cout << "e: " << e <<endl;
            updateJ=true;
        }
        else
        {
            updateJ=false;
            lambda*=10.0;
        }

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

    Mat K(4,4,CV_64F);

	params.at<double>(0,0)=1;
	params.at<double>(0,1)=2;
	params.at<double>(0,2)=3;
	params.at<double>(0,3)=4;
	params.at<double>(0,4)=5;
	params.at<double>(0,5)=6;


	lm(100,K,params);

	cout << t0 << endl << R0 << endl;

	return 0;

}
