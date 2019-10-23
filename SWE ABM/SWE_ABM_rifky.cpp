//Linear shallow water solver using Adam-Bashforth-Moulton (ABM) predictor-corrector method
//this code is written by rifkyfauzi9@gmail.com
//pleas let me know if there is an error

#include<iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include <algorithm>
#include <complex>
#include <time.h>
#include <sstream>
#include<string>

//# define M_PI           3.14159265358979323846  /* pi */
using namespace std;



float fun4( float k ) {
    return( 0.25 * exp(-500*(k-1)*(k-1)) );
}

double fh(double h[], double u[], double dx){
	double a = - 0.5 * (u[2] - u[0]) / (2 * dx);
	return a; //<h3>
}
double fu(double h[], double u[],  double g, double dx){
	double b =  - g * (h[2] - h[0]) / (2 * dx)   ;
	return b; //<u3>
}


void SWE() {
//int main() {

    double L = 2, dx = 0.001, dt = dx*0.01;

    long int X = round(L/dx);
    long int T = round(3.5/dt);
    double x[X-1], h[X - 1],  u[X - 1], hi[X - 1], ui[X - 1], hh[X - 1], uu[X - 1], hn[X - 1], un[X - 1], H1[X - 1], H2[X - 1], H3[X - 1], HH[X - 1], U1[X - 1], U2[X - 1], U3[X - 1], UU[X - 1], Hs[4][X - 1], Us[4][X - 1], umax;
    double ht[23][X-1],  ta, tb, tc, td;
		float g = 9.8;
    int n,j, vv=1;
    int pembagi = 100, maks = T/100;
    int v = 1, skip = round( T / 1000 );


    if (skip == 0){
        skip = 1;
    }

    // space & time discretization
    for (j = 0; j <= X ; j++) {
        x[j + 1] = x[j] + dx;
    }


    printf("Please wait \n");


    /// initial condition
    x[0] = dt;
    for (j = 0; j <= X - 1; j++) {
        h[j] = fun4( x[j] );
        //u[j] = fun2(x[j], h0,  k,  co,  F);
        u[j] = 0;
    }
    // initial wave profile
    for (j = 0; j <= X - 1; j++) {
        hi[j] = fun4( x[j] );
        ui[j] = 0;
    }

    for (j = 0; j <= X - 1; j++) {
        H1[j] = 0;
        U1[j] = 0;
        H2[j] = 0;
        U2[j] = 0;
    }

		ofstream data_plot, u_plot;
    data_plot.open ("elevation.csv");
    u_plot.open ("velocity.csv");

    for (int j = 0; j <= X - 2; j++){
        data_plot << x[j] <<"," ;
    }
    data_plot << x[X - 1] <<"\n";

    for (int j = 0; j <= X - 2; j++){
        data_plot << hi[j] <<"," ;
				u_plot << ui[j] <<"," ;
    }
    data_plot << hi[X - 1] <<"\n";
    u_plot << ui[X - 1] <<"\n";

    /// main iteration starts
    for (n = 0; n <= T; n++) {

        // Predictor Operator
        for (j = 1; j <= X - 2 ; j++){
			double ha[] = {h[j-1], h[j], h[j+1]}, ua[] = {u[j-1], u[j], u[j+1]};
			H3[j] = fh(ha, ua, dx);
			U3[j] = fu(ha, ua,  g, dx);
        }

        // Periodic Boundary Condition
        /*
         H3[0] = H3[X - 2];
         H3[X - 1] = H3[1];
         U3[0] = U3[X - 2];
         U3[X - 1] = U3[1];
         */

        // Predictor steps
        for (j = 1; j <= X - 2; j++){
            hh[j] = h[j] + dt * (23 * H3[j] - 16 * H2[j] + 5 * H1[j])/12;
            uu[j] = u[j] + dt * (23 * U3[j] - 16 * U2[j] + 5 * U1[j])/12;
        }

        hh[X - 1] = hh[1];
        hh[0] = hh[X - 2];

        uu[X - 1] = uu[1];
        uu[0] = uu[X - 2];


        // Corrector Operator
        for (j = 1; j <= X - 2; j++) {
			double ha[] = {hh[j-1], hh[j], hh[j+1]}, ua[] = {uu[j-1], uu[j], uu[j+1]};
			HH[j] = fh(ha, ua, dx);
			UU[j] = fu(ha, ua,  g, dx);
    	}

        // Periodic Boundary Condition
        /*
         HH[0] = HH[X - 2];
         HH[X - 1] = HH[1];
         UU[0] = UU[X - 2];
         UU[X - 1] = UU[1];
        */

        // Corrector steps
        for (j = 1; j <= X - 2; j++) {
            hn[j] = h[j] + dt * (5 * HH[j] + 8 * H3[j] - H2[j])/12;
            un[j] = u[j] + dt * (5 * UU[j] + 8 * U3[j] - U2[j])/12;
        }

        hn[X - 1] = hn[1];
        hn[0] = hn[X - 2];

        un[X - 1] = un[1];
        un[0] = un[X - 2];



        //printf("hhX-1 =  %f \n", hh[X-1]);
        // update operator

        for (j = 0; j <= X - 1; j++) {
            H1[j] = H2[j];
            H2[j] = H3[j];
            U1[j] = U2[j];
            U2[j] = U3[j];
        }


        // update
        for (int i = 0; i <= X - 1; i++) {
            h[i] = hn[i];
            u[i] = un[i];
        }

        //plot
        if ( n == skip * v ) {
        	// fluid depth data save
            for (j = 0; j <= X - 2; j++){
                data_plot << hn[j] <<"," ;
            }
            data_plot << hn[X - 1];
            data_plot <<"\n";
            // fluid velocity data save
            for (j = 0; j <= X - 2; j++){
                u_plot << un[j] <<"," ;
            }
            u_plot << un[X - 1];
            u_plot <<"\n";

            v++;
        }


        //umax = *std::max_element(un, un + X - 1);
        //umaxx << umax << "\n";
    } //main iteartion stops
    cout<< "\n" <<endl;


    // Print Results

		if (hn[1] != hn[1]){
			printf(" blow up\n");
		}


		data_plot.close();
		u_plot.close();



    cout << "Iterasi Berhasil \n";
}

int main() {

	int ind = 1;
	for (int j = 0; j < 1; j++){
		SWE( );
		ind = ind + 1;
	}


}
