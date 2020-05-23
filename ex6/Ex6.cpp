#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;
double pi = 3.1415926;
//############################ define function #######################
int make_x(double* x_axis,double A,int N);
int make_real_ft(double* in_arr,double* axis,int N);
int make_image_ft(double* in_arr,double* axis,int N);
int DFT(double* image,double* real,double* gj,int N,double A);
//###############################################################
int main(){
    //################ constant #######################
    int N = 2048;
    double A = 7;
    //################ make x_axis ##################
    double x_axis[N];

    make_x(x_axis, A,N); // call func
    //################ fk #####################
    double real_ft[N];
    double image_ft[N];

    make_real_ft(real_ft,x_axis,N); // call func
    make_image_ft(image_ft,x_axis,N); // call func
    //############## DFT ###################
    double gj[N];

    DFT(real_ft,image_ft,gj,N,A); // call func
    //############### write file ################
    ofstream file;
    file.open("output.txt");
    for (int i=0;i<N;i++){
        file << x_axis[i] <<"  "<<gj[i] <<"\n";
    }
    file.close();
     //########################################

    system("draw.py");
    return 0;}
    
int make_x(double* x_axis,double A,int N){
    double step = 2*A/(N-1);
    double now = -A;
    for (int i=0;i<N;i++){
        x_axis[i] = now;
        now = now + step;
    }}
int make_real_ft(double* in_arr,double* axis,int N){
    double t;
    for(int i=0;i<N;i++){
        t = axis[i];
        in_arr[i] = exp(-t*t)*sin(2*t);
    }}
int make_image_ft(double* in_arr,double* axis,int N){
    double t;
    for(int i=0;i<N;i++){
        t = axis[i];
        in_arr[i] = exp(-t*t)*cos(2*t);
    }}
int DFT(double* image,double* real,double* gj,int N,double A){
    double shift_j,exp_it;
    double d_N = N;
    double tt = 0;
    double t2 = 0;
    double step = 2*A/(d_N-1);
    for (int j=0;j<N;j++){
        tt = 0.0;
        t2 = 0.0;
        shift_j = j-d_N/2;
        for (int k=0;k<N;k++){
            
            tt += (cos(2*pi*shift_j*k/d_N)*real[k] + sin(2*pi*shift_j*k/d_N)*image[k]);

            t2 += (cos(2*pi*shift_j*k/d_N)*image[k] - sin(2*pi*shift_j*k/d_N)*real[k]);
        }
        tt = tt/sqrt(d_N);
        t2 = t2/sqrt(d_N);
        exp_it = (cos(2*A*pi*shift_j/d_N/step)*tt - sin(2*A*pi*shift_j/d_N/step)*t2);

        gj[j] = exp_it/(sqrt(2*pi)/sqrt(d_N)/step) ;
    }}