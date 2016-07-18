// matrices.cpp
//
#include <iomanip>
#include <iostream>
#include <cmath>
#include <stdio.h>
#include "matrices.hpp"

using namespace std;

extern bool multithread;

void libereD (int n, double * * A)
                               {
    for (int i=0;i<n;i++) delete[] A[i];
    delete [] A;
}
void libereL (int n, long double * * A)
                                    {
    for (int i=0;i<n;i++) delete[] A[i];
    delete [] A;
}

long double** copymatL(int n, int m, long double **A) {
	long double **AA;
	AA = new long double*[n];
	for (int i=0;i<n;i++){ AA[i] = new long double[m];}
	for (int i=0;i<n;i++){
		for (int j=0;j<m;j++) {
			AA[i][j] = A[i][j];
		}
	}
	return AA;
}
double * * transpose (int n, int m, double * * A)
                                             {
    double **C;
    C = new double*[m];
    for (int i=0;i<m;i++) C[i] = new double[n];
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) C[i][j] = A[j][i];
    }
    return C;
}
double * * prodM (int n, int m, int p, double * * A, double * * B)
                                                           {
    double **C;
    C = new double*[n];
    for (int i=0;i<n;i++) {
        C[i] = new double[p];
        for(int j=0;j<p;j++) C[i][j]=0.0;
    }
    for (int i=0;i<n;i++) {
        for (int j=0;j<p;j++) {
            for (int k=0;k<m;k++) C[i][j] += A[i][k]*B[k][j];
        }  
    }
    return C;
}
double * * prodMs (int n, int m, double * * A, double b)
                                                   {
    double **C;
    C = new double*[n];
    for (int i=0;i<n;i++) {
        C[i] = new double[m];
        for(int j=0;j<m;j++) C[i][j]=A[i][j]*b;
    }
    return C;
}
long double * * transposeL (int n, int m, long double * * A)
                                                        {
    long double **C;
    C = new long double*[m];
    for (int i=0;i<m;i++) C[i] = new long double[n];
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) C[i][j] = A[j][i];
    }
    return C;
}
long double * * prodML (int n, int m, int p, long double * * A, long double * * B)
                                                                           {
    long double **C;
	int j,k;
    C = new long double*[n];
    for (int i=0;i<n;i++) {
        C[i] = new long double[p];
        for(j=0;j<p;j++) C[i][j]=0.0;
    }
#pragma omp parallel for private(k,j) if(multithread)
    for (int i=0;i<n;i++) {
        for (j=0;j<p;j++) {
            for (k=0;k<m;k++) {
				C[i][j] += A[i][k]*B[k][j];
				/*if (C[i][j]!=C[i][j]) {
					cout<<"dans prodML\n";
					cout<<"A["<<i<<"]["<<k<<"] = "<<A[i][k]<<"\n";
					cout<<"B["<<k<<"]["<<j<<"] = "<<B[k][j]<<"\n";
				}*/
			}
        }
    }
    return C;
}
long double * * prodMsL (int n, int m, long double * * A, long double b)
                                                                   {
    long double **C;
    C = new long double*[n];
    for (int i=0;i<n;i++) {
        C[i] = new long double[m];
        for(int j=0;j<m;j++) C[i][j]=A[i][j]*b;
    }
    return C;
}
double * * invM (int n, double * * A)
{
    int i,j,k,l,err=0;
    double max,pivot,coef, **T, **C;
    T = new double*[n];for (i=0;i<n;i++) T[i]= new double[2*n]; 
    C = new double*[n];for (i=0;i<n;i++) C[i]= new double[n];
    for (i=0;i<n;i++) 
        {for (j=0;j<n;j++) 
            {T[i][j]=A[i][j];
             if(i==j) T[i][j+n]=1.0; else T[i][j+n]=0.0;
            }
        }
    k=0;
    while ((err==0)&&(k<n))
        {max=fabs(T[k][k]);
         l=k;
         for (i=k+1;i<n;i++)
            {if (max<fabs(T[i][k])) {max=fabs(T[i][k]);l=i;}
            }
         if (max!=0.0)
            {for (j=k;j<2*n;j++)
                {pivot=T[k][j];
                 T[k][j]=T[l][j];
                 T[l][j]=pivot;
                }
             pivot=T[k][k];
             if (pivot!=0.0)
                {
                 for (j=k+1;j<2*n;j++) T[k][j]/=pivot;
#pragma omp parallel for private(coef,j) if(multithread)
                 for (i=0;i<n;i++)
                    {if (i!=k) 
                        {coef=T[i][k];
                         for (j=k+1;j<2*n;j++) T[i][j]-=coef*T[k][j];
                        }
                    }
                } else err=3;
            } else err=4;
         k++;   
        }
    if (err==0)
        {for (i=0;i<n;i++)
            {for (j=0;j<n;j++) C[i][j]=T[i][j+n];}
        }
    for (i=0;i<n;i++) delete[] T[i];delete[] T;
    return C; 
}
int inverse (int n, long double * * A, long double * * C)
{
    int i,j,k,l,err=0;
    long double max,pivot,coef, **T;
    T = new long double*[n];for (i=0;i<n;i++) T[i]= new long double[2*n]; 
//cout<<"début de inverse   n="<<n<<"\n";
    for (i=0;i<n;i++) 
        {for (j=0;j<n;j++) 
            {T[i][j]=A[i][j];
             if(i==j) T[i][j+n]=1.0; else T[i][j+n]=0.0;
            }
        }
//cout<<"fin de l'initialisation de T\n";
	k=0;
    while ((err==0)&&(k<n))
        {max=fabs(T[k][k]);
         l=k;
         for (i=k+1;i<n;i++)
            {if (max<fabs(T[i][k])) {max=fabs(T[i][k]);l=i;}
            }
         if (fabs(max)>EPSILON)//(max!=0.0)
            {for (j=k;j<2*n;j++)
                {pivot=T[k][j];
                 T[k][j]=T[l][j];
                 T[l][j]=pivot;
                }
             pivot=T[k][k];
             if (fabs(pivot)>EPSILON)   //pivot!=0.0)
                {
                 for (j=k+1;j<2*n;j++) T[k][j]/=pivot;
#pragma omp parallel for private(coef,j)
                 for (i=0;i<n;i++)
                    {if (i!=k) 
                        {coef=T[i][k];
                         for (j=k+1;j<2*n;j++) T[i][j]-=coef*T[k][j];
                        }
                    }
                } else err=3;
            } else err=4;
         k++;   
        }
//cout<<"fin du while   err="<<err<<"   n="<<n<<"\n";
    if (err==0)
        {for (i=0;i<n;i++)
            {for (j=0;j<n;j++) C[i][j]=T[i][j+n];}
        }
//cout<<"apres l'ecriture de C\n";
    if (err==4) {
	    std::cout<<"err="<<err<<"\n";
		std::cout<<"k+1="<<k+1<<"      n="<<n<<"\n";
		for  (i=k+1;i<n;i++) std::cout<<"T["<<i<<"]["<<k<<"]="<<T[i][k]<<"\n";
		std::cout <<"\n";
		/*for (i=0;i<n;i++) {
			for (j=0;j<n;j++) std::cout<<A[i][j]<<"  ";
		    std::cout<<"\n";
		}*/
	}
    for (i=0;i<n;i++) delete[] T[i];delete[] T;
//cout<<"fin de inverse\n";
    return err; 
}
int jacobi (int n, double * * A, double * D, double * * V)
                                                         {
        int nrot=0,nm;
        double *b,*z,tresh,theta,tau,t,sm,s,h,g,cc;
        b = new double[n];
        z = new double[n];
        for (int ip=0;ip<n;ip++) {for (int iq=0;iq<n;iq++) V[ip][iq]=0.0; V[ip][ip]=1.0;}
        for (int ip=0;ip<n;ip++) {b[ip]=A[ip][ip];D[ip]=b[ip];z[ip]=0.0;}
        std::cout<<"matrice A dans jacobi  n="<<n<<"\n";
		if (n<10) nm=n; else nm=10;
		for (int i=0;i<nm;i++) {
		  for (int j=0;j<nm;j++) std::cout<< A[i][j]<<"  ";
		  std::cout<<"\n";
		}
		std::cout<<"\n";
        
        for (int i=0;i<51;i++) {
            sm=0.0;
            for (int ip=0;ip<n-1;ip++) {for (int iq=ip+1;iq<n;iq++) sm += fabs(A[ip][iq]);}
            std::cout<<"jacobi sm="<<sm<<"\n";
            if (sm==0.0) {delete []b;delete []z;return nrot;}
            if (i<4) tresh=0.2*sm/(double)n/(double)n; else tresh=0.0;
            for (int ip=0;ip<n-1;ip++) {
                for (int iq=ip+1;iq<n;iq++) {
                    g=100.0*fabs(A[ip][iq]);
                    if (((i>4)and(fabs(D[ip]+g)==fabs(D[ip])))and(fabs(D[iq]+g)==fabs(D[iq]))) A[ip][iq]=0.0;
                    else if (fabs(A[ip][iq])>tresh) {
                        h=D[iq]-D[ip];
                        if (fabs(h)+g == fabs(h)) t=A[ip][iq]/h;
                        else  {
                            theta=0.5*h/A[ip][iq];
                            t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                            if (theta<0.0) t=-t;
                        }
                        cc=1.0/sqrt(1+t*t);
                        s=t*cc;
                        tau=s/(1.0+cc);
                        h=t*A[ip][iq];
                        z[ip] -=h;
                        z[iq] +=h;
                        D[ip] -=h;
                        D[iq] +=h;
                        A[ip][iq]=0.0;
                        for (int j=0;j<ip;j++) {
                            g=A[j][ip];
                            h=A[j][iq];
                            A[j][ip] = g-s*(h+g*tau);
                            A[j][iq] = h+s*(g-h*tau);
                        }
                        for (int j=ip+1;j<iq;j++) {
                            g=A[ip][j];
                            h=A[j][iq];
                            A[ip][j] = g-s*(h+g*tau);
                            A[j][iq] = h+s*(g-h*tau);
                        }
                       for (int j=iq+1;j<n;j++) {
                            g=A[ip][j];  
                            h=A[iq][j];  
                            A[ip][j] = g-s*(h+g*tau);   
                            A[iq][j] = h+s*(g-h*tau);   
                        }
                        for (int j=0;j<n;j++) {
                            g=V[j][ip];
                            h=V[j][iq];
                            V[j][ip] = g-s*(h+g*tau);
                            V[j][iq] = h+s*(g-h*tau);
                        }
                    }
                }  
            } 
            nrot++;
        }
        for (int ip=0;ip<n-1;ip++) {
            b[ip]=b[ip]+z[ip];
            D[ip]=b[ip];
            z[ip]=0.0;
        }
        delete []b;delete []z;
        return nrot;    
    }
int jacobiL (int n, long double * * A, long double * D, long double * * V)
                                                                             {
        int nrot=0,iq;
        long double *b,*z,tresh,theta,tau,t,sm,s,h,g,cc;
        b = new long double[n];
        z = new long double[n];
        for (int ip=0;ip<n;ip++) {for (int iq=0;iq<n;iq++) V[ip][iq]=0.0; V[ip][ip]=1.0;}
        for (int ip=0;ip<n;ip++) {b[ip]=A[ip][ip];D[ip]=b[ip];z[ip]=0.0;}
        //std::cout<<"matrice A dans jacobi  n="<<n<<"\n";
		/*if (n<10) nm=n; else nm=10;
		for (int i=0;i<nm;i++) {
		  for (int j=0;j<nm;j++) std::cout<< A[i][j]<<"  ";
		  std::cout<<"\n";
		}
		std::cout<<"\n";*/
        
        for (int i=0;i<51;i++) {
            sm=0.0;
            for (int ip=0;ip<n-1;ip++) {for (int iq=ip+1;iq<n;iq++) sm += fabs(A[ip][iq]);}
            //std::cout<<"jacobi sm="<<sm<<"\n";
            if (sm==0.0) {delete []b;delete []z;return nrot;}
            if (i<4) tresh=0.2*sm/(long double)n/(long double)n; else tresh=0.0;
//#pragma omp parallel for shared (i,n,A,D,z,V,tresh) private(ip,iq,g,theta,t,cc,s,tau,h)
            for (int ip=0;ip<n-1;ip++) {
                for (iq=ip+1;iq<n;iq++) {
                    g=100.0*fabs(A[ip][iq]);
                    if (((i>4)and(fabs(D[ip]+g)==fabs(D[ip])))and(fabs(D[iq]+g)==fabs(D[iq]))) A[ip][iq]=0.0;
                    else if (fabs(A[ip][iq])>tresh) {
                        h=D[iq]-D[ip];
                        if (fabs(h)+g == fabs(h)) t=A[ip][iq]/h;
                        else  {
                            theta=0.5*h/A[ip][iq];
                            t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                            if (theta<0.0) t=-t;
                        }
                        cc=1.0/sqrt(1+t*t);
                        s=t*cc;
                        tau=s/(1.0+cc);
                        h=t*A[ip][iq];
                        z[ip] -=h;
                        z[iq] +=h;
                        D[ip] -=h;
                        D[iq] +=h;
                        A[ip][iq]=0.0;
                        for (int j=0;j<ip;j++) {
                            g=A[j][ip];
                            h=A[j][iq];
                            A[j][ip] = g-s*(h+g*tau);
                            A[j][iq] = h+s*(g-h*tau);
                        }
                        for (int j=ip+1;j<iq;j++) {
                            g=A[ip][j];
                            h=A[j][iq];
                            A[ip][j] = g-s*(h+g*tau);
                            A[j][iq] = h+s*(g-h*tau);
                        }
                       for (int j=iq+1;j<n;j++) {
                            g=A[ip][j];  
                            h=A[iq][j];  
                            A[ip][j] = g-s*(h+g*tau);   
                            A[iq][j] = h+s*(g-h*tau);   
                        }
                        for (int j=0;j<n;j++) {
                            g=V[j][ip];
                            h=V[j][iq];
                            V[j][ip] = g-s*(h+g*tau);
                            V[j][iq] = h+s*(g-h*tau);
                        }
                    }
                }  
            } 
            nrot++;
        }
        for (int ip=0;ip<n-1;ip++) {
            b[ip]=b[ip]+z[ip];
            D[ip]=b[ip];
            z[ip]=0.0;
        }
        delete []b;delete []z;
        return nrot;    
    }
double kappa (int n, long double * * A) {
		long double **vcprop, *valprop,piv,**B, vpn1,vp0;
        vcprop  = new long double*[n];for (int i=0;i<n;i++) vcprop[i]=new long double [n];
        valprop = new long double[n];
        B  = new long double*[n];for (int i=0;i<n;i++) B[i]=new long double [n];
		for (int i=0;i<n;i++) {for (int j=0;j<n;j++) B[i][j]=A[i][j];}
		int nrot=jacobiL(n,B,valprop,vcprop); if (nrot==0) cout<<"kappa nrot=0\n";
        for (int i=0;i<n-1;i++) {
            for (int j=i+1;j<n;j++) {
                if (valprop[i]<valprop[j]) {piv=valprop[i];valprop[i]=valprop[j];valprop[j]=piv;}
            }
        }
		vpn1=valprop[n-1];vp0=valprop[0];
		for (int i=0;i<n;i++) delete [] B[i];delete [] B;       
		for (int i=0;i<n;i++) delete [] vcprop[i];delete [] vcprop;
		delete [] valprop;
	    if (vpn1>0) return vp0/vpn1; else return 1E100;
	}
	
int inverse_Tik (int n, long double * * A, long double * * C)
                                                                 {
		long double t,coeff,**AA,seuil_kappa=1E8;
		double kap;
		int err = 0;
		AA = new long double*[n]; for (int i=0;i<n;i++) AA[i] = new long double [n];
//		std::cout<<"debut d'inverse_Tik  n="<<n<<"\n";
		coeff=0;
		
		kap=kappa(n,A);
		//std::cout<<"coeff = "<<coeff<<"   err="<<err<<"   kappa="<<kap<<"\n";
		/*for (int i=0;i<n;i++) {
			for (int j=0;j<n;j++) printf(" %Le8",A[i][j]);printf("\n");
		}*/
		err=inverse(n,A,C);
		//std::cout<<"err="<<err<<"\n";
		if ((err==0)and(kap<seuil_kappa)) return 0;
		
		coeff=1E-15;
		for (int i=0;i<n;i++) {for (int j=0;j<n;j++) AA[i][j] = A[i][j];}
		t=0.0;
		for (int i=0;i<n;i++) t = t + fabs(A[i][i]);
		t /= n;
		
		for (int i=0;i<n;i++)  AA[i][i] = A[i][i] + coeff*t;
		
		while (((err!=0)or(kap>seuil_kappa))and(coeff<0.1)) {
		    kap = kappa(n,AA);
			err=inverse(n,AA,C);
			//std::cout<<"coeff = "<<coeff<<"   err="<<err<<"   kappa="<<kap<<"\n";
			if ((err!=0)or(kap>seuil_kappa)) {
			  coeff *=10.0;
			  for (int i=0;i<n;i++)  AA[i][i] = A[i][i] + coeff*t;
			}
		}
		return 1;  
	}
int inverse_Tik2 (int n, long double * * A, long double * * C, long double coeff)
                                                                                     {
		long double t,**AA;
		int err;
		AA = new long double*[n]; for (int i=0;i<n;i++) AA[i] = new long double [n];
		for (int i=0;i<n;i++) {for (int j=0;j<n;j++) AA[i][j] = A[i][j];}
		t=0.0;
		for (int i=0;i<n;i++) t = t + fabs(A[i][i]);
		t /= n;		
		for (int i=0;i<n;i++)  AA[i][i] = A[i][i] + coeff*t;
		err=inverse(n,AA,C);
		for (int i=0;i<n;i++) delete [] AA[i];delete [] AA;
		return err;  
	}
