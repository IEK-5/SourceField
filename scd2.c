#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define eps 1e-15
double ell_int(double k)
{
	double K, a0, b0, c, a1, b1;
	a0=1;
	b0=sqrt(1-k*k);
	c=1;
	while(fabs(c)>eps)
	{
		a1=(a0+b0)/2;
		b1=sqrt(a0*b0);
		c=(a1-b1)/(a1+b1);
		a0=a1;
		b0=b1;
	}
	K=M_PI/(b1+a1);
	return K;
}
#undef eps


# define eps 1e-15
#define itermax 50
double ell_mod (double Kratio)
{
	double q, nom=1, denom=0;
	int m=1;
	
	q=exp(-M_PI*Kratio);
	while (exp(-M_PI*Kratio*(double)(m*m-2*m+1))>eps && m < itermax)
	{
		denom+=2*exp(-M_PI*Kratio*((double)m-0.5)*((double)m-0.5));
		nom+=2*exp(-M_PI*Kratio*(double)(m*m));
		m++;
	}
	q=(denom/nom);
	return q*q;
}
#undef eps
#undef intermax

/* Fallback Test routine using bisection to find the elliptic modulus *//*
# define eps 1e-15
double ell_mod (double Kratio)
{
	double kmin=1e-8, kmax=1-1e-8,k, Emin, Emax, E;
	int i=0;
	
	Emin=Kratio-ell_int(sqrt(1-kmin*kmin))/ell_int(kmin);
	Emax=Kratio-ell_int(sqrt(1-kmax*kmax))/ell_int(kmax);
	do
	{
		k=(kmin+kmax)/2;
		E=Kratio-ell_int(sqrt(1-k*k))/ell_int(k);
		if (Emin*E>0)
		{
			Emin=E;
			kmin=k;
		}
		else
		{
			Emax=E;
			kmax=k;
		}
		
		i++;
	} while ((kmax-kmin)>eps && (fabs(E)>eps));
	return k;
}
#undef eps */

int scd2(double u,double mc, double * s ,double * c, double * d)
{
	double B10=1.0/24.0, B11=1.0/6.0, B20=1.0/720.0;
	double B21=11.0/180.0,B22=1.0/45.0;
	double m,uA,uT,u0,v,a,b,y,z,my,mc2,m2,x,xz,w;
	int n,i,j;
	m=1.0-mc;
	uA=1.76269+mc*1.16357;
	uT=5.217e-3-m*2.143e-3;
	u0=u;
	for (n=0;n<=20;n++)
	{
		if (u0<uT)
			break;
		u0/=2.0;
	}
	if (n==20)
	{
		fprintf(stderr,"error: argument u too large in scd2: u=%e",u);
		return 1;
	}
	v=u0*u0;
	a=1.0;
	b=v*(0.5-v*(B10+m*B11-v*(B20+m*(B21+m*B22))));
	if (u<uA)
	{
		for (j=0;j<n;j++)
		{
			y=b*(a*2.0-b);
			z=a*a;
			my=m*y;
			b=y*2.0*(z-my);
			a=z*z-my*y;
		}
	}
	else
	{
		for (j=0;j<n;j++)
		{
			y=b*(a*2.0-b);
			z=a*a;
			my=m*y;
			if (z<my*2.0)
			{
				*c=a-b;
				mc2=mc*2.0;
				m2=m*2.0;
				for (i=j;i<=n;i++)
				{
					x=(*c)*(*c);
					z=a*a;
					w=m*x*x-mc*z*z;
					xz=x*z;
					*c=mc2*xz+w;
					a=m2*xz-w;
					
				}
				*c/=a;
				x=(*c)*(*c);
				*s=sqrt(1.0-x);
				*d=sqrt(mc+m*x);
				return 0;
			}
			b=y*2.0*(z-my);
			a=z*z-my*y;
		}
	
	}
	b/=a;
	y=b*(2.0-b);
	*c=1.0-b;
	*s=sqrt(y);
	*d=sqrt(1.0-m*y);
	return 0;	

}

int scd(double u, double K, double k, double *s, double *c, double *d)
{
	double m;
	int Ndiv2=0;
	int N4K=0;
	int K2=0;
	int NEG=0;
	m=k*k;
	
	/* gsl_sf_elljac_e (u, m, s, c, d); */
	if (u<0)
	{
		NEG=1;
		u=fabs(u);
	}
	N4K=(int)floor(u/(4*K));
	u-=4*K*((double)N4K);
	
	if(u>2.0*K)
	{
		u-=2.0*K;
		K2=1;
	}
	
	while (u>=K/2.0)
	{
		u/=2;
		Ndiv2++;
	}
	
	if(scd2(u,1-m, s ,c, d))
		return 1;
	else
	{
		if (Ndiv2)
		{
			double ss,cc,dd;
			double nom;			
			while (Ndiv2)
			{
				nom=(1-m*(*s)*(*s)*(*s)*(*s));		
				ss=2*(*s)*(*c)*(*d)/nom;
				cc=(1-2*(*s)*(*s)+m*(*s)*(*s)*(*s)*(*s))/nom;
				dd=(1-2*m*(*s)*(*s)+m*(*s)*(*s)*(*s)*(*s))/nom;
				(*s)=ss;
				(*c)=cc;
				(*d)=dd;
				Ndiv2--;		
			}
		}
		if (K2)
		{
			(*s)=-(*s);
			(*c)=-(*c);
		}
		if (NEG)
		{
			(*s)=-(*s);
			(*c)=(*c);
		}
	}
	return 0;
}

int sn_i(double u, double v, double K, double k, double *sn_re, double *sn_im, double *dsn_redu, double *dsn_redv,double *dsn_imdu,double *dsn_imdv)
{
	double s,c,d,s1,c1,d1;
	double kk, KK, nom;
	if(scd(u, K, k, &s, &c, &d))
		return 1;
	kk=sqrt(1-k*k);
	KK=ell_int(kk);
	if (scd(v, KK, kk, &s1, &c1, &d1))
		return 1;
	nom=1-d*d*s1*s1;
	
	(*sn_re)=s*d1/nom;
	(*sn_im)=c*d*s1*c1/nom;
	
	if (dsn_redu!=NULL)
		(*dsn_redu)=c*d*d1/nom-2*(k*k)*c*d*(s*s)*d1*(s1*s1)/(nom*nom);
	if (dsn_imdu!=NULL)
		(*dsn_imdu)=-(d*d)*s*c1*s1/nom-(k*k)*(c*c)*s*c1*s1/nom-2*(k*k)*(c*c)*(d*d)*s*c1*(s1*s1*s1)/(nom*nom);
	if (dsn_redv!=NULL)
		(*dsn_redv)=2*(d*d)*s*c1*(d1*d1)*s1/(nom*nom)-(kk*kk)*s*c1*s1/nom; 
	if (dsn_imdv!=NULL)
		(*dsn_imdv)=-c*d*d1*(s1*s1)/nom+c*d*(c1*c1)*d1/nom+2*c*(d*d*d)*(c1*c1)*d1*(s1*s1)/(nom*nom);
	return 0;
}
