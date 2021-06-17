/*This is unsed to simulate the flash spray of LOx at 200mbar...
The gas and liquid droplet properties are fitted with NIST database...
gas properties is unsed the isoproperities table with constant pressure
DPM properties is used the saturation table*/


#include "udf.h"
#include "dpm.h"
#include "mem.h"
#include "math.h"

/* The following law determines how the droplet loses its mass and enegry*/
/* This law determine the properties of a droplet e.g diameter,density,mass,temperature etc for a single droplet!*/

DEFINE_DPM_LAW(Flashing_Evaporation_Law,p,ci) 
{
Material*m = P_MATERIAL(p); 
cell_t c = P_CELL(p);                                                                
Thread *t = P_CELL_THREAD(p);     
real alpha;                                                                                                                                                                                                                                                         
real m_dot_flash;
real m_dot_heat;
real m_dot_radi;
real Tc=C_T(c,t); 
real Rho_c= C_R(c,t);  
real Tp=P_T(p);
real mp=P_MASS(p);
real Dp=P_DIAM(p);
real Ap=DPM_AREA(Dp);                                                                                                                                                                                                                                                                                                                                                                                                               
real Delta_L=DPM_LATENT_HEAT(p,m);   
real Tb=DPM_BOILING_TEMPERATURE(p,m);   
real cp_d = DPM_SPECIFIC_HEAT(p,Tp);     
real Delta_T=Tp-Tb;  
real u_p= P_VEL(p)[0];
real v_p= P_VEL(p)[1]; 
real w_p= P_VEL(p)[2];
real u_g=C_U(c,t);
real v_g=C_V(c,t);
real w_g=C_W(c,t);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
real xn, xn1, xneu, fx, fxl, fc, dfx, tol, abserror, guess, cc, Er;
real xl=1e-15;                                                                                                                                                                                                                        
real xu=1e-10; 
int n;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
int k=1;                                                                                                                                          		                        				 																	                                              
int l=1; 
real T_wall=230;

/* Calculation of gas phase properties at specific refer temperature, using rule one-third*/
real T_spe= (1./3.)*Tc+(2./3.)*Tb;                                                                                                                                                                                   
real visc_mu_spe= 1.8542e-16*pow(T_spe,4)-6.8643e-14*pow(T_spe,3)-5.6412e-11*pow(T_spe,2)+8.8871e-08*T_spe-5.6932e-07;                                                                                                                                                 
real cond_lambda_spe=4.1868e-13*pow(T_spe,4)-2.6237e-10*pow(T_spe,3)+8.0071e-09*pow(T_spe,2)+0.00010188*T_spe-0.0010885;                                                    
real cp_spe= 5.5376e-07*pow(T_spe,4)-0.00037957*pow(T_spe,3)+0.094696*pow(T_spe,2)-10.153*T_spe+1304.2;   
real h_amb=4.0786e-07*pow(T_spe,4)-0.0001704*pow(T_spe,3)+0.012296*pow(T_spe,2)+913.19*T_spe-829.35;
real h_b=9.6769e-06*pow(Tb,4)-0.082121*pow(Tb,3)+15.491*pow(Tb,2)-90.222*Tb+21421;                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
real Re=Rho_c*sqrt((u_p-u_g)*(u_p-u_g)+(v_p-v_g)*(v_p-v_g)+(w_p-w_g)*(w_p-w_g))*Dp/visc_mu_spe; 
real Pr=visc_mu_spe*cp_spe/cond_lambda_spe;                                                                                                                                                                                               
real BT=cp_spe*(Tc-Tb)/Delta_L;                                                                                                                                                                                                                                                           
real FT=pow((1+BT),0.7)*log(1+BT)/BT;   
real Nu=2.0+0.552*sqrt(Re)*pow(Pr,1./3.)/FT; 
real a=M_PI*cond_lambda_spe/cp_spe*Dp*Nu;
real b= (h_amb-h_b)/Delta_L;

/*Superheat coffecient according to Adachi's correlation*/             
if (Delta_T >= 25)  

	{alpha=13.8*pow(Delta_T, 0.39); }                                                                                                                                        

else if (5 < Delta_T && Delta_T < 25) 

	{alpha= 0.027*pow(Delta_T, 2.33);}                                                                                                                                                                             

else if (0 < Delta_T && Delta_T <=5) 

	{alpha= 0.76*pow(Delta_T, 0.26); }  

/* Calculation of mass vaporization due to internal superheat->Flashing (m_dot_flash) */
m_dot_flash=Ap*alpha*Delta_T/(Delta_L/1000);                                                                                                                                                                           


/* Calculation of mass vaporization due to external heat conduction and convection->(m_dot_heat) */
/* First step-> Bisection method for finding the initial guess for fixed point iteration method*/
for (k;k<5;k++)                                                                                                              
{
	cc=(xl+xu)/2;                                                                                                             
	fxl=-xl+a/(1+m_dot_flash/xl)*log(1+(1+m_dot_flash/xl)*b);
	fc=-cc+a/(1+m_dot_flash/cc)*log(1+(1+m_dot_flash/cc)*b);   
 
	if ( fc == 0 )                                                                                                                                                                                    
		{break;}                                                                                               
	else if (fxl*fc < 0 )                                                                                  
		{xu = cc;}                                                                                         
	else                                                                                                  
		{xl = cc;}  
}  

guess=cc;                             

/* Second Step-> Now use the guess from bisection method as the initial value for the fixed point iteration method*/
do                                                                               
{  
	xneu=a/(1+m_dot_flash/guess)*log(1+(1+m_dot_flash/guess)*b);   
    abserror = fabs(xneu-guess);                                                                                                                
    guess=xneu; 
} while (abserror>1e-17);

if(Tc>Tb)
	{m_dot_heat=xneu;}
else
	{m_dot_heat=0;}


/* Calculation of mass vaporization due to raditon (real m_dot_radi)*/
m_dot_radi=0.96*5.670367e-08*Ap*(1*pow(T_wall,4)-1*pow(Tb,4))/Delta_L;


/* Determine the new particle mass, diameter and temperature*/
if ( P_MASS(p)> 0 && P_T(p)>Tb)                                                     
{  
    P_MASS(p)-=(m_dot_heat+m_dot_flash+ m_dot_radi)*P_DT(p); 
    P_T(p)-=m_dot_flash*Delta_L/(cp_d*P_MASS(p))*P_DT(p);
    P_RHO(p)=-1.2945e-05*pow(Tp,4)+0.0048875*pow(Tp,3)-0.69914*pow(Tp,2)+39.693*Tp+518.69;
    P_DIAM(p)=pow(6.0*P_MASS(p)/ (P_RHO(p)*M_PI ), 1./3.);   
}
//Message( " %e %e %e %e %e %e %e\n",u_p,v_p, w_p, u_g, v_g, w_g, Re); 
}


/*Calculation of the temperature dependent properties of the droplet*/
/* Determine the temperature dependent density of the droplet*/
DEFINE_DPM_PROPERTY(droplet_density,c,t,p,T)     
{ 
real rho;
rho=-1.2945e-05*pow(T,4)+0.0048875*pow(T,3)-0.69914*pow(T,2)+39.693*T+518.69;
return rho;
}

/*Determine the temperature dependent specific heat capacity of the droplet*/
DEFINE_DPM_PROPERTY(droplet_cp,c,t,p,T)                           
{
real mp0= P_INIT_MASS(p);
real mp = P_MASS(p);
real D_cp=0.00060785*pow(T,4)-0.24247*pow(T,3)+35.77*pow(T,2)-2305.9*T+56387;
p->enthalpy = D_cp*(T-T_REF);
return D_cp;
}

/*Determine the temperature dependent surface tension of the droplet*/
DEFINE_DPM_PROPERTY(droplet_surftension,c,t,p,T)                           
{
real D_surf_sigma;
D_surf_sigma=4.4486e-11*pow(T,4)-1.4839e-08*pow(T,3)+2.2499e-06*pow(T,2)-4.2403e-04*T+0.041034;
return D_surf_sigma;
}

/*Determine the temperature dependent viscosity of the droplet*/
DEFINE_DPM_PROPERTY(droplet_viscosity,c,t,p,T)                           
{
real D_visc_mu;
D_visc_mu=2.4046e-11*pow(T,4)-1.1587e-08*pow(T,3)+2.0953e-06*pow(T,2)-1.7086e-04*T+ 0.0054666;
return D_visc_mu;
}

/*Determine the temperature dependent properties of the gas phase*/ 
/*Calculation of gas phase viscosity */
DEFINE_PROPERTY(cell_viscosity,c,t)                                                                                                                    
{                                                                                                                                                                                                                                                        
real visc_mu;                                                                                            
real Tc=C_T(c,t);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
visc_mu= 1.8542e-16*pow(Tc,4)-6.8643e-14*pow(Tc,3)-5.6412e-11*pow(Tc,2)+8.8871e-08*Tc-5.6932e-07;                  

return visc_mu;                                                                                                                                  
}   

/* Calculation of gas phase thermal conductivity */
DEFINE_PROPERTY(cell_thermalconductiviy,c,t)                                                                                                        
{                                                                                                     
real cond_lambda;                                                                                                                                                                                                                                        
real Tc=C_T(c,t);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
cond_lambda=4.1868e-13*pow(Tc,4)-2.6237e-10*pow(Tc,3)+8.0071e-09*pow(Tc,2)+0.00010188*Tc-0.0010885;         
return cond_lambda;  
} 

/*Calculation of gas phase specific heat*/
DEFINE_SPECIFIC_HEAT(cp_gas,T,Tref,h,yi)                                                                                                                
{                                                                                            
real cp;   
cp=5.5376e-07*pow(T,4)-0.00037957*pow(T,3)+0.094696*pow(T,2)-10.153*T+1304.2;                        			
*h= cp*(T-Tref);                                                                                                                                                                 
return cp;                                                                                                                      
}

                                               
DEFINE_DPM_DRAG(DPM_dragforce,Re,p)
{
 Material*m = P_MATERIAL(p); 
 cell_t c = P_CELL(p);                                                                
 Thread *t = P_CELL_THREAD(p);  
 real dragforce;
 real cd0;
 real cd;
 real Tc=C_T(c,t);   
 real Tb=DPM_BOILING_TEMPERATURE(p,m); 
 real Delta_L=DPM_LATENT_HEAT(p,m);
 real T_spe= (1./3.)*Tc+(2./3.)*Tb;                                                                                                                                                                                                                                   
 real cp_spe = 5.5376e-07*pow(T_spe,4)-0.00037957*pow(T_spe,3)+0.094696*pow(T_spe,2)-10.153*T_spe+1304.2; 
 //real cp_spe=C_CP(c,t); 
 real Bm=cp_spe*(Tc-Tb)/Delta_L;  
 
 if (Re< 1000.0)
{
	cd0=24/Re*(1+0.15*pow(Re, 0.687));
	cd=cd0/(1+Bm);
	dragforce = 18*cd*Re/24;
return (dragforce);
}

else
{
	cd0 = 0.44 ;
	cd=cd0/(1+Bm);
	dragforce = 18*cd*Re/24;
return (dragforce);
}

}

/*The following UDF defines the source terms to the gas phase*/
DEFINE_DPM_SOURCE(dpm_source,c,t,S,strength,p)                               
{
 Material *m=P_MATERIAL(p);   
 int n2=0;                                                                                                                                                                                        
 real m_dot_total, heatsource;                                                                                                                                                                                                          
 real Dp=P_DIAM(p);
 real P_rho=P_RHO(p);                                                                                                                                                                                                                                                                                                    
 real Tp_entr = P_T0(p);
 real Tp=P_T(p);
 real Tb = DPM_BOILING_TEMPERATURE(p,m);   
 real Delta_L=DPM_LATENT_HEAT(p,m);                                                                       
 real Delta_T = Tp_entr-Tb;                         
 real cp_d = DPM_SPECIFIC_HEAT(p,Tp);

 /**********************/
 real cd, cd0, Re;
 real Tc=C_T(c,t); 
 real rho_c= C_R(c,t); 
 real u_p= P_VEL(p)[0];
 real v_p= P_VEL(p)[1]; 
 real w_p= P_VEL(p)[2];
  real u_g=C_U(c,t);
 real v_g=C_V(c,t);
 real w_g=C_W(c,t); 
 real T_spe= (1./3.)*Tc+(2./3.)*Tb;                                                                                                                                                                                                                                   
 real cp_spe = 5.5376e-07*pow(T_spe,4)-0.00037957*pow(T_spe,3)+0.094696*pow(T_spe,2)-10.153*T_spe+1304.2;     
 real Bm=cp_spe*(Tc-Tb)/Delta_L; 
 Re=rho_c*sqrt((u_p-u_g)*(u_p-u_g)+(v_p-v_g)*(v_p-v_g)+(w_p-w_g)*(w_p-w_g))*Dp/C_MU_T(c,t);

 /* Mass source to the gas phase*/
/*(The mass source of a single dropelt is the difference in mass between the entry  of droplet to a cell and exit from this cell)*/
/*Total mass source from a cell is mass source of single dropelt multiplied by strength,which is number of particles/s in stream*/
m_dot_total=(P_MASS0(p) - P_MASS(p))*strength;                                                                                                      

/*Energy source term to the gas phase*/ 
heatsource=(P_MASS0(p)*cp_d*(P_T0(p)-T_REF)-cp_d*P_MASS(p)*(P_T(p)-T_REF)-(P_MASS0(p)-P_MASS(p))*P_INJECTION(p)->latent_heat_ref)*strength;
 
/*Switching between the boiling law and Flash_Evaporation_Law*/
if (P_CURRENT_LAW(p) == DPM_LAW_USER_1)
{
	S->species [n2]+= m_dot_total;
	S->mass = m_dot_total;
	S->energy = heatsource;

	if (Re< 1000.0)
	{	cd0=24/Re*(1+0.15*pow(Re, 0.687));
		cd=cd0/(1+Bm);
	}
	else
	{	cd0 = 0.44 ;
		cd=cd0/(1+Bm);
	}

	S->momentum_s[0]+=18/24*C_MU_T(c,t)*cd*Re/P_rho/Dp/Dp*(u_p-u_g)*P_DT(p)*P_FLOW_RATE(p);
	S->momentum_s[1]+=18/24*C_MU_T(c,t)*cd*Re/P_rho/Dp/Dp*(v_p-v_g)*P_DT(p)*P_FLOW_RATE(p);
	S->momentum_s[2]+=18/24*C_MU_T(c,t)*cd*Re/P_rho/Dp/Dp*(w_p-w_g)*P_DT(p)*P_FLOW_RATE(p);
	
}
//Message( " %e %e  %e %e %e %e %e\n",heatsource,S->energy, S->mass,S->momentum_s[0],S->momentum_s[1],S->momentum_s[2], Re); 

}


/*If particle temperature is higher than the boiling temperature-> switch to Flashing law*/
/*if not then use depending on the temperature the boiling law or vaporitzation law*/
DEFINE_DPM_SWITCH(dpm_switch, p, coupled)
{
cell_t c = P_CELL(p);
Thread *t = P_CELL_THREAD(p);
Material *m = P_MATERIAL(p);

if ( P_T(p) > DPM_BOILING_TEMPERATURE(p, m))
	P_CURRENT_LAW(p) = DPM_LAW_USER_1;
else
	P_CURRENT_LAW(p) = DPM_LAW_BOILING;
}