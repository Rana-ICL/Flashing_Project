 clc;
 clear all;
 close all;
  
%initial values boundary condition% 
                    
Tb=73;                           
deltahv=217000;                  
Tinj=97;                         
Tamb=100;                         
v_d=45;
p_amb=0.110;
droplet_d=30e-06;               

%Reference Temperaure%

Tref=(2/3)*Tb+(1/3)*Tamb;       

%Calculation of gas properties at reference Temperature%

hb=9.6769e-06*Tb^4-0.082121*Tb^3+15.491*Tb^2-90.222*Tb+21421;
cp= 5.5376e-07*Tref^4-0.00037957*Tref^3+0.094696*Tref^2-10.153*Tref+1304.2;                    
ktc=4.1868e-13*Tref^4-2.6237e-10*Tref^3+8.0071e-09*Tref^2+0.00010188*Tref-0.0010885;                
mu= 1.8542e-16*Tref^4-6.8643e-14*Tref^3-5.6412e-11*Tref^2+8.8871e-08*Tref-5.6932e-07;                                                                                                                                         
hamb=4.0786e-07*Tref^4-0.0001704*Tref^3+0.012296*Tref^2+913.19*Tref-829.35;              
rho =9.4899e-101*Tamb^4-6.8501e-07*Tamb^3+0.00019183*Tamb^2-0.02597*Tamb+1.6936;
Pr=(mu*cp)/ktc;                                                                       


% simulating the reduction of superheat with intial delta T = 29 K%
%Assumption the superheat reduces 1K every time step%


 for k=1:24

  Tinjneu(k)=97-k;                                     
  deltat(k)=(Tinjneu(k)-Tb);                           

  if (deltat(k)>= 25)                                   
      alpha(k)=13.8*((deltat(k))^0.39);                
  
  elseif  (5<=deltat(k)<25)
      alpha(k)=0.027*((deltat(k))^2.33);                
          
  else  
     alpha(k)=0.76*((deltat(k))^0.26);                  
         
  end
    
 alphaneu=[50.6140067450742,49.9011959749860,49.1720927476472,48.4256757547041,44.3870167765003,40.1966279448231,36.2417041235893,32.5188157159495,29.0244269598294,25.7548873340318,22.7064217946134,19.8751196107489,17.2569215080159,14.8476047447673,12.6427656350376,10.6377988754845,8.82787281268690,7.20789946594158,5.77249764165603,4.51594673828569,3.43212766423040,2.51444532955883,1.75572372411015,1.14805811908165,1.089805988,1.011265331,0.9100838155,0.76,0];

  
 %Radius reduction of the dropet every time step%
  
  red=((droplet_d)/2)/29;
  r(k)=(droplet_d/2)-(k*red);                          
  
 %Area of the droplets at every time step %
  A(k)=4*pi*r(k)^2;                                     
 
 %mass vaporization due to flashing%
  mflash(k)=(A(k)*alphaneu(k)*1000*(deltat(k)))/deltahv;     
   
% Calculating reynolds number, nusselt number, transfer number and FT %
  
Re(k)=(rho*v_d*(2*r(k)))/mu;                            
 
%Standard nusselt number%  
Nuo(k)=2+0.552*((Re(k)^1/2)*Pr^1/3);   
  
Bt=(cp*(Tamb-Tb))/deltahv;              
  
Ft=((1+Bt)^0.7)*(log(1+Bt)/Bt);         

%Modified Nusselt number% 

Nu(k)=2+((Nuo(k)-2)/Ft);           


% Guess for a good initial value for the Bisection method%
%xl=1e-15;   %upper guess%    % change it according to the boundary conditions  %very important%
%xu=1e-10;   %lower guess%    % change it according to the boundary condtions  % %very important%

%up(k)=-xu+((2*pi*(ktc/cp)*r(k))*(Nu(k)/(1+(mflash(k)/xu)))*(log(1+(1+(mflash(k)/xu))*(((hamb-hb)/deltahv)))));
%low(k)=-xl+((2*pi*(ktc/cp)*r(k))*(Nu(k)/(1+(mflash(k)/xl)))*(log(1+(1+(mflash(k)/xl))*(((hamb-hb)/deltahv)))));


 end
 
 
% Bisection method for finding the initial valure for fixed point iteration%

for l=1:24;

  xl=1e-15;  % upper guess
  xu=1e-10;  % lower guess
    
  r1=r(l);
  Nu1=Nu(l);
  mflash1=mflash(l);

for u=1:20;

c=(xl+xu)/2;

fxl=-xl+((2*pi*(ktc/cp)*r(k))*(Nu1/(1+(mflash(k)/xl)))*(log(1+(1+(mflash(k)/xl))*(((hamb-hb)/deltahv))))); 
fc=-c+((2*pi*(ktc/cp)*r1)*(Nu1/(1+(mflash1/c)))*(log(1+(1+(mflash1/c))*(((hamb-hb)/deltahv)))));           
 
if ( fc == 0 )                                                                                             
       break;                                                                                              
  elseif (fxl*fc < 0 )
       xu = c;                                                                                             
    else
       xl = c;                                                                                             
end                                                                                                        
 
end

mitt(l)=c;   % initial guesses for the fixed iteration method%                                             

end

% Fixed point iteration method%                                
  
for we=1:24

    r12=r(we);
    Nu12=Nu(we);
    mflash12=mflash(we);
    mitt2=mitt(we);

for er=1:80
    tol=1e-60;
    x1=((2*pi*(ktc/cp)*r12)*(Nu12/(1+(mflash12/mitt2)))*(log(1+(1+(mflash12/mitt2))*(((hamb-hb)/deltahv)))));
    
    abserror=abs(x1-mitt2);   

    if (abserror<tol)
             break  

end
    mitt2=x1;
end
 % Mass vaporization due to heat conduction%
mittneu(we)=mitt2;  

end 

% check, whether the found values for the root are good  best case if the variable 'zero' is 0 for every time step%

for  tr=1:24
  
    r13=r(tr);
    Nu13=Nu(tr);
    mflash13=mflash(tr);
    mitt3=mittneu(tr);
    
    zero(tr)=-mitt3+((2*pi*(ktc/cp)*r13)*(Nu13/(1+(mflash13/mitt3)))*(log(1+(1+(mflash13/mitt3))*(((hamb-hb)/deltahv)))));
    
end

% plot of the droplet mass evaporation due to flashing%

figure
plot(deltat,mflash,'black')
title('Droplet mass evaporation due to flashing')
xlabel('{\Delta} T [K]')
ylabel('mdot-flash [kg/s]')

% plot of the droplet mass evaporation due to heat conduction%
figure
plot(deltat,mittneu,'black')
title('Droplet mass evaporation due to heat conduction')
xlabel('{\Delta} T [K]')
ylabel('mdot-con [kg/s]')

%total mass evaporation including flashing and heat conduction %
mres=mittneu+mflash;

% plot of total droplet mass evaporation%
figure
plot(deltat,mres,'black');
title('total droplet mass evaporation');
xlabel('{\Delta} T [K]');
ylabel('mdot-result [kg/s]');