%%%Resonance driving term calculation program
%%%read magnets list file generated from drivingTermsCal 
%%%in nonlinear directory in linux
%%%author: T. Pulampong 13/11/12
%%%--> with shorter formula (summation notation)
%%% use k2l from processed file instead of k*l
%%% improved calculation with data output with more precision*
%%% *used script txtGen generate *.txt instead of normal sddsprintout
%%%--------
%%% edited 1 by T.P.23/7/13
%%% add 1st order 
%%% all optics and kl processed from drivingtersCal script
%%% (split elements and averaged optics over each element)
%%% ***Automatically generated out-*.sumMagnets and out-*.sumMagnets2 by 
%%% drivingtermsCal script and working Tunes are required***
%%%--------
%%% edited 2 by T.P. 16/8/13
%%% add the higher order chromatic terms
%%%--------

%%%select data text file
 filename = uigetfile('out*.sumMagnets','Select the data file')
 %read quadrupole and sextupole parameters
 data   = textread(filename);
 s1     = data(:,1);
 betax1 = data(:,2);
 betay1 = data(:,3);
 etax1  = data(:,4);
 mux1   = data(:,5);
 muy1   = data(:,6);
 kl     = data(:,7); 
 typeNo = data(:,8);  %typeNo=1 for quadrupole,2 for sextupole,3 for dipole
 %read sextupole parameters from filename2
 data2 = textread(sprintf('%s2',filename));
 s     = data2(:,1);
 betax = data2(:,2);
 betay = data2(:,3);
 etax  = data2(:,4);
 mux   = data2(:,5);
 muy   = data2(:,6);
 k2l   = data2(:,7); %for sextupoles only

%%%plese put the right tunes for your lattice
disp('input tunes:')
Nux = input('TuneX:');
Nuy = input('TuneY:');


%%%firt order chromatic terms
%%%use sextupole and quadrupole in differrent itteration 

H11001sum=0; 
H00111sum=0; 
H20001sum=0; 
H00201sum=0;
H10002sum=0;

for iqs=1:length(typeNo)
    if typeNo(iqs)==1||typeNo(iqs)==3
       H11001sum=-0.25*kl(iqs)*betax1(iqs) +H11001sum; 
       H00111sum= 0.25*kl(iqs)*betay1(iqs) +H00111sum ; 
       H20001sum= 0.125*kl(iqs)*betax1(iqs)*exp(1i*2*mux1(iqs)) +H20001sum; 
       H00201sum=-0.125*kl(iqs)*betay1(iqs)*exp(1i*2*muy1(iqs)) +H00201sum;
       H10002sum= 0.5*kl(iqs)*etax1(iqs)*sqrt(betax1(iqs))*exp(1i*mux1(iqs)) +H10002sum;        
    else
       H11001sum= 0.25*kl(iqs)*betax1(iqs)*etax1(iqs) +H11001sum; 
       H00111sum=-0.25*kl(iqs)*betay1(iqs)*etax1(iqs) +H00111sum ; 
       H20001sum=-0.125*kl(iqs)*betax1(iqs)*etax1(iqs)*exp(1i*2*mux1(iqs)) +H20001sum; 
       H00201sum= 0.125*kl(iqs)*betay1(iqs)*etax1(iqs)*exp(1i*2*muy1(iqs)) +H00201sum;
       H10002sum=-0.25*kl(iqs)*etax1(iqs)^2*sqrt(betax1(iqs))*exp(1i*mux1(iqs)) +H10002sum;  
    end
 h11001sum(iqs)=H11001sum;
 h00111sum(iqs)=H00111sum;
 h20001sum(iqs)=H20001sum;
 h00201sum(iqs)=H00201sum;
 h10002sum(iqs)=H10002sum;
end


%%%higher order off-momentum terms
%%%use sextupole and quadrupole strength at the same time
disp('Calculating higher order off-momentum Driving terms...')
 h21001sumre=0; h21001sumim=0;
 h21001sum=0;
 h30001sum=0;
 h10021sum=0;
 h10111sum=0;
 h10201sum=0;
 h11002sum=0;
 h20002sum=0;
 h00112sum=0;
 h00202sum=0;
 h10003sum=0;
 h00004sum=0;
for jj=1:length(typeNo)    
    if typeNo(jj)==1||typeNo(jj)==3
       klqj=kl(jj); 
       klsj=0;
    else
       klqj=0;
       klsj=kl(jj)/2;
    end
        sumH21001re=0;     sumH21001im=0;               
    sumH21001=0;    
    sumH30001=0;
    sumH10021=0;
    sumH10111=0;
    sumH10201=0;
    sumH11002=0;
    sumH20002=0;
    sumH00112=0;
    sumH00202=0;
    sumH10003=0;
    sumH00004=0;
    for ii=1:length(typeNo)
        if typeNo(ii)==1||typeNo(ii)==3
           klqi=kl(ii); 
           klsi=0;
        else
           klqi=0;
           klsi=kl(ii)/2;
        end
                
        termSign=sign(s1(jj)-s1(ii)); %important bit controling position i<j or i>j
      
        H21001=termSign*( (-1i/32)*klsi*klqj*betax1(ii)^1.5*betax1(jj)*(exp(1i*mux1(ii))+exp(1i*(3*mux1(ii)-2*mux1(jj)))-2*exp(-1i*(mux1(ii)-2*mux1(jj))))...
               -(1i/16)*klsi*klsj*betax1(ii)*betax1(jj)^1.5*etax1(ii)*(exp(1i*mux1(jj))-2*exp(1i*(2*mux1(ii)-mux1(jj)))+exp(-1i*(2*mux1(ii)-3*mux1(jj)))) );
        H21001im=termSign*( (-1/32)*klsi*klqj*betax1(ii)^1.5*betax1(jj)*(cos(mux1(ii))+cos(3*mux1(ii)-2*mux1(jj))-2*cos(mux1(ii)-2*mux1(jj)))...
               -(1/16)*klsi*klsj*betax1(ii)*betax1(jj)^1.5*etax1(ii)*(cos(mux1(jj))-2*cos((2*mux1(ii)-mux1(jj)))+cos((2*mux1(ii)-3*mux1(jj)))) );
       H21001re=termSign*( (1/32)*klsi*klqj*betax1(ii)^1.5*betax1(jj)*(sin(mux1(ii))+sin(3*mux1(ii)-2*mux1(jj))+2*sin(mux1(ii)-2*mux1(jj)))...
               +(1/16)*klsi*klsj*betax1(ii)*betax1(jj)^1.5*etax1(ii)*(sin(mux1(jj))-2*sin((2*mux1(ii)-mux1(jj)))-sin((2*mux1(ii)-3*mux1(jj)))) );
       
        H30001=termSign*( (-1i/32)*klsi*klqj*betax1(ii)^1.5*betax1(jj)*(exp(1i*3*mux1(ii))-exp(1i*(mux1(ii)+2*mux1(jj))))...
               -(1i/16)*klsi*klsj*betax1(ii)*betax1(jj)^1.5*etax1(ii)*(exp(1i*3*mux1(jj))-exp(1i*(2*mux1(ii)+mux1(jj)))) );
        H10021=termSign*( (1i/32)*klsi*klqj*sqrt(betax1(ii))*betax1(jj)*betay1(ii)*(exp(1i*(mux1(ii)-muy1(ii)))-exp(-1i*(mux1(ii)-2*mux1(jj)+2*muy1(ii))))...
               +(1i/16)*klsi*klqj*sqrt(betax1(ii))*betay1(jj)*betay1(ii)*(exp(1i*(mux1(ii)-2*muy1(ii)))-exp(1i*(mux1(ii)-2*muy1(jj))))...
               +(1i/16)*klsi*klsj*betax1(ii)*sqrt(betax1(jj))*betay1(jj)*etax1(ii)*(exp(1i*(mux1(jj)-2*muy1(jj)))-exp(1i*(2*mux1(ii)-mux1(jj)-2*muy1(jj))))...        
               -(1i/8) *klsi*klsj*sqrt(betax1(jj))*betay1(ii)*betay1(jj)*etax1(ii)*(exp(1i*(mux1(jj)-2*muy1(ii)))-exp(1i*(mux1(jj)-2*muy1(jj)))) );        
        H10111=termSign*( (1i/16)*klsi*klqj*sqrt(betax1(ii))*betax1(jj)*betay1(ii)*(exp(1i*(mux1(ii)))-exp(-1i*(mux1(ii)-2*mux1(jj))))...
               +(1i/16)*klsi*klqj*sqrt(betax1(ii))*betay1(jj)*betay1(ii)*(exp(1i*(mux1(ii)-2*muy1(ii)+2*muy1(jj)))-exp(1i*(mux1(ii)+2*muy1(ii)-2*muy1(jj))))...
               +(1i/8) *klsi*klsj*betax1(ii)*sqrt(betax1(jj))*betay1(jj)*etax1(ii)*(exp(1i*(mux1(jj)))-exp(1i*(2*mux1(ii)-mux1(jj))))...        
               -(1i/8) *klsi*klsj*sqrt(betax1(jj))*betay1(ii)*betay1(jj)*etax1(ii)*(exp(1i*(mux1(jj)-2*muy1(ii)+2*muy1(jj)))-exp(1i*(mux1(jj)+2*muy1(ii)-2*muy1(jj)))) );     
        H10201=termSign*( (1i/32)*klsi*klqj*sqrt(betax1(ii))*betax1(jj)*betay1(ii)*(exp(1i*(mux1(ii)+2*muy1(ii)))-exp(-1i*(mux1(ii)-2*mux1(jj)-2*muy1(ii))))...
               -(1i/16)*klsi*klqj*sqrt(betax1(ii))*betay1(jj)*betay1(ii)*(exp(1i*(mux1(ii)+2*muy1(ii)))-exp(1i*(mux1(ii)+2*muy1(jj))))...
               +(1i/16)*klsi*klsj*betax1(ii)*sqrt(betax1(jj))*betay1(jj)*etax1(ii)*(exp(1i*(mux1(jj)+2*muy1(jj)))-exp(1i*(2*mux1(ii)-mux1(jj)+2*muy1(jj))))...        
               +(1i/8) *klsi*klsj*sqrt(betax1(jj))*betay1(ii)*betay1(jj)*etax1(ii)*(exp(1i*(mux1(jj)+2*muy1(ii)))-exp(1i*(mux1(jj)+2*muy1(jj)))) );   
        H11002=termSign*( (1i/16)*betax1(ii)*betax1(jj)*( (klqi*klqj-2*klsi*klqj*etax1(ii)+4*klsi*klsj*etax1(ii)*etax1(jj))*exp(1i*2*(mux1(ii)-mux1(jj)))+2*klsi*klqj*etax1(ii)*exp(-1i*2*(mux1(ii)-mux1(jj))) )...
               +(1i/8)*sqrt(betax1(ii))*betax1(jj)^1.5*etax1(ii)*(klsi*etax1(ii)-klqi)*klsj*(exp(1i*(mux1(ii)-mux1(jj)))-exp(-1i*(mux1(ii)-mux1(jj)))) ); 
        H20002=termSign*( (1i/16)*betax1(ii)*betax1(jj)*( (klqi*klqj-2*klsi*klqj*etax1(ii)+4*klsi*klsj*etax1(ii)*etax1(jj))*exp(1i*2*(mux1(ii)))+2*klsi*klqj*etax1(ii)*exp(1i*2*(mux1(jj))) )...
               +(1i/16)*sqrt(betax1(ii))*betax1(jj)^1.5*etax1(ii)*(klsi*etax1(ii)-klqi)*klsj*(exp(1i*(mux1(ii)+mux1(jj)))-exp(-1i*(mux1(ii)-3*mux1(jj)))) ); 
        H00112=termSign*( (1i/16)*betax1(ii)*betax1(jj)*( (klqi*klqj-2*klsi*klqj*etax1(ii)+4*klsi*klsj*etax1(ii)*etax1(jj))*exp(1i*2*(muy1(ii)-muy1(jj)))+2*klsi*klqj*etax1(ii)*exp(-1i*2*(muy1(ii)-muy1(jj))) )...
               -(1i/8)*sqrt(betax1(ii)*betax1(jj))*betay1(jj)*etax1(ii)*(klsi*etax1(ii)-klqi)*klsj*(exp(1i*(mux1(ii)-mux1(jj)))-exp(-1i*(mux1(ii)-mux1(jj)))) ); 
        H00202=termSign*( (1i/16)*betax1(ii)*betax1(jj)*( (klqi*klqj-2*klsi*klqj*etax1(ii)+4*klsi*klsj*etax1(ii)*etax1(jj))*exp(1i*2*(muy1(ii)))+2*klsi*klqj*etax1(ii)*exp(1i*2*(muy1(jj))) )...
               -(1i/16)*sqrt(betax1(ii)*betax1(jj))*betay1(jj)*etax1(ii)*(klsi*etax1(ii)-klqi)*klsj*(exp(1i*(mux1(ii)-mux1(jj)+2*muy1(jj)))-exp(-1i*(mux1(ii)-mux1(jj)-2*muy1(jj)))) );            
        H10003=termSign*( (1i/4)*klsi*klqj*betax1(ii)*sqrt(betax1(jj))*etax1(ii)*etax1(jj)*(exp(1i*(mux1(jj)))-exp(1i*(2*mux1(ii)-mux1(jj))))...
               +(1i/8)*sqrt(betax1(ii))*betax1(jj)*etax1(ii)*(klqi*klqj-klsi*etax1(ii)*(klqj-2*klsj*etax1(jj)))*(exp(1i*(mux1(ii)))-exp(-1i*(mux1(ii)-2*mux1(jj)))) );       
        H00004=termSign*( (1i/4)*sqrt(betax1(ii)*betax1(jj))*etax1(ii)*etax1(jj)*((klqi*klqj-klsi*etax1(ii)*(klqj-klsj*etax1(jj)))*exp(1i*(mux1(ii)-mux1(jj)))+klsi*klqj*etax1(ii)*exp(-1i*(mux1(ii)-mux1(jj)))) );   
           
        sumH21001re=sumH21001re+H21001re;
        sumH21001im=sumH21001im+H21001im; 
        sumH21001=sumH21001+H21001;    
        sumH30001=sumH30001+H30001;  
        sumH10021=sumH10021+H10021;  
        sumH10111=sumH10111+H10111;  
        sumH10201=sumH10201+H10201;  
        sumH11002=sumH11002+H11002; 
        sumH20002=sumH20002+H20002;        
        sumH00112=sumH00112+H00112;          
        sumH00202=sumH00202+H00202;          
        sumH10003=sumH10003+H10003;         
        sumH00004=sumH00004+H00004;        
    end
    h21001sumre=h21001sumre+sumH21001re; h21001SUMre(jj)=h21001sumre;
    h21001sumim=h21001sumim+sumH21001im; h21001SUMim(jj)=h21001sumim;
    h21001sum=h21001sum+sumH21001; h21001SUM(jj)=h21001sum;
    h30001sum=h30001sum+sumH30001; h30001SUM(jj)=h30001sum;
    h10021sum=h10021sum+sumH10021; h10021SUM(jj)=h10021sum;
    h10111sum=h10111sum+sumH10111; h10111SUM(jj)=h10111sum;
    h10201sum=h10201sum+sumH10201; h10201SUM(jj)=h10201sum;
    h11002sum=h11002sum+sumH11002; h11002SUM(jj)=h11002sum;
    h20002sum=h20002sum+sumH20002; h20002SUM(jj)=h20002sum;    
    h00112sum=h00112sum+sumH00112; h00112SUM(jj)=h00112sum;    
    h00202sum=h00202sum+sumH00202; h00202SUM(jj)=h00202sum;    
    h10003sum=h10003sum+sumH10003; h10003SUM(jj)=h10003sum;    
    h00004sum=h00004sum+sumH00004; h00004SUM(jj)=h00004sum;
end




%%%use sextupole strength only
disp('Calculating 2nd order Driving terms...')
 H21000sum=0;
 H30000sum=0;
 H10110sum=0;
 H10020sum=0;
 H10200sum=0;

h31000sum=0;
h22000sum=0;
h11110sum=0;
h11200sum=0;
h40000sum=0;
h20020sum=0;
h20110sum=0;
h20200sum=0;
h00220sum=0;
h00310sum=0;
h00400sum=0;
dnux_dJxsum=0;
dnux_dJysum=0;
dnuy_dJysum=0;
for jj=1:length(k2l)    
    H21000sum=H21000sum-(1/16)*k2l(jj)*(betax(jj)^1.5)*exp(1i*mux(jj));
    H30000sum=H30000sum-(1/48)*k2l(jj)*(betax(jj)^1.5)*exp(1i*3*mux(jj));
    H10110sum=H10110sum+(1/8)*k2l(jj)*sqrt(betax(jj))*betay(jj)*exp(1i*mux(jj));
    H10020sum=H10020sum+(1/16)*k2l(jj)*sqrt(betax(jj))*betay(jj)*exp(1i*(mux(jj)-2*muy(jj)));
    H10200sum=H10200sum+(1/16)*k2l(jj)*sqrt(betax(jj))*betay(jj)*exp(1i*(mux(jj)+2*muy(jj)));
    h21000sum(jj)=H21000sum;
    h30000sum(jj)=H30000sum;
    h10110sum(jj)=H10110sum;
    h10020sum(jj)=H10020sum;
    h10200sum(jj)=H10200sum;
    
    sumH31000=0;    
    sumH22000=0;
    sumH11110=0;
    sumH11200=0;
    sumH40000=0;
    sumH20020=0;
    sumH20110=0;
    sumH20200=0; 
    sumH00220=0;
    sumH00310=0;
    sumH00400=0;
    sumdnux_dJx=0;
    sumdnux_dJy=0;
    sumdnuy_dJy=0;
  %  disp(ii)
    for ii=1:length(k2l)
        termSign=sign(s(jj)-s(ii)); %important bit controling position i<j or i>j
    H31000=termSign*((1i/128)*k2l(ii)*k2l(jj)*(betax(ii)*betax(jj))^1.5)*(exp(1i*(3*mux(ii)-mux(jj))));     
    H22000=termSign*((1i/256)*k2l(ii)*k2l(jj)*(betax(ii)*betax(jj))^1.5)*(exp(1i*(3*(mux(ii)-mux(jj))))...
           +3*(exp(1i*(mux(ii)-mux(jj)))))  ;  
    H11110=termSign*((1i/64)*k2l(ii)*k2l(jj)*sqrt(betax(ii)*betax(jj)))...
           *( betay(ii)*betax(jj)*(exp(-1i*(mux(ii)-mux(jj)))-exp(1i*(mux(ii)-mux(jj))))...
             +betay(ii)*betay(jj)*(exp(1i*(mux(ii)-mux(jj)+2*(muy(ii)-muy(jj))))+exp(-1i*(mux(ii)-mux(jj)-2*(muy(ii)-muy(jj))))));
    H11200=termSign*((1i/128)*k2l(ii)*k2l(jj))*sqrt(betax(ii)*betax(jj))...
            *( betay(ii)*betax(jj)*(exp(-1i*(mux(ii)-mux(jj)-2*muy(ii)))-exp(1i*(mux(ii)-mux(jj)+2*muy(ii))))...
             +2*betay(ii)*betay(jj)*(exp(1i*(mux(ii)-mux(jj)+2*muy(ii)))+exp(-1i*(mux(ii)-mux(jj)-2*muy(ii)))));
    H40000=termSign*((1i/256)*k2l(ii)*k2l(jj)*(betax(ii)*betax(jj))^1.5)*(exp(1i*(3*mux(ii)+mux(jj))));   
    H20020=termSign*((1i/256)*k2l(ii)*k2l(jj)*sqrt(betax(ii)*betax(jj)))...
           *( betay(ii)*betax(jj)*(exp(-1i*(mux(ii)-3*mux(jj)+2*muy(ii))))...
             -betay(ii)*(betax(jj)+4*betay(jj))*(exp(1i*(mux(ii)+mux(jj)-2*muy(ii)))));     
    H20110=termSign*((1i/128)*k2l(ii)*k2l(jj)*sqrt(betax(ii)*betax(jj)))...
           *( betay(ii)*betax(jj)*(exp(-1i*(mux(ii)-3*mux(jj)))-exp(1i*(mux(ii)+mux(jj))))...
             +2*betay(ii)*betay(jj)*exp(1i*(mux(ii)+mux(jj)+2*muy(ii)-2*muy(jj))));       
    H20200=termSign*((1i/256)*k2l(ii)*k2l(jj)*sqrt(betax(ii)*betax(jj)))...
           *( betay(ii)*betax(jj)*(exp(-1i*(mux(ii)-3*mux(jj)-2*muy(ii))))...
             -betay(ii)*(betax(jj)-4*betay(jj))*(exp(1i*(mux(ii)+mux(jj)+2*muy(ii)))));    
    H00220=termSign*((1i/256)*k2l(ii)*k2l(jj)*sqrt(betax(ii)*betax(jj))*betay(ii)*betay(jj))...
           *(exp(1i*(mux(ii)-mux(jj)+2*muy(ii)-2*muy(jj)))+4*exp(1i*(mux(ii)-mux(jj)))-exp(-1i*(mux(ii)-mux(jj)-2*muy(ii)+2*muy(jj)))); 
    H00310=termSign*((1i/128)*k2l(ii)*k2l(jj)*sqrt(betax(ii)*betax(jj))*betay(ii)*betay(jj))...
           *(exp(1i*(mux(ii)-mux(jj)+2*muy(ii)))-exp(-1i*(mux(ii)-mux(jj)-2*muy(ii))));
    H00400=termSign*((1i/256)*k2l(ii)*k2l(jj)*sqrt(betax(ii)*betax(jj))*betay(ii)*betay(jj))...
           *(exp(1i*(mux(ii)-mux(jj)+2*muy(ii)+2*muy(jj))));   
    
    dnux_dJx= (-1/(4*16*pi))*k2l(jj)*k2l(ii)*(betax(ii)*betax(jj))^1.5...
              *(3*(cos(abs(mux(ii)-mux(jj))-pi*Nux))/sin(pi*Nux) +  (cos(3*abs(mux(ii)-mux(jj))-3*pi*Nux))/sin(3*pi*Nux) );          

    dnux_dJy= (1/(4*8*pi))*k2l(jj)*k2l(ii)*sqrt(betax(ii)*betax(jj))*betay(jj)...
              *(2*betax(ii)*cos(abs(mux(ii)-mux(jj))-pi*Nux)/sin(pi*Nux)...
              -(betay(ii)*cos(abs(mux(ii)-mux(jj))+2*abs(muy(ii)-muy(jj))-pi*(Nux+2*Nuy))/sin(pi*(Nux+2*Nuy)))...
              +(betay(ii)*cos(abs(mux(ii)-mux(jj))-2*abs(muy(ii)-muy(jj))-pi*(Nux-2*Nuy))/sin(pi*(Nux-2*Nuy))));
    
    dnuy_dJy= (-1/(4*16*pi))*k2l(jj)*k2l(ii)*sqrt(betax(ii)*betax(jj))*betay(jj)*betay(ii)...
              *(4*cos(abs(mux(ii)-mux(jj))-pi*Nux)/sin(pi*Nux)...
              +(cos(abs(mux(ii)-mux(jj))+2*abs(muy(ii)-muy(jj))-pi*(Nux+2*Nuy))/sin(pi*(Nux+2*Nuy)))...
              +(cos(abs(mux(ii)-mux(jj))-2*abs(muy(ii)-muy(jj))-pi*(Nux-2*Nuy))/sin(pi*(Nux-2*Nuy))));      
               
    sumH31000=sumH31000+H31000;       
    sumH22000=sumH22000+H22000;     
    sumH11110=sumH11110+H11110;  
    sumH11200=sumH11200+H11200; 
    sumH40000=sumH40000+H40000; 
    sumH20020=sumH20020+H20020;    
    sumH20110=sumH20110+H20110;
    sumH20200=sumH20200+H20200;
    sumH00220=sumH00220+H00220;
    sumH00310=sumH00310+H00310;
    sumH00400=sumH00400+H00400;
    sumdnux_dJx=sumdnux_dJx+dnux_dJx;
    sumdnux_dJy=sumdnux_dJy+dnux_dJy;
    sumdnuy_dJy=sumdnuy_dJy+dnuy_dJy;
    end
    h31000sum=h31000sum+sumH31000;
    h31000SUM(jj)=h31000sum;
    h22000sum=h22000sum+sumH22000;
    h22000SUM(jj)=h22000sum;
    h11110sum=h11110sum+sumH11110;
    h11110SUM(jj)=h11110sum;
    h11200sum=h11200sum+sumH11200;
    h11200SUM(jj)=h11200sum;   
    h40000sum=h40000sum+sumH40000;
    h40000SUM(jj)=h40000sum;  
    h20020sum=h20020sum+sumH20020;
    h20020SUM(jj)=h20020sum;    
    h20110sum=h20110sum+sumH20110;
    h20110SUM(jj)=h20110sum;     
    h20200sum=h20200sum+sumH20200;
    h20200SUM(jj)=h20200sum;     
    h00220sum=h00220sum+sumH00220;
    h00220SUM(jj)=h00220sum;     
    h00310sum=h00310sum+sumH00310;
    h00310SUM(jj)=h00310sum;  
    h00400sum=h00400sum+sumH00400;
    h00400SUM(jj)=h00400sum; 
    dnux_dJxsum=dnux_dJxsum+sumdnux_dJx;
    dnux_dJxSUM(jj)=dnux_dJxsum;
    dnux_dJysum=dnux_dJysum+sumdnux_dJy;
    dnux_dJySUM(jj)=dnux_dJysum;
    dnuy_dJysum=dnuy_dJysum+sumdnuy_dJy;
    dnuy_dJySUM(jj)=dnuy_dJysum;
end
disp('1stOrder Chromatic terms--------')
cmap=jet(5);
figure;
 h11001=h11001sum;
 fprintf('h11001:%f\n',h11001(end));
 plot(s1,h11001,'-','Color',cmap(1,:));hold on;
 h00111=h00111sum; 
 fprintf('h00111:%f\n',h00111(end));
 plot(s1,h00111,'-','Color',cmap(2,:));hold on;
 h20001=abs(h20001sum); 
 fprintf('h20001:%f\n',h20001(end));
 plot(s1,h20001,'-','Color',cmap(3,:));hold on;
 h00201=abs(h00201sum);
 fprintf('h00201:%f\n',h00201(end));
 plot(s1,h00201,'-','Color',cmap(4,:));hold on;
 h10002=abs(h10002sum); 
 fprintf('h10002:%f\n',h10002(end));
 plot(s1,h10002,'-','Color',cmap(5,:));hold on;
xlabel('S[m]');ylabel('Driving terms')
legend('h11001','h00111','h20001','h00201','h10002','Location','NorthWest');
 
disp('1st Order Geometric terms-------')
cmap=jet(5);
figure;
 h21000=abs(h21000sum);
 fprintf('h21000:%f\n',h21000(end));
 plot(s,h21000,'-','Color',cmap(1,:));hold on;
 h30000=abs(h30000sum);
 fprintf('h30000:%f\n',h30000(end));
 plot(s,h30000,'-','Color',cmap(2,:));hold on;
 h10110=abs(h10110sum);
 fprintf('h10110:%f\n',h10110(end));
 plot(s,h10110,'-','Color',cmap(3,:));hold on;
 h10020=abs(h10020sum);
 fprintf('h10020:%f\n',h10020(end));
 plot(s,h10020,'-','Color',cmap(4,:));hold on;
 h10200=abs(h10200sum);
 fprintf('h10200:%f\n',h10200(end)); 
 plot(s,h10200,'-','Color',cmap(5,:));hold on;

xlabel('S[m]');ylabel('Driving terms')
legend('h21000','h30000','h10110','h10020','h10200','Location','NorthWest');

disp('2nd Order Geometric terms-------')

cmap=jet(11);
figure;
h31000=abs(h31000SUM);
fprintf('h31000:%f\n',h31000(end));
plot(s,h31000,'-','Color',cmap(1,:));hold on;

h22000=abs(h22000SUM);
fprintf('h22000:%f\n',h22000(end));
plot(s,h22000,'-','Color',cmap(2,:));

h11110=abs(h11110SUM);
fprintf('h11110:%f\n',h11110(end));
plot(s,h11110,'-','Color',cmap(3,:));

h11200=abs(h11200SUM);
fprintf('h11200:%f\n',h11200(end));
plot(s,h11200,'-','Color',cmap(4,:));

h40000=abs(h40000SUM);
fprintf('h40000:%f\n',h40000(end));
plot(s,h40000,'-','Color',cmap(5,:));

h20020=abs(h20020SUM);
fprintf('h20020:%f\n',h20020(end));
plot(s,h20020,'-','Color',cmap(6,:));

h20110=abs(h20110SUM);
fprintf('h20110:%f\n',h20110(end));
plot(s,h20110,'-','Color',cmap(7,:));

h20200=abs(h20200SUM);
fprintf('h20200:%f\n',h20200(end));
plot(s,h20200,'-','Color',cmap(8,:));

h00220=abs(h00220SUM);
fprintf('h00220:%f\n',h00220(end));
plot(s,h00220,'-','Color',cmap(9,:));

h00310=abs(h00310SUM);
fprintf('h00310:%f\n',h00310(end));
plot(s,h00310,'-','Color',cmap(10,:));

h00400=abs(h00400SUM);
fprintf('h00400:%f\n',h00400(end));
plot(s,h00400,'-','Color',cmap(11,:));

xlabel('S[m]');ylabel('Driving terms')
legend('h31000','h22000','h11110','h11200','h40000','h20020','h20110','h20200','h00220','h00310','h00400','Location','NorthWest');

%%%higer order off-momentum terms -----------
disp('higher order off-momentum terms--------')
cmap=jet(11);
figure;
h21001=abs(h21001SUM);
fprintf('h21001:%f\n',h21001(end));
plot(s1,h21001,'-','Color',cmap(1,:)); hold on;

h30001=abs(h30001SUM);
fprintf('h30001:%f\n',h30001(end));
plot(s1,h30001,'-','Color',cmap(2,:));

h10021=abs(h10021SUM);
fprintf('h10021:%f\n',h10021(end));
plot(s1,h10021,'-','Color',cmap(3,:));

h10111=abs(h10111SUM);
fprintf('h10111:%f\n',h10111(end));
plot(s1,h10111,'-','Color',cmap(4,:));

h10201=abs(h10201SUM);
fprintf('h10201:%f\n',h10201(end));
plot(s1,h10201,'-','Color',cmap(5,:));

h11002=abs(h11002SUM);
fprintf('h11002:%f\n',h11002(end));
plot(s1,h11002,'-','Color',cmap(6,:));

h20002=abs(h20002SUM);
fprintf('h20002:%f\n',h20002(end));
plot(s1,h20002,'-','Color',cmap(7,:));

h00112=abs(h00112SUM);
fprintf('h00112:%f\n',h00112(end));
plot(s1,h00112,'-','Color',cmap(8,:));

h00202=abs(h00202SUM);
fprintf('h00202:%f\n',h00202(end));
plot(s1,h00202,'-','Color',cmap(9,:));

h10003=abs(h10003SUM);
fprintf('h10003:%f\n',h10003(end));
plot(s1,h10003,'-','Color',cmap(10,:));

h00004=abs(h00004SUM);
fprintf('h00004:%f\n',h00004(end));
plot(s1,h00004,'-','Color',cmap(11,:));

xlabel('S[m]');ylabel('Driving terms')
legend('h21001','h30001','h10021','h10111','h10201','h11002','h20002','h00112','h00202','h10003','h00004','Location','NorthWest');

fprintf('dnux/dJx:%f\n',dnux_dJxSUM(end));

fprintf('dnux/dJy:%f\n',dnux_dJySUM(end));

fprintf('dnuy/dJy:%f\n',dnuy_dJySUM(end));
figure;
plot(s,dnux_dJxSUM,'-r');hold on;
plot(s,dnux_dJySUM,'-g');
plot(s,dnuy_dJySUM,'-b');

hLegend=legend('dnux/dJx','dnux/dJy','dnuy/dJy','Location','NorthWest');
% hTitle=title('Detuning with amplitude');