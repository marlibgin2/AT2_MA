function [LAT_tune,Ks, rpara]=ExGl(ScanResults,i, plotf, splitn, fitf, tunes)
%Examines result of Glass Scan
% matches disersion and achromat tunes 
% A global variable ACHRORBSSIM containing the lattice (simplified version
% only one slice per dipole) must exist before the function is called
%
% inputs:
%           ScanResults: structure containig glass results (see glass_scan
%           i : linde in "Values" in ScanResults structure from whihc data will be retrived
%           plotf : if "Y", plots lattice functions for the scaned lattice
%                   as well for the descendent version with dspersion and tne fits
%           splitn : breaks the elements 2^splitn (if splitn is larger than
%           fitf : if "Y" fits lattice to match dispersion and fit tunes
%           tunes : 1x2 matrix of derired achromat fractional tunes 
%
global ACHRORBSSIM

scan_fams=ScanResults.scan_fams;
nvars = size(scan_fams,1);
Vals =ScanResults.Values;
Ks=Vals(i,1:nvars)';
props=atCheckRingProperties(ACHRORBSSIM);

LAT=ACHRORBSSIM;
if (splitn>1)
    for i=1:splitn
        LAT=atinsertelems(LAT,[1:length(LAT)],0.5,[]);
    end
end

for j=1:nvars
    I_fams{j} =  find(atgetcells(LAT,'FamName',scan_fams{j}));
    LAT   =  atsetfieldvalues(LAT,I_fams{j},'K', Vals(i,j));
    LAT   =  atsetfieldvalues(LAT,I_fams{j},'PolynomB',{1,2}, Vals(i,j));
end

%twissdata=findm44_fast(LAT,length(LAT),props);
%TRx = twissdata(1,1)+twissdata(2,2);

%TRy = twissdata(3,3)+twissdata(4,4);

%evalin('caller','global LAT');
%tic;
%rpara=atsummary(LAT);
%toc;
%tic;
%rpara=ringpara(LAT);
%toc;
if (strcmp(plotf,'Y'))
    lindata=atlinopt(LAT,0.0,1:length(LAT)+1);
    PlotBetaDisp(lindata,'GLASS Scan extracted');
end

if (strcmp(fitf,'Y'))
    Variab1 = atVariableBuilder(LAT,{'reversebend_sim','reversebendmc1_sim',...
             'reversebendmc2_sim','ssdipm','sshdip'},...
             {{'PolynomB',{1,2}},{'PolynomB',{1,2}},...
              {'PolynomB',{1,2}},{'PolynomB',{1,2}},...
              {'PolynomB',{1,2}}},...
              {0,0.0,0.0,-1.5,-1.5},{5.0,5.0,5.0,-0.5,-0.5});
         
    Constr1  = atlinconstraint(1,...
                          {{'Dispersion',{1}},...
                           {'Dispersion',{2}},...
                           {'beta',{1}},{'beta',{2}}},...
                           [0 0 0 0], [0 0 20 15], [1 1 1 1]);

    [LAT_md, penalty_md, dmin_md] = atmatch(LAT,Variab1,Constr1,...
                         1E-8,1000,1,@fminsearch); % @fmincon
    
    fprintf('Dispersion matched with penalty = %6.2e \n', sum(penalty_md.^2));
  

    for j=1:nvars
        K_fams    =  atgetfieldvalues(LAT_md,I_fams{j},'PolynomB',{1,2});
        LAT_md    =  atsetfieldvalues(LAT_md,I_fams{j},'K', K_fams);
        Ks(j)   = K_fams(1);
    end
    
    LAT_tune = atfittune(LAT_md, tunes, 'qfend_sim', 'qdend_sim');
    LAT_tune = atfittune(LAT_tune, tunes, 'qfend_sim', 'qdend_sim');
    
    if (splitn>1)
        for i=1:splitn
            LAT_tune=atinsertelems(LAT_tune,[1:length(LAT_tune)],0.5,[]);
        end
    end

    rpara=ringpara(LAT_tune);
    if (strcmp(plotf,'Y'))
        lindata=atlinopt(LAT_tune,0.0,1:length(LAT_tune)+1);
        PlotBetaDisp(lindata,'Glass Scan Extracted');
    end

end
end

