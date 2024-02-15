function rpar = ExTuneScanDet(TuneScanDet,index,...
                              plotf, splitn,saveOPAf)
%Analyzes lattices generated from Tune Scan
%   
LatticeOptData=TuneScanDet{1}.inputs.LatticeOptData;
ACHRO       = LatticeOptData.ACHRO;
isdipole    = LatticeOptData.isdipole;
nvars       = LatticeOptData.nvars;
Trb         = LatticeOptData.Trb;

TuneScan    = TuneScanDet{index};
DVs         = TuneScan.outputs.DVs;

ACHRO = setDVs(2, ACHRO,LatticeOptData, DVs);

try
    rpar=atsummary_fast(ACHRO,isdipole);
catch
    fprintf('Problems at summary \n');
    rpar=NaN;
    return
end

if (strcmp(plotf,'Y'))
    ACHRO_SP=ACHRO;
    if (splitn>1)
       for i=1:splitn
           ACHRO_SP=atinsertelems(ACHRO_SP,1:length(ACHRO_SP),0.5,[]);
        end
    end
    try
        lindata=atlinopt(ACHRO_SP,0.0,1:length(ACHRO_SP)+1);
        PlotBetaDisp(lindata,'Tune scan');
    catch
        fprintf('Error in atlinopt for Tune scan lattice \n');
    end
end

if(strcmp(saveOPAf,'Y'))
   DVs=getDVs(2,ACHRO,LatticeOptData); 
   filein   = 'm4_Studies_InputfromMOGA_AT_Template_COMP.opa';
   fileout  = 'm4_Studies_InputfromTuneScan_AT_COMP.opa';
   fileIDin = fopen(filein,'r');
   pathnameout = 'C:\Users\pedtav\Documents\Accelerators\NewMachineStudies\Lattices\OrbitShift\OPA\';
   fileIDout = fopen([pathnameout fileout],'w');
   fprintf(fileIDout, '{ Decision Variables from Matched Lattice }\n');
   for i=1:nvars
       fprintf(fileIDout, 'DV%1d = %8.5f ; \n', i, DVs(i));
   end
   fprintf(fileIDout, 'RBK_MOGA = %8.5f ;\n', Trb);
    
   while 1
        tline = fgetl(fileIDin);
        if ~ischar(tline), break, end
        fprintf(fileIDout,'%s \n', tline);
    end
     fclose(fileIDin);
     fclose(fileIDout);
     fprintf('saved OPA file %s \n',fileout);
end       


