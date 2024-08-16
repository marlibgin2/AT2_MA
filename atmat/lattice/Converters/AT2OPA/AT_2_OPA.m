function AT_2_OPA(AT_ring,linename)
%% function AT_2_OPA(AT_ring,linename)
% Converts the AT2.0 lattice AT_ring in OPA format.
% 
% Elements in the input structure are expected to have an "OPAType" field
% 
% file ['' linename '_lattice.opa'] is generated containing the lattice
% element definitions and the LINE. no other commands introduced
% 
%
% OPA may be found here: http://people.web.psi.ch/streun/opa/
% OPA Types are
% 'Dipole'
% 'Quadrupole'
% 'Sextupole'
% 'Octupole'
% 'Multipole'

%% History
% PFT 2024/06/05 :updated to handle mutipoles in OPA. Added OPAType
% PFT 2024/06/08 :updated to read in the periodicity and export to OPA

outfile=['' linename '_lattice.opa'];

%outfile='madXelemdef.elem';
%elelat=['{com   madX lattice elements: ' linename ' com}\n{com   Created: ' datestr(now) ' com}\n'];

%% get family names for definitions
%[families,ind_first_oc_ring]=...
%    unique(getcellstruct(AT_ring,'FamName',1:length(AT_ring)),'first');

[~,nbper]=atenergy(AT_ring);
for i=1:numel(AT_ring)
    if (isfield(AT_ring{i},'OPAFam'))
        AT_ring{i}.FamName=AT_ring{i}.OPAFam;
    end
end

[families,ind_first_oc_ring]=...
    unique(getcellstruct(AT_ring,'FamName',1:length(AT_ring)),'first');

%elelat=[elelat '{com DEFINITIONS com}\n\n'];

elelat=['Energy = ' num2str(AT_ring{1}.Energy*1e-9) ';\n\n'];

format='%8.10f';

%% loop families for definitions
for i=1:length(families)
   el= AT_ring{ind_first_oc_ring(i)};
   if isfield(el,'BetaCode')
       type=el.BetaCode;
   elseif isfield(el,'OPAType')
       type=el.OPAType;
   elseif isfield(el,'Class')
       type=el.Class;
   else
       type='Marker';
   end
      
    switch type
        case {'DI','Dipole','Bend'} % dipole
            di=[' ' el.('FamName')   ': '...
                ' Combined, L = ' num2str(el.('Length'),format) ', '...
                'T = ' num2str(180/pi*el.('BendingAngle'),format) ', '...
                'K = ' num2str(el.('PolynomB')(2),format) ', '...
                'T1 = ' num2str(180/pi*el.('EntranceAngle'),format) ', '...
                'T2 = ' num2str(180/pi*el.('ExitAngle'),format) ', '...
                'K1in = ' num2str(0,format) ', '...
                'k1ex = ' num2str(0,format) ', '...
                'k2in = ' num2str(0,format) ', '...
                'k2ex = ' num2str(0,format) ', '...
                'gap = ' num2str(0,format) ', '...
                'ax = ' num2str(32,format) ', '...
                'ay = ' num2str(20,format) '; '...
               ];
           
            elelat=[elelat di '\n\n']; %#ok<*AGROW>
        case {'QP','Quadrupole'} % quadrupole
            
            if el.('MaxOrder')==3
            if el.('PolynomB')(4)==0
            qp=[' ' el.('FamName') ': '...
                ' quadrupole,  L = ' num2str(el.('Length'),format)  ', '...
                ...'K = ' num2str(el.('PolynomB')(2)/el.('Length'),format) ', '...
                'K = ' num2str(el.('PolynomB')(2),format) ', '...
                 'ax = ' num2str(32,format) ', '...
                'ay = ' num2str(20,format) '; '...
               ];
            else % if octupole in quadrupole, split quadrupole as in mad8 q oc q q oc q
              
              qp=[' ' el.('FamName') '_sl: '...
                ' quadrupole,  L = ' num2str(el.('Length')/4,format)  ', '...
                'K = ' num2str(el.('PolynomB')(2),format) ', '...
                 'ax = ' num2str(32,format) ', '...
                'ay = ' num2str(20,format) '; '...
               '\n'...
               ]; 
              qp=[qp ' ' el.('FamName')   '_oc:Multipole,'...
                ' N= ' num2str(4,format)  ','...
                ...' K= ' num2str(el.('PolynomB')(4)*6*el.('Length')/2,format)  ','...
                ' K= ' num2str(el.('PolynomB')(4)*el.('Length')/2,format)  ','...
                'ax = ' num2str(32,format) ', '...
                'ay = ' num2str(20,format) '; '...
               '\n'...
               ];
              qp=[qp el.('FamName') ':'...
                   el.('FamName') '_sl,'...
                   el.('FamName') '_oc,'...
                   el.('FamName') '_sl,'...
                   el.('FamName') '_sl,'...
                   el.('FamName') '_oc,'...
                   el.('FamName') '_sl,'...
               '\n'...
               ];
              qp(end-2:end)=[];
              qp=[qp ';\n'];
            end
            else
             qp=[' ' el.('FamName') ': '...
                ' quadrupole,  L = ' num2str(el.('Length'),format)  ', '...
                 'K = ' num2str(el.('PolynomB')(2),format) ', '...
                 'ax = ' num2str(32,format) ', '...
                'ay = ' num2str(20,format) '; '...
               ];
               
            end
            
            
            elelat=[elelat qp '\n\n'];
        case {'SX','Sextupole'} % sextupole
            sx=[' ' el.('FamName')   ': '...
                ' sextupole,  L = ' num2str(el.('Length'),format)  ', '...
                ...' K = ' num2str(el.('PolynomB')(3)*2,format) ', n = 4 , '...
                ' K = ' num2str(el.('PolynomB')(3),format) ', n = 6 , '...
                'ax = ' num2str(32,format) ', '...
                'ay = ' num2str(20,format) '; '...
               ];
            elelat=[elelat sx '\n\n'];
        case {'OC','Octupole'} % octupole
            sx=[];
            nslice=2;
            sx=[sx ' ' el.('FamName')   '_dr: drift, L= ' num2str(el.('Length')/2/nslice,format)  ';' '\n'];
            sx=[sx ' ' el.('FamName')   '_sl:Multipole,'...
                ' N= ' num2str(el.('MaxOrder')+1,format)  ','...
                ...' K= ' num2str(el.('PolynomB')(el.('MaxOrder')+1)*factorial(el.('MaxOrder')-1)*el.('Length')/nslice,format)  ','...
                ' K= ' num2str(el.('PolynomB')(el.('MaxOrder')+1)*el.('Length')/nslice,format)  ','...
                'ax = ' num2str(32,format) ', '...
                'ay = ' num2str(20,format) '; '...
               '\n'...
               ];
            sx=[sx el.('FamName') ':'];
            for islice=1:nslice
              sx=[sx el.('FamName') '_dr, ' el.('FamName') '_sl, ' el.('FamName') '_dr, '];
            end
            sx(end-1:end)=[];
            sx=[sx ';\n'];
            
            elelat=[elelat sx '\n\n'];
        case {'MP','Multipole'} % multipole
           mp=[' ' el.('FamName')   ':Multipole,'...
                ' N= ' num2str(el.('MaxOrder')+1,format)  ','...
                ...' K= ' num2str(el.('PolynomB')(el.('MaxOrder')+1)*factorial(el.('MaxOrder')),format)  ','...
                ' K= ' num2str(el.('PolynomB')(el.('MaxOrder')+1),format)  ','...
                ' ax = ' num2str(32,format) ', '...
                'ay = ' num2str(20,format) '; '...
               ];
            elelat=[elelat mp '\n\n'];
        case {'ThinMultipole'} % multipole
        warning('still to be defined!')
        case {'PU','Monitor'} % bpm
            pu=[' ' el.('FamName') ': monitor' ', '...
                  'ax = ' num2str(32,format) ', '...
                'ay = ' num2str(20,format) '; '...
              ];
            elelat=[elelat pu '\n\n'];
        case {'MRK','Marker'}% marker
            mrk=[' ' el.('FamName') ': marker' '; '...
              ];
            elelat=[elelat mrk '\n\n'];
        case {'KI','Corrector'} % kicker
            ki=[' ' el.('FamName') ': marker' '; '...
               ];
            
             elelat=[elelat ki '\n\n'];
        case {'SD','DR','Drift'} % drift
            
            dr=[' ' el.('FamName') ': drift, L = ' num2str(el.('Length'),format) ', '...
                'ax = ' num2str(32,format) ', '...
                'ay = ' num2str(20,format) '; '...
               ];
           
            elelat=[elelat dr '\n\n']; 	
            
         case {'CAV','RFCavity'} % drift
            %rfc=[' ' el.('FamName') ' : RFCavity, L=' num2str(el.('Length'),format)...
                %',VOLT='  num2str(el.('Voltage'),format) ' ;'...
                %', freq=' num2str(el.('Frequency'),format) '' ' ;'];
             rfc=[' ' el.('FamName') ':marker;'];
             elelat=[elelat rfc '\n\n'];
           
        otherwise
            warning(['Element: ' el.('FamName') ' was not converted, since it does not match any Class.'])
            mrk=[' ' el.('FamName') ': marker' '; '...
                ];
            elelat=[elelat mrk '\n\n'];
    end
  
end


%elelat=[elelat '\n\n {com ----- table of segments -------------- com}\n\n'];

elelat=[elelat linename ':'];

%% define lattice line
% loop all elements
for i=1:length(AT_ring)
    if i~=1
        elelat=[elelat ',' AT_ring{i}.('FamName') '\n'];
    else
        elelat=[elelat '' AT_ring{i}.('FamName') '\n'];
    end
end

elelat=[elelat ', nper = ' num2str(nbper) '\n'];

elelat=[elelat ';'];

elelat=strrep(elelat,'RFC','RFCav');
%% print to file

of=fopen(outfile,'w');
fprintf(of,elelat);

fclose('all');




return