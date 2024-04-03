function RDT=computeRDT_MS(ring, index, varargin)
%COMPUTERDT Computes Hamiltonian resonance driving terms (RDTs)
%   This function calls RDTElegantAT mex function and returns the
%   hamiltonian resonance driving terms, using the elegant c++ 
%   function computeDrivingTerms().
%   
%   RDT=computeRDT(ring, index, varargin)
%
%   ring is the AT lattice
%   index is the vector of indexes where one wants to compute RDTs
%   The additional arguments specify what set of driving terms should be
%   computed;
%       'chromatic'
%       'coupling'
%       'geometric1' 
%       'geometric2'
%       'tuneshifts'
%   If the linear lattice has not changed it is possible to re-use data
%   from the previous call, by using the string input
%       'reuselinear'
%   This is particularly useful if trying to minimize RDTs via optimization
%   procedures, using sextupoles and octupoles as knobs
%
%   Note that if a lattice contains an atringparam element specifying a
%   periodicity greater than 1, the RDT calculation will be done for the
%   full lattice.
%
%   example:
%   RDT=computeRDT(ring, indexBPM, 'geometric1', 'tuneshifts');
%   creates an array of structs (the length of the array is the number of 
%   indexes where you want to compute driving terms) with first order
%   geometric driving terms and tune shifts with amplitude.
%   The driving terms are complex numbers, the tune shifts are real.


% Declare persistent variables, i.e. the data that may be reused if linear
% optics have not changed.
persistent indDQSO Lin AVEBETA AVEDISP AVEMU elemL
persistent betax betay etax phix phiy a2L b2L Mux Muy Tunex Tuney nElem indDQSOe s sIndex sEnd

% Default values
reuse_linear = 0;

% INPUT HANDLING
chromatic=0;
coupling=0;
geometric1=0;
geometric2=0;
tuneshifts=0;
for ii=numel(varargin):-1:1
%     if iscell(varargin{ii})
%         if numel(varargin{ii}) == 3 && all(cellfun(@isnumeric,varargin{ii}))
%             optFieldInput = varargin{ii};
%         end
%         varargin{ii} = [];
%     end
    if ischar(varargin{ii})
        switch varargin{ii}
        case 'chromatic'
            chromatic=1;
            varargin(ii) = [];
        case 'coupling'
            coupling=1;            
            varargin(ii) = [];
        case 'geometric1'
            geometric1=1;            
            varargin(ii) = [];
        case 'geometric2'
            geometric2=1;            
            varargin(ii) = [];
        case 'tuneshifts'
            tuneshifts=1;
            varargin(ii) = [];
        case 'reuselinear'
            reuse_linear = 1;
            varargin(ii) = [];
        otherwise
            disp(['The input number ' num2str(ii+2) ' must be one of these:']);
            disp('''chromatic'', ''coupling'', ''geometric1'',''geometric2'', ''tuneshifts'',''reuselinear''');
            disp('your input will be ignored');
        end
    end
end

% If no RDT output specified, calculate all supported by the mex file
if ~any([chromatic, coupling, geometric1, geometric2, tuneshifts])
    chromatic=1;
    coupling=1;
    geometric1=1;
    geometric2=1;
    tuneshifts=1;
end

% Check periodicity of the lattice
if isfield(ring{1},'Periodicity')
    nPeriods = ring{1}.Periodicity;
else
    nPeriods = 1;
end

% If first run, or the linear optics should not be reused from previous
% function calls, run through the calculations.
if isempty(indDQSO) || ~reuse_linear 
    % Check that the last element is a L=0 IdentityPass, otherwise add it.
    % Useful to get all needed data from atavedata and avoid an additional call
    % to linopt...
    if ring{end}.Length ~= 0 || ~strcmpi(ring{end}.PassMethod,'IdentityPass')
        ring{end+1} = atmarker('End');
    end

    % Calculate all the expensive things that may be kept static during
    % sextupole/octupole optimizations, i.e. the linear machine functions
    % along with extracting the quadrupole strengths.

    % Note that the below does not take thin multipoles into account, i.e.
    % elements with ThinMPolePass are ignored currently.
    indDQSO=findcells(ring,'Class','Bend','Quadrupole','Sextupole','Octupole','Multipole');
    [Lin,AVEBETA,AVEMU,AVEDISP,~,~]=atavedata_MS(ring,0,1:length(ring));
    elemL = getcellstruct(ring,'Length',indDQSO);
    a2 = getcellstruct(ring,'PolynomA',indDQSO,1,2);
    b2 = getcellstruct(ring,'PolynomB',indDQSO,1,2);

    % Create input arguments for the mex function
    s = cat(2,Lin.SPos);
    sIndex = s(indDQSO);
    sEnd = s(end);
    betax=AVEBETA(indDQSO,1);
    betay=AVEBETA(indDQSO,2);
    etax=AVEDISP(indDQSO,1);
    phix=AVEMU(indDQSO,1);
    phiy=AVEMU(indDQSO,2);

    a2L=a2.*elemL;
    a2L(isnan(a2L))=0;
    b2L=b2.*elemL;
    b2L(isnan(b2L))=0;

    Mux=Lin(end).mu(1);    % Hor. phase adv.
    Tunex=Mux/2/pi;
    Muy=Lin(end).mu(2);    % Ver. phase adv.
    Tuney=Muy/2/pi;
    nElem=length(indDQSO);

    % If nPeriods > 1, expand the machine functions and input data
    if nPeriods > 1
        betax   = repmat(betax,nPeriods,1);
        betay   = repmat(betay,nPeriods,1);
        etax    = repmat(etax,nPeriods,1);
        tmpx    = ones(numel(indDQSO),1) * ((1:nPeriods) - 1)*Mux;
        phix    = repmat(phix,nPeriods,1) + tmpx(:);
        tmpy    = ones(numel(indDQSO),1) * ((1:nPeriods) - 1)*Muy;
        phiy    = repmat(phiy,nPeriods,1) + tmpy(:);
        Mux     = Mux * nPeriods;
        Muy     = Muy * nPeriods;
        a2L     = repmat(a2L,nPeriods,1);
        b2L     = repmat(b2L,nPeriods,1);
        Tunex   = nPeriods*Tunex;
        Tuney   = nPeriods*Tuney;

        tmpi    = ones(numel(indDQSO),1) * ((1:nPeriods) - 1)*numel(ring);
        indDQSOe = repmat(indDQSO,1,nPeriods) + tmpi(:)';
        tmps    = ones(numel(s),1) * ((1:nPeriods) - 1)*s(end);
        s       = repmat(s,1,nPeriods) + tmps(:)';
        sIndex  = s(indDQSOe);
        sEnd    = s(end);
        nElem   = length(indDQSOe);
    else
        indDQSOe = indDQSO;
    end
end


% Extract higher order fields (sextupole, octupole), which may have changed
% while linear lattice kept static
B = getcellstruct(ring,'PolynomB',indDQSO); % <-- 2nd most significant expense next to running the RDTelegantAT mex-function
% b3L = cellfun(@(x) x(3*double(3 <= numel(x))), b) .* elemL;
% b4L = cellfun(@(x) x(4*double(4 <= numel(x))), b) .* elemL;
b = zeros(size(B,1),4); for n = 1:size(B,1), b(n,1:numel(B{n})) = B{n}; end
% b = atgetfieldvalues(ring(indDQSO),'PolynomB',3:4);
% b = cat(1,b{:});
% 
% b3L = cellfun(@(x,n) x(n*(n <= numel(x))), b);

b3L=b(:,3).*elemL;
b3L(isnan(b3L))=0;
b4L=b(:,4).*elemL;
b4L(isnan(b4L))=0;

% If nPeriods > 1, expand the sextupole/octupole vectors
if nPeriods > 1
    b3L     = repmat(b3L,nPeriods,1);
    b4L     = repmat(b4L,nPeriods,1);
end



% signMat = zeros(numel(b3L));
% for m = 1:numel(b3L)
%     for n = 1:numel(b3L)
%         signMat(m,n) = sign(s(m) - s(n));
%     end
% end
% 
% h22000 = (b3L.*betax.^(3/2).*exp(1i*3*(phix-phix(1)))) * (b3L.*betax.^(3/2).*exp(1i*-3*(phix-phix(1))))' + ...
%     3*(b3L.*betax.^(3/2).*exp(1i*(phix-phix(1)))) * (b3L.*betax.^(3/2).*exp(1i*-1*(phix-phix(1))))';
% h22000 = h22000 .* signMat .* 1i ./ 64; 
% h22000 = abs(sum(sum(h22000)))
% 



for ii=1:length(index)
    % Adjust mex-function input for which element in which the RDTs are to
    % be computed
    FromindexDQSO=sum(indDQSOe<index(ii))+1;
    betax_Fromindex=[betax(FromindexDQSO:end);betax(1:FromindexDQSO-1)];
    betay_Fromindex=[betay(FromindexDQSO:end);betay(1:FromindexDQSO-1)];
    etax_Fromindex=[etax(FromindexDQSO:end);etax(1:FromindexDQSO-1)];
    phix_Fromindex=[phix(FromindexDQSO:end)-AVEMU(index(ii),1);phix(1:FromindexDQSO-1)+Mux-AVEMU(index(ii),1)];
    phiy_Fromindex=[phiy(FromindexDQSO:end)-AVEMU(index(ii),2);phiy(1:FromindexDQSO-1)+Muy-AVEMU(index(ii),2)];
    s_Fromindex=[sIndex(FromindexDQSO:end)-s(index(ii)),sIndex(1:FromindexDQSO-1)+sEnd-s(index(ii))];
    a2L_Fromindex=[a2L(FromindexDQSO:end);a2L(1:FromindexDQSO-1)];
    b2L_Fromindex=[b2L(FromindexDQSO:end);b2L(1:FromindexDQSO-1)];
    b3L_Fromindex=[b3L(FromindexDQSO:end);b3L(1:FromindexDQSO-1)];
    b4L_Fromindex=[b4L(FromindexDQSO:end);b4L(1:FromindexDQSO-1)];

    % Call the mex-function
    [ReRDT, ImRDT, TSwA]=RDTelegantAT(s_Fromindex,betax_Fromindex,betay_Fromindex,...
        etax_Fromindex,phix_Fromindex,phiy_Fromindex,a2L_Fromindex,b2L_Fromindex,...
        b3L_Fromindex,b4L_Fromindex,Tunex,Tuney,nElem,...
        chromatic,coupling,geometric1,geometric2,tuneshifts);

        
    %% Assign outputs
    RDT = struct([]);

    %chromatic
    if(chromatic)
        RDT(ii).h11001=ReRDT(6)+1i*ImRDT(6);
        RDT(ii).h00111=ReRDT(7)+1i*ImRDT(7);
        RDT(ii).h20001=ReRDT(8)+1i*ImRDT(8);
        RDT(ii).h00201=ReRDT(9)+1i*ImRDT(9);
        RDT(ii).h10002=ReRDT(10)+1i*ImRDT(10);
    end
    %coupling
    if(coupling)
        RDT(ii).h10010=ReRDT(11)+1i*ImRDT(11);
    	RDT(ii).h10100=ReRDT(12)+1i*ImRDT(12);
    end
    %geometric1
    if(geometric1)
    	RDT(ii).h21000=ReRDT(1)+1i*ImRDT(1);
    	RDT(ii).h30000=ReRDT(2)+1i*ImRDT(2);
    	RDT(ii).h10110=ReRDT(3)+1i*ImRDT(3);
    	RDT(ii).h10020=ReRDT(4)+1i*ImRDT(4);
    	RDT(ii).h10200=ReRDT(5)+1i*ImRDT(5);
    end
    %geometric2
    if(geometric2)
        RDT(ii).h22000=ReRDT(13)+1i*ImRDT(13);
    	RDT(ii).h11110=ReRDT(14)+1i*ImRDT(14);
    	RDT(ii).h00220=ReRDT(15)+1i*ImRDT(15);
    	RDT(ii).h31000=ReRDT(16)+1i*ImRDT(16);
    	RDT(ii).h40000=ReRDT(17)+1i*ImRDT(17);
    	RDT(ii).h20110=ReRDT(18)+1i*ImRDT(18);
    	RDT(ii).h11200=ReRDT(19)+1i*ImRDT(19);
    	RDT(ii).h20020=ReRDT(20)+1i*ImRDT(20);
    	RDT(ii).h20200=ReRDT(21)+1i*ImRDT(21);
    	RDT(ii).h00310=ReRDT(22)+1i*ImRDT(22);
    	RDT(ii).h00400=ReRDT(23)+1i*ImRDT(23);
    end
    %tuneshifts
    if(tuneshifts)
        RDT(ii).dnux_dJx=TSwA(1);
    	RDT(ii).dnux_dJy=TSwA(2);
    	RDT(ii).dnuy_dJy=TSwA(3);
    end

end

