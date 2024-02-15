function RDT_Penalty = RDT_Penalty(RDTs)
%
% Calculates a penalty function from given RDTs
% Action amplitudes, dispersion and Weights are fixed parameters defined in
% the function itself
%
Js = [30 30]*1E-3/2;
DP = 0.03;
Ws = 2*ones(1,23);
WR = 0;

fields=fieldnames(RDTs);
Indices=NaN(5);
RDT_Penalty=0.0;
  for i=1:21
     for s=1:5
      Indices(s)=str2double(fields{i+2}(s+1));
     end
     j=Indices(1);
     k=Indices(2);
     l=Indices(3);
     m=Indices(4);
     p=Indices(5);
     rdt = abs(RDTs.(fields{i+2}));
     if ((j~=k)||(l~=m))
        RDT_Penalty=RDT_Penalty+(rdt*((2*Js(1))^((j+k)/2))*...
                                         ((2*Js(2))^((l+m)/2))*...
                                         DP^p*...
                                        (2^Ws(i)-1)*10^(WR))^2;
     else
        RDT_Penalty=RDT_Penalty+(rdt*((2*Js(1))^((j+k)/2))*...
                                         ((2*Js(2))^((l+m)/2))*...
                                         DP^p*...
                                        (2^Ws(i)-1))^2;
%        fprintf('jklm = %1d %1d %1d %1d %1d \n', j,k,l,m,p);
     end
  end

end

