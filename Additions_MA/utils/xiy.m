function Xiy=xiy(ring, periodicity)
% get value of horizontal dispersion for  Seq(indx)
% % global THERING
% % THERING = ring;

% % AAA=atsummary;
% % Qx=mod(AAA.tunes(1) * periodicity,1);
% % %Qy=AAA.tunes(2) * periodicity;

if nargin<2
    periodicity =1;
end

[~, tunes, chrom] = twissring(ring, 0, 1:length(ring)+1, 'chrom', 1e-8);

Xiy = chrom(2) * periodicity; % 
