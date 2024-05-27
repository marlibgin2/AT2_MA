function [LAT_C,Penalty,its] = fitchroit(LAT, chrom_fams, chroms0, Nitchro, TolChrom)
% Fits the chromaticty returning a penalty function and the updated lattice
% by iterative calls to atfitchrome  
%% inputs
%      LAT: Lattice cell array
%      chrom_fams: (1X2) cell array of strings with  sextupole families to
%                  be used
%      chrosm0: (1X2) array of target chromatricities
%      Nitchro: # max of iterations
%      TolChrom: (1X2) array of tolerances
%% outputs
% LAT_C : lattice with corrected
% Penalty: penalty achieved
% its: iterations performed
%
%% Usage Examples
%[ACHR_ZC,Penalty,its] = fitchroit(ACHR, {'S3';'S5'}, [0 0], 10, [1E-6 1E-6})
% 


  try
        k=0; its=0.0;   
        LAT_C=LAT;
        [~,chr]=tunechrom(LAT);
        Penalty = (chr(1)-chroms0(1))^2 +(chr(2)-chroms0(2))^2;
        while ( ((abs(chr(1)-chroms0(1))>TolChrom(1))||...
                (abs(chr(2)-chroms0(2))>TolChrom(2))) &&(k<Nitchro))
           LAT_C=atfitchrom(LAT_C,chroms0,chrom_fams{1},chrom_fams{2});
           [~,chr]=tunechrom(LAT_C);
           Penalty = (chr(1)-chroms0(1))^2 +(chr(2)-chroms0(2))^2;
           k=k+1;
        end
        its=k;
  catch ME
        LAT_C = LAT;
        Penalty = NaN;
        rethrow(ME);
  end
