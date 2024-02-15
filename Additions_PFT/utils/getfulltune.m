function tun = getfulltune(ring)
%UNTITLED Full tunes with integer part
%   
   [ringdata,TD] = atlinopt4(ring,1:length(ring)+1,'get_chrom','coupled',false); 
   Qx = TD(length(TD)).mu(1)/(2*pi);
   Qy = TD(length(TD)).mu(2)/(2*pi);
   tun=[Qx Qy];
end

