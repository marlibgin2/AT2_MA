function emitx = emitx(ACHRO,isdipole)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    rpara = atsummary_fast(ACHRO,isdipole);
    emitx = rpara.naturalEmittance*1E12;
end