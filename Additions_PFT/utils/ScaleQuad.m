function LATTICE_new = ScaleQuad(LATTICE,Kval,fam)
%ScaleQuad(LATTICE,Kval,fam) changes all quad strengths in "LATTICE" for family "fam" by the same scaling
%factor. Kval is the desired maximum absolute value of the quad strengths 
%in the family.
%   
    I_fam = find(atgetcells(LATTICE, 'FamName', fam));
    K_fam = atgetfieldvalues(LATTICE, I_fam, 'PolynomB', {1,2});
    factor = Kval/max(abs(K_fam));
    LATTICE_new = atsetfieldvalues(LATTICE, I_fam, 'PolynomB',{1,2}, K_fam*factor); 
end