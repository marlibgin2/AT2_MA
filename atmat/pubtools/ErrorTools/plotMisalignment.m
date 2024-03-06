function varargout=plotMisalignment(varargin) 
%PLOTMISALIGNMENT Plots element misalignments in the lattice
%
%Helper function for atplot: plot the misalignments
%
%  USAGE:
% >> atbaseplot(ring,@plotAperture);
% >> atplot(ring,@plotAperture);        (obsolete)
%
%See also atplot atbaseplot

if nargout == 1 % From atplot
    ring=varargin{2};
%     index=findcells(ring,'T1');
    index = cellfun(@(x) isfield(x,'T1') && any(isfield(x,{'PolynomA','PolynomB','K','BendingAngle'})),ring);
    indexB = circshift(index,1);

    % Establish output vectors
    Xe=[nan(size(ring)); nan];
    Ye=Xe;
    Pe=Xe;

    % Retrieve coord. transformation at the entrance and exit
    xe   = [-getcellstruct(ring,'T1',index,1), getcellstruct(ring,'T2',index,1) ];
    ye   = [-getcellstruct(ring,'T1',index,3), getcellstruct(ring,'T2',index,3) ];
    phie = [-asin(getcellstruct(ring,'R1',index,1,3)), asin(getcellstruct(ring,'R2',index,1,3)) ];
    

    Xe(indexB)=xe(:,2);
    Ye(indexB)=ye(:,2);
    Pe(indexB)=phie(:,2);

    Xe(index)=xe(:,1);
    Ye(index)=ye(:,1);
    Pe(index)=phie(:,1);



    plotdata(1).values=[Xe Ye; nan, nan]*1e6;
    plotdata(1).labels={'\Deltax','\Deltay'};
    plotdata(1).axislabel='\Deltax,y [µm]';
    
    plotdata(2).values=[Pe; nan]*1e3;%
    plotdata(2).labels={'\Delta\phi'};
    plotdata(2).axislabel='Roll \phi [mrad]';
    varargout={plotdata};
else                % From atbaseplot
    s=findspos(varargin{1},1:length(varargin{1})+1);
    varargout={s,plotMisalignment([],varargin{:})};
end

end 
