function RING = calculateGirderMaps(RING)
% CALCULATEGIRDERMAPS adds maps for girder misalignments to the lattice
%  NEWRING = calculateGirderMaps(RING)
%
% NOTES
% 1. Relies on girder markers being present in the lattice. See
%    getMagGroupsFromGirderIndex for details.
% 2. The added linear maps describe the relation between shifts in the
%    girder Cartesian coordinate system to shifts in the element
%    Frenet-Serret coordinate system. The primary use is to, based on
%    girder error seeds and assuming the girder is a stiff body, accurately
%    calculate the T1, T2, R1 and R2 misalignment fields in AT2.
% 3. Once the lattice ideal orbit and girder boundaries are defined, these
%    maps are static. It is recommended to run this function once
%    immediately after lattice definition.
% 4. The girder coordinate system choice can be easily altered by replacing
%    the 'centerGirder' nested function with another custom function. It is
%    expected this may be needed based on alignment team preferences.
%    The default choice girder coordinate system takes x as a straight line
%    from the entry point to the exit point, z is the upward direction, and
%    y is the cross-product of z and x (giving a vector pointing outward
%    from the ring centre. The origo and hence the rotation point is taken
%    as the center point on the x-axis, i.e. halfway between the entrance
%    and exit points.
% 6. Currently assumes a flat accelerator, i.e. the function cannot cope
%    with vertical bends.
%
% See also getMagGroupsFromGirderIndex

% Get the positions of each element in a centered Cartesian system. This is
% done once for the entire ring, as the call is somewhat expensive.
POSDATA_RING=atgeometry3(RING,1:numel(RING)+1);

% Go through the lattice girder by girder to compute and store the
% transformation maps for entry and exit to each girder element
girderElemIndices = getMagGroupsFromGirderIndex(RING);
for girderNbr = 1:numel(girderElemIndices)
    calculatemap(girderElemIndices{girderNbr});
end


    function calculatemap(fullGirderIndex)

        % Clear out zero-length elements (required for establishing the basis)
        I = cellfun(@(x) x.Length ~= 0, RING(fullGirderIndex));
        girderIndexThick = fullGirderIndex(I);
        girderIndexThin = fullGirderIndex(~I);

        % Add the next thick element to the girderIndex, since it will be required to
        % compute the last exit map. Use mod to cover for the case when girderIndex
        % covers the last element of the lattice.

        girderIndexThick(end+1) = getnextthick(girderIndexThick(end));

        function k = getnextthick(x)
            k = mod(x,numel(RING))+1;
            L = 0;
            while L == 0
                L = RING{k}.Length;
                if L == 0, k = mod(k,numel(RING))+1; end
            end
        end


        % Get the positions of each element in a centered Cartesian system
        POSDATA=POSDATA_RING(girderIndexThick);

        % Define the girder origo and axes. This is a somewhat arbitrary choice
        % based on convenience and can be altered if required.
        centerGirder_origoGirderMidpoint;



        %% CALCULATE MAP TO LOCAL FRENET-SERRET
        % NB! Only valid for flat accelerators, i.e. no vertical bends.

        % Get any bending angles
        %         deflection = cellfun(@getBendingAngle,RING(girderIndexThick));
        %         deflectionA = atgetfieldvalues(RING(girderIndexThick),1:numel(girderIndexThick),'BendingAngle','Default',0);

        deflection = cellfun(@getDeflectionAngle,RING(girderIndexThick));

        function y = getDeflectionAngle(x)
            if isfield(x,'BendingAngle')
                y = x.BendingAngle; 
            else 
                y = 0; 
            end
        end

        for n = 1:numel(girderIndexThick)
            % Establish the local basis describing the Frenet-Serret (ex,ey,es) system in
            % the local cartesian system (as defined in @centerGirder). This is done by...
            %
            % ...first identifying the s-vector.
            if n < numel(girderIndexThick)
                es = [  POSDATA(n+1).x - POSDATA(n).x; ...
                    POSDATA(n+1).y - POSDATA(n).y; ...
                    0];
                % ...compensate for any bending angle (reverse the rotation around the z-axis)
                es = [cos(-deflection(n)/2) sin(-deflection(n)/2) 0; -sin(-deflection(n)/2) cos(-deflection(n)/2) 0; 0 0 1]*es;
            else
                es = [  POSDATA(n).x - POSDATA(n-1).x; ...
                    POSDATA(n).y - POSDATA(n-1).y; ...
                    0];
                % ...compensate for any bending angle (reverse the rotation around the z-axis)
                es = [cos(deflection(n-1)/2) sin(deflection(n-1)/2) 0; -sin(deflection(n-1)/2) cos(deflection(n-1)/2) 0; 0 0 1]*es;
            end

            % ...and normalize
            es = es./norm(es);

            % ... then identify ey as Ez
            ey = [0; 0; 1];

            % ... and take ex as the cross product ey x es
            ex = cross(ey,es);

            % Also produce the translation vector from the Cartesian origo to the
            % lattice element entrance
            et = [POSDATA(n).x; POSDATA(n).y; POSDATA(n).z];

            % Assemble the linear affine maps for element entry and exit
            % and store them in the lattice elements for future use; it's
            % static once the girder coordinate system has been defined.
            % NB! Note that thin elements should get the same map as any
            % succeding element.
            M = [ex, ey, es, et; 0 0 0 1];

            if n <= numel(girderIndexThick) - 1
                % Apply the calculated map as the entry map for the current
                % element
                RING{girderIndexThick(n)}.misalignmentMapGirderEntry = M;
            end
            if n > 1
                % Apply the calculated map as the exit map for the previous
                % element
                RING{girderIndexThick(n-1)}.misalignmentMapGirderExit = M;
            end
        end

        % Do a second pass for the full girder to assign maps to all thin elements
        for n = 1:numel(girderIndexThin)
            thisElem = girderIndexThin(n);
            prevElem = mod(thisElem-2,numel(RING))+1;
            nextElem = mod(thisElem,numel(RING))+1;

            if isfield(RING{prevElem},'misalignmentMapGirderExit')
                RING{thisElem}.misalignmentMapGirderEntry = RING{prevElem}.misalignmentMapGirderExit;
                RING{thisElem}.misalignmentMapGirderExit  = RING{prevElem}.misalignmentMapGirderExit;
            elseif isfield(RING{nextElem},'misalignmentMapGirderEntry')
                RING{thisElem}.misalignmentMapGirderEntry  = RING{nextElem}.misalignmentMapGirderEntry;
                RING{thisElem}.misalignmentMapGirderExit  = RING{nextElem}.misalignmentMapGirderEntry;
            end
        end


        function centerGirder_origoGirderMidpoint
            % The girder coordinate system is defined based on the entry
            % and exit markers. A straight line is drawn between them and is
            % used as the x-axis. z-axis is up, y-axis orthonormal to the x and
            % z-axes. The origo is taken as the half-way point.
            xm = (POSDATA(end).x + POSDATA(1).x)/2;
            ym = (POSDATA(end).y + POSDATA(1).y)/2;
            zm = (POSDATA(end).z + POSDATA(1).z)/2;

            xa = (POSDATA(end).x - POSDATA(1).x);
            ya = (POSDATA(end).y - POSDATA(1).y);
            za = (POSDATA(end).z - POSDATA(1).z);

            gx = [xa; ya; za]/norm([xa; ya; za]);
            gz = [0; 0; 1];
            gy = cross(gz,gx);

            for p = 1:numel(POSDATA)
                POSDATA(p).x = POSDATA(p).x - xm; POSDATA(p).y = POSDATA(p).y - ym; POSDATA(p).z = POSDATA(p).z - zm;
                RES = pinv([gx, gy, gz]) * [POSDATA(p).x; POSDATA(p).y; POSDATA(p).z];
                POSDATA(p).x = RES(1); POSDATA(p).y = RES(2); POSDATA(p).z = RES(3);
            end
        end


        function centerGirder_origoTrajectoryMidpoint
            % The girder coordinate system is defined based on the entry
            % and exit markers. A straight line is drawn between them and is
            % used as the x-axis. z-axis is up, y-axis orthonormal to the x and
            % z-axes. The origo is located on the reference particle
            % trajectory midpoint.
            xm = (POSDATA(end).x + POSDATA(1).x)/2;
            ym = (POSDATA(end).y + POSDATA(1).y)/2;
            zm = (POSDATA(end).z + POSDATA(1).z)/2;

            xa = (POSDATA(end).x - POSDATA(1).x);
            ya = (POSDATA(end).y - POSDATA(1).y);
            za = (POSDATA(end).z - POSDATA(1).z);

            gx = [xa; ya; za]/norm([xa; ya; za]);
            gz = [0; 0; 1];
            gy = cross(gz,gx);

            for p = 1:numel(POSDATA)
                POSDATA(p).x = POSDATA(p).x - xm; POSDATA(p).y = POSDATA(p).y - ym; POSDATA(p).z = POSDATA(p).z - zm;
                RES = pinv([gx, gy, gz]) * [POSDATA(p).x; POSDATA(p).y; POSDATA(p).z];
                POSDATA(p).x = RES(1); POSDATA(p).y = RES(2); POSDATA(p).z = RES(3);
            end
        end
        
    end
end