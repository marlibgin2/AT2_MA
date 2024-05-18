function RING = applyGirderError(RING, gI, surge, sway, heave, yaw, pitch, roll, varargin)
% APPLYGIRDERERROR applies element misalignments based on a rigid body girder model
%
%  NEWRING = applyGirderError(RING, girderIndex, surge, sway, heave, yaw, pitch, roll, varargin)
%
% NOTES
% 1. The function relies on pre-computed linear maps going from the
%    Cartesian girder system to the local Frenet-Serret coordinate system
%    at the entrance and exit of each element. These maps are static as
%    long as the reference particle trajectory for the lattice doesn't
%    change, which is usually not the case in most common scenarios. To
%    generate them, call calculateGirderMaps. 
% 2. Affine matrix transformations are used as the linear maps. These can
%    handle both translations, rotations, and scaling. The latter is for
%    obvious reasons disabled by enforcing the matrix norm to be 1.
% 3. Rotations are applied according to the Tate convention, i.e. yaw
%    first, then pitch, and roll last.
% 
% See also calculateGirderMaps, applyErrorModel


% This step should be possible to avoid with a slight rewrite below.
% RING = RING;

% Generate the girder affine transformation matrix
R = eye(4);
R(1:3,1:3) = genGirderRotMat3(yaw,pitch,roll);
T = [surge, sway, heave]';
RT = R; RT(1:3,4) = T;

% Apply the transformation to the lattice girder elements
applyRTmat(gI);


    function applyRTmat(index)

        for n = 1:numel(index)
            % Assemble the affine matrix for the coordinate errors
            % NB! Energy deviation and time lag are ignored, as the former
            % isn't affected by misalignment and the latter is irrelevant
            % for non-cavities.
            a1 = eye(4); p1 = eye(4);
            a2 = eye(4); p2 = eye(4);

            % Note that the 'exit fields' are inverted when assembling the
            % affine matrix, before applying the transformation.
            if isfield(RING{index(n)},'R1')
                a1(1:2,1:2) =     RING{index(n)}.R1([1 3],[1 3]);
                a2(1:2,1:2) = inv(RING{index(n)}.R2([1 3],[1 3]));
                p1(1:2,1:2) =     RING{index(n)}.R1([2 4],[2 4]);
                p2(1:2,1:2) = inv(RING{index(n)}.R2([2 4],[2 4]));
            else
                % If the R1 and R2 fields didn't exist, create them in the
                % output lattice
                RING{index(n)}.R1 = eye(6);
                RING{index(n)}.R2 = eye(6);
            end
            if isfield(RING{index(n)},'T1')
                a1(1:2,4) = +RING{index(n)}.T1([1 3]);
                a2(1:2,4) = -RING{index(n)}.T2([1 3]);
                p1(1:2,4) = +RING{index(n)}.T1([2 4]);
                p2(1:2,4) = -RING{index(n)}.T2([2 4]);
            else
                % If the T1 and T2 fields didn't exist, create them in the
                % output lattice
                RING{index(n)}.T1 = zeros(6,1);
                RING{index(n)}.T2 = zeros(6,1);
            end

            % Apply the girder error for the spatial coordinates
            % RT is inverted as it represents the girder shift and a1, p1
            % etc. represent the change in particle coordinates.
            A1 = RING{index(n)}.misalignmentMapGirderEntry \ inv(RT) * RING{index(n)}.misalignmentMapGirderEntry * a1;
            A2 = RING{index(n)}.misalignmentMapGirderExit \ inv(RT) * RING{index(n)}.misalignmentMapGirderExit * a2;

            % Apply the girder error for the momentum coordinates; here
            % the translations should not be applied so only use R..
            P1 = affine2rot(RING{index(n)}.misalignmentMapGirderEntry) \ inv(R) * affine2rot(RING{index(n)}.misalignmentMapGirderEntry) * p1;
            P2 = affine2rot(RING{index(n)}.misalignmentMapGirderExit) \ inv(R) * affine2rot(RING{index(n)}.misalignmentMapGirderExit) * p2;

            % Build the new T1, T2, R1, R2
            % NB! cT (time lag) may be treated as just a position. dP
            % (particle energy deviation) is not affected by any girder
            % rotation or shift however (assuming only magnetic elements);
            % hence the difference when assembling the T1/R1/R2/T2.
            RING{index(n)}.T1([1 3 6])           = A1(1:3,4);
            RING{index(n)}.T1([2 4])             = P1(1:2,4);

            RING{index(n)}.T2([1 3 6])           = -A2(1:3,4);
            RING{index(n)}.T2([2 4])             = -P2(1:2,4);

            RING{index(n)}.R1([1 3 6],[1 3 6])   = A1(1:3,1:3);
            RING{index(n)}.R1([2 4],[2 4])       = P1(1:2,1:2);

            RING{index(n)}.R2([1 3 6],[1 3 6])   = inv(A2(1:3,1:3));
            RING{index(n)}.R2([2 4],[2 4])       = inv(P2(1:2,1:2));

            % To be safe, enforce the R1, R2 norm to be identically 1, in
            % order to avoid spurious growth or damping.
            RING{index(n)}.R1 = RING{index(n)}.R1 ./ norm(RING{index(n)}.R1);
            RING{index(n)}.R2 = RING{index(n)}.R2 ./ norm(RING{index(n)}.R2);
        end
    end


    
    function Rout = affine2rot(Rin)
    % AFFINE2ROT strip out the translation part of the affine transformation
        Rout = Rin;
        Rout(1:3,4) = zeros(3,1);
    end


    function R = genGirderRotMat3(yaw, pitch, roll)
        % GENGIRDERROTMAT generates a 3D rotation matrix
        % Order of angle application follows Tate convention, i.e. the sequence is yaw,
        % pitch and roll. This is done in the girder system, i.e. Cartesian x,y,z.
        %
        % NB! The rotation matrix assumes x axis is along the girder, z axis is
        % down-up, and y is side-to-side of the girder.
        Ry = rotz(yaw*180/pi);
        Rp = roty(pitch*180/pi);
        Rr = rotx(roll*180/pi);

        R = Rr*Rp*Ry;
    end



    function R = genGirderRotMat6(yaw, pitch, roll)
        % GENGIRDERROTMAT generates a 6D rotation matrix
        % Angle order of application follows Tate convention, i.e. sequence is yaw,
        % pitch and roll. This is done in the girder system, i.e. Cartesian x,y,z.
        %
        % NB! The rotation matrix assumes x axis is along the girder, z axis is
        % down-up, and y is side-to-side of the girder.

        Ry = eye(6);
        Ry(1:2:end,1:2:end) = rotz(yaw);
        Ry(2:2:end,2:2:end) = rotz(yaw);

        Rp = eye(6);
        Rp(1:2:end,1:2:end) = roty(pitch);
        Rp(2:2:end,2:2:end) = roty(pitch);

        Rr = eye(6);
        Rr(1:2:end,1:2:end) = rotx(roll);
        Rr(2:2:end,2:2:end) = rotx(roll);


        R = Rr*Rp*Ry;

    end

end
