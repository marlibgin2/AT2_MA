function rotate_floor_coordinates(theta)
Routb = maxivU_b1_1_1; Routa = maxivU_a1_1_1; 
[s2d,x2d, y2d, a2d,baa, ban] = Survey2D(Routa,0e-3);
Va = [x2d; y2d];
Va = Va - [0; -167.6291964/2]; % [0; -528/pi/2]

[s2d,x2d, y2d, a2d,baa, ban] = Survey2D(Routb,0e-3);
Vb = [x2d; y2d];
Vb = Vb - [0; -167.6291964/2]; % [0; -528/pi/2]

figure(123); hold on
plot(Va(1,:), Va(2,:),'k.')
plot(Vb(1,:), Vb(2,:),'r.')

VRa = RM(theta)*Va; VRb = RM(theta)*Vb; 

figure(321); hold on; grid on
plot(VRa(1,:), VRa(2,:),'k.')
plot(VRb(1,:), VRb(2,:),'r.')
axis([-3 3 84.08 84.18])
end


function RoMa = RM(theta)
    RoMa = [cos(theta) sin(theta); -sin(theta) cos(theta)];
end
