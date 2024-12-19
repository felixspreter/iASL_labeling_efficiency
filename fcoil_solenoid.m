% Function creates coil geometrie of solenoid coil and the orientation of
% the current I at each position
% returns Nx6 array where the first three entries describe the position of
% each N coil elemt and the other three entries describe the current in
% this segment
function coil = fcoil_solenoid(I)

% define geometrie
radius = 10.75;
windings = 7;
height = 50;

% number of segments
N = 610;

coil = zeros([N,6]);                % (x,y,z,jx,jy,jz) for each element


% the position of each coil segment is manually defined
coil(1:30,1) = linspace(-height/2-30,-height/2,30);
coil(1:30,2) = sin(pi/4)*(radius+0.75);
coil(1:30,3) = cos(pi/4)*(radius+0.75);

coil(31:530,2) = sin(linspace(pi/4, windings*2*pi, 500))*radius;
coil(31:530,3) = cos(linspace(pi/4, windings*2*pi, 500))*radius;
coil(31:530,1) = linspace(-height/2,height/2,500);

coil(531:610,1) = linspace(height/2,-height/2-30,80);
coil(531:610,2) = 0;
coil(531:610,3) = radius+0.75;

% the current in each segment is defined from the vector between subsequent
% segments
for i = (1:N-1)
    dx = (coil(i+1,1)-coil(i,1));
    dy = (coil(i+1,2)-coil(i,2));
    dz = (coil(i+1,3)-coil(i,3));
    coil(i,4) = dx*I;
    coil(i,5) = dy*I;
    coil(i,6) = dz*I;
end
end