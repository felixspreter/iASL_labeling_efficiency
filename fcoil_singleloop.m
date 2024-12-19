% Function creates coil geometrie of single loop saddle coil and the orientation of
% the current I at each position
% returns Nx6 array where the first three entries describe the position of
% each N coil elemt and the other three entries describe the current in
% this segment

function coil = fcoil_singleloop(I)

% define geometrie
radius = 10.75;
length = 50;
width = 21.5;

% number of segments
N = 170;

coil = zeros([N,6]);                % (x,y,z,jx,jy,jz) for each element

% the position of each coil segment is manually defined
coil(1:50,1) = linspace(-length/2,length/2,50);
coil(1:50,2) = -width/2;
coil(1:50,3) = sqrt(radius^2-(width/2)^2);

coil(51:85,1) = length/2;
coil(51:85,2) = cos(linspace(pi,0,35))*radius;
coil(51:85,3) = sin(linspace(pi,0,35))*radius;

coil(86:135,1) = linspace(length/2,-length/2,50);
coil(86:135,2) = width/2;
coil(86:135,3) = sqrt(radius^2-(width/2)^2);

coil(136:170,1) = -length/2;
coil(136:170,2) = cos(linspace(0,pi,35))*radius;
coil(136:170,3) = sin(linspace(0,pi,35))*radius;

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