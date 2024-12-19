% Calculates the B1+field crated by the coil using the Biot_Savart Law
% input: voxelsize(intager value, resolution of simulation in 0.1mm)
% coil(Nx6 double, (x,y,z,Ix,Iy,Iz)) 
% volume( boundaries of simulation field, (x_min x_max y_min y_max z_min
% z_max)
%return: (4-D double (x-position,y-position,z-position,B-field component)
function Bfield = fbiotsavart(voxelsize,coil,volume)

mu_0 = 1.2564*10^-6;                % in N/A^2
diameter = 1.5;                     % wire diameter in 0.1 mm
x_min = volume(1);   
x_max = volume(2); 
y_min = volume(3); 
y_max = volume(4); 
z_min = volume(5); 
z_max = volume(6); 

x = x_min:voxelsize:x_max;
y = y_min:voxelsize:y_max;
z = z_min:voxelsize:z_max;

% create array of position of each voxel
positions = zeros([length(x),length(y),length(z),3]);
for i = (1:length(x))
    for j = (1:length(y))
        for k = (1:length(z))
            positions(i,j,k,:) = [x(i),y(j),z(k)];
        end
    end
end

Bfield = zeros([length(x),length(y),length(z),3]);
A = zeros(size(positions));                             % distances r-r'
Wire = ones([length(x),length(y),length(z)]);           % coil wire
for i = (1:size(coil,1)) 
    A(:,:,:,1) = positions(:,:,:,1)-coil(i,1);
    A(:,:,:,2) = positions(:,:,:,2)-coil(i,2);
    A(:,:,:,3) = positions(:,:,:,3)-coil(i,3);

    Aabs = sqrt(A(:,:,:,1).^2+A(:,:,:,2).^2+A(:,:,:,3).^2);
    Bfield(:,:,:,1) = Bfield(:,:,:,1) + mu_0/(4*pi)*(coil(i,5).*A(:,:,:,3)-coil(i,6).*A(:,:,:,2)).*Aabs(:,:,:).^(-3).*10^4;
    Bfield(:,:,:,2) = Bfield(:,:,:,2) + mu_0/(4*pi)*(coil(i,6).*A(:,:,:,1)-coil(i,4).*A(:,:,:,3)).*Aabs(:,:,:).^(-3).*10^4;
    Bfield(:,:,:,3) = Bfield(:,:,:,3) + mu_0/(4*pi)*(coil(i,4).*A(:,:,:,2)-coil(i,5).*A(:,:,:,1)).*Aabs(:,:,:).^(-3).*10^4;
    Wire = Wire.*Aabs>diameter;  
end

Bfield(:,:,:,1) = Bfield(:,:,:,1).*Wire;
Bfield(:,:,:,2) = Bfield(:,:,:,2).*Wire;
Bfield(:,:,:,3) = Bfield(:,:,:,3).*Wire;

Bfield = Bfield./2;                     	            %only considering B1+ field
end