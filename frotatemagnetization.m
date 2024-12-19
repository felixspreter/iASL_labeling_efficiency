% Calculates flipangle along streamline, using stepwise rotation around
% arbitrary axis using rodrigues rotation algorithm resulting in
% flipanglemap, modulo 180°
% input: flow profile
% voxelsize(intager,determines resoluntion of simulation in 0.1mm)
% Bfield(4-D double, matrix of Blief component in total simulation volume,
% (x-position,y-position,z-postion, Bfield-component)
% alpha(angle between B0 and main coil/catheter axis in degree)
function flipanglemap = frotatemagnetization(flowprofile,voxelsize,Bfield,theta)

gamma = 42.577*10^6;     % MHz/T

flipanglemap = zeros([size(Bfield,2), size(Bfield,3)]);

for k = (1:size(Bfield,2))
    for l = (1:size(Bfield,3))
        if flowprofile(k,l) > 0
            v = flowprofile(k,l);
            dt = voxelsize*0.0001/(v*0.01);         % time to travese one voxel

            % determine all B1s along one streamline
            Bxline = Bfield(:,k,l,1);
            Byline = Bfield(:,k,l,2);
            Bzline = Bfield(:,k,l,3);
            B1line = sqrt((Bxline.*cos(theta)).^2+Byline.^2+(Bzline.*sin(theta)).^2);
            
            % calculate flip angle in each voxel along one streamline
            phi = -gamma*dt*abs(B1line)*2*pi;
            
            % determine rotation axis in each voxel along one streamline
            kxline = Bxline.*cos(theta)./B1line;
            kxline(isnan(kxline)) = 0;
            kyline = Byline./B1line;
            kyline(isnan(kyline)) = 0;
            kzline = Bzline.*sin(theta)./B1line;
            kzline(isnan(kzline)) = 0;

            % define rotation matrices
            R = zeros([length(Bxline),3,3]);
            R(:,1,1) = cos(phi) + kxline.^2.*(1-cos(phi));
            R(:,1,2) = kxline.*kyline.*(1-cos(phi))-kzline.*sin(phi);
            R(:,1,3) = kyline.*sin(phi)+kxline.*kzline.*(1-cos(phi));
            R(:,2,1) = kzline.*sin(phi)+kxline.*kyline.*(1-cos(phi));
            R(:,2,2) = cos(phi) + kyline.^2.*(1-cos(phi));
            R(:,2,3) = -kxline.*sin(phi)+kyline.*kzline.*(1-cos(phi));
            R(:,3,1) = -kyline.*sin(phi)+kxline.*kzline.*(1-cos(phi));
            R(:,3,2) = kxline.*sin(phi)+kyline.*kzline.*(1-cos(phi));
            R(:,3,3) = cos(phi)+kzline.^2.*(1-cos(phi));
            
            % apply all rotations to the inital magnetization
            M = [sin(theta),0,cos(theta)];
            for m = (1:length(Byline))
                M = M*squeeze(R(m,:,:));
            end

            % calculate flip angle from the final orientation of M
            flipanglemap(k,l) = acos(M(1)*sin(theta)+M(3)*cos(theta))*360/(2*pi);

        end
    end
end
end