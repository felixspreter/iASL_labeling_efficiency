%   Main simulation file for labeling efficiency simulation
%   returns mean longitudinal magnetization of labeled blood as function of
%   labeling pulse amplitude I
%
%   Felix Spreter
%   18.11.2024


%% settings

% variable parameters
theta = 70/360*2*pi;                        % angle between catheter axis and B0 in rad
position = [0, 0]; 	                        % position of catheter axis relative to vessel axis in 0.1 mm
v = 10;                                     % flow velocity in cm/s
coil_geometry = @fcoil_singleloop;          % handle for function determining coil geometry
IlistdBm = -20:0.5:19;                        % simulated labeling pulse amplitudes in dBm
Ilist = sqrt(10.^((IlistdBm-30)./10)./50);  % pulse amplitude in A

% simulation parameter
catheter_radius = 11.5;                     % radius of catheter a coil position in 0.1 mm 
vessel_radius = 30;                         % radius of vessel at coil position in 0.1 mm
voxelsize = 2;                              % isotrop size of simulated voxels in 0.1 mm
x_min = -75;                                % volume boundaries
x_max = 75;
y_min = -50;
y_max = 50;
z_min = -50;
z_max = 50;
volume = [x_min x_max y_min y_max z_min z_max];



%% initialize

% create flow profile
x = x_min:voxelsize:x_max;
y = y_min:voxelsize:y_max;
z = z_min:voxelsize:z_max;
catheterposition = position./voxelsize;
flowprofile = zeros(length(y),length(z));

for i = (1:length(y))
    for j = (1:length(z))
        r1 = sqrt((i-length(y)/2-0.5-catheterposition(1))^2+(j-length(z)/2-0.5-catheterposition(2))^2)*voxelsize;
        r2 = sqrt((i-length(y)/2-0.5)^2+(j-length(z)/2-0.5)^2)*voxelsize;
        if r1 < vessel_radius
            if r2 > catheter_radius
                flowprofile(i,j) = v;
            end
        end
    end
end


%% Main simulation
tic

MzList = zeros(length(Ilist),1);

for step = (1:length(Ilist))
    I = Ilist(step);

    %calculate B-field of coil geometry with pulse amplitude I
    Bfield = fbiotsavart(voxelsize,coil_geometry(I),volume);
    
    % calculate flip angle map by rotating magnetization for each streamline
    flipanglemap = frotatemagnetization(flowprofile,voxelsize,Bfield,theta);

    % calculate mean z-magnetization from the flip angle map
    MzList(step) = fmeanmz(flipanglemap,voxelsize,catheterposition);
end

toc

Elist = (1-MzList)./2;
%% display results
powerlist = Ilist.^2*50;

figure
plot(powerlist*1000, Elist)

ylim([-1 1])
xlabel('P [mW]')
ylabel('E')
