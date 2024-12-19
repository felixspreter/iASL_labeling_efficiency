% Calculates mean flip angle of magnetization flowing through a pipe
% cross section downstream of coil, by weighing flipanglemap with flow
% profile and caqlculating mean
% input: voxelsize(intager,determines resoluntion of simulation in 0.1mm)
% flipanglemap(2-D double, flipanglemap in ° in crosssection of pipe)function meanfa = fmeanflipangleflowcompensated(flow,voxelsize,flipanglemap)
% return: meanfa(floating point value, mean flip angle through crosssection)
function meanmz = fmeanmz(flipanglemap,voxelsize,cathposition)

cradius = 11.5;
pradius = 30;
totalMz = 0;
n = 0;

% loop over all pixels in the flip angle map
for i = (1:size(flipanglemap,1))
    for j = (1:size(flipanglemap,2))
        % calculate distance to simulation/vessel center
        r1 = sqrt((i-size(flipanglemap,1)/2-0.5)^2+(j-size(flipanglemap,2)/2-0.5)^2)*voxelsize;
        % calculate distance to catheter center
        r2 = sqrt((i-size(flipanglemap,1)/2-0.5-cathposition(1))^2+(j-size(flipanglemap,2)/2-0.5-cathposition(2))^2)*voxelsize;

        % if pixel is inside vessel and outside catheter consider it in the
        % mean
        if r1>cradius
            if r2 < pradius
                if ~isnan(flipanglemap(i,j))
                    % add the cosine of the flip angle to the total Mz
                    totalMz = totalMz + cos(flipanglemap(i,j)/360*2*pi);
                    % count the number of relevant pixels
                    n = n + 1;
                end
            end
        end
    end
end
meanmz = totalMz/n;
end