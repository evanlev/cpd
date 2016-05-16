% mask = removeCorners(maskin)
% 
% maskin = input mask
%
% mask   = input mask with corners cut
function mask = removeCorners(maskin)
    [ny nz] = size(maskin);
    % Elliptical sampling for Random and Uniform
    [x, y] = meshgrid(1:ny, 1:nz);
    f = sqrt((x-ny/2).^2/(ny/2)^2 + (y-nz/2).^2/(nz/2)^2)  < 1;
    mask = maskin.*(f');

end
