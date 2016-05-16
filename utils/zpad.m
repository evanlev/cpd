% yz = zpad(y, dimsNew);
% 
% zero padding 
%
% y       = data (kx,ky,kz,coils or ky,kz,coils)
% dimsNew = new zero pad dimensions
function yNew = zpad(y, dimsNew)
dims = size(y);
if length(dims) ~= length(dimsNew)
    dims
    dimsNew
    error('dims (above) and dimsNew (second above) are not equal');
end

yNew = zeros(dimsNew(:)');

for k = 1:length(dims)
    middle{k} = ((-dims(k)/2 + 1):dims(k)/2) + dimsNew(k)/2; 
end
switch( length(dimsNew) )
    case 2
        yNew(middle{1}, middle{2}) = y;
    case 3
        yNew(middle{1}, middle{2}, middle{3}) = y;
    case 4
        yNew(middle{1}, middle{2}, middle{3}, middle{4}) = y;
    case 5
        yNew(middle{1}, middle{2}, middle{3}, middle{4}, middle{5}) = y;
    case 6
        yNew(middle{1}, middle{2}, middle{3}, middle{4}, middle{5}, middle{6}) = y;
    otherwise
        length(dims)
        error('length(dims) (above) is wrong');
end

end
