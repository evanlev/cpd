%% Complementary Poisson-Disc Sampling
%
% (c) Evan G. Levine (egl@stanford.edu) 2014.
%
% <http://med.stanford.edu/bmrgroup.html BMR page>.
%
% <https://github.com/evanlev/cpd.git Download> the demos and sampling functions
%%

%% References
% 1) EG Levine, M Saranathan, and B Hargreaves. “Complementary Poisson-Disc Sampling” Proceedings of ISMRM, Milan, Italy, May 2014
% 
% 2) EG Levine, B Quist, B Daniel, B Hargreaves, and M Saranathan. “View-sharing and Compressed Sensing in Two-Point Dixon-based DCE-MRI” Proceedings of ISMRM, Milan, Italy, May 2014
%
% 3) EG Levine, B Hargreaves, B Daniel, S Vasanawala, and M Saranathan. "3D Cartesian MRI with Compressed Sensing and Variable View Sharing Using Complementary Poisson-disc Sampling"  Magnic Resonance in Medicine 2016
%%

%% Setup
% For Examples 1-3, uniform density sampling patterns will contain a fully-sampled 
% central region and annular region with 2x1 regular under-sampling.
%%
clear, close all;
cd('./src');
% makeMex;
cd('..');
addpath('utils/');
addpath('src/');

FOVRatio = 0.4; % FOVz / FOVy
nt = 4;         % # temporal phases
ny = 180;       % y-dimension
nz = 60;        % z-dimension
% create a fully-sampled central region
Areg = zpad(removeCorners(ones(24,24)), [ny nz]); 
% Annular region with 2x1 regular under-sampling
F = removeCorners(ones(ny, nz) - Areg);
F(1:2:end,:) = 0;

%% Example 1: Variable density CPD sampling by region-wise UD-CPD
%%
Ry = sqrt(ny*nz / sum(sum(F/nt)));
Rz = Ry;
R = [Ry Rz];             % Ry = Rz
shapeOpt = 'cross';      % Union of a line along k = 0 and an ellipse at t = 0
distRelaxationOpt = 'k'; % Only relax k-space min. distance constraint
vd_exp = 1;              % 1/kr^1 density
Rmax = 22;               % reduction factor at kr = kmax
verbose = 0;
C = 1;

tic;
M1 = genVDCPD(nt,Rmax,vd_exp,FOVRatio,F,shapeOpt,verbose,C);
toc;

figure(1);
imshow( [reshape(M1, [ny nz*nt]), sum(M1,3)/nt]);
title('Uniform density CPD sampling patterns from "genVDCPD.m" and sum');


%% Example 2: Variable density random sampling by region-wise UD random
%%
Ry = sqrt(ny*nz / sum(sum(F/nt)));
Rz = Ry;
R = [Ry Rz];             % Ry = Rz
shapeOpt = 'cross';      % Union of a line along k = 0 and an ellipse at t = 0
distRelaxationOpt = 'k'; % Only relax k-space min. distance constraint
vd_exp = 1;              % 1/kr^1 density
Rmax = 22;               % reduction factor at kr = kmax
tic;
M1r = genVDCPD(nt,Rmax,vd_exp,FOVRatio,F,shapeOpt,0,1,0);
toc;

figure(2);
imshow( [reshape(M1r, [ny nz*nt]), sum(M1r,3)/nt]);
title('Uniform density random sampling patterns from "genVDCPD.m" and sum');

%% Example 3: Efficient uniform density CPD implementation
% Uniform density sampling is efficiently implemented with a
% random queue and requires little tuning.
%%
tic;
C = 0.5;
M2 = genUDCPD(nt,FOVRatio,F,shapeOpt) + repmat(Areg, [1 1 nt]);
toc;

figure(3);
imshow( [reshape(M2, [ny nz*nt]), sum(M2,3)/nt]);
title('Uniform density CPD sampling patterns from "genUDCPD.m" and sum');

%% Example 4: Use of min. distance criterion for temporal phase combination
% Combining multiple phases should preserve the uniformity and incoherence of 
% the Poisson-disc sample distribution. Here different min. distance criteria are
% compared. The 'cross' parameter generates an ellipse in k-space around any sample
% drawn and a line through the ellipse along the time axis, fixed at the sample k-space.
% Subsequent samples are drawn from the region of k-t space outside of the ellipse and
% line. Parameters 'ellipsoid', 'l1 ball', and 'cones' correspond to other shapes
% that better preserve the sampling properties after combining temporal phases, as
% shown in this example.
%%
figure(4);
nt = 20; 
F = removeCorners(ones(100,100));
M3= genUDCPD(nt,FOVRatio,F,'cross'); %+ repmat(Areg, [1 1 nt]);
M6= genUDCPD(nt,FOVRatio,F,'cones'); % + repmat(Areg, [1 1 nt]);
nVS = 4; % View-share this many phases
imshow( [M3(:,:,1), sum(M3(:,:,1:2),3), sum(M3(:,:,1:4),3); M6(:,:,1) sum(M6(:,:,1:2),3), sum(M6(:,:,1:4),3)]);
title({'Logical OR of 1,2,and 4 sequential patterns with ', '''cross'' (top),''cones'' (bottom) shapes (left to right)'});


