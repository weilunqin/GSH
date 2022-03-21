% this file initiates the global spherical harmonics analysis and synthesis 
% for general purporse | inserted functions are written by Wouter van der
% Val, Bart Root, and Weilun Qin from Delft University of Technology

clear all
close all
clc


%% important inputs
% 1- topography
data = load('MOLA_topo_180*360.mat'); 
topo = data.h*1e3;
root = rho_c/d_rho*topo; % calculate the root from Airy's isostasy

% 2- Effective Elastic thickness of lithosphere[km]
Te = 120e3;      % influence the effects of the flexural repsonse function
% 3- Crustal thickness with zero topography [m] 
T = 50e3;        
% 4-Max D/O in GSH transformation
N = 179;        


%% define parameters and constants for Mars
E = 1e11;       % Young's modulus 
v = 0.25;       % Poisson's ratio    
rho_m = 3500;   % mantle density [kg/m^3]
rho_c = 2900;   % crustal density [kg/m^3]
d_rho = rho_m - rho_c; % density difference [kg/m^3]
g = 3.711;      % m/s^2
R = 3389.5e3;   % radius in [m]
D = E*Te^3/(12*(1-v^2)); % Flexural rigidity



%% Run GSHA (Analysis)
tic;   
% Clm & Slm in |C\S| format [matrix]
    cs = GSHA(root,N); 
% converts spherical harmonics coefficients in |C\S| storage format into 
% a rectangular (L+1)x(2L+1) matrix in  /S|C\format.
    sc = cs2sc(cs);
% converts into SH degree,order and coefficients in matrix
    [Clm,Slm,llvec,mmvec] = sc2vecml(sc,N);
    V = [llvec' mmvec' Clm Slm];    
toc


%% load Airy Moho depth profile in Spherical Harmonics format
Moho = V;

% the flexural models 
a = size(Moho); b = a(1);

for l = 1:b
    % degree
    n = Moho(l,1); 
    % Infinite plate model
    FRS_plate(l) = 1/(1+D/(d_rho)/g*((2*n+1)/(2*R))^4);
	% Thin shell model
    AA = D/((d_rho)*g);
    BB = (n^3*(n+1)^3-4*n^2*(n+1)^2+4*n*(n+1))/(n*(n+1)-(1-v))/R^4;
    CC = 12*(1-v^2)/(R^2*Te^2);
    DD = (n*(n+1)-2)/(n*(n+1)-(1-v));
    FRS_shell(l) = 1/(1+AA*(BB+CC*DD));

    new_Moho_plate(l,1) = Moho(l,1);             new_Moho_shell(l,1) = Moho(l,1);
    new_Moho_plate(l,2) = Moho(l,2);             new_Moho_shell(l,2) = Moho(l,2);
    new_Moho_plate(l,3) = Moho(l,3).*FRS_plate(l);    new_Moho_shell(l,3) = Moho(l,3).*FRS_shell(l); 
    new_Moho_plate(l,4) = Moho(l,4).*FRS_plate(l);    new_Moho_shell(l,4) = Moho(l,4).*FRS_shell(l);
end


%% Run GSHS (Synthesis)
tic;
    lmax=N;
    Clm_plate = new_Moho_plate(:,3); Slm1 = new_Moho_plate(:,4);
    Clm_shell = new_Moho_shell(:,3); Slm2 = new_Moho_shell(:,4);
    sc_plate = vecml2sc(Clm_plate,Slm1,lmax);
    sc_shell = vecml2sc(Clm_shell,Slm2,lmax);

% VISU2GSHSAG calculates a global spherical harmonic synthesis for any grid
% defined by lam and phi (each vectors). The radius must be scalar. The output
% is the distrubing potential and any derivative up to the fourth.

    lam = [0.5:1:359.5];        % lam   [n x 1]   longitude [deg]
    th = [0.5:1:179.5];         % th    [m x 1]   co-latitude [deg]      【90-lat】
    h = 0;                      % h     [1 x 1]   height [m] (optional)  【0】
    ldesired = N;              % ldesired  [1 x 1]   maximum degree (optional)    
    cap = 0;                    % cap   [1 x 1]   radius of smoothing cap [degree] (default: 0) 
    jflag = 0;                  % quant [string]  optional argument, defining the field quantity:    'none'
    
    % call the function
    r_plate = visu2gshsag_ww(sc_plate,lam,th,h,ldesired,cap,'none',jflag);
    r_shell = visu2gshsag_ww(sc_shell,lam,th,h,ldesired,cap,'none',jflag);
    
toc


T_plate = r_plate + T*ones(size(r_plate));
T_shell = r_shell + T*ones(size(r_shell));

%% visulize the data
figure(1);clf(1);
axesm ('mollweid', 'Frame', 'on', 'Grid', 'on'); axis off;
meshm(T_plate/1e3,[1 90 180],size(T_plate));
load('roma50.mat');
colormap(roma50)
colorbar('southoutside')
ax = gca; ax.YDir = 'reverse';
c = colorbar('southoutside'); c.Label.String = 'Thickness [km] ';
caxis([0 100]);
title('Crustal Thickness (Infinite Plate - Airy)')
set(gca,'FontSize',18)
set(gcf, 'Position',  [500, 200, 900, 600])
txt1 = { ['T_{e}= ' num2str(Te/1e3) ' km']};
    text(-3,1.2,txt1,'Color','black','FontSize',11)   

output = ['Output for the infinite plate model:'];
disp(output)   
output = ['Max Tc = ' , num2str(max(max(T_plate/1e3))) 'km'];
disp(output) 
output = ['Min Tc = ' , num2str(min(min(T_plate/1e3))) 'km'];
disp(output) 
output = ['Mean Tc = ' , num2str(mean(mean(T_plate/1e3))) 'km'];
disp(output)

%%
figure(2);clf(2);
axesm ('mollweid', 'Frame', 'on', 'Grid', 'on'); axis off;
meshm(T_shell/1e3,[1 90 180],size(T_shell));
load('roma50.mat');
colormap(roma50)
colorbar('southoutside')
ax = gca; ax.YDir = 'reverse';
c = colorbar('southoutside'); c.Label.String = 'Thickness [km] ';
caxis([0 100]);
title('Crustal Thickness (Thin Shell - Airy)')
set(gca,'FontSize',18)
set(gcf, 'Position',  [500, 200, 900, 600])
txt1 = { ['T_{e}= ' num2str(Te/1e3) ' km']};
    text(-3,1.2,txt1,'Color','black','FontSize',11)   

output = ['Output for the thin shell model:'];
disp(output)
output = ['Max Tc = ' , num2str(max(max(T_shell/1e3))) 'km'];
disp(output) 
output = ['Min Tc = ' , num2str(min(min(T_shell/1e3))) 'km'];
disp(output) 
output = ['Mean Tc = ' , num2str(mean(mean(T_shell/1e3))) 'km'];
disp(output)


