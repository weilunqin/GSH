function tf = visu2isotf_ww(l,fstr,h,cap)

% VISU2ISOTF(l,fstr,h,cap) returns the isotropic spectral transfer
%      (or: eigenvalues) of several gravity related quantities. 
%      Upward continuation and cap-smoothing may be included.
%    
% IN  l    - spherical harmonic degree. 		Vector.
%     fstr - string, denoting the functional under consideration:
%               'none', 'geoid', 'dg' or 'gravity' (grav. anomaly),
%               'tr' (grav. disturbance), 'trr' (d^2/dr^2),  
%               'tew' (d^2/dlambda^2),'tns' (d^2/dtheta^2)
%               'trt' (d^2/drdtheta)
%               'potential','slope' (size of surface gradient).
%            When FSTR is left empty, a menu pops up.
%     h    - height [m] above equator. 		Scalar. Default 0.
%     cap  - radius [deg] of smoothing cap(s).          Default 0.
%            If scalar: low-pass filter. If two-vector: band-pass.
% OUT tf   - transfer. Size and shape equal to l. Units are respectively
%               [none], [m], [mGal], [mGal], [E], [m^2/s^2], [rad]
%
% NB .In case fstr='none' ISOTF represents pure upward continuation and/or
%     cap-smoothing.
%
% See also UPWCON, PELLINEN

% Nico Sneeuw                        Munich                           03/05/95
%
% uses ISSCALAR, ISINTEGER, ISVECTOR, UPWCON, PELLINEN
%      constants.m 

% rev.1  27/9/96 NS:
%    - menu-driven input to get FSTR.
%    - slope option added.
%    - band smoothing, when CAP is two-element vector. 

% march 9 2006 - Wouter
% pellinen filter changed to gaussian filter
% filter coefficients loaded, zeros added from degree 60 to 120
%
% May 7 2006
% if length(cap) == 1 CHANGED TO
% if cap > 0
% so that no gaussian filter is applied for cap = 0
%
% November, 21, 2007
% Use the routine gaussian from Matthias to compute the gaussian weights
% difference with previously compute weights is in the order of 10^-13
%
% November 27, 2007
%   tf = -(l+1) * GM/r/r * 1e5;			    % [mGal]
% changed to: 
%   tf = (l+1) * GM/r/r * 1e5;			    % [mGal]
%
% December 23, 2008
% added: transfer functions for horizontal gravity
%
% January 30, 2009
% added: transfer functions for uplift rate (Wahr et al.  2000)
%
% June 14, 2010
% added: transfer functions for east-west and north-south component of the gravity tensor
% (T_{lambda lambda} and T_{theta theta})
%
% defaults
% ----------------------------------
if nargin < 4,   cap  = 0;       end
if nargin < 3,   h    = 0;       end
if nargin < 2,   fstr = 'none';  end
if isempty(cap), cap  = 0;       end
if isempty(h),   h    = 0;       end

% general checks
% --------------------------------------------------------
if nargin < 1,         error('Gimme more'),            end
if prod(size(cap))>2,  error('CAP should be scalar or 2-vector.'), end
if ~(length(h) == 1),  error('H should be scalar.'),   end
if ~isstr(fstr),       error('FSTR must be string.'),  end
if ~any(~rem(l,1)),    error('L must be integer.'),    end
if ~(min(size(l))==1), error('L must be vector.'),     end
if min(l) < 0,         error('Negative L occurs.'),    end

% treat eigenvalue part and dimensioning
% --------------------------------------------------------
% visu2constants_champ					  % load constants
 constants_Grace % (Wouter dec. 6)
%constants_geophys;
r = ae + h;
if strcmp(fstr,'none')
   tf = ones(size(l));				        % []
elseif strcmp(fstr,'geoid')
   tf =  ones(size(l)) * r;			        % [m]
elseif strcmp(fstr,'potential')
   tf = ones(size(l)) * GM/r;			    % [m^2/s^2]
elseif strcmp(fstr,'horns')
   tf = ones(size(l)) * GM/r/r * 1e5;   	% [mGal]
elseif strcmp(fstr,'dispns')
   tf = ones(size(l));                     % [m]
elseif strcmp(fstr,'horew')
   tf = ones(size(l)) * GM/r/r * 1e5;   	% [mGal]
elseif strcmp(fstr,'dispew')
   tf = ones(size(l));                     	% [m]
elseif strcmp(fstr,'gravity') | strcmp(fstr,'dg')
   tf = (l-1) * GM/r/r * 1e5;			    % [mGal]
elseif strcmp(fstr,'tr')
   tf = (l+1) * GM/r/r * 1e5;			    % [mGal]
elseif strcmp(fstr,'dfl_nth')
   tf = -ones(size(l)) * 180 / pi * 3600;   % [arcsec]
elseif strcmp(fstr,'dfl_east')
   tf = -ones(size(l)) * 180 / pi * 3600;   % [arcsec]
elseif strcmp(fstr,'trr')
   tf = (l+1).*(l+2) * GM/r/r/r * 1e9;		% [E]
elseif strcmp(fstr,'tew')
   tf = ones(size(l)) * GM/r/r/r * 1e9;                     % [E]
elseif strcmp(fstr,'tns')
   tf = ones(size(l)) * GM/r/r/r * 1e9;                     % [E]
elseif strcmp(fstr,'trt')
   tf = -(l+1) * GM/r/r/r * 1e9;                     % [E]
elseif strcmp(fstr,'slope') % (delta T_n / delta x, delta T_n  delta y; see Rummel and van Gelderen 1995).
   tf = sqrt(l.*(l+1)) * 180 / pi * 3600;	% [arcsec]
elseif strcmp(fstr,'upl')
   tf = (2*l+1)/2 * r;                       % [m]
%    tf = (1.15*l +0.5)* r;                       % [m]
else
   error('Requested functional FSTR not available.')
end

% treat upward continuation part: when r=R+h is used in the dimensioning
% it is equal for all functionals
if h > 0, tf = visu2upwcon(l,h) .* tf; end   	% add upw. cont. to transfer

% treat cap-smoothing part
% if length(cap) == 1                             	% low-pass only
%      'low-pass filter'
if cap == 0
else
%     Gauss = gaussian_holger(length(l)-1,cap);
    Gauss = gaussian(length(l)-1,cap);
%     save d:/data/Benchmark/checkgauss_400km.dat Gauss -ascii
%     save d:/data/Benchmark/transf_before.dat tf -ascii
    tf = Gauss.* tf;
%     save d:/data/Benchmark/transf_after.dat tf -ascii
% elseif cap == 50                             	
% %     'apply gaussian filter 50 km smoothing radius'
%     load 'D:/codes/gauss_50km.dat'
%     Gauss = gauss_50km(1:size(l))';
%     tf = Gauss.* tf;
%     tf = Gauss.* tf;
% else
%     error('Filter does not exist for specified radius')
end
    
% PELLINEN FILTER
% if length(cap)==1
%     tf = visu2pellinen(l,cap) .* tf;
% else                                            % band pass filter
%     'band-pass pellinen filter'
%     tf = visu2pellinen(l,min(cap)).*(1-visu2pellinen(l,max(cap))).*tf;
end   
