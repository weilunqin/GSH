function f = visu2gshsag(field,lam,th,h,ldesired,cap,quant,jflag)

% VISU2GSHSAG calculates a global spherical harmonic synthesis for any grid
% defined by lam and phi (each vectors). The radius must be scalar. The output
% is the distrubing potential and any derivative up to the fourth.
%
% HOW:     f = visu2gshsag(field,lam,th,h,ldesired,cap,quant)
%
% Input: field [c x s]   gravity field in cs or sc format
%        lam   [n x 1]   longitude [deg]
%        th    [m x 1]   co-latitude [deg]      【90-latitude】
%        h     [1 x 1]   height [m] (optional)  【0】
%        ldesired  [1 x 1]   maximum degree (optional)     【0】
%        cap   [1 x 1]   radius of smoothing cap [degree] (default: 0)  【0】
%        quant [string]  optional argument, defining the field quantity:    【none】
%                         - 'none'      = coefficients define the output
%                         - 'potential' = potential [m^2/s^2]
%                         - 'horns'     = derivative in theta direction [m/s^2]
%                         - 'horew'     = derivative in lambda direction [m/s^2]
%                         - 'geoid'     = geoid height [m]
%                         - 'tr'        = gravity disturbance [mGal]
%                         - 'dg'        = gravity anomaly [mGal]
%                         - 'dfl_nth'   = north comp. of deflection [arcsec]
%                         - 'dfl_east'  = east comp. of deflection [arcsec]
%                         - 'trr'       = 2nd rad. derivative [E]
%                         - 'tew'       = 2nd derivative (lambda lambda) [E]
%                         - 'tns'       = 2nd derivative (theta theta) [E]
%                         - 'slope'     = slope [arcsec]
%                         - 'upl'       = uplift rate according to Wahr et
%                         al. (2000)
%                        default: 'none'
%       jflag            determines whether GRS80 is subtracted.	【0】	- def: true (1)
%
% Output: f    [n x m]   field quantity
%
% Note:    - h must be scalar!
%
% Weigelt, DoGE @ UofC                                             25.10.05

%----------------------------------------------------------------------------
% uses
% m-files: visu2normalklm.m, cs2sc.m, visu2isotf_ww.m
%----------------------------------------------------------------------------
% revision history
%
% Wouter
% ldesired, synthesis up to this degree
% constants_Grace
%
% January 30, 2009
% Add conversion to uplift rate
%
% June 14, 2010
% added: east-west component of the gravity tensor (T_{lambda lambda})
% see Novak and Grafarend (Stud Geophys Geod 2006) equation 55
%
% January 14, 2010
% added: horizontal displacement (same as horizontal gravity but with unit
% zero).
%
% March 10, 2011
% Corrected: horizontal gravity, derivative of Plm(costheta) does not give
% a sine(theta). This term is now removed.
%
%----------------------------------------------------------------------------
% remarks
%----------------------------------------------------------------------------
%

%----------------------------------------------------------------------------
% INPUT CHECK and PREPARATION
%----------------------------------------------------------------------------
error(nargchk(4,8,nargin));
if nargin < 8 || isempty(jflag), jflag = 1; end
if nargin < 7 || isempty(quant), quant = 'none'; end
if nargin < 6 || isempty(cap),   cap   = [];      end

% load necessary constants
% constants_Grace;

% Size determination for field
[row col]=size(field);
if (row~=col) && (col~=2*row-1)
    error('Input ''field'' not in cs or sc format');
elseif col~=row
    % if data is in cs-format we transfer it to sc-format
    field=visu2sc2cs(field);
    [row,col]=size(field);
end
lmax = row-1;

% check if r, lam and th are vectors
if ~any(size(lam)==1)
    error('''lam'' must be scalar or vectorial.');
elseif ~any(size(th)==1)
    error('''phi'' must be scalar or vectorial.');
elseif prod(size(h))~=1
    error('''h'' must be scalar.');
end

% prepare coordinates
th   = sort(th(:));
lam  = lam(:).*pi/180;

% use the desired degree
if nargin < 5 || isempty(ldesired), ldesired = lmax; end
if ldesired < lmax,
    field = field(1:ldesired+1,1:ldesired+1);
    lmax  = ldesired;
end

% rearrange field and substract reference field
field     = visu2cs2sc(field);
[row,col] = size(field);
if jflag, field = field - visu2cs2sc(full(visu2normalklm(lmax,'wgs84'))); end

% prepare cosine and sine --> cos(m*lam) and sin(m*lam)
m       = [0:lmax];
l       = m';
mlam    = [lam*m]';
cosmlam = cos(mlam);
sinmlam = sin(mlam);

% apply transfer function
transf  = visu2isotf_ww(l,quant,h,cap);
field   = field .* (transf(:) * ones(1,2*lmax+1));


%----------------------------------------------------------------------------
% CALCULATION
%----------------------------------------------------------------------------
for m = 0:row-1
    if m==0
        if strcmp(quant,'dfl_east') || strcmp(quant,'horew') || strcmp(quant,'dispew') % derivative to lambda cos <--> sine. this is taken care of by interchanging Clm and Slm
            Cnm = zeros(row,1);             % there are no Sn0 coefficients
            Snm = -m.*field(:,row+m);       % get column with order 0
        elseif strcmp(quant,'tew') % for a double derivative coefficients should not be interchanged.
            Cnm = -m.^2.*field(:,row+m);    % get column with order 0
            Snm = zeros(row,1);             % there are no Sn0 coefficients
        else
            Cnm = field(:,row+m);           % get column with order 0
            Snm = zeros(row,1);             % there are no Sn0 coefficients
        end
        if strcmp(quant,'dfl_nth') || strcmp(quant,'horns') || strcmp(quant,'trt') || strcmp(quant,'dispns')
            [dummy,P] = visu2plm_ww(l,m,th);   % calc derivative of fully normalized Legendre Polynoms for deflection of vertical
        elseif strcmp(quant,'tns') 
            [dummy,dummy2,P] = visu2plm_ww(l,m,th);   % calc second derivative of fully normalized Legendre Polynoms
        else
            P = visu2plm_ww(l,m,th);         % calc fully normalized Legendre Polynoms
        end
        % ------------------------------------------------------------------
        TA(:,1) = P*Cnm;
        TB(:,1) = P*Snm;
    else
        if strcmp(quant,'dfl_east') || strcmp(quant,'horew') || strcmp(quant,'dispew')
            Cnm = m.*field(:,row-m);        % get Cnm coefficients for order m
            Snm = -m.*field(:,row+m);        % get Snm coefficients for order m
        elseif strcmp(quant,'tew') % for a double derivative coefficients should not be interchanged.
            Cnm = -m.^2.*field(:,row+m);        % get Cnm coefficients for order m
            Snm = -m.^2.*field(:,row-m);        % get Snm coefficients for order m
        else
            Cnm = field(:,row+m);           % get Cnm coefficients for order m
            Snm = field(:,row-m);           % get Snm coefficients for order m
        end
        if strcmp(quant,'dfl_nth') | strcmp(quant,'horns') | strcmp(quant,'trt') || strcmp(quant,'dispns')
            [dummy,P] = visu2plm_ww(l,m,th);   % calc derivative of fully normalized Legendre Polynoms for deflection of vertical
        elseif strcmp(quant,'tns') 
            [dummy,dummy2,P] = visu2plm_ww(l,m,th);   % calc second derivative of fully normalized Legendre Polynoms
        else
            P = visu2plm_ww(l,m,th);         % calc fully normalized Legendre Polynoms
        end
        % ------------------------------------------------------------------
        TA(:,m+1) = P*Cnm;
        TB(:,m+1) = P*Snm;
    end
end

% now do the final summation
f = TA*cosmlam+TB*sinmlam;

if strcmp(quant,'horew') || strcmp(quant,'dispew')
    f = repmat(1./sind(th),1,length(lam)).*f;
elseif strcmp(quant,'tew')
    f = repmat(1./sind(th).^2,1,length(lam)).*f;
end


