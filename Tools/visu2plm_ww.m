function [p, dp, ddp] = visu2plm_ww(l,m,th)

% PLM Fully normalized associated Legendre functions for a selected order M
%
% HOW p       = plm(l,th)			- assumes M=0
%     p       = plm(l,m,th)
%     [p,dp]  = plm(l,m,th)
%     [p,ddp] = plm(l,m,th)
%
%
% IN  l  - degree (vector). Integer, but not necessarily monotonic.
%          For l < m a vector of zeros will be returned.
%     m  - order (scalar). If absent, m=0 is assumed.
%     th - co-latitude [deg] (vector)
% OUT p  - Matrix with Legendre functions. The matrix has length(TH) rows
%          and length(L) columns, unless L or TH is scalar. Then the output
%          vector follows the shape of respectively L or TH. 
%    dp  - Matrix with first derivative of Legendre functions. The matrix 
%          has length(TH) rows and length(L) columns, unless L or TH is 
%          scalar. Then the output vector follows the shape of respectively 
%          L or TH. 
%    ddp - Matrix with second (colatitude) derivative of Legendre functions. The matrix 
%          has length(TH) rows and length(L) columns, unless L or TH is 
%          scalar. Then the output vector follows the shape of respectively 
%          L or TH. 
% 
% See also LEGPOL, YLM, IPLM

%-----------------------------------------------------------------------------
% Nico Sneeuw, IAPG, TU-Munich                                       08/08/94
%-----------------------------------------------------------------------------
% Uses none
%-----------------------------------------------------------------------------
% Revision history:
%  - NS09/06/97: help text brushed up
%  - NS13/07/98: Pmm non-recursive anymore 
%  - NS0299:     further help text brush-up
%  - MW13/08/04: extension for first derivative
%  - MW24/11/04: speed up calculation
%  - Wouter van der Wal, June 15, 2006: extension for second derivative
%    using recursion formulas of Novak and Grafarend (2006)
%    tested with analytical second derivatives (coefficient pairs 1,1; 2,1; 2,2; 3,1; 3,3; 4,1; 4,2; 4,4) 
%    see: 2ndDerivativeLegendre.pdf
%-----------------------------------------------------------------------------


% Some input checking.
if nargin == 2
   th = m;
   m  = 0;
end
if min(size(l)) ~= 1,  error('Degree l must be vector (or scalar)'), end
if any(rem(l,1) ~= 0), error('Vector l contains non-integers.'), end
if max(size(m)) ~= 1,  error('Order m must be scalar.'), end
if rem(m,1) ~=0,       error('Order m must be integer.'), end


% Preliminaries.
[lrow,lcol] = size(l);
[trow,tcol] = size(th);
lmax = max(l);
if lmax < m, error('Largest degree still smaller than order m.'), end
n    = length(th);				% number of latitudes
t    = th(:)*pi/180;
x    = cos(t);
y    = sin(t);
lvec = l(:)';					% l can be used now as running index.

% Recursive computation of the temporary matrix ptmp, containing the Legendre
% functions in its columns, with progressing degree l. The last column of
% ptmp will contain zeros, which is useful for assignments when l < m.
ptmp  = zeros(n,lmax-m+2);
if nargout >= 2, dptmp = zeros(n,lmax-m+2); end
if nargout == 3, ddptmp = zeros(n,lmax-m+2); end

%--------------------------------------------------------------------
% sectorial recursion: PM (non-recursive, though)
%--------------------------------------------------------------------
%WW: produces sqrt( (2n+1)/2n )
% Novak and Grafarend (2006), eq 64
if m == 0
   fac = 1;
else
   mm  = 2*(1:m);
   fac = sqrt(2*prod((mm+1)./mm));   % extra sqrt(2) because summation not over negative orders
end

ptmp(:,1) = fac*y.^m;                                      % The 1st column of ptmp.
% TEST: for m = 1, theta = 10, fac = sqrt(6/2)
% ptmp(1,1) = sqrt(6/2)*sind(10)*1
if nargout >= 2
    dptmp(:,1) = m*fac*(y.^(m-1).*x);
end     % The 1st column of dptmp.
% recursion is beta_n,n * beta_n-1,n-1 * ...
% sin(theta) * sin(theta) * ...
% * cos(theta) * P_0,0 (which is 1, dP_0,0 is 0)
% note that the term beta_n,n*cos(theta)*P_n-1,n-1 of Novak and
% Grafarend(2006) equation 72 is taken care of by the m.
if nargout == 3
    ddptmp(:,1) = -m*fac*(y.^m) + m*(m-1)*fac*(y.^(m-2).*x.^2);
end     % The 1st column of ddptmp.


%--------------------------------------------------------------------
% l-recursion: P
%--------------------------------------------------------------------
for l = m+1:lmax
   col   = l - m + 1;			% points to the next column of ptmp
   root1 = sqrt( (2*l+1)*(2*l-1)/((l-m)*(l+m)) ) ;                      % beta_n,m (65) 
   root2 = sqrt( (2*l+1)*(l+m-1)*(l-m-1) / ( (2*l-3)*(l-m)*(l+m) ) );   % beta_n,m (65) * gamma_n,m (66)

   % recursion
   if l == m+1
       ptmp(:,col) = root1 *x.*ptmp(:,col-1);
   else
       ptmp(:,col) = root1 *x.*ptmp(:,col-1) - root2 *ptmp(:,col-2);
   end
       
   if nargout >= 2, 
       if l == m+1
           dptmp(:,col) = root1 *(x.*dptmp(:,col-1)-y.*ptmp(:,col-1)); 
       else
           dptmp(:,col) = root1 *(x.*dptmp(:,col-1)-y.*ptmp(:,col-1)) - root2 *dptmp(:,col-2); 
       end
   end

   if nargout == 3,
       if l == m+1
           ddptmp(:,col) = root1 *(-x.*ptmp(:,col-1) -2*y.*dptmp(:,col-1) + x.*ddptmp(:,col-1) );
       else
           ddptmp(:,col) = root1 *(-x.*ptmp(:,col-1) -2*y.*dptmp(:,col-1) + x.*ddptmp(:,col-1) ) - root2*ddptmp(:,col-2); 
       end
   end

end


% The Legendre functions have been computed. What remains to be done, is to
% extract the proper columns from ptmp, corresponding to the vector lvec. 
% If l or theta is scalar the output matrix p reduces to a vector. It should
% have the shape of respectively theta or l in that case.

p          = zeros(n,length(lvec));         % size declaration.
lind       = find(lvec < m);			    % index into l < m
pcol       = lvec - m + 1;			        % index into columns of ptmp
pcol(lind) = (lmax-m+2)*ones(size(lind));	% Now l < m points to last col.
p          = ptmp(:,pcol);			        % proper column extraction 
if nargout >= 2, dp = dptmp(:,pcol); end    % proper column extraction
if nargout == 3, ddp = ddptmp(:,pcol); end  % proper column extraction

if max(size(lvec))==1  & min(size(th))==1 & (trow == 1), 
    p = p'; 
    if nargout >= 2, dp = dp'; end
    if nargout == 3, ddp = ddp'; end
end
if max(size(th))==1 & min(size(lvec))==1  & (lcol == 1), 
    p = p'; 
    if nargout >= 2, dp = dp'; end
    if nargout == 3, ddp = ddp'; end
end