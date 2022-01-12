function ff =fluxswRoe_templ(v,vp) 
%
% function for Final Application Lesson MF206 
% solution of 1D Shallow Water Equations via Rusanov's method (or Roe's)
% TEMPLATE
% Numerical flux F_{i+1/2} = F(u_i,u_{i+1}) 
% v = u_i, vp = u_{i+1}
%
% initialization
grav = 1.; % gravity constant
%-------------------------------------------------------------------
%initialization
flul = zeros(1,2);   % f(u_i)      % SWE flux evaluated at u_i
flur = zeros(1,2);   % f(u_{i+1})  % SWE flux evaluated at u_{i+1}

ff = zeros(1,2);     % numerical flux (output)

%-----------------------------------------------------------------------
% set f(u_i) 
rl = v(1);   % water height
rul = v(2);  % momentum

if(rl>0)
 ul = rul/rl;   % velocity
else
 ul=0;  % set zero velocity if dry state
end


flul(1) = rul
flul(2) = (rul^2)/rl + grav * (rl^2)/2

%------------------------------------------------------------------------
% set f(u_{i+1}) 
rr = vp(1);   % water height
rur = vp(2);  % momentum

if(rr>0)
 ur = rur/rr;   % velocity
else
 ur=0;  % set zero velocity if dry state
end
    
flur(1) = rur
flur(2) = (rur^2)/rr + grav * (rr^2)/2


if ((rl ==0)&&(rr==0))  % flux is zero if both left and right states are dry
return                 
end

%-----------------------------------------------------------------
% for Roe flux need to define Roe averages,
% e.g. hroe = 0.5*(rl+rr); % water height,

% Initialization of eigenvalues, eigenvectors and alpha coefficients : 
lambdav = zeros(2); % Roe eigenvalues
Rmatv = zeros(2, 2); % Roe eigenvectors 
alphav = zeros(2);  % Roe alpha coefficients, alpha = R^{-1}\Deltaq

delta = vp - v

% Roe averages :
hroe = 0.5*(rl+rr);
uroe = (sqrt(rl)*rul/rl + sqrt(rr)*rur/rr)/(sqrt(rr) + sqrt(rl));

% Roe eigenvalues :
lambdav(1)= uroe - sqrt(grav * hroe);
lambdav(2)= uroe + sqrt(grav * hroe);

% Roe eigenvectors :
Rmatv(1,1) = 1.;
Rmatv(1,2) = uroe - grav * sqrt(hroe);

Rmatv(2,1) = 1.;
Rmatv(2,2) = uroe + grav * sqrt(hroe);  

% Roe alpha coefficients :
alphav(1)= 1/(2*sqrt( grav * hroe)) * ((uroe + sqrt(grav * hroe)) * delta(1) - delta(2));
alphav(2)= 1/(2*sqrt( grav * hroe)) * (-(uroe - sqrt(grav * hroe)) * delta(1) + delta(2));

%------------------------------------------------------------------
% Rusanov (or Roe) numerical flux

ff = (0.5*(flul+flur) - 0.5 * (abs(lambdav(1)) * alphav(1) * Rmatv(1,:)  +  abs(lambdav(2)) * alphav(2) * Rmatv(2,:)));


 