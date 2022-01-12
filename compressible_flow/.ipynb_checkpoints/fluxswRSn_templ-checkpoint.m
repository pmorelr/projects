function ff =fluxswRSn_templ(v,vp)

% function for Rusanov's flux
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
flul(2) = ul * rul + grav * (rl^2)/2

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
flur(2) = ur * rur + grav * (rr^2)/2


if ((rl ==0)&&(rr==0))  % flux is zero if both left and right states are dry
return                 
end

% Rusanov numerical flux

% amax = Local maximum absolute value of characteristic speeds
amax = max( max(abs(ur + sqrt(grav*rr)), abs(ur - sqrt(grav*rr))),  max(abs(ul + sqrt(grav*rl)), abs(ul - sqrt(grav*rl))) )

% Finally, we compute :
ff = 0.5*((flul +flur) - amax*(vp-v));
 