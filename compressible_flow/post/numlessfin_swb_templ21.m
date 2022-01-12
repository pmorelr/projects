% MF206 Introduction to CFD 2021
% Instructor: Marica Pelanti
% Final Application Lesson - FINAL EXAM 
% ---------------------------------------------------------
% solution of the 1D shallow water equations via FV conservative schemes
% SWE with bottom topography b
% Hydrostatic Reconstruction Method 
% Uses Rusanov for homogeneous system
% TEMPLATE TO FILL
%-------------------------------------------------------------------------------
%
% qv = matrix of conserved variables   Nct \times 2
% columns = conserved variables (2)
% rows = cell values (Nct)
% qv(:,1) = density
% qv(:,2) = momentum 
%
%for final_time = [0.1 0.4 0.7 2.0]
grav = 1.; % gravity constant
%---------------------------------------------------------------------------
% grid and time step
%
Nc = 200;        % number of grid cells over the interval [xl xr]
Nct = Nc+4;      % total number of grid cells (including 4 ghost cells) 
rat = 0.8;       % fixed ratio dt/h % chosen so that CFL = max(lambda) dt/h \leq 1.
xl= 0.;          % left limit of space interval  [xl xr]      
xr = 1.;         % right limit of space interval [xl xr]
tf= 19.8 ;   %0.7 % final time %
%
% Type of BC (= 1 zero-order extrapolation, = 2 solid fixed wall)
bcl=1; % Left BC
bcr=1; % Right BC

%-------------------------------------------------------------
% Set method ('Rusanov' or other implemented)
method = 'Rusanov';
%-------------------------------------------------------------
% Finite Volume discretization 
% assume uniform grid with grid spacing = h
%
% 	
%              xl    					xr
%  X   |   X   |---X---|---X---|---X---  ...  --|---X---|   X   |   X
%  1       2       3       4       5              Nc+2    Nc+3    Nc+4=Nct
%  
%  1,2 : ghost cells on the left of xl
%  Nc+3,Nc+4 : ghost cells on the right of xr
%-------------------------------------------------------------
% mesh width
h = (xr-xl)/Nc;

% time step (here fixed)
dt = rat*h;

x  = [xl+h/2:h:xr-h/2]';        % coordinates of cell centers  
xt = [xl-3*h/2:h:xr+3*h/2]';    % coordinates of cell centers including ghost cells


MaxStep = ceil(tf/dt);          % number of time steps
%-----------------------------------------------------------------------

% Initialization

t=0.0 ;                         % initial time 

qv  = zeros(Nct,2);

bot = zeros(Nct,1);   % bottom topography

%------------------------------------------------------------------
% Set initial conditions 
% iproblb=1 : small perturbation problem; iproblb=2 : oscillating lake
iproblb=2;

%---------------------------------------------------------------------------

if (iproblb==1)
% Initial conditions [LeVeque test J. Comput. Phys., 146:346, 1998]
  hfix = 1;
  pert = 0.2;

  for k = 1:Nct		  % water height  perturbation
     qv(k,1) =hfix;     
   if ((xt(k)> 0.1)&&(xt(k) < 0.2)) % better && than & - 15/03/2020
     qv(k,1) =hfix+pert;	
   end 
  end 

  for k = 1:Nct	          % set topography
    htot = qv(k,1);
   if ((xt(k)> 0.4)&&(xt(k) < 0.6)) % better && than & - 15/03/2020
      bot(k) = 0.25*(cos(pi*(xt(k)-0.5)/0.1)+1.0);	
      qv(k,1)= htot - bot(k);
   end 
  end
  
else % iproblb=2 (oscillating lake) 
% Initial conditions [Audusse et al. test SIAM J. Sci. Comput., 25:2050, 2004]

 for k = 1:Nct	          % set topography 
      bot(k) = 0.5*(1.-0.5*(cos(pi*(xt(k)-0.5)/0.5)+1.0));	
 end

 for k = 1:Nct		  % water height
  htemp(k) = .4-bot(k)+0.04*sin((xt(k)-0.5)/.25);
  htottemp(k) = bot(k) + htemp(k);
  htotv(k) = max (bot(k), htottemp(k));
  qv(k,1) = htotv(k)-bot(k);
 end 

end
%-------------------------------------------------------------------
Ncp2 = Nc+2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for  istep=1:MaxStep, % starts time integration loop
     w = qv;  % store old qv 
   
% update solution qv
   
     
     for i=3:Ncp2

     if(w(i-1,1)<0)
     ur = rur/rr;   % velocity
     else
     ur=0;  % set zero velocity if dry state
     end

     hl1 = max(w(i-1,1) + bot(i-1) - max(bot(i-1), bot(i)), 0);
     hr1 = max(w(i,1) + bot(i) - max(bot(i-1), bot(i)), 0);

     if(w(i-1,1)>0)
     Ul1 = [hl1 hl1*w(i-1,2)/w(i-1,1)];   % velocity
     else
     Ul1 = [hl1 0];  % set zero velocity if dry state
     end

     if(w(i,1)>0)
     Ur1 = [hr1 hr1*w(i,2)/w(i,1)]   % velocity
     else
     Ur1 = [hr1 0]  % set zero velocity if dry state
     end

     hl2 = max(w(i,1) + bot(i) - max(bot(i), bot(i+1)), 0);
     hr2 = max(w(i+1,1) + bot(i+1) - max(bot(i), bot(i+1)), 0);

     if(w(i,1)>0)
     Ul2 = [hl2 hl2*w(i,2)/w(i,1)];   % velocity
     else
     Ul2 = [hl2 0];  % set zero velocity if dry state
     end

     if(w(i+1,1)>0)
     Ur2 = [hr2 hr2*w(i+1,2)/w(i+1,1)]   % velocity
     else
     Ur2 = [hr2 0]  % set zero velocity if dry state
     end

     flui = fluxswRSn_templ(Ul2, Ur2) + [0 0.5*(grav*(w(i,1)**2 - hl2**2))]
     fluip = fluxswRSn_templ(Ul1, Ur1) + [0 0.5*(grav*(w(i,1)**2 - hr1**2))]

     qv(i,:) = w(i,:) - dt/h*(flui-fluip);
          
    end 
     
%------------------------------------------------------------- ---       
%
% set Bondary conditions
% left ghost cells
     for k = 1:2      
       qv(2,k) =  qv(3,k); % zero-order extrapolation
       qv(1,k) =  qv(3,k);	
     end  
     if (bcl==2)  % solid wall     
       qv(2,1) =  qv(3,1); % symmetry
       qv(1,1) =  qv(4,1);
     
       qv(2,2) = -qv(3,2); % negate momentum
       qv(1,2) = -qv(4,2);
     end
% right ghost cells      
     for k = 1:2      
       qv(Nct-1,k) =  qv(Ncp2,k); % zero-order extrapolation
       qv(Nct,k) =  qv(Ncp2,k);	
     end         
      if (bcr==2)  % solid wall
       qv(Nct-1,1) = qv(Ncp2,1); % symmetry
       qv(Nct,1) =  qv(Nc+1,1);
      
       qv(Nct-1,2) = -qv(Ncp2,2); % negate momentum
       qv(Nct,2) = -qv(Nc+1,2);
     end       
%                
%-------------------------------------    
    t=t+dt;  % time increment 
%------------------------------------------------------------------------
end  % end time loop
%------------------------------------------------------------------------
%
%-- Plot numerical solution at final time in terms of  variables h+b,u
%
%  variables at interior cells 
rhosol = qv(3:Ncp2,1);   % height
rhousol = qv(3:Ncp2,2);  % momentum
bottom = bot(3:Ncp2,1);  % bottom

usol = zeros(size(rhosol));

htotsol = bottom +rhosol;  % water level h+b

for ii=1:Nc
 if (rhosol(ii)>0)
  usol(ii) = rhousol(ii)./rhosol(ii);  % velocity
 end
end

%  Plot results 

      figure(1) % (total) water height
      plot(x,htotsol,'bo')
      axis([xl xr  0  1.22])
      if (iproblb==2)
        axis([xl xr  0  0.5])
      end
      hold on
      plot(x,bottom,'--k')
      grid     
      title(sprintf('Height at t = %d, Nc = %d',istep*dt,Nc))
      set(gca,'FontSize',20)
      hold off
      
      figure(2) % velocity
      plot(x,usol,'bo')
%      axis([xl xr  0  1.1])
      grid     
      title(sprintf('Velocity at t = %d, Nc = %d',istep*dt,Nc))
      set(gca,'FontSize',20)
      
      figure(3) % water height - zoom
      plot(x,htotsol,'bo')
      axis([0 1 0.98 1.12])    
      if (iproblb==2)
        axis([xl xr  0  1.])
      end
      grid     
      title(sprintf('Height at t = %d, Nc = %d',istep*dt,Nc))
      set(gca,'FontSize',20)

      csvwrite('hdry.csv',[x htotsol])
      csvwrite('udry.csv',[x usol])
      csvwrite('bdry.csv',[x bottom])


      if final_time == 0.1
      csvwrite('h01.csv',[x htotsol])
      csvwrite('u01.csv',[x usol])
      end

      if final_time == 0.4
      csvwrite('h04.csv',[x htotsol])
      csvwrite('u04.csv',[x usol])
      end

      if final_time == 0.7
      csvwrite('h07.csv',[x htotsol])
      csvwrite('u07.csv',[x usol])
      end

      if final_time == 2.0
      csvwrite('h20.csv',[x htotsol])
      csvwrite('u20.csv',[x usol])
      end
      

%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% depends on functions: ....m (numerical fluxes)
%

 