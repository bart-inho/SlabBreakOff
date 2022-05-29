% Final project - Slab Break off

% Clear memory, close figures and set-up
clear
close all
warning off
format long
colormap jet
fs = 8 ;

% Define numerical model--------------------------------------------------%

    % Define initial grid size and parameters

xsize = 5e5 ; % x-model size [m] 
ysize = 4e5 ; % y-model size [m
Nx = 101 ; % x-grid resolution 
Ny = 81 ; % y-grid resolution
n = 3 ; % number of equations
N = (Nx+1)*(Ny+1)*n ; % total number of operations
NT = N/n ; % number of unknowns temperature
dx = xsize / (Nx - 1) ; % x-grid step
dy = ysize / (Ny - 1) ; % y-grid step
angle = deg2rad(45) ;

    % coordinate vectors for staggered grid

x = 0 : dx : xsize + dx ; % x-coordinate vector
y = 0 : dy : ysize + dy ; % y-coordinate vector
xvx = x ;
yvx = -dy / 2 : dy : ysize + dy / 2 ; % half step back and forward
xvy = -dx / 2 : dx : xsize + dx / 2 ; % half step back and forward
yvy = y ;
xp = -dx / 2 : dx : xsize + dx / 2 ; % half step back and forward
yp = -dy / 2 : dy : ysize + dy / 2 ; % half step back and forward

    % define material properties and matrixes

g = 10 ; % gravity acceleration

crust = [3400 %1 RHO0 [kg/m^3]
    3400 %2 RHOini
    2e-8 %3 Hr
    1e+8 %4 strength [Pa]
    273 %5 T [K] top
    1573 %6 T [K] bot
    1e+23 %7 ETAini [Pa s]
    nan %8 thermal conductivity
    1e3 %9 thermal capacity
    3e-5 %10 alpha
    1e-11 %11 beta
    ]; 
weak = [3400 %1 RHO0 [kg/m^3]
    3400 %2 RHOini
    2e-8 %3 Hr
    2e+7 %4 strength [Pa]
    273 %5 T [K] top
    1573 %6 T [K] bot
    1e+23 %7 ETAini [Pa s]
    nan %8 thermal conductivity
    1e3 %9 thermal capacity
    3e-5 %10 alpha
    1e-11 %11 beta
    ]; 
mantl = [3300 %1 RHO0 [kg/m^3]
    3250 %2 RHOini
    3e-8 %3 Hr
    5e+7 %4 strength [Pa]
    1573 %5 T [K] top
    1573 %6 T [K] bot
    1e+20 %7 ETAini [Pa s]
    nan %8 thermal conductivity
    1e3 %9 thermal capacity
    3e-5 %10 alpha
    1e-11 %11 beta
    ]; 
air = [1 %1 RHO0 [kg/m^3]
    1 %2 RHOini
    0 %3 Hr
    0 %4 strength [Pa]
    273 %5 T [K] top
    273 %6 T [K] bot
    1e+18 %7 ETAini [Pa s]
    3000 %8 thermal conductivity
    3.3e6 %9 thermal capacity
    0 %10 alpha
    0 %11 beta
    ]; 

ETAB = ones(Ny + 1, Nx + 1)*mantl(7) ; % viscosity matrix for basic nodal points
ETAP = ones(Ny + 1, Nx + 1)*mantl(7) ; % viscosity matrix for pressure points
RHOvy = ones(Ny + 1, Nx + 1)* mantl(1) ; % density matrix 
ALPHA = zeros(Ny+1, Nx+1) ;
HR = zeros(Ny+1, Nx+1) ;
HA = zeros(Ny+1, Nx+1) ;

ADdry = 2.5e-17 ;
ndry = 3.5 ;
Vadry = 8e-6 ;
Eadry = 532e3 ;

ADwet = 2e-21 ;
nwet = 4 ;
Vawet = 4e-6 ;
Eawet = 471e3 ;

ETAmin=1e+18 ; 
ETAmax=1e+24 ;

% sum martices -----------------------------------------------------------%

T0 = zeros(Ny+1, Nx+1) ;
RHOvy_SUM = zeros(Ny+1, Nx+1) ;
RHOCP_SUM = zeros(Ny+1, Nx+1) ;
ETAP_SUM  = zeros(Ny+1, Nx+1) ;
ETAB_SUM  = zeros(Ny+1, Nx+1) ;
WTVX_SUM  = zeros(Ny+1, Nx+1) ;
WTVY_SUM  = zeros(Ny+1, Nx+1) ;
WTB_SUM   = zeros(Ny+1, Nx+1) ;
WTP_SUM   = zeros(Ny+1, Nx+1) ;
WTT_SUM   = zeros(Ny+1, Nx+1) ;
Ky_SUM    = zeros(Ny+1, Nx+1) ;
Kx_SUM    = zeros(Ny+1, Nx+1) ;
T_SUM     = zeros(Ny+1, Nx+1) ;
HR_SUM    = zeros(Ny+1, Nx+1) ;
AL_SUM    = zeros(Ny+1, Nx+1) ;
DTSUM     = zeros(Ny+1, Nx+1) ;
RHOCPSUMsub = zeros(Ny+1, Nx+1) ;

% Lagragian - define marker arrays ---------------------------------------%

Nxm  = (Nx-1)*4 ; % number of markers in horizontal direction
Nym  = (Ny-1)*4 ; % number of markers in vertical direction
Nm   = Nxm * Nym     ; % total number or markers
xm   = zeros(Nm, 1)  ; % x - coordinate of markers 
ym   = zeros(Nm, 1)  ; % y - coordinate of markers
vxm  = zeros(1, 4) ; % x - velocity marker vector
vym  = zeros(1, 4) ; % y - velocity marker vector
Tm   = zeros(Nm, 1) ; % temperature [K]
ETAm = zeros(Nm, 1) ; % viscosity [Pa s]
RHOm = zeros(Nm, 1) ; % denstity [kg/m^3]
Km   = zeros(Nm, 1) ; % conductivity [W/(m K)]
CPm  = zeros(Nm, 1) ; % thermal capacity [J/K]
DTm  = zeros(Nm, 1) ; % temperature difference [K]
HRm  = zeros(Nm, 1) ; % radiogenic heat [W/m^3]
ALm  = zeros(Nm, 1) ; % alpha
T0m  = zeros(Nm, 1) ;
DTm0 = zeros(Nm, 1) ;
bigDTm = zeros(Nm, 1) ;
Pm   = zeros(Nm, 1) ;
EXXm = zeros(Nm, 1) ;
EYYm = zeros(Nm, 1) ;
EXYm = zeros(Nm, 1) ;
DTrm = zeros(Nm, 1) ;
type = zeros(Nm, 1) ;

% we could introduce marker types RHOm = [3200 3300 3350] -> class ? You 
% limit the memory used and visualise all the rock types ! -> final project

% Lagragian material conditions - positions of markers + RHO & ETA values %

dxm = xsize/Nxm ; % average horizontal distance [m] btw markers
dym = ysize/Nym ; % average vertical distance [m] btw markers

% initialize vx, vy and p matrix -----------------------------------------%

p     = zeros(Ny+1, Nx+1) ;
vx    = zeros(Ny+1, Nx+1) ;
vy    = zeros(Ny+1, Nx+1) ;
ky    = zeros(Ny+1, Nx+1) ;
kx    = zeros(Ny+1, Nx+1) ;
RHOCP = zeros(Ny+1, Nx+1) ;
DTsub = zeros(Ny+1, Nx+1) ;

SXX   = zeros(Ny+1, Nx+1) ;
SXY   = zeros(Ny+1, Nx+1) ;
EXX   = zeros(Ny+1, Nx+1) ;
EXY   = zeros(Ny+1, Nx+1) ;
EYY   = zeros(Ny+1, Nx+1) ;

SEAV = zeros(Ny+1, Nx+1) ;

% initialize matrixes solver matrices ------------------------------------%

LT = sparse(NT, NT) ; % thermal solver left hand side
RT = zeros(NT, 1) ; % thermal solver right hand side
L  = sparse(N, N) ; % mechanical solver left hand side
R  = zeros(N, 1) ; % mechanical solver right hand side

% Boundary condition type (free slip / no slip) --------------------------%

bctype = -1 ; % free slip : -1, no slip : 1

% dt and max displacement ------------------------------------------------%

dt = 6e+10 ;
% dtthermal=min(dx,dy)^2/(4*...
%     max([k_cru/RHO_cru/cp_cru k_in/RHO_in/cp_in  k_ou/RHO_ou/cp_ou]));
dismax = 0.5 ; % number of grid steps per time step
DTmax = 50 ; % max temperature difference / time step [K]
timesum = 0 ;

m = 1 ;

x1 = (500 - (200 + 70))*1e3 ;
x2 = (500 - (200 + 70 + 60))*1e3 ;
b = 50e3;
slope = (50-100)/(220-230) ;

for im = 1:1:Nym % design slab breakoff model 
    for jm = 1:1:Nxm

        % coordinates
        xm(m) = dxm/2 + (jm-1)*dxm + (rand - 0.5)*dxm ; % horizontal
        ym(m) = dym/2 + (im-1)*dym + (rand - 0.5)*dym ; % vertical

        % sticky air
        if ym(m) <= 50e3
            type(m) = 4 ;
            RHOm(m) = air(2) ;
            HRm(m) = air(3) ;
            Tm(m) = air(5) ;
            ETAm(m) = air(7) ;
            Km(m) = air(8) ;
            CPm(m) = air(9) ;
            ALm(m) = air(10) ;
        elseif ym(m) <= 100e3 && ym(m) > 50e3 % crust
            type(m) = 1 ;
            RHOm(m) = crust(2) ;
            HRm(m) = crust(3) ;
            Tm(m) = crust(5) + (crust(6)-crust(5))/(100e3 - 50e3)*(ym(m)-50e3) ;
            ETAm(m) = crust(7) ;
            Km(m) = .73 + 1293/(Tm(m)+77) ;
            CPm(m) = crust(9) ;
            ALm(m) = crust(10) ;
        elseif ym(m) > 100e3 ... % weak part
                && ym(m) <= 150e3 ...
                && ym(m) >= tan(angle)*(xm(m)-x1)+b...
                && ym(m) <= tan(angle)*(xm(m)-x2)+b
            type(m) = 3 ;
            RHOm(m) = weak(2) ;
            HRm(m) = weak(3) ;
            Tm(m) = weak(6) - (xm(m)+(150e3-ym(m))/tan(angle)-220e3-50e3)*((weak(6)-weak(5))/60e3) ;
            ETAm(m) = weak(7) ;
            Km(m) = .73 + 1293/(Tm(m)+77) ;
            CPm(m) = weak(9) ;
            ALm(m) = weak(10) ;
        elseif ym(m) > 150e3... % slab's tail
                && ym(m) <= 250e3 ...
                && ym(m) >= tan(angle)*(xm(m)-x1)+b...
                && ym(m) <= tan(angle)*(xm(m)-x2)+b ...
            type(m) = 1 ;
            RHOm(m) = crust(2) ;
            HRm(m) = crust(3) ;
            Tm(m) =  crust(6) - (xm(m)+(250e3-ym(m))/tan(angle)-220e3-150e3)*((crust(6)-crust(5))/60e3) ;
            ETAm(m) = crust(7) ;
            Km(m) = .73 + 1293/(Tm(m)+77) ;
            CPm(m) = crust(9) ;
            ALm(m) = crust(10) ;
        else % mantle
            type(m) = 2 ;
            RHOm(m) = mantl(2) ;
            HRm(m) = mantl(3) ;
            Tm(m) = mantl(5) ;
            ETAm(m) = mantl(7) ;
            Km(m) = .73 + 1293/(Tm(m)+77) ;
            CPm(m) = mantl(9) ;
            ALm(m) = mantl(10) ;
        end

        if ym(m) > 50e3... % inside part
                && ym(m) <= 100e3...
                && xm(m) >= 220e3 + (100e3-ym(m))/slope...
                && xm(m) <= 220e3+60e3-(100e3-ym(m))/tan(angle)
            type(m)=1;
            RHOm(m) = crust(1) ;
            HRm(m) = crust(3) ;
            Tm(m) = crust(6) - (xm(m)+(100e3-ym(m))/tan(angle)-220e3)*((crust(6)-crust(5))/60e3) ;
            ETAm(m) = crust(7) ;
            Km(m) = .73 + 1293/(Tm(m)+77) ;
            CPm(m) = crust(9) ;
            ALm(m) = crust(10) ;
        end

        % update m
        m = m + 1 ;
    end
end

% compute RHOCP
RHOCPm = RHOm .* CPm ;

% time loop --------------------------------------------------------------%

nstep = 200 ; % number of timestep

for ts = 1:1:nstep

    % initialize matrixes      
    RHOvy_SUM(:, :) = 0 ;
    ETAP_SUM(:, :)  = 0 ;
    ETAB_SUM(:, :)  = 0 ;
    Kx_SUM(:, :)    = 0 ;
    Ky_SUM(:, :)    = 0 ;
    RHOCP_SUM(:, :) = 0 ;
    T_SUM(:, :)     = 0 ;
    WTB_SUM(:, :)   = 0 ;
    WTP_SUM(:, :)   = 0 ;
    WTT_SUM(:, :)   = 0 ;
    WTVY_SUM(:, :)  = 0 ;
    WTVX_SUM(:, :)  = 0 ;
    HR_SUM(:, :)    = 0 ;
    AL_SUM(:, :)    = 0 ;
    DTSUM(:,:)      = 0 ;
    RHOCPSUMsub(:,:) = 0 ;
    
    for m = 1:1:Nm % FROM LAGRANGIAN POINTS TO EULERIAN GRID
        
        % make sure to take into account the markers that are located in
        % the defined grid. The fact that marker leave the model is very
        % unlikely as the timestep is sufficiently low !

        if (xm(m)>0 && xm(m)<xsize && ym(m)>0 && ym(m)<ysize)

            % Compute for each unknowns sums and weights
            
            % Basic Nodes ------------------------------------------------%
            j = fix((xm(m) - x(1))/dx) + 1 ;
            i = fix((ym(m) - y(1))/dy) + 1 ;
            
            % weight
            wtmij   = (1-(xm(m)-x(j))/dx) * (1-(ym(m)-y(i))/dy);
            wtmi1j  = (1-(xm(m)-x(j))/dx) * ((ym(m)-y(i))/dy) ;
            wtmij1  = ((xm(m)-x(j))/dx) * (1-(ym(m)-y(i))/dy) ;
            wtmi1j1 = ((xm(m)-x(j))/dx) * ((ym(m)-y(i))/dy) ;
            
            % update ETAB sum
            ETAB_SUM(i, j)     = ETAB_SUM(i, j) + ETAm(m)*wtmij ;
            ETAB_SUM(i+1, j)   = ETAB_SUM(i+1, j) + ETAm(m)*wtmi1j ;
            ETAB_SUM(i, j+1)   = ETAB_SUM(i, j+1) + ETAm(m)*wtmij1 ;
            ETAB_SUM(i+1, j+1) = ETAB_SUM(i+1, j+1) + ETAm(m)*wtmi1j1 ;
            
            % update weight sum
            WTB_SUM(i, j)     = WTB_SUM(i, j) + wtmij ;
            WTB_SUM(i+1, j)   = WTB_SUM(i+1, j) + wtmi1j ;
            WTB_SUM(i, j+1)   = WTB_SUM(i, j+1) + wtmij1 ;
            WTB_SUM(i+1, j+1) = WTB_SUM(i+1, j+1) + wtmi1j1 ;
    
            % Pressure Nodes ---------------------------------------------%
            j = fix((xm(m) - xp(1))/dx) + 1 ;
            i = fix((ym(m) - yp(1))/dy) + 1 ;
            
            % weight
            wtmij   = (1-(xm(m)-xp(j))/dx) * (1-(ym(m)-yp(i))/dy);
            wtmi1j  = (1-(xm(m)-xp(j))/dx) * ((ym(m)-yp(i))/dy) ;
            wtmij1  = ((xm(m)-xp(j))/dx) * (1-(ym(m)-yp(i))/dy) ;
            wtmi1j1 = ((xm(m)-xp(j))/dx) * ((ym(m)-yp(i))/dy) ;
            
            % update ETAP sum
            ETAP_SUM(i, j)     = ETAP_SUM(i, j) + ETAm(m)*wtmij ;
            ETAP_SUM(i+1, j)   = ETAP_SUM(i+1, j) + ETAm(m)*wtmi1j ;
            ETAP_SUM(i, j+1)   = ETAP_SUM(i, j+1) + ETAm(m)*wtmij1 ;
            ETAP_SUM(i+1, j+1) = ETAP_SUM(i+1, j+1) + ETAm(m)*wtmi1j1 ;
            
            % update RHOCP sum
            RHOCP_SUM(i,j)     = RHOCP_SUM(i,j)+RHOCPm(m)*wtmij ;
            RHOCP_SUM(i+1,j)   = RHOCP_SUM(i+1,j)+RHOCPm(m)*wtmi1j ;
            RHOCP_SUM(i,j+1)   = RHOCP_SUM(i,j+1)+RHOCPm(m)*wtmij1 ;
            RHOCP_SUM(i+1,j+1) = RHOCP_SUM(i+1,j+1)+RHOCPm(m)*wtmi1j1 ;
    
            % update T sum
            T_SUM(i,j)     = T_SUM(i,j) + Tm(m)*RHOCPm(m)*wtmij ;
            T_SUM(i+1,j)   = T_SUM(i+1,j) + Tm(m)*RHOCPm(m)*wtmi1j ;
            T_SUM(i,j+1)   = T_SUM(i,j+1) + Tm(m)*RHOCPm(m)*wtmij1 ;
            T_SUM(i+1,j+1) = T_SUM(i+1,j+1) + Tm(m)*RHOCPm(m)*wtmi1j1 ;

            % update Hr sum 
            HR_SUM(i,j) = HR_SUM(i,j) + HRm(m)*wtmij ;
            HR_SUM(i+1,j) = HR_SUM(i+1,j) + HRm(m)*wtmi1j ;
            HR_SUM(i,j+1) = HR_SUM(i,j+1) + HRm(m)*wtmij1 ;
            HR_SUM(i+1,j+1) = HR_SUM(i+1,j+1) + HRm(m)*wtmi1j1 ;
            
            % update alpha sum
            AL_SUM(i,j) = AL_SUM(i,j) + ALm(m)*wtmij ;
            AL_SUM(i+1,j) = AL_SUM(i+1,j) + ALm(m)*wtmi1j ;
            AL_SUM(i,j+1) = AL_SUM(i,j+1) + ALm(m)*wtmij1 ;
            AL_SUM(i+1,j+1) = AL_SUM(i+1,j+1) + ALm(m)*wtmi1j1 ;
                        
            % update weight sum for pressure
            WTP_SUM(i, j) = WTP_SUM(i, j) + wtmij ;
            WTP_SUM(i+1, j) = WTP_SUM(i+1, j) + wtmi1j ;
            WTP_SUM(i, j+1) = WTP_SUM(i, j+1) + wtmij1 ;
            WTP_SUM(i+1, j+1) = WTP_SUM(i+1, j+1) + wtmi1j1 ;     
    
            % update weight sum for temperature (including RHOCP)
            WTT_SUM(i,j) = WTT_SUM(i,j)+RHOCPm(m)*wtmij ;
            WTT_SUM(i+1,j) = WTT_SUM(i+1,j)+RHOCPm(m)*wtmi1j ;
            WTT_SUM(i,j+1) = WTT_SUM(i,j+1)+RHOCPm(m)*wtmij1 ;
            WTT_SUM(i+1,j+1) = WTT_SUM(i+1,j+1)+RHOCPm(m)*wtmi1j1 ;
    
            % Vx nodes ---------------------------------------------------%
            j = fix((xm(m)-xvx(1))/dx)+1 ;
            i = fix((ym(m)-yvx(1))/dy)+1 ;
    
            wtmij = (1- (xm(m) - xvx(j))/dx)*(1-(ym(m)-yvx(i))/dy) ;
            wtmi1j = (1- (xm(m) - xvx(j))/dx)*((ym(m)-yvx(i))/dy) ;
            wtmij1 = ((xm(m) - xvx(j))/dx)*(1-(ym(m)-yvx(i))/dy) ;
            wtmi1j1 = ((xm(m) - xvx(j))/dx)*((ym(m)-yvx(i))/dy) ;
            
            % update k on vx nodes sum
            Kx_SUM(i,j) = Kx_SUM(i,j) + Km(m)*wtmij ;
            Kx_SUM(i+1,j) = Kx_SUM(i+1,j) + Km(m)*wtmi1j ;
            Kx_SUM(i,j+1) = Kx_SUM(i,j+1) + Km(m)*wtmij1 ;
            Kx_SUM(i+1,j+1) = Kx_SUM(i+1,j+1) + Km(m)*wtmi1j1 ;
            
            % update weights on vy nodes
            WTVX_SUM(i,j) = WTVX_SUM(i,j)+wtmij ;
            WTVX_SUM(i+1,j) = WTVX_SUM(i+1,j)+wtmi1j ;
            WTVX_SUM(i,j+1) = WTVX_SUM(i,j+1)+wtmij1 ;
            WTVX_SUM(i+1,j+1) = WTVX_SUM(i+1,j+1)+wtmi1j1 ;
                    
            % Vy nodes ---------------------------------------------------%
            j = fix((xm(m) - xvy(1))/dx) + 1 ;
            i = fix((ym(m) - yvy(1))/dy) + 1 ;
    
            % weight
            wtmij =  (1-(xm(m)-xvy(j))/dx) * (1-(ym(m)-yvy(i))/dy);
            wtmi1j = (1-(xm(m)-xvy(j))/dx) * ((ym(m)-yvy(i))/dy) ;
            wtmij1 = ((xm(m)-xvy(j))/dx) * (1-(ym(m)-yvy(i))/dy) ;
            wtmi1j1 = ((xm(m)-xvy(j))/dx) * ((ym(m)-yvy(i))/dy) ;
            
            % update rho sum
            RHOvy_SUM(i, j) = RHOvy_SUM(i, j) + RHOm(m)*wtmij ;
            RHOvy_SUM(i+1, j) = RHOvy_SUM(i+1, j) + RHOm(m)*wtmi1j ;
            RHOvy_SUM(i, j+1) = RHOvy_SUM(i, j+1) + RHOm(m)*wtmij1 ;
            RHOvy_SUM(i+1, j+1) = RHOvy_SUM(i+1, j+1) + RHOm(m)*wtmi1j1 ;
    
            % update k on vy nodes sum
            Ky_SUM(i,j) = Ky_SUM(i,j) + Km(m)*wtmij ;
            Ky_SUM(i+1,j) = Ky_SUM(i+1,j) + Km(m)*wtmi1j ;
            Ky_SUM(i,j+1) = Ky_SUM(i,j+1) + Km(m)*wtmij1 ;
            Ky_SUM(i+1,j+1) = Ky_SUM(i+1,j+1) + Km(m)*wtmi1j1 ;
    
            % update weight sum
            WTVY_SUM(i, j) = WTVY_SUM(i, j) + wtmij ;
            WTVY_SUM(i+1, j) = WTVY_SUM(i+1, j) + wtmi1j ;
            WTVY_SUM(i, j+1) = WTVY_SUM(i, j+1) + wtmij1 ;
            WTVY_SUM(i+1, j+1) = WTVY_SUM(i+1, j+1) + wtmi1j1 ;
        end
    end
    
    for i = 1:1:Ny+1 % COMPUTE VALUES FOR EULERIAN GRID
        for j = 1:1:Nx+1
            if WTP_SUM(i,j) > 0
                RHOCP(i,j) = RHOCP_SUM(i,j)/WTP_SUM(i,j) ;
                ETAP(i,j) = ETAP_SUM(i,j)/WTP_SUM(i,j) ;
                T0(i,j) = T_SUM(i,j)/WTT_SUM(i,j) ;
                HR(i,j) = HR_SUM(i,j)/WTP_SUM(i,j) ;
                ALPHA(i,j) = AL_SUM(i,j)/WTP_SUM(i,j) ;
            end
            if WTB_SUM(i,j) > 0
                ETAB(i,j) = ETAB_SUM(i,j)/WTB_SUM(i,j) ;
            end
            if WTVX_SUM(i,j) > 0
                kx(i,j) = Kx_SUM(i,j)/WTVX_SUM(i,j) ;
            end
            if WTVY_SUM(i,j) > 0
                RHOvy(i,j) = RHOvy_SUM(i,j)/WTVY_SUM(i,j) ;
                ky(i,j) = Ky_SUM(i,j)/WTVY_SUM(i,j) ;
            end
        end
    end

    T0(:,Nx+1) = T0(:,Nx) ;
    T0(:,1) = T0(:,2);
    T0(1,2:Nx) = 2*air(5)-T0(2,2:Nx);
    T0(Ny+1,2:Nx) = 2*mantl(6) - T0(Ny,2:Nx);

    % increase current timestep size dt by 20% ---------------------------%

    dt = dt*1.2 ;

    % introduce dt adjustment --------------------------------------------%

    ndtiter = 3 ;
    for dtiter = 1:1:ndtiter % MECHANICAL AND TEMPERATURE SOLVER
    
        % mechanical solver ----------------------------------------------%
        
        for i = 1:1:Ny+1 % MECHANICAL SOLVER
            for j = 1:1:Nx+1 
                
                % general indexes
                gp = ((j - 1) * (Ny + 1) + (i - 1)) * n + 1 ;
                gvx = gp + 1 ;
                gvy = gp + 2 ;
                
                % x - stokes
                if (i == 1 || j == 1 || i == Ny + 1 || j == Nx) % general BC
                    L(gvx, gvx) = 1 ; 
                    R(gvx,   1) = 0 ;
        
                    if (i == 1 && j>1 && j<Nx) % top BC vx 
                        L(gvx, gvx + n) = bctype ;
        
                    elseif (i == Ny + 1 && j>1 && j<Nx) % bottom BC vx
                        L(gvx, gvx - n) = bctype ;
        
                    end
        
                elseif (j == Nx + 1) % FC
                    L(gvx, gvx) = 1 ;
                    R(gvx, 1) = 0 ;
                    
                else % x-stokes equations coefficients
                    L(gvx, gvx - n*(Ny+1)) =  2*ETAP(i, j)/dx^2 ; % vx1
                    L(gvx, gvx - n)        =  ETAB(i-1, j)/dy^2 ; % vx2
                    L(gvx, gvx)            = -2*ETAP(i, j)/dx^2 - 2*ETAP(i, j+1)/dx^2 - ETAB(i-1,j)/dy^2 - ETAB(i,j)/dy^2; % vx3
                    L(gvx, gvx + n)        =  ETAB(i, j)/dy^2 ; % vx4
                    L(gvx, gvx + n*(Ny+1)) =  2*ETAP(i,j+1)/dx^2 ; % vx5
        
                    L(gvx, gvy - n)              =  ETAB(i-1, j)/(dx*dy) ; % vy1
                    L(gvx, gvy)                  = -ETAB(i, j)/(dx*dy) ; % vy2
                    L(gvx, gvy - n + (Ny + 1)*n) = -ETAB(i-1, j)/(dx*dy) ; % vy3
                    L(gvx, gvy + (Ny+1)*n)       =  ETAB(i, j)/(dx*dy) ; %vy 4
        
                    L(gvx, gp) = 1/dx ; % P1
                    L(gvx, gp + (Ny+1)*n) = -1/dx ; % P2
        
                    R(gvx, 1) = 0 ; % RHS
                end
                
                % y - stokes 
                if (i == 1 || j == 1 || i == Ny || j == Nx + 1) % general BC
                    L(gvy, gvy) = 1 ;
                    R(gvy,   1) = 0 ;
        
                    if (j == 1 && i > 1 && i < Ny) % left BC vy
                        L(gvy, gvy + n * (Ny + 1)) = bctype ;
        
                    elseif (j == Nx + 1 && i > 1 && i < Ny) % right BC vy
                        L(gvy, gvy - n * (Ny + 1)) = bctype ;
                    end
                    
                elseif (i == Ny + 1) % FC
                    L(gvy, gvy) = 1 ;
                    R(gvy, 1) = 0 ;
        
                else % y-stokes equations coefficients
                    dRHOdx = (RHOvy(i,j+1)-RHOvy(i,j-1))/(2*dx) ;
                    dRHOdy = (RHOvy(i+1,j)-RHOvy(i-1,j))/(2*dy) ;
    
                    L(gvy, gvy - n * (Ny + 1))   =  ETAB(i, j-1)/dx^2 ; % vy1
                    L(gvy, gvy - n)              =  2*ETAP(i, j)/dy^2 ; % vy2
                    L(gvy, gvy)                  = -2*ETAP(i+1, j)/dy^2 - 2*ETAP(i, j)/dy^2 - ETAB(i, j)/dx^2 - ETAB(i, j-1)/dx^2 - dt*dRHOdy*g ; % vy3
                    L(gvy, gvy + n)              =  2*ETAP(i+1, j)/dy^2; % vy4
                    L(gvy, gvy + n * (Ny + 1))   =  ETAB(i, j)/dx^2; % vy5
                    
                    L(gvy, gvx - n * (Ny + 1))   =  ETAB(i, j-1)/(dx*dy)-.25*dt*dRHOdx*g ; %vx1
                    L(gvy, gvx - (Ny + 1)*n + n) = -ETAB(i, j-1)/(dx*dy)-.25*dt*dRHOdx*g ; % vx2
                    L(gvy, gvx)                  = -ETAB(i, j)/(dx*dy)-.25*dt*dRHOdx*g  ; % vx3
                    L(gvy, gvx + n)              =  ETAB(i, j)/(dx*dy)-.25*dt*dRHOdx*g  ; % vx4
        
                    L(gvy, gp) = 1/dy ; % P1
                    L(gvy, gp + n) = -1/dy ; % P2
        
                    R(gvy, 1) = -g * RHOvy(i, j) ; % RHS
                end
                
                % continuity 
                if (j == 1 || i == 1 || j == Nx + 1 || i == Ny+1) % BC
                    L(gp, gp) = 1 ;
                    R(gp,  1) = 0 ;
    
                elseif (i == 2 && j == 2) % FC
                    L(gp, gp) = 1 ;
                    R(gp,  1) = 1e+5 ;
    
                else 
                    L(gp, gvx) = 1/dx ; % vx2
                    L(gp, gvy) = 1/dy ; % vy2
                    L(gp, gvx - n * (Ny + 1)) = -1/dx ; % vx1
                    L(gp, gvy - n) = -1/dy ; % vy1
                    R(gp, 1) = 0 ;
                end
            end
        end
        
        S = L\R ; % compute solution
        
        for i = 1:1:Ny + 1 % STORE DATA IN EULERIAN GRID
            for j = 1:1:Nx + 1  
                gp  = ((j - 1) * (Ny + 1) + (i - 1)) * n + 1 ;
                gvx = gp + 1 ;
                gvy = gp + 2 ;
                vx(i, j) = S(gvx, 1) ;
                vy(i, j) = S(gvy, 1) ;
                p (i, j) = S(gp, 1) ;
            end
        end
    
        % dt (defined) ---------------------------------------------------%

        if (dtiter<ndtiter)
            vxmax = max(max(abs(vx))) ;
            if (vxmax*dt > dx*dismax)
                dt = dx*dismax/vxmax ;
            end
            vymax = max(max(abs(vy))) ;
            if (vymax*dt > dy*dismax)
                dt = dy*dismax/vymax ;
            end
        end

        % Hr, Hs, Ha -----------------------------------------------------%
        
        % compute shear heat Hs ------------------------------------------%

        % sigma and epsilon
        for i = 1:1:Ny+1 % epsilon and sigma XX
            for j = 2:1:Nx+1
                EXX(i,j) = (vx(i,j) - vx(i,j-1))/dx ;
                SXX(i,j) = 2*ETAP(i,j) * EXX(i,j) ;
            end
        end

        for i = 1:1:Ny % epsilon and sigma XY
            for j = 2:1:Nx+1
                EXY(i,j) = .5*((vx(i+1,j)-vx(i,j))/dy + (vy(i,j) - vy(i,j-1))/dx) ;
                SXY(i,j) = 2*ETAB(i,j)*EXY(i,j) ;
            end
        end

        for i = 2:1:Ny+1 % mean of the 4 nodal points
            for j = 2:1:Nx+1
                SEAV(i,j) = (SXY(i,j)*EXY(i,j)+SXY(i-1,j)*EXY(i-1,j)+...
                    SXY(i,j-1)*EXY(i,j-1)+SXY(i-1,j-1)*EXY(i-1,j-1))/4 ;
            end
        end
        
        HS = EXX(i,j)*SXX(i,j) + 2*SEAV ; % given formula for Hs

        % compute adiabatic heat Ha --------------------------------------%

        for i = 2:1:Ny+1 % Ha
            for j = 1:1:Nx+1
                if dtiter == 1 % needed for the first iteration
                    TDT = T0 ;
                end
                HA(i,j) = ALPHA(i,j)*(T0(i,j)+TDT(i,j))/2*...
                    (RHOvy(i-1,j)*vy(i-1,j)*g/2 + RHOvy(i,j)*vy(i,j)*g/2) ;
            end
        end
    
        % thermal solver -------------------------------------------------%
    
        for i = 2:1:Ny % THERMAL SOLVER
            for j = 2:1:Nx
                gt = ((j-1)*(Ny+1) + (i-1)) + 1 ;
                LT(gt,gt-(Ny+1)) = -kx(i,j-1)/dx^2;
                LT(gt,gt-1)      = -ky(i-1,j)/dy^2;
                LT(gt,gt)        =  RHOCP(i,j)/dt + ...
                    (kx(i,j)+kx(i,j-1))/dx^2 + (ky(i,j)+ky(i-1,j))/dy^2;
                LT(gt,gt+1)      = -ky(i,j)/dy^2;
                LT(gt,gt+(Ny+1)) = -kx(i,j)/dx^2;
                RT(gt)          =  RHOCP(i,j)*T0(i,j)/dt + HS(i,j) + HR(i,j) + HA(i,j);
            end
        end
    
        for i = 1:1:Ny+1 % THERMAL SOLVER BC
            for j = 1:1:Nx+1
                gT = ((j-1)*(Ny+1) + (i-1)) + 1 ;
                if i == 1
                    LT(gT, gT + 1) = 1 ;
                    LT(gT,gT) = 1 ;
                    RT(gT) = 2*air(5) ;
                elseif i == Ny+1
                    LT(gT, gT - 1) = 1 ;
                    LT(gT,gT) = 1 ;
                    RT(gT) = 2*mantl(5) ;
                elseif j == 1 
                    LT(gT,gT)=1;
                    LT(gT,gT+(Ny+1))=-1;
                    RT(gT)=0;
                elseif j == Nx+1
                    LT(gT,gT)=1;
                    LT(gT,gT-(Ny+1))=-1;
                    RT(gT)=0;
                end
            end
        end
    
        ST = LT\RT ; % compute solution
        TDT = reshape(ST, [Ny+1 Nx+1]) ; % STORE DATA IN EULERIAN GRID
    
        DT = TDT - T0 ; % compute temperature difference
    
        if dtiter<ndtiter % condition for DT max
            maxcurrentDT = max(max(abs(DT))) ;
            if maxcurrentDT>DTmax
                dt = .7*dt*DTmax/maxcurrentDT ;
            end
        end

    % end of the adjustment
    end

    % plot the datas -----------------------------------------------------%

    if (ts == 1 || mod(ts, 20) == 0)

        f = figure(1) ;     
        f.Position = [200 200 1100 700]; % change figure size if needed
        % plot the results

        subplot(3,4,1)
        h = pcolor(x, y, RHOvy) ;
%         hold on
%         quiver(x, y, vx, vy, 'k')
        c = colorbar ;
        title('Density')
        ylabel(c, '\rho [kg/m^3]')
        set(h, 'EdgeColor', 'none')
        axis ij

        subplot(3,4,2)
        h = pcolor(x, y, log10(ETAB)) ;
        c = colorbar ;
        title('Viscosity B')
        ylabel(c, '\eta_b [Pa s]')
        set(h, 'EdgeColor', 'none')
        axis ij
        shading interp
        
        subplot(3,4,3)
        pcolor(xp, yp, TDT)
        c = colorbar;
        shading interp
        title('Temperature')
        ylabel(c, 'T [K]')
        axis ij

        subplot(3,4,4)
        pcolor(xp, yp, p)
        c = colorbar ;
        title('Pressure')
        ylabel(c, 'p [Pa]')
        shading interp
        axis ij

        subplot(3,4,5)
        pcolor(xvx, yvx, vx)
        c = colorbar ;
        title('x-velocity')
        ylabel(c, 'v_x [m/s]')
        shading interp
        axis ij

        subplot(3,4,6)
        pcolor(xvy, yvy, vy)
        c = colorbar;
        shading interp
        title('y-velocity')
        ylabel(c, 'v_y [m/s]')
        axis ij

        subplot(3,4,7)
        pcolor(xp, yp, RHOCP)
        c = colorbar;
        shading interp
        title('RHOCP')
        ylabel(c, '\rho C_p [J/m^3/K]')
        axis ij

        subplot(3,4,8)
        pcolor(xvx, yvx, log10(kx))
        c = colorbar ;
        title('Conudctivity')
        ylabel(c, 'log10(K_x) [W/m/K]')
        shading interp
        axis ij

        subplot(3,4,9)
        pcolor(xp, yp, ALPHA)
        c = colorbar ;
        title('Alpha')
        ylabel(c, '\alpha [ ]')
        shading interp
        axis ij

        subplot(3,4,10)
        pcolor(xp, yp, HR)
        c = colorbar ;
        title('Radiogenic')
        ylabel(c, 'H_r [W/m^3]')
        shading interp
        axis ij

        subplot(3,4,11)
        pcolor(xp, yp, HS)
        c = colorbar ;
        title('Shear')
        ylabel(c, 'H_s [W/m^3]')
        shading interp
        axis ij

        subplot(3,4,12)
        pcolor(xp, yp, HA)
        c = colorbar ;
        title('Adiabatic')
        ylabel(c, 'H_a [W/m^3]')
        shading interp
        axis ij

        set(findall(gcf,'-property','FontSize'),'FontSize',fs)
        sgtitle(['Slab break-off modelling, timestep = ', num2str(ts)])
    end

    % interpolation ------------------------------------------------------%
    
    for m = 1:1:Nm % managing dt
        i = fix((ym(m) - yp(1))/dy) + 1 ;
        j = fix((xm(m) - xp(1))/dx) + 1 ;
        dxmj = (xm(m) - xp(j))/dx ;
        dymi = (ym(m) - yp(i))/dy ;

        wtmij = (1-dxmj)*(1-dymi) ;
        wtmi1j = (1-dxmj)*dymi ;
        wtmij1 = dxmj*(1-dymi) ;
        wtmi1j1 = dxmj*dymi ;

        T0m(m) = T0(i,j)*wtmij + T0(i+1,j)*wtmi1j + T0(i,j+1)*wtmij1 + T0(i+1,j+1)*wtmi1j1 ;
        DTm0(m) = Tm(m)-T0m(m) ;
        DTm(m) = DTm0(m)*exp((-2*dt*Km(m))*(1/dx^2 + 1/dy^2)/RHOCPm(m)) ;
        bigDTm(m) = DTm(m) - DTm0(m) ;

        DTSUM(i,j) = DTSUM(i,j)+bigDTm(m)*RHOCPm(m)*wtmij ;
        DTSUM(i+1,j) = DTSUM(i+1,j)+bigDTm(m)*RHOCPm(m)*wtmi1j ;
        DTSUM(i,j+1) = DTSUM(i,j+1)+bigDTm(m)*RHOCPm(m)*wtmij1 ;
        DTSUM(i+1,j+1) = DTSUM(i+1,j+1)+bigDTm(m)*RHOCPm(m)*wtmi1j1;

        RHOCPSUMsub(i,j) = RHOCPSUMsub(i,j)+RHOCPm(m)*wtmij ; 
        RHOCPSUMsub(i+1,j) = RHOCPSUMsub(i+1,j)+RHOCPm(m)*wtmi1j ;
        RHOCPSUMsub(i,j+1) = RHOCPSUMsub(i,j+1)+RHOCPm(m)*wtmij1 ;
        RHOCPSUMsub(i+1,j+1) = RHOCPSUMsub(i+1,j+1)+RHOCPm(m)*wtmi1j1 ;
    end

    if RHOCPSUMsub(i,j)>0 % calculate DTsub
        DTsub(i,j) = DTSUM(i,j)/RHOCPSUMsub(i,j) ;
    end


    DTr = DT - DTsub ; % dt remaining

    for m = 1:1:Nm  % THERMAL INTERPOLATION
        % same condition as in the beggining
        if (xm(m)>0 && xm(m)<xsize && ym(m)>0 && ym(m)<ysize)
            i = fix((ym(m) - yp(1))/dy) + 1 ;
            j = fix((xm(m) - xp(1))/dx) + 1 ;
    
            dxmj = (xm(m) - xp(j))/dx ;
            dymi = (ym(m) - yp(i))/dy ;
    
            wtmij = (1-dxmj)*(1-dymi) ;
            wtmi1j = (1-dxmj)*dymi ;
            wtmij1 = dxmj*(1-dymi) ;
            wtmi1j1 = dxmj*dymi ;
    
            if ts == 1 % don't use DT for the first interpolation
                Tm(m) = TDT(i,j)*wtmij+TDT(i+1,j)*wtmi1j+TDT(i,j+1)*wtmij1+TDT(i+1,j+1)*wtmi1j1 ;
            else
                DTrm(m) = DTr(i,j)*wtmij + DTr(i+1,j)*wtmi1j + DTr(i,j+1)*wtmij1 + DTr(i+1,j+1)*wtmi1j1 ;
                Tm(m) = Tm(m) + bigDTm(m) + DTrm(m) ;
            end

            % we must constrain the model in the grid
            if i<2
                i=2 ;
            elseif i>Ny-1
                i=Ny-1;
            end
            if j<2
                j=2;
            elseif j>Nx-1
                j = Nx-1 ;
            end

            % pressure and XX epsilon marker
            Pm(m) = p(i,j)*wtmij + p(i+1,j)*wtmi1j + p(i,j+1)*wtmij1 + p(i+1,j+1)*wtmi1j1 ;
            EXXm(m) = EXX(i,j)*wtmij+EXX(i+1,j)*wtmi1j+EXX(i,j+1)*wtmij1+EXX(i+1,j+1)*wtmi1j1 ;

            j = fix((xm(m)-x(1))/dx)+1;
            i = fix((ym(m)-y(1))/dy)+1;

            dxmj = (xm(m) - x(j))/dx ;
            dymi = (ym(m) - y(i))/dy ;

            wtmij = (1-dxmj)*(1-dymi) ;
            wtmi1j = (1-dxmj)*dymi ;
            wtmij1 = dxmj*(1-dymi) ;
            wtmi1j1 = dxmj*dymi ;

            EXYm(m) = EXY(i,j)*wtmij+EXY(i+1,j)*wtmi1j+EXY(i,j+1)*wtmij1+EXY(i+1,j+1)*wtmi1j1 ;

            if type(m) == 1 % crust parameters
                Km(m) = .73 + 1293/(Tm(m)+77) ;
                RHOm(m) = crust(1)*(1+crust(11)*(Pm(m)-1e+5))...
                    /(1+crust(10)*(Tm(m)-273)) ;
                RHOCPm(m) = RHOm(m)*crust(9) ;
                EII = sqrt(EXXm(m)^2 + EXYm(m)^2) ;
                ETAm(m) = 1/2/ADdry^(1/ndry)*EII^(1/ndry-1)...
                    *exp((Eadry+Pm(m)*Vadry)/8.314/Tm(m)/ndry) ;
                
                if ETAm(m)>ETAmax
                    ETAm(m)=ETAmax ;
                elseif ETAm(m)<ETAmin
                    ETAm(m) = ETAmin ;
                end
                SII = 2*ETAm(m)*EII ;
                if SII > crust(4)
                    ETAm(m)=crust(4)/2/EII ;
                end
            end

            if type(m) == 2 % mantle parameters
                Km(m) = .73 + 1293/(Tm(m)+77) ;
                RHOm(m) = mantl(1)*(1+mantl(11)*(Pm(m)-1e+5))...
                    /(1+mantl(10)*(Tm(m)-273)) ;
                RHOCPm(m) = RHOm(m)*mantl(9) ;
                EII = sqrt((EXXm(m)^2+EYYm(m)^2)/2 + EXYm(m)^2) ;
                ETAm(m) = 1/2/ADwet^(1/nwet)*EII^(1/nwet-1)...
                    *exp((Eawet+Pm(m)*Vawet)/8.314/Tm(m)/nwet) ;
                
                if ETAm(m)>ETAmax
                    ETAm(m)=ETAmax ;
                elseif ETAm(m)<ETAmin
                    ETAm(m) = ETAmin ;
                end
                SII = 2*ETAm(m)*EII ;
                if SII > mantl(4)
                    ETAm(m)=mantl(4)/2/EII ;
                end

            end

            if type(m) == 3 % weak part parameters
                Km(m) = .73 + 1293/(Tm(m)+77) ;
                RHOm(m) = weak(1)*(1+weak(11)*(Pm(m)-1e+5))...
                    /(1+weak(10)*(Tm(m)-273)) ;
                RHOCPm(m) = RHOm(m)*weak(9) ;
                EII = sqrt((EXXm(m)^2+EYYm(m)^2)/2 + EXYm(m)^2) ;
                ETAm(m) = 1/2/ADdry^(1/ndry)*EII^(1/ndry-1)...
                    *exp((Eadry+Pm(m)*Vadry)/8.314/Tm(m)/ndry) ;
                if ETAm(m)>ETAmax
                    ETAm(m)=ETAmax ;
                elseif ETAm(m)<ETAmin
                    ETAm(m) = ETAmin ;
                end
                SII = 2*ETAm(m)*EII ;
                if SII > weak(4)
                    ETAm(m)=weak(4)/2/EII ;
                end
            end
        end
    end

    for m = 1:1:Nm % MECHANICAL INTERPOLATION    
        % same condition as in the beggining
        if (xm(m)>0 && xm(m)<xsize && ym(m)>0 && ym(m)<ysize)
            xa = xm(m) ;
            ya = ym(m) ;
            xrk = xa ;
            yrk = ya ;
            for rk = 1:1:4
                % vx :
                i = fix((yrk - yvx(1))/dy) + 1 ;
                j = fix((xrk - xvx(1))/dx) + 1 ;
            
                dxmj = xrk - xvx(j) ; 
                dymi = yrk - yvx(i) ;
                vx130 = (1-dxmj/dx) * vx(i,j) + (dxmj/dx) * vx(i,j+1) ;
                vx240 = (1-dxmj/dx) * vx(i+1,j) + (dxmj/dx) * vx(i+1,j+1) ;
                vx13cor = 0 ; % initialize 
                vx24cor = 0 ; % initialize  
    
                if (i>=1 && j>1 && i<=Ny+1 && j<Nx-1)
                    if (dxmj > dx/2)
                        % vx correction on point 1 and 3 (right hand)
                        vx13cor = .5*(dxmj/dx-0.5)^2*(vx(i,j)-2*vx(i,j+1)+vx(i,j+2)) ;
                        % vx correction on point 2 and 4 (right hand)
                        vx24cor = .5*(dxmj/dx-0.5)^2*(vx(i+1,j)-2*vx(i+1,j+1)+vx(i+1,j+2)) ; 
                    elseif (dxmj < dx/2)
                        % vx correction on point 1 and 3 (left hand)
                        vx13cor = .5*(dxmj/dx-0.5)^2*(vx(i,j-1)-2*vx(i,j)+vx(i,j+1)) ; 
                        % vx correction on point 2 and 4 (left hand)
                        vx24cor = .5*(dxmj/dx-0.5)^2*(vx(i+1,j-1)-2*vx(i+1,j)+vx(i+1,j+1)) ;
                    end
                end
                
                vxm(rk) = (vx130+vx13cor)*(1-dymi/dy) + (vx240+vx24cor)*(dymi/dy) ;
    
                % vy :
                i = fix((yrk - yvy(1))/dy) + 1 ;
                j = fix((xrk - xvy(1))/dx) + 1 ;
    
                dxmj = xrk - xvy(j) ;
                dymi = yrk - yvy(i) ;
                vy130 = (1-dymi/dy) * vy(i,j) + (dymi/dy) * vy(i+1,j) ;
                vy240 = (1-dymi/dy) * vy(i,j+1) + (dymi/dy) * vy(i+1,j+1) ;
                vy13cor = 0 ; % initialization
                vy24cor = 0 ; % initialization
    
                if (i>1 && j>=1 && i<Ny-1 && j<=Nx+1)
                    if (dymi < dy/2)
                        % vy correction on point 1 and 3 (right hand)
                        vy13cor = .5*(dymi/dy-0.5)^2*(vy(i-1,j)-2*vy(i,j)+vy(i+1,j)) ;
                        % vy correction on point 2 and 4 (right hand)
                        vy24cor = .5*(dymi/dy-0.5)^2*(vy(i-1,j+1)-2*vy(i,j+1)+vy(i+1,j+1)) ;
                    elseif (dymi > dy/2)
                        % vy correction on point 1 and 3 (left hand)
                        vy13cor = .5*(dymi/dy-0.5)^2*(vy(i,j)-2*vy(i+1,j)+vy(i+2,j)) ;
                        % vy correction on point 2 and 4 (left hand)
                        vy24cor = .5*(dymi/dy-0.5)^2*(vy(i,j+1)-2*vy(i+1,j+1)+vy(i+2,j+1)) ;
                    end
                end
    
                vym(rk) = (vy130+vy13cor) * (1-dxmj/dx) + (vy240+vy24cor) * (dxmj/dx) ;            
    
                % update variables
                if rk <= 2
                    xrk = xa + vxm(rk)*dt/2 ;
                    yrk = ya + vym(rk)*dt/2 ;
                elseif rk == 3
                    xrk = xa + vxm(rk)*dt ;
                    yrk = ya + vym(rk)*dt ;
                end
            end
            
            vxmeff = 1/6*(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4)) ;
            vymeff = 1/6*(vym(1)+2*vym(2)+2*vym(3)+vym(4)) ;
    
            % movement of the marker
            xm(m) = xm(m) + vxmeff*dt ;
            ym(m) = ym(m) + vymeff*dt ;
        end
    end

    timesum = timesum + dt ;

    % references values --------------------------------------------------%

    if ts == 1 || ts == 10
        p (27,12) 
        vx(27,12)  
        vy(27,12) 
        TDT(27,12) 
        dt 
        timesum 
    end
end