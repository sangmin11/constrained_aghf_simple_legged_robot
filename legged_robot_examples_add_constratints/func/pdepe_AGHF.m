function varargout = pdepe_AGHF(xmesh,t,X,params,options)
%customized PDEPE for AGHF (modified from original pdepe)

m=0;

nt = length(t);
if nt < 3
  error(message('MATLAB:pdepe:TSPANnotEnoughPts'))
end
if any(diff(t) <= 0)
  error(message('MATLAB:pdepe:TSPANnotIncreasing'))
end

xmesh = xmesh(:);
if m > 0 && xmesh(1) < 0
  error(message('MATLAB:pdepe:NegXMESHwithPosM'))
end
nx = length(xmesh);
if nx < 3
  error(message('MATLAB:pdepe:XMESHnotEnoughPts'))
end
if any(diff(xmesh) <= 0)
  error(message('MATLAB:pdepe:XMESHnotIncreasing'))
end

% Initialize the nx-1 points xi where functions will be evaluated
% and as many coefficients in the difference formulas as possible.
% For problems with coordinate singularity, use 'singular' formula
% on all subintervals.
singular = (xmesh(1) == 0 && m > 0);
xL = xmesh(1:end-1);
xR = xmesh(2:end);
xM = xL + 0.5*(xR-xL);
switch m
 case 0
  xi = xM;
  zeta = xi;
 case 1
  if singular
    xi = (2/3)*(xL.^2 + xL.*xR + xR.^2) ./ (xL+xR);
  else
    xi = (xR-xL) ./ log(xR./xL);
  end
  zeta = (xi .* xM).^(1/2);
 case 2
  if singular
    xi = (2/3)*(xL.^2 + xL.*xR + xR.^2) ./ (xL+xR);
  else
    xi = xL .* xR .* log(xR./xL) ./ (xR-xL);
  end
  zeta = (xL .* xR .* xM).^(1/3);
end
if singular
  xim = (zeta .^ (m+1))./xi;
else
  xim = xi .^ m;
end
zxmp1 = (zeta.^(m+1) - xL.^(m+1)) / (m+1);
xzmp1 = (xR.^(m+1) - zeta.^(m+1)) / (m+1);

% Form the initial values with a column of unknowns at each
% mesh point and then reshape into a column vector for dae.
temp = mypdexic_generated(xmesh(1),X);
if size(temp,1) ~= length(temp)
  error(message('MATLAB:pdepe:InvalidOutputICFUN'))
end

npde = length(temp);
y0 = zeros(npde,nx);
y0(:,1) = temp;
for j = 2:nx
  y0(:,j) = mypdexic_generated(xmesh(j),X);
end

% Classify the equations so that a constant, diagonal mass matrix
% can be formed. The entries of c are to be identically zero or
% vanish only at the entries of xmesh, so need be checked only at one
% point not in the mesh.
[U,Ux] = pdentrp_AGHF(singular,m,xmesh(1),y0(:,1),xmesh(2),y0(:,2),xi(1));
[c,f,s] = mypdexpde_generated(xi(1),t(1),U,Ux,params);
if any([size(c,1),size(f,1),size(s,1)]~=npde)
  error(message('MATLAB:pdepe:UnexpectedOutputPDEFUN',sprintf('%d',npde)))
end
[pL,qL,pR,qR] = mypdexbc_generated(xmesh(1),y0(:,1),xmesh(nx),y0(:,nx),t(1),X);
if any([size(pL,1),size(qL,1),size(pR,1),size(qR,1)]~=npde)
  error(message('MATLAB:pdepe:UnexpectedOutputBCFUN',sprintf('%d',npde)))
end

D = ones(npde,nx);
D( c == 0, 2:nx-1) = 0;
if ~singular
  D( qL == 0, 1) = 0;
end
D( qR == 0, nx) = 0;
M = spdiags( D(:), 0, npde*nx, npde*nx);

% Construct block-diagonal pattern of Jacobian matrix
S = kron( spdiags(ones(nx,3),-1:1,nx,nx), ones(npde,npde));

% Extract relevant options and augment with new ones
reltol = odeget(options,'RelTol',[],'fast');
abstol = odeget(options,'AbsTol',[],'fast');
normcontrol = odeget(options,'NormControl',[],'fast');
initialstep = odeget(options,'InitialStep',[],'fast');
maxstep = odeget(options,'MaxStep',[],'fast');
events = odeget(options,'Events',[],'fast');  % events(m,t,xmesh,umesh)
hasEvents = ~isempty(events);
if hasEvents
  eventfcn = @(t,y) events(m,t,xmesh,y);
  %eventfcn = @(t,y) events(m,t,xmesh,y,varargin{:});
else
  eventfcn = [];
end
opts = odeset('RelTol',reltol,'AbsTol',abstol,'NormControl',normcontrol,...
              'InitialStep',initialstep,'MaxStep',maxstep,'Events',eventfcn,...
              'Mass',M,'JPattern',S);

% Call DAE solver
tfinal = t(end);
try
  if hasEvents
    [t,y,te,ye,ie] = ode15s(@pdeodes,t,y0(:),opts);
  else
    [t,y] = ode15s(@pdeodes,t,y0(:),opts);
  end
catch ME
  if strcmp(ME.identifier,'MATLAB:daeic12:IndexGTOne')
    error(message('MATLAB:pdepe:SpatialDiscretizationFailed'))
  else
    rethrow(ME);
  end
end

% Verify and process the solution
if t(end) ~= tfinal
  nt = length(t);
  if ~hasEvents || (te(end) ~= t(end))  % did not stop on a terminal event
    warning(message('MATLAB:pdepe:TimeIntegrationFailed',sprintf('%e',t(end))));
  end
end

sol = zeros(nt,nx,npde);
for i = 1:npde
  sol(:,:,i) = y(:,i:npde:end);
end
varargout{1} = sol;
if hasEvents % [sol,t,sole,te,ie] = pdepe(...)
  varargout{2} = t(:);
  sole = zeros(length(te),nx,npde);
  for i = 1:npde
    sole(:,:,i) = ye(:,i:npde:end);
  end
  varargout{3} = sole;
  varargout{4} = te(:);
  varargout{5} = ie(:);
end

%---------------------------------------------------------------------------
% Nested functions
%---------------------------------------------------------------------------

  function dudt = pdeodes(tnow,y)
  %PDEODES  Assemble the difference equations and evaluate the time derivative
  %   for the ODE system.

    u = reshape(y,npde,nx);
    up = zeros(npde,nx);
    [U,Ux] = pdentrp_AGHF(singular,m,xmesh(1),u(:,1),xmesh(2),u(:,2),xi(1));
    [cL,fL,sL] = mypdexpde_generated(xi(1),tnow,U,Ux,params);

    %  Evaluate the boundary conditions
    [pL,qL,pR,qR] = mypdexbc_generated(xmesh(1),u(:,1),xmesh(nx),u(:,nx),tnow,X);

    %  Left boundary
    if singular
      denom = cL;
      denom(denom == 0) = 1;
      up(:,1) = (sL + (m+1) * fL / xi(1)) ./ denom;
    else
      up(:,1) = pL;
      idx = (qL ~= 0);
      denom = (qL(idx)/xmesh(1)^m) .* (zxmp1(1)*cL(idx));
      denom(denom == 0) = 1;
      up(idx,1) = ( pL(idx) + (qL(idx)/xmesh(1)^m) .* ...
                    (xim(1)*fL(idx) + zxmp1(1)*sL(idx))) ./ denom;
    end
    %  Interior points
    for ii = 2:nx-1
      [U,Ux] = pdentrp_AGHF(singular,m,xmesh(ii),u(:,ii),xmesh(ii+1),u(:,ii+1),xi(ii));
      [cR,fR,sR] = mypdexpde_generated(xi(ii),tnow,U,Ux,params);

      denom = zxmp1(ii) * cR + xzmp1(ii-1) * cL;
      denom(denom == 0) = 1;
      up(:,ii) = ((xim(ii) * fR - xim(ii-1) * fL) + ...
                  (zxmp1(ii) * sR + xzmp1(ii-1) * sL)) ./ denom;

      cL = cR;
      fL = fR;
      sL = sR;
    end
    %  Right boundary
    up(:,nx) = pR;
    idx = (qR ~= 0);
    denom = -(qR(idx)/xmesh(nx)^m) .* (xzmp1(nx-1)*cL(idx));
    denom(denom == 0) = 1;
    up(idx,nx) = ( pR(idx) + (qR(idx)/xmesh(nx)^m) .* ...
                   (xim(nx-1)*fL(idx) - xzmp1(nx-1)*sL(idx))) ./ denom;

    dudt = up(:);
  end  % pdeodes

% --------------------------------------------------------------------------

end  % pdepe

