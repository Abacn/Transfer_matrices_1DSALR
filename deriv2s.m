% Second order differential
function [dfddx, dfddy, dfdxy, dfdx]=deriv2s(f,x,y,h0,method)
  if nargin<4
      absx = abs(x);
      if absx<1e-10
          h = 1e-14;
      elseif absx<10
          h = absx*1e-4;
      else
          h = 1e-3;
      end
      absy = abs(y);
      if absy<1e-10
          k = 1e-14;
      elseif absy<10
          k = absy*1e-4;
      else
          k = 1e-3;
      end
  else
      h = h0(1);
      k = h0(2);
  end
  if nargin<5
      method = 'mid';
  end
  fneg2x = f(x-2*h, y);
  fnegx = f(x-h, y);
  fme = f(x, y);
  fposx = f(x+h, y);
  fpos2x = f(x+2*h, y);
  fneg2y = f(x, y-2*k);
  fnegy = f(x, y-k);
  fposy = f(x, y+k);
  fpos2y = f(x, y+2*k);
  
  fposxy = f(x+h, y+k);
  fnegxy = f(x-h, y-k);
  
  dfddx = (16*(fnegx+fposx)-30*fme-(fneg2x+fpos2x))/(12*h^2);
  dfddy = (16*(fnegy+fposy)-30*fme-(fneg2y+fpos2y))/(12*k^2);
  dfdxy = (fposxy-fposx-fposy+2*fme-fnegx-fnegy+fnegxy)/(2*h*k);
  dfdx = (fneg2x-fpos2x+8*(fposx-fnegx))/(12*h);

end