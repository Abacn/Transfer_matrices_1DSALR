% Second order differential
function [dfddx, dfdx]=deriv2(f,x,h,method)
  if nargin<3
      absx = abs(x);
      if absx<1e-10
          h = 1e-14;
      elseif absx<10
          h = absx*1e-4;
      else
          h = 1e-3;
      end
  end
  if nargin<4
      method = 'mid';
  end
  fneg2 = f(x-2*h);
  fneg = f(x-h);
  fme = f(x);
  fpos = f(x+h);
  fpos2 = f(x+2*h);
  
  dfddx = (16*(fneg+fpos)-30*fme-(fneg2+fpos2))/(12*h^2);
  dfdx = (fneg2-fpos2+8*(fpos-fneg))/(12*h);

end