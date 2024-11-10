
function [s_rule_val,v] = f_gss(w,sval,sgrid,r,beta,V,snum)

eps = 1e-8;

aa = min(sgrid);
bb = (1+r)*sval + w;

rr = (3-sqrt(5))/2;

cc = (1-rr)*aa + rr*bb;
dd = rr*aa + (1-rr)*bb;

fc = find_value(cc,w,sval,sgrid,r,eps,beta,V,snum);
fd = find_value(dd,w,sval,sgrid,r,eps,beta,V,snum);

while abs(dd - cc) >= eps
    
    if fc >= fd
        aa = cc;
        cc = dd;
        fc = fd;
        dd = rr*aa + (1-rr)*bb;
        fd = find_value(dd,w,sval,sgrid,r,eps,beta,V,snum);
    else
        bb = dd;
        dd = cc;
        fd = fc;
        cc = (1-rr)*aa + rr*bb;
        fc = find_value(cc,w,sval,sgrid,r,eps,beta,V,snum);
    end
end

s_rule_val = (aa+bb)/2;
v = -find_value(s_rule_val,w,sval,sgrid,r,eps,beta,V,snum);

end

function [f] = find_value(x,w,sval,sgrid,r,eps,beta,V,snum)

consump = (1+r)*sval + w - x;
if consump <= 0 
    consump=eps;
    x = (1+r)*sval + w - eps;
end
util = log(consump);

    %%%% importantly, interpolate future value associated with s_rule and the grid sgrid: 
    
    unit = ones(size(sgrid));
    resource = (x+0.00000001)*unit - sgrid;
    sindex = sum(resource > 0);
    if sindex <= 0
      sindex = 1;
    end    
    
    if(sindex < snum)
      weight = sgrid(sindex+1) - x;
      weight = weight/(sgrid(sindex+1) - sgrid(sindex));
    else
      weight = 0.0;
      sindex = sindex - 1;
    end    
   
    valf = V(sindex)*weight + V(sindex+1)*(1.0 - weight);
    
fv = util + beta*valf;
f = -fv;

end

