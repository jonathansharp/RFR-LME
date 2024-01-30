function zz = regrid(x,y,z,xq,yq)

xx = repmat(x',length(y),1);
yy = repmat(y,1,length(x));

xxq = repmat(xq',length(yq),1);
yyq = repmat(yq,1,length(xq));

idx = ~isnan(z);
interpolant = scatteredInterpolant(xx(idx),yy(idx),z(idx));
zz = interpolant(xxq,yyq)';