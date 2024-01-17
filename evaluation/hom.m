function [HOM]=hom(a,b,na,nb)

HOM=(a/na-b/nb)^2/(a+b);
end