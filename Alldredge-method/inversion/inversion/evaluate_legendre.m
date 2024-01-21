function legendre_polynomial = evaluate_legendre(degree, order, param)

legendre_polynomial = legendre(degree, cos(param), 'sch');
%legendre_polynomial = legendre_polynomial(order);
if order == 1
    legendre_polynomial = legendre_polynomial(2);
else
    legendre_polynomial = legendre_polynomial(order+1);

end