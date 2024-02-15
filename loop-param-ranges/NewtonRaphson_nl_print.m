function [v1 , no_itr, norm1] = NewtonRaphson_nl_print(v,fn,jacob_fn,no_itr,error,lam)
    v1 = v;
    vsize = size(v);
    disp(vsize)
    fnv1 = feval(fn,v1);
    II = eye(vsize(1));
    i = 0;
    fprintf('      Iteration|    Ilmin     |     Ilmax     |     Rlmin     |     Rlmax     |     nmin     |     nmax     | Error      | \n')
    while true
        norm1 = norm(fnv1);
        fprintf('%10d     | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4d |\n',i,v1(1),v1(2),v1(3),v1(4),v1(5),v1(6),norm1)
        jacob_fnv1 = feval(jacob_fn,v1);
        H = (jacob_fnv1'*jacob_fnv1+lam*II)\(jacob_fnv1'*fnv1);
        v1 = v1 - H;
        fnv1 = feval(fn,v1);
        i = i + 1 ;
        norm1 = norm(fnv1);
        if i > no_itr && norm1 < error, break , end %
        %if norm(fnv1) < error , break , end
        
    end
end
