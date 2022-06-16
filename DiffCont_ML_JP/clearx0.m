function x0po = clearx0(x0po_iter, mu, OPTIONS)

len = size(x0po_iter, 1) ;

x0po = zeros(1,7) ;
j=1 ;
for k = 1:len
    x0 = x0po_iter(k,1:6) ;
    tf = x0po_iter(k, end) ;
    [x,t,phi_t1,PHI] = stateTransMat3BP3d(x0, tf, mu, OPTIONS) ;
    eigens = eig(phi_t1) 
    eigens_ = round(eigens, 2) ;
    idx1 = imag(eigens_) == 0 ;
    eigens = round(eigens(idx1)) ;
    idx2 = real(eigens) == 1 ; 
    if sum( idx2) >= 2
        x0po(j,:) = x0po_iter(k,:) ;
        j=j+1 ;
    end
end


end