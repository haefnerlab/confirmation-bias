function pyes = compute_sigmoid(x,k,m,a)
z = (x-m).*k;
pyes = a+( (1.0-a)./ (1+ exp( - z ) ));
end