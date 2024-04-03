function OBJ = OBJfcn_TEST(X)

o1    = X(1)+X(2)^2-X(3)*X(4)+X(5)*X(6)^2-X(7)*X(8);
o2    = X(9)-X(10)^2+X(11)*X(12)-X(13)*sin(X(14)*X(15)) + X(16)*cos(X(17)*X(18));

OBJ = [o1, o2];

end


