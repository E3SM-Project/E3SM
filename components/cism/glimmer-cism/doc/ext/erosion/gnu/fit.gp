set angles deg

c0  = 3.75
phi = 32

t0(N) = 3.75+N*0.625

f1(x,y) =  B1* x**n1 *y**m1
f2(x,y) = B2*(x-t0(y))**n2 *y**m2
f3(x,y) = B3*(x-t0(y))**1*y**(-2)

fit f1(x,y) "stress-strain.data" u 2:3:4:(1) via B1,n1,m1
fit f2(x,y) "stress-strain.data" u 2:3:4:(1) via B2,n2,m2
fit f3(x,y) "stress-strain.data" u 2:3:4:(1) via B3

i1(N,e) = (e/(B1*N**m1))**(1/n1)
i2(N,e) = (e/(B2*N**m2))**(1/n2)+t0(N)
i3(N,e) = (e/(B3*N**(-2)))+t0(N)
