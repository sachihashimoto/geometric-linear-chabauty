
#Sage code for Example 5.8
#Written by Sachi Hashimoto and Pim Spelier
#The goal of the sage code is to rule out the extra points R1 and R2 which are found by the Coleman code by running Chabauty-Coleman

R.<x>=QQ['x']
f=4*x^7 - 12*x^6 + 20*x^5 - 20*x^4 + 12*x^3 - 4*x^2 + 1  #eqn

f4 = f/4 #make monic, since sage only does coleman ints on monic polys
p=3 #prime
H=HyperellipticCurve(f/4)
Hpad=H.change_ring(Qp(p,15))
H.change_ring(GF(p)).points() #residue disks
R = Qp(3,15)

#compute M0bar
G1 = Hpad(0 , 1/2 ,1) #assume that generator for MW group is inf - P, then P is this point
G2 = Hpad(1,  1/2, 1)


logG1=Hpad.coleman_integrals_on_basis(Hpad(0,1,0), G1) #compute the integral from inf to P over omega0, omega1
logG2=Hpad.coleman_integrals_on_basis(Hpad(0,1,0), G2)
print("logs of G1, G2")
print(logG1)
print(logG2)

orderG1 = [-6,-4] #generators for the kernel ofreduction
orderG2 = [-10,11]
logG1tilde = [orderG1[0] * logG1[i]/p+ orderG1[1]*logG2[i]/p for i in range(3)]
logG2tilde = [orderG2[0] * logG1[i]/p+ orderG2[1]*logG2[i]/p for i in range(3)]
print("logG1tilde")
print(logG1tilde)
print("logG2tilde")
print(logG2tilde)

#Compute DQ
#let xcoord be the x-coordinate of the residue disk Q you are considering
xcoord = 2
Qlift = Hpad.lift_x(xcoord) #lift this point
C = Hpad.coleman_integrals_on_basis(Hpad(0,1,0), Qlift)
print("this is log v")
print(C[0]/p) #this gives you constants for Qmu - inf
print(C[1]/p)
print(C[2]/p)

xres, yres = Hpad.local_coordinates_at_nonweierstrass(Hpad.lift_x(xcoord)) #local coords at Q the residue disk
#this compute tiny integrals for DQ and also for Qmu
I0 = (xres.derivative() / (2*yres)).integral() #this sets up the integral of d2x(t)/2y(t) dt
I1=(xres *(xres.derivative() )/ (2*yres)).integral() # this sets up x(t) *dx(t) / 2y(t) dt
I2=(xres^2 *(xres.derivative() )/ (2*yres)).integral() # this sets up x(t) *dx(t) / 2y(t) dt
L= I0.parent()
t=L.gen()
integral0 = I0(t=p*t).add_bigoh(2)/p
integral1 = I1(t=p*t).add_bigoh(2)/p
integral2 = I2(t=p*t).add_bigoh(2)/p
print("this is DQ")
print(integral0)
print(integral1)
print(integral2)

R1 = Hpad(R(3757457), R(52603757)/2) #values from Magma Coleman code for R1, R2 the extra p-adic points found in Chabauty-Coleman
R2 =  Hpad(R(3757457), R(-52603757)/2)

logR1 = Hpad.coleman_integrals_on_basis(Hpad(0,1,0), R1)
logR2 = Hpad.coleman_integrals_on_basis(Hpad(0,1,0), R1)
logR1 = (logR1/p)[0:3]
logR2 = (logR2/p)[0:3]
print("this is logR1, logR2")
print(logR1)
print(logR2)
v = vector(logR1)
logGimat = Matrix([logG1tilde, logG2tilde])
print(logGimat.solve_left(v))
