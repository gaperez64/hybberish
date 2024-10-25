using TaylorModels

# Example 3.3.6
domain = IntervalBox([-1..1,-0.5..0.5,0..0.02])
zd = zero(domain)  # creates a hyperrectangle of 0's with same dimension
x, y, t = set_variables("x y t", order=3)
# focusing on the second dimension of the RHS of dynamics
f(a) = -a^2
# focusing on the second dimension of polynomial approx
p(a,b,c) = a + c + b * c
tm3 = TaylorModelN(p(x,y,t), -0.1..0.1, zd, domain)
# FIXME: this is needed as a hack to allow TM multiplication
# we update the order two twice the original one (or more)
x, y, t = set_variables("x y t", order=8)
ftm3 = f(tm3)
# I happen to know the total degree of the result is 3
# so we just need to get rid of that one, i.e. it is p_e
(polynomial(ftm3)[3](-1..1,-0.5..0.5,0..0.02) + remainder(ftm3)) * (0..0.02)
