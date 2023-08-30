x = I(3)
hx = homomatrix(x)
@test dehomomatrix(hx) == x
x = rand(1:3,3,4)
hx = homomatrix(x)
@test dehomomatrix(hx) == x
x = rand(3,3)
hx = homomatrix(x)
@test dehomomatrix(hx) == x

x = I(3)
hx = homovector(x)
@test dehomovector(hx) == x
x = rand(1:3,3,6)
hx = homovector(x)
@test dehomovector(hx) == x
x = rand(3,6)
hx = homovector(x)
@test dehomovector(hx) == x
x = rand(3)
hx = homovector(x)
@test dehomovector(hx) == x

