module NewtonOptimisation where

import Vectors

f :: Vec -> Double
f xs = (x**3) - (2*x*y) - (y**6)
    where x = xs!!0
          y = xs!!1

dfdx :: Vec -> Vec
dfdx xs = [(3*(x**2)) - (2*y), -(2*x)-(6*(y**5))]
    where x = xs!!0
          y = xs!!1

hess :: Vec -> Matrix
hess xs = [[6*x, -2], [-2, -30*(y**4)]]
    where x = xs!!0
          y = xs!!1

