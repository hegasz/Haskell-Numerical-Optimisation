module NewtonOptimisationSpec where
import SpecHelper

import NewtonOptimisation
import VectorsInternal
import Data.Maybe

f :: Vec -> Double
f xs = (x-1)**2 + (y+4)**2 + 3
    where x = xs!!0
          y = xs!!1

dfdx :: Vec -> Vec
dfdx xs = [(x-1)*2, (y+4)*2]
    where x = xs!!0
          y = xs!!1

hess :: Vec -> Matrix
hess xs = [[2,0],[0,2]]
    where x = xs!!0
          y = xs!!1

bananaF :: Vec -> Double
bananaF xs = (1-x)**2 + 100*((y-(x**2))**2)
    where x = xs!!0
          y = xs!!1

bananaDfdx :: Vec -> Vec
bananaDfdx xs = [-2*(1-x) - 400*(y-(x**2))*x, 200*(y-(x**2))]
    where x = xs!!0
          y = xs!!1

bananaHess :: Vec -> Matrix
bananaHess xs = [[2-400*(y-(3*x**2)),-400*x],[-400*x,200]]
    where x = xs!!0
          y = xs!!1

f2 :: Vec -> Double
f2 xs = (x**3) - (2*x*y) - (y**6)
    where x = xs!!0
          y = xs!!1

dfdx2 :: Vec -> Vec
dfdx2 xs = [(3*(x**2)) - (2*y), -(2*x)-(6*(y**5))]
    where x = xs!!0
          y = xs!!1


hess2 :: Vec -> Matrix
hess2 xs = [[6*x, -2], [-2, -30*(y**4)]]
    where x = xs!!0
          y = xs!!1

spec :: Spec
spec = describe "All" $ do
        describe "NewtonOptimise" $ do
            it "works on a standard polynomial function" $ do
                let result = fromJust $ newtonOptimisation f2 dfdx2 hess2 [-0.7,0.7] 0.01 100000 0.000001
                (abs (result!!0 + 0.706575) < 1e-5) && (abs (result!!1 - 0.748872) < 1e-5) `shouldBe` True
            it "works on a quadratic function" $ do
                let result = fromJust $ newtonOptimisation f dfdx hess [100,214] 0.01 20 0.000000001
                (abs (result!!0 - 1.003021240234375) < 1e-10) && (abs (result!!1 + 3.99334716796875) < 1e-10) `shouldBe` True
            it "works on Rosenbrock banana function" $ do
                let result = fromJust $ newtonOptimisation bananaF bananaDfdx bananaHess [2,2] 0.01 20000000000 0.000000000001
                (abs (result!!0 - 1.0003898717649136) < 1e-10) && (abs (result!!1 - 1.0007778479460452) < 1e-10) `shouldBe` True

