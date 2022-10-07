module GradientDescentSpec where
import SpecHelper

import GradientDescent
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

bananaF :: Vec -> Double
bananaF xs = (1-x)**2 + 100*((y-(x**2))**2)
    where x = xs!!0
          y = xs!!1

bananaDfdx :: Vec -> Vec
bananaDfdx xs = [-2*(1-x) - 400*(y-(x**2))*x, 200*(y-(x**2))]
    where x = xs!!0
          y = xs!!1

spec :: Spec
spec = describe "All" $ do
        describe "Armijo" $ do
            -- satisfiesArmijo f dfdx [3,4] alpha sigma should
            -- be true only if alpha<=1-sigma
            it "Armijo 1" $ do
                satisfiesArmijo f dfdx [3,4] 0.79 0.2 `shouldBe` True
            it "Armijo 2" $ do
                satisfiesArmijo f dfdx [3,4] 0.81 0.2 `shouldBe` False
            it "Armijo 3" $ do
                satisfiesArmijo f dfdx [3,4] 0.09 0.9 `shouldBe` True
            it "Armijo 4" $ do
                satisfiesArmijo f dfdx [3,4] 0.11 0.9 `shouldBe` False

        describe "NextStepGuess" $ do
            -- nextStepGuess dfdx [3,4] [5,6] prevStep should
            -- return 272/464 * prevStep
            it "NextStepGuess 1" $ do
                abs (nextStepGuess dfdx [3,4] [5,6] 7 - (272/464) * 7) < 1e-15 `shouldBe` True
            it "NextStepGuess 2" $ do
                abs (nextStepGuess dfdx [3,4] [5,6] 1 - (272/464)) < 1e-15 `shouldBe` True
        
        describe "Backtrack" $ do
            -- test: backtrack f dfdx [3,4] sigma initialGuess should
            -- halve initialGuess until it is less than 1-sigma
            it "Backtrack 1" $ do
                abs (backtrack f dfdx [3,4] 0.2 7 - 0.4375) < 1e-15 `shouldBe` True
            it "Backtrack 2" $ do
                abs (backtrack f dfdx [3,4] 0.5 0.49 - 0.49) < 1e-15 `shouldBe` True
            it "Backtrack 3" $ do
                abs (backtrack f dfdx [3,4] 0.9 2 - 0.0625) < 1e-15 `shouldBe` True

        describe "NextStep" $ do
            -- test: nextStep f dfdx [3,4] [2,1] 3.5 sigma should
            -- return (3-4alpha, 4-16alpha) where alpha is the first halving of
            -- 3.5*104/272 = 1.3382... (the initial step guess) that falls under 1-sigma
            -- so nextStep f dfdx [3,4] [2,1] 3.5 0.9 should give x_n+1 = [2.6654411764705883,2.6617647058823533]
            it "Next Step 1" $ do
                let result = (\(x,_,_)->x) (nextStep f dfdx [3,4] [2,1] 3.5 0.9)
                (abs (result!!0 - 2.6654411764705883) < 1e-15) && (abs (result!!1 - 2.6617647058823533) < 1e-15) `shouldBe` True

        describe "GradientDescent" $ do
            it "Descent 1" $ do
                let result = fromJust (gradientDescent f dfdx [3,4] 2 0.001 1000000 1e-7)
                (abs (result!!0 - 1) < 1e-5) && (abs (result!!1 + 4) < 1e-5) `shouldBe` True
            it "Rosenbrock banana function descent" $ do
                let result = fromJust (gradientDescent bananaF bananaDfdx [3,4] 2 0.1 1000000 1e-14)
                (abs (result!!0 - 1) < 0.01) && (abs (result!!1 - 1) < 0.01) `shouldBe` True


main :: IO ()
main = hspec spec