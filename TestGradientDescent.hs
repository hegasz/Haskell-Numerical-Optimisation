module TestGradientDescent where

import Test.HUnit
import GradientDescent
import Vectors
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

-- satisfiesArmijo f dfdx [3,4] alpha sigma should
-- be true only if alpha<=1-sigma
testArmijo1 = TestCase (assertEqual "Armijo 1" True (satisfiesArmijo f dfdx [3,4] 0.79 0.2))
testArmijo2 = TestCase (assertEqual "Armijo 2" False (satisfiesArmijo f dfdx [3,4] 0.81 0.2))
testArmijo3 = TestCase (assertEqual "Armijo 3" True (satisfiesArmijo f dfdx [3,4] 0.09 0.9))
testArmijo4 = TestCase (assertEqual "Armijo 4" False (satisfiesArmijo f dfdx [3,4] 0.11 0.9))
armijoTests = TestList [testArmijo1,
                        testArmijo2,
                        testArmijo3,
                        testArmijo4]

-- nextStepGuess dfdx [3,4] [5,6] prevStep should
-- return 272/464 * prevStep
testNextStepGuess1 = TestCase (assertBool "NextStepGuess 1"
    (abs (nextStepGuess dfdx [3,4] [5,6] 7 - (272/464) * 7) < 1e-15))
testNextStepGuess2 = TestCase (assertBool "NextStepGuess 2"
    (abs (nextStepGuess dfdx [3,4] [5,6] 1 - (272/464)) < 1e-15))
nextStepGuessTests = TestList [testNextStepGuess1,
                               testNextStepGuess2]

-- test: backtrack f dfdx [3,4] sigma initialGuess should
-- halve initialGuess until it is less than 1-sigma
testBacktrack1 = TestCase (assertBool "Backtrack 1"
    (abs (backtrack f dfdx [3,4] 0.2 7 - 0.4375) < 1e-15))
testBacktrack2 = TestCase (assertBool "Backtrack 2"
    (abs (backtrack f dfdx [3,4] 0.5 0.49 - 0.49) < 1e-15))
testBacktrack3 = TestCase (assertBool "Backtrack 3"
    (abs (backtrack f dfdx [3,4] 0.9 2 - 0.0625) < 1e-15))
backtrackTests = TestList [testBacktrack1, testBacktrack2, testBacktrack3]

-- test: nextStep f dfdx [3,4] [2,1] 3.5 sigma should
-- return (3-4alpha, 4-16alpha) where alpha is the first halving of
-- 3.5*104/272 = 1.3382... (the initial step guess) that falls under 1-sigma
-- so nextStep f dfdx [3,4] [2,1] 3.5 0.9 should give x_n+1 = [2.6654411764705883,2.6617647058823533]
testNextStep1 = TestCase (assertBool "Next Step 1"
    ((abs (result!!0 - 2.6654411764705883) < 1e-15) && (abs (result!!1 - 2.6617647058823533) < 1e-15)))
        where result = first (nextStep f dfdx [3,4] [2,1] 3.5 0.9)
              first (x,_,_) = x


testDescent1 = TestCase (assertBool "Descent 1"
    ((abs (result!!0 - 1) < 1e-5) && (abs (result!!1 + 4) < 1e-5)))
        where result = fromJust (gradientDescent f dfdx [3,4] 2 0.001 1000000 1e-7)
testDescentBanana = TestCase (assertBool "Rosenbrock banana function descent"
    ((abs (result!!0 - 1) < 0.01) && (abs (result!!1 - 1) < 0.01)))
        where result = fromJust (gradientDescent bananaF bananaDfdx [3,4] 2 0.1 1000000 1e-14)
descentTests = TestList [testDescent1, testDescentBanana]

tests = TestList [armijoTests, nextStepGuessTests, backtrackTests, testNextStep1, descentTests]

main :: IO Counts
main = runTestTT tests
