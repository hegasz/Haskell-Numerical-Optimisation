module GradientDescent(
    satisfiesArmijo,
    nextStepGuess,
    backtrack,
    nextStep,
    gradSmallEnough,
    first,
    gradientDescent) where


import Vectors
import Data.Maybe


-- returns true iff this step size satisfies armijo rule
satisfiesArmijo :: (Vec -> Double) ->
                   (Vec -> Vec) ->
                   Vec ->
                   Double ->
                   Double ->
                   Bool
satisfiesArmijo f dfdx x_n alpha sigma = f(fromJust(x_n-^(alpha*^g_n)))<=f x_n-alpha*sigma*fromJust(dot g_n g_n)
    where g_n = dfdx x_n

-- returns an initial guess for the step size for the next x value
nextStepGuess :: (Vec -> Vec) ->
                 Vec ->
                 Vec ->
                 Double ->
                 Double
nextStepGuess dfdx x_n_min1 x_n prevStep = prevStep * fromJust(dot g_n_min1 g_n_min1) / fromJust(dot g_n g_n)
    where g_n = dfdx x_n
          g_n_min1 = dfdx x_n_min1

-- returns a backtracked acceptable step size for x_n
backtrack :: (Vec -> Double) ->
             (Vec -> Vec) ->
             Vec ->
             Double ->
             Double ->
             Double
backtrack f dfdx x_n sigma = head .
            dropWhile (\alpha -> not (satisfiesArmijo f dfdx x_n alpha sigma)) .
            iterate (/2)


-- DEBUGGING: try replacing first line with:
-- *function declaration* ... = (initialGuess, alpha){-x_n -^ (alpha *^ (dfdx x_n))-}
-- calculates x_n+1 from x_n, using x_n-1 and alpha_n-1 to guide first step size guess
-- returns (x_n+1, x_n, alpha_n)
nextStep :: (Vec -> Double) ->
            (Vec -> Vec) ->
            Vec ->
            Vec ->
            Double ->
            Double ->
            (Vec, Vec, Double)
nextStep f dfdx x_n x_n_min1 prevStep sigma = (fromJust(x_n -^ (alpha *^ dfdx x_n)), x_n, alpha)
    where alpha = backtrack f dfdx x_n sigma initialGuess
          initialGuess = nextStepGuess dfdx x_n_min1 x_n prevStep

-- returns whether gradient is within tolerance (absolutely or relatively to start)
gradSmallEnough :: (Vec -> Vec) ->
                   Vec ->
                   Double ->
                   Vec ->
                   Bool
gradSmallEnough dfdx firstGrad tol x = fromJust(dot grad grad) < tol * (1 + fromJust(dot firstGrad firstGrad))
    where grad = dfdx x

first :: (a,b,c) -> a
first (x,_,_) = x


-- Takes in a function, its derivative, an initial vector, a first step size,
-- a sigma value for the Armijo rule, a max number of iterations allowed, and
-- a tolerance for an acceptable final gradient - and returns a vector that
-- is at a stationary point of the function to within this tolerance.
-- pre: alpha_0 is large enough
gradientDescent :: (Vec -> Double) ->
                   (Vec -> Vec) ->
                   Vec ->
                   Double ->
                   Double ->
                   Int ->
                   Double ->
                   Maybe Vec
gradientDescent f dfdx x_0 alpha_0 sigma n tol | alpha_0 <= 0 || sigma <= 0 || sigma >= 1 || n <= 0 || tol <= 0 = Nothing
                                               | otherwise = Just(
    (first .
    head .
    dropWhile (\(x,_,_) -> not (gradSmallEnough dfdx (dfdx x_0) tol x)) .
    take n . -- only allowed max n iterations
    iterate nextStepFunc) (x_1, x_0, alpha_0))
      where nextStepFunc = \(x_n, x_n_min1, prevStep) -> nextStep f dfdx x_n x_n_min1 prevStep sigma
            x_1 = fromJust(x_0 -^ (backtrack f dfdx x_0 sigma alpha_0 *^ dfdx x_0))

-- nextStep has constant parameters f, dfdx and sigma
-- but churns through new values of x_n, x_n_min1 and prevStep