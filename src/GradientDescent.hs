module GradientDescent(
    satisfiesArmijo,
    nextStepGuess,
    backtrack,
    nextStep,
    gradSmallEnough,
    first,
    gradientDescent) where


import VectorsInternal
import Data.Maybe
import OptimisationUtils



-- Calculates x_n+1 from x_n, using x_n-1 and alpha_n-1 to guide first step size guess
-- Returns (x_n+1, x_n, alpha_n)
-- nextStep has constant parameters f, dfdx and sigma
-- but churns through new values of x_n, x_n_min1 and prevStep
nextStep :: (Vec -> Double) ->
            (Vec -> Vec) ->
            Vec ->
            Vec ->
            Double ->
            Double ->
            (Vec, Vec, Double)
nextStep f dfdx x_n x_n_min1 prevStep sigma = (x_n -^ (alpha *^ dfdx x_n), x_n, alpha)
    where alpha = backtrack f dfdx x_n sigma initialGuess
          initialGuess = nextStepGuess dfdx x_n_min1 x_n prevStep


-- Takes in a function, its derivative, an initial vector, a first step size,
-- a sigma value for the Armijo rule, a max number of iterations allowed, and
-- a tolerance for an acceptable final gradient - and returns a vector that
-- is at a stationary point of the function to within this tolerance.
-- Uses classic gradient descent.
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
                                               | otherwise =
    let result = dropWhile (\(x,_,_) -> not (gradSmallEnough dfdx (dfdx x_0) tol x)) .
                 take n . -- only allowed max n iterations
                 iterate nextStepFunc $ (x_1, x_0, alpha_0)
    in case result of
        [] -> Nothing -- not enough iterations to get within tolerance
        xs -> Just ((\(x,_,_) -> x) . head $ result)
      where nextStepFunc = \(x_n, x_n_min1, prevStep) -> nextStep f dfdx x_n x_n_min1 prevStep sigma
            x_1 = x_0 -^ (backtrack f dfdx x_0 sigma alpha_0 *^ dfdx x_0)