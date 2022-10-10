module NewtonOptimisation where

import VectorsInternal
import LinearEquationSolver (solve)
import Data.Either
import OptimisationUtils


-- Performs one step of Newton optimisation
nextStep :: (Vec -> Entry) ->
            (Vec -> Vec) ->
            (Vec -> Matrix) ->
            Double ->
            Vec ->
            Vec
nextStep f dfdx hess sigma x_n = x_n +^ (alpha *^ dir_n)
    where dir_n = fromRight ([] :: Vec) (solve (hess x_n) ((-1) *^ (dfdx x_n)))
          alpha = backtrack f dfdx x_n sigma 1

-- Takes in a function, its derivative, its Hessian, an initial vector,
-- a sigma value for the Armijo rule, a max number of iterations allowed, and
-- a tolerance for an acceptable final gradient - and returns a vector that
-- is at a stationary point of the function to within this tolerance.
-- Uses Newton's method for optimisation.
newtonOptimisation :: (Vec -> Entry) ->
                   (Vec -> Vec) ->
                   (Vec -> Matrix) ->
                   Vec ->
                   Double ->
                   Int ->
                   Double ->
                   Maybe Vec
newtonOptimisation f dfdx hess x_0 sigma n tol | sigma <= 0 || sigma >= 1 || n <= 0 || tol <= 0 = Nothing
                                               | otherwise =
    let result = dropWhile (\x -> not (gradSmallEnough dfdx (dfdx x_0) tol x)) .
             take n . -- only allowed max n iterations
             iterate (nextStep f dfdx hess sigma) $ x_0                           
    in case result of
        [] -> Nothing -- not enough iterations to get within tolerance
        xs -> Just (head result)
    
