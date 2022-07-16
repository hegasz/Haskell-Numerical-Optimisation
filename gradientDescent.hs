
f :: (Double, Double) -> Double
f (x,y) = (x-1)**2 + (y+4)**2 + 3

dfdx :: (Double, Double) -> (Double, Double)
dfdx (x,y) = ((x-1)*2, (y+4)*2)

norm :: (Double, Double) -> (Double, Double) -> Double
norm (a,b) (c,d) = a*c + b*d

(-^) :: (Double, Double) -> (Double, Double) -> (Double, Double)
(-^) (a,b) (c,d) = (a-c, b-d)


(*^) :: Double -> (Double, Double) -> (Double, Double)
(*^) a (c,d) = (a*c, a*d)

-- test: satisfiesArmijo f dfdx (3,4) alpha sigma should
-- be true only if alpha<=1-sigma
-- returns true iff this step size satisfies armijo rule
satisfiesArmijo :: ((Double, Double) -> Double) ->
                   ((Double, Double) -> (Double, Double)) ->
                   (Double, Double) ->
                   Double ->
                   Double ->
                   Bool
satisfiesArmijo f dfdx x_n alpha sigma = f(x_n-^(alpha*^g_n))<=f(x_n)-alpha*sigma*(norm g_n g_n)
    where g_n = dfdx(x_n)

-- test: nextStepGuess dfdx (3,4) (5,6) prevStep should
-- return 272/464 * prevStep
-- returns an initial guess for the step size for the next x value
nextStepGuess :: ((Double, Double) -> (Double, Double)) ->
                 (Double, Double) ->
                 (Double, Double) ->
                 Double ->
                 Double
nextStepGuess dfdx x_n_min1 x_n prevStep = prevStep * (norm g_n_min1 g_n_min1) / (norm g_n g_n)
    where g_n = dfdx x_n
          g_n_min1 = dfdx x_n_min1

-- test: backtrack f dfdx (3,4) sigma initialGuess should
-- halve initialGuess until it is less than 1-sigma
-- returns a backtracked acceptable step size for x_n
backtrack :: ((Double, Double) -> Double) ->
             ((Double, Double) -> (Double, Double)) ->
             (Double, Double) ->
             Double ->
             Double ->
             Double
backtrack f dfdx x_n sigma initialGuess =
            (head .
            dropWhile (\alpha -> not (satisfiesArmijo f dfdx x_n alpha sigma)) .
            iterate (/2)) initialGuess

-- test: nextStep f dfdx (3,4) (2,1) 3.5 sigma should
-- return (3-4alpha, 4-16alpha) where alpha is the first halving of
-- 3.5*104/272 = 1.3382... (the initial step guess) that falls under 1-sigma
-- so nextStep f dfdx (3,4) (2,1) 3.5 0.9 should be (2.6654411764705883,2.6617647058823533)
-- DEBUGGING: try replacing first line with:
-- *function declaration* ... = (initialGuess, alpha){-x_n -^ (alpha *^ (dfdx x_n))-}
-- calculates x_n+1 from x_n, using x_n-1 and alpha_n-1 to guide first step size guess
-- returns (x_n+1, x_n, alpha_n)
nextStep :: ((Double, Double) -> Double) ->
            ((Double, Double) -> (Double, Double)) ->
            (Double, Double) ->
            (Double, Double) ->
            Double ->
            Double ->
            ((Double, Double), (Double, Double), Double)
nextStep f dfdx x_n x_n_min1 prevStep sigma = (x_n -^ (alpha *^ (dfdx x_n)), x_n, alpha)
    where alpha = backtrack f dfdx x_n sigma initialGuess
          initialGuess = nextStepGuess dfdx x_n_min1 x_n prevStep

gradSmallEnough :: ((Double, Double) -> (Double, Double)) ->
                   (Double, Double) ->
                   Double ->
                   (Double, Double) ->
                   Bool
gradSmallEnough dfdx firstGrad tol x = norm grad grad < tol * (1 + (norm firstGrad firstGrad))
    where grad = dfdx x

first :: (a,b,c) -> a
first (x,_,_) = x


-- Takes in a function, its derivative, an initial vector, a first step size,
-- a sigma value for the Armijo rule, a max number of iterations allowed, and
-- a tolerance for an acceptable final gradient - and returns a vector that
-- is at a stationary point of the function to within this tolerance.
gradientDescent :: ((Double, Double) -> Double) ->
                   ((Double, Double) -> (Double, Double)) ->
                   (Double, Double) ->
                   Double ->
                   Double ->
                   Int ->
                   Double ->
                   (Double, Double)
gradientDescent f dfdx x_0 alpha_0 sigma n tol =
    (first .
    head .
    dropWhile (\(x,_,_) -> not (gradSmallEnough dfdx (dfdx x_0) tol x)) .
    take n . -- only allowed max n iterations
    iterate nextStepFunc) (x_1, x_0, alpha_0)
      where nextStepFunc = (\(x_n, x_n_min1, prevStep) -> nextStep f dfdx x_n x_n_min1 prevStep sigma)
            x_1 = x_0 -^ ((backtrack f dfdx x_0 sigma alpha_0) *^ (dfdx x_0))

-- nextStep has constant parameters f, dfdx and sigma
-- but churns through new values of x_n, x_n_min1 and prevStep