module OptimisationUtils where


-- Returns true iff this step size satisfies Armijo rule
-- with sigma value given.
satisfiesArmijo :: (Vec -> Entry) ->
                   (Vec -> Vec) ->
                   Vec ->
                   Double ->
                   Double ->
                   Bool
satisfiesArmijo f dfdx x_n alpha sigma = f(x_n-^(alpha*^g_n))<=f x_n-alpha*sigma*(dot g_n g_n)
    where g_n = dfdx x_n


-- Returns a backtracked acceptable step size for x_n
-- according to the Armijo rule.
-- @param f the function being optimised
--        df_dx its gradient function
--        x_n the current guess
--        sigma the sigma value for armijo rule
--        an initial guess for the step size (eta reduced parameter)
backtrack :: (Vec -> Entry) ->
             (Vec -> Vec) ->
             Vec ->
             Double ->
             Double ->
             Double
backtrack f dfdx x_n sigma = head .
            dropWhile (\alpha -> not (satisfiesArmijo f dfdx x_n alpha sigma)) .
            iterate (/2)


-- Returns an initial guess for the step size for the next x value
-- in such a way as to minimise amount of backtracking needed
nextStepGuess :: (Vec -> Vec) ->
                 Vec ->
                 Vec ->
                 Double ->
                 Double
nextStepGuess dfdx x_n_min1 x_n prevStep = prevStep * (dot g_n_min1 g_n_min1) / (dot g_n g_n)
    where g_n = dfdx x_n
          g_n_min1 = dfdx x_n_min1



-- Returns whether gradient is within tolerance (absolutely or relatively to start)
gradSmallEnough :: (Vec -> Vec) ->
                   Vec ->
                   Double ->
                   Vec ->
                   Bool
gradSmallEnough dfdx firstGrad tol x = (dot grad grad) < tol * (1 + (dot firstGrad firstGrad))
    where grad = dfdx x