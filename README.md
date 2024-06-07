# -utl-r-generating-independent-random-variables-from-a-bivariate-normal-distribution
R-generating independent random variables from a bivariate normal distribution
    %let pgm=utl-r-generating-independent-random-variables-from-a-bivariate-normal-distribution;

     R-generating independent random variables from a bivariate normal distribution

     Hope I got this right?

       Agenda
          1 theory
          2 application
          3 Symbolic math

        In the case of a bivariate normal distribution, if the correlation coef is zero
        then the independently generated random variables are independent.

        github
        https://tinyurl.com/2hfzxbfh
        https://github.com/rogerjdeangelis/utl-r-generating-independent-random-variables-from-a-bivariate-normal-distribution


        if X and Y are drawn independently from a bivariate normal density, they will be uncorrelated.
        Lets make sure we have independence. To theorectically check that any multivariate
        distribution has idependent variables, the joint
        mutiivariate density must equal the product of marginal densities.

    /*   _
    / | | |_| |__   ___  ___  _ __ _   _
    | | | __| `_ \ / _ \/ _ \| `__| | | |
    | | | |_| | | |  __/ (_) | |  | |_| |
    |_|  \__|_| |_|\___|\___/|_|   \__, |
                                   |___/
    */

         For independence the joint density must equal the product of the marginals.
         Independence implies uncorrelated but not the reverse.

         pdf(x,y)=pdf(x)*pdf(y)

         BIVARIATE NORMAL WITH CORRELATION=0 (independent random samples)

                      _                                    _
                     |                 2                 2  |
                     |     (-mu2 + y)        (-mu1 + x)     |
                     | - --------------- - ---------------  |
                     |             2                 2      |
                     |_    2*sigma2          2*sigma1      _|

         pdf(x,y)  = ----------------------------------------
                             2*pi*sigma1*sigma2

         NORMAL X:
                             _               _
                            |             2   |
                            |  -(-mu1 + x)    |
                            |  -------------  |
                            |            2    |
                       ___  |_   2*sigma1    _|
                     \/ 2 *e
         pdf(x)    = ---------------------------
                               ____
                           2*\/ pi *sigma1

                     NORMAL Y:
                             _              _
                            |            2   |
                            | -(-mu2 + y)    |
                            | -------------  |
                            |           2    |
                       ___  |_  2*sigma2    _|
                     \/ 2 *e
         pdf(y  )  =  ------------------------
                           ____
                       2*\/ pi *sigma2

                                                                              2              2
                         PRODUCT EQUALS THE BIVARIATE PDF     Note: (-mu2 + y)   =  (mu2 - y)
                           _                          _        _                                    _
                          |            2            2  |      |                 2                 2  |
                          |   (mu2 - y)    (mu1 - x)   |      |     (-mu2 + y)        (-mu1 + x)     |
                          | - ---------- - ----------  |      | - --------------- - ---------------  |
                          |           2            2   |      |             2                 2      |
                          |_  2*sigma2     2*sigma1   _|      |_    2*sigma2          2*sigma1      _|
                         e                                  e
         pdf(x)*pdf(y) =  ------------------------------  = ------------------------------------------
                             2*pi*sigma1*sigma2                       2*pi*sigma1*sigma2

    /*                   _
    (_)_ __  _ __  _   _| |_
    | | `_ \| `_ \| | | | __|
    | | | | | |_) | |_| | |_
    |_|_| |_| .__/ \__,_|\__|
            |_|
    */

     means

      mu1=0
      mu2=0

    covariance matrix

       V1 V2
        1  0
        0  1
    For generality we wiil use symbos mu1, mu2 and sigma1,sigma4
    but subsitute 0 or correlation
                    _                                    _
                   |                 2                 2  |
                   |     (-mu2 + y)        (-mu1 + x)     |
                   | - --------------- - ---------------  |
                   |             2                 2      |
                   |_    2*sigma2          2*sigma1      _|
      pdf(x,y) =  e
                  -----------------------------------------
                            2*pi*sigma1*sigma2

    /*___                      _ _           _   _
    |___ \    __ _ _ __  _ __ | (_) ___ __ _| |_(_) ___  _ __
      __) |  / _` | `_ \| `_ \| | |/ __/ _` | __| |/ _ \| `_ \
     / __/  | (_| | |_) | |_) | | | (_| (_| | |_| | (_) | | | |
    |_____|  \__,_| .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
                  |_|   |_|
     _ __  _ __ ___   ___ ___  ___ ___
    | `_ \| `__/ _ \ / __/ _ \/ __/ __|
    | |_) | | | (_) | (_|  __/\__ \__ \
    | .__/|_|  \___/ \___\___||___/___/
    |_|
    */

    %utl_rbegin;
    parmcards4;
    library(mvtnorm)
    source("c:/temp/fn_tosas9.R")
    mean_vec <- c(0,0)
    cov_mat=as.matrix(read.table(header=FALSE,text = "
    1 0
    0 1
    "))
    cov_mat
    want <- rmvnorm(10000000, mean_vec, cov_mat)
    cor(want)
    want<-want[1:5,]
    fn_tosas9(dataf=want)
    ;;;;
    %utl_rend;

    libname tmp "c:/temp";
    proc print data=tmp.want;
    run;quit;

    /*           _               _
      ___  _   _| |_ _ __  _   _| |_
     / _ \| | | | __| `_ \| | | | __|
    | (_) | |_| | |_| |_) | |_| | |_
     \___/ \__,_|\__| .__/ \__,_|\__|
                    |_|
    */

    /**************************************************************************************************************************/
    /*                                                                                                                        */
    /* covariance(x,y)=0 and corr(x,y)=0                                                                                      */
    /*                                                                                                                        */
    /* > cor(want) (correlations)                                                                                             */
    /*             [,1]        [,2]                                                                                           */
    /* [1,] 1.000000000 0.000153453                                                                                           */
    /* [2,] 0.000153453 1.000000000                                                                                           */
    /*                                                                                                                        */
    /* random sample                                                                                                          */
    /*                                                                                                                        */
    /*    COL_0       COL_1                                                                                                   */
    /*                                                                                                                        */
    /*   0.66522    -1.07444                                                                                                  */
    /*  -1.37979    -0.52272                                                                                                  */
    /*  -1.19242    -0.50048                                                                                                  */
    /*   0.56343     0.27889                                                                                                  */
    /*   0.67474    -0.54999                                                                                                  */
    /*                                                                                                                        */
    /**************************************************************************************************************************/
    /*____                       _           _ _                        _   _
    |___ /   ___ _   _ _ __ ___ | |__   ___ | (_) ___   _ __ ___   __ _| |_| |__
      |_ \  / __| | | | `_ ` _ \| `_ \ / _ \| | |/ __| | `_ ` _ \ / _` | __| `_ \
     ___) | \__ \ |_| | | | | | | |_) | (_) | | | (__  | | | | | | (_| | |_| | | |
    |____/  |___/\__, |_| |_| |_|_.__/ \___/|_|_|\___| |_| |_| |_|\__,_|\__|_| |_|
                 |___/
    */
    options ls=255 ps=255;
    %utl_pybegin;
    parmcards4;

    import sympy as sp
    from sympy.stats import Normal, density

    # Define symbols
    x, y = sp.symbols('x y')
    mu1, mu2 = sp.symbols('mu1 mu2')
    sigma1, sigma2 = sp.symbols('sigma1 sigma2')

    # Define the univariate normal distributions
    X = Normal('X', mu1, sigma1)
    Y = Normal('Y', mu2, sigma2)

    # Get the density functions of the univariate normal distributions
    density_X = density(X)(x)
    density_Y = density(Y)(y)

    # Define the bivariate normal distribution
    rho = sp.Symbol('rho')
    cov_matrix = sp.Matrix([[sigma1**2, rho*sigma1*sigma2], [rho*sigma1*sigma2, sigma2**2]])
    mean_vector = sp.Matrix([mu1, mu2])
    xy_vector = sp.Matrix([x, y])

    rho=0;

    # Bivariate normal density function
    bivariate_density = (1 / (2 * sp.pi * sigma1 * sigma2 * sp.sqrt(1 - rho**2))) * \
                        sp.exp(-1 / (2 * (1 - rho**2)) * \
                               ((x - mu1)**2 / sigma1**2 + (y - mu2)**2 / sigma2**2 - \
                                2 * rho * (x - mu1) * (y - mu2) / (sigma1 * sigma2)))

    # Verify the product of the univariate densities
    product_density = density_X * density_Y

    print("Intermediate Product Density:")
    sp.pprint(product_density)

    # Simplify and compare
    simplified_bivariate_density = sp.simplify(bivariate_density)
    simplified_product_density = sp.simplify(product_density)

    # Print the results
    print("Normal X:")
    sp.pprint(density_X)
    print("Normal Y")
    sp.pprint(density_Y)

    print("Bivariate Normal Density Function:")
    sp.pprint(simplified_bivariate_density)

    print("\nProduct of Univariate Normal Density Functions:")
    sp.pprint(simplified_product_density)

    # Check if they are equal
    are_equal = sp.simplify(simplified_bivariate_density - simplified_product_density) == 0
    print("\nAre the bivariate density and the product of univariate densities equal? ", are_equal)
    ;;;;
    %utl_pyend;

    /*              _
      ___ _ __   __| |
     / _ \ `_ \ / _` |
    |  __/ | | | (_| |
     \___|_| |_|\__,_|

    */
