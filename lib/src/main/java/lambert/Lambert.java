package lambert;

/**
 * Implementation of an algorithm for the Lambert <font face="serif"><i>W<sub>0</sub></i></font> and
 * <font face="serif"><i>W<sub>-1</i></sub></font> real-only functions.
 * <p>The Lambert W function is defined as: <font face="serif"><i>W(z ⋅ e<sup>z</sup>) = z</i></font>
 * </p>
 * <p>
 *     This code is based on the FORTRAN algorithm by Toshio Fukushima,
 *     "Precise and fast computation of Lambert W-functions without transcendental function evaluations",
 *     Journal of Computational and Applied Mathematics, 244, 77-89.
 * </p>
 */
public final class Lambert {
    private Lambert() {}

    /**
     * Compute W(z) by a truncated Maclaurin series around the branch point W=-1
     *
     * <p>Reference: Toshio Fukushima (2013) J.Comp.Appl.Math., 244, 77-89
     *                "Precise and fast computation of Lambert W-functions
     *                 without transcendental function evaluations"
     */
    private static double lambertWSeries(final double p) {
        final double ap = Math.abs(p);

        return  -1 + p * (+1 +
                     p * (-0.33333333333333333333 +
                     p * (+0.15277777777777777778 +
                     p * (-0.079629629629629629630 +
                     p * (+0.044502314814814814815 +
                     p * (-0.025984714873603762493 + (ap < 0.01159 ? 0
                           : p * (+0.015635632532333921223 +
                             p * (-0.0096168920242994317068 +
                             p * (+0.0060145432529561178610 +
                             p * (-0.0038112980348919992267 + (ap < 0.0766 ? 0
                                   : p * (+0.0024408779911439826659 +
                                     p * (-0.0015769303446867842539 +
                                     p * (+0.0010262633205076071544 +
                                     p * (-0.00067206163115613620400 +
                                     p * (+0.00044247306181462090993 +
                                     p * (-0.00029267722472962744485 +
                                     p * (+0.00019438727605453931782 +
                                     p * (-0.00012957426685274881888 +
                                     p * (+0.000086650358052081271660 +
                                     p * (-0.000058113607504413816772))))))))))))))))))))));
    }

    private static double finalResult(double w, double y) {
        final double f0 = w - y;
        final double f1 = 1 + y;
        final double f00 = f0 * f0;
        final double f11 = f1 * f1;
        final double f0y = f0 * y;
        return w - 4 * f0 * (6 * f1 * (f11 + f0y) + f00 * y) / (f11 * (24 * f11 + 36 * f0y) + f00 * (6 * y * y + 8 * f1 * y + f0y));
    }

    /**
     * Compute <font face="serif"><i>W<sub>0</sub>(z)</i></font>, Lambert W-function of the branch 0
     */
    public static double lambertW0(double z) {
        return W0.lambert(z);
    }

    /**
     * Compute <font face="serif"><i>W<sub>-1</sub>(z)</i></font>, Lambert W-function of the branch -1
     */
    public static double lambertWm1(double z) {
        return Wm1.lambert(z);
    }

    /**
     * Halley iterate to refine Lambert_w estimate, taking at least one {@link #halley_step}.
     * <p>Repeat Halley steps until the <strong>last step</strong> had fewer than half the digits wrong,
     * the step we've just taken should have been sufficient to have completed the iteration.
     *
     * @param z Argument z for Lambert_w function.
     * @param w_est Lambert w estimate.
     * @return  New w estimate
     */
    private static double halley_iterate(double w_est, double z)
    {
        final double max_diff = Math.sqrt(Math.ulp(w_est));

        double diff = 2*max_diff;
        while (true) {
            double w_new = halley_step(w_est, z);
            if (diff <= max_diff) return w_new;
            diff = Math.abs(w_est - w_new);
            w_est = w_new;
        }
  }

    private static double halley_step(double w_est, double z) {
        double exp_west = Math.exp(w_est);
        double z_est = w_est * exp_west;
        double z_diff = z_est - z;
        return w_est - z_diff / (z_est + exp_west - z_diff * (1 + 1/(w_est+1)) / 2);
    }

    private static final class W0 {
        private static final double[] E = new double[66];
        private static final double[] G = new double[65];
        private static final double[] A = new double[12];
        private static final double[] B = new double[12];

        static {
            final double Em1 = 1.0 / Math.E;
            E[0] = Math.E;
            double ej = 1.0;
            E[1] = 1.0;
            G[0] = 0.0;
            for (int j = 1; j <= 64; j++) {
                ej *= Math.E;
                E[j + 1] = E[j] * Em1;
                G[j] = j * ej;
            }
            A[0] = Math.sqrt(Em1);
            B[0] = 0.5;
            for (int j = 1; j < 12; j++) {
                A[j] = Math.sqrt(A[j - 1]);
                B[j] = B[j - 1] * 0.5;
            }
        }

        private static double lambertW0ZeroSeries(final double z) {
            return z * (1 -
                   z * (1 -
                   z * (1.5 -
                   z * (2.6666666666666666667 -
                   z * (5.2083333333333333333 -
                   z * (10.8 -
                   z * (23.343055555555555556 -
                   z * (52.012698412698412698 -
                   z * (118.62522321428571429 -
                   z * (275.57319223985890653 -
                   z * (649.78717234347442681 -
                   z * (1551.1605194805194805 -
                   z * (3741.4497029592385495 -
                   z * (9104.5002411580189358 -
                   z * (22324.308512706601434 -
                   z * (55103.621972903835338 -
                   z *  136808.86090394293563
                   ))))))))))))))));
        }

        private static double lambert(final double z) {
            if (Math.abs(z) < 0.05) {
                return lambertW0ZeroSeries(z);
            }
            if (z < -0.35) {
                final double p2 = 2 * (Math.E * z + 1);
                if (p2 > 0) {
                    return lambertWSeries(Math.sqrt(p2));
                }
                if (p2 == 0) {
                    return -1;
                }
                // Argument z out of range
                return Double.NaN;
            }
            int n;
            for (n = 0; n <= 2; ++n) {
                if (G[n] > z) {
                    break;
                }
            }

            if (n > 2) {
                for (n = 4; n <= 64; n *= 2) {
                    if (G[n] > z) {
                        break;
                    }
                }
                if (n > 64) {
                    if (z == Double.POSITIVE_INFINITY) return z;
                    // Approximate lambert_w0 for z values that are outside range of lookup tables
                    double lz = Math.log(z);
                    double llz = Math.log(lz);
                    double w = lz - llz + (llz / lz); // Corless equation 4.19, page 349, and Chapeau-Blondeau equation 20, page 2162.
                    return halley_iterate(w, z);
                }

                int nh = n / 2;
                for (int j = 1; j <= 5; ++j) {
                    nh /= 2;
                    if (nh == 0) {
                        break;
                    }
                    if (G[n - nh] > z) {
                        n -= nh;
                    }
                }
            }
            --n;
            int jmax = 8;
            if (z <= -0.3) {
                jmax = 11;
            } else if (n <= 0) {
                jmax = 10;
            } else if (n == 1) {
                jmax = 9;
            }
            double y = z * E[n + 1];
            double w = n;
            for (int j = 0; j < jmax; ++j) {
                final double wj = w + B[j];
                final double yj = y * A[j];
                if (wj < yj) {
                    w = wj;
                    y = yj;
                }
            }
            return finalResult(w, y);
        }
    }


    private static final class Wm1 {
        private static final double[] E = new double[64];
        private static final double[] G = new double[64];
        private static final double[] A = new double[12];
        private static final double[] B = new double[12];

        static {
            final double Em1 = 1 / Math.E;
            double emj = Em1;
            E[0] = Math.E;
            G[0] = -Em1;
            for (int j = 1; j < 64; j++) {
                emj *= Em1;
                E[j] = E[j - 1] * Math.E;
                G[j] = -(j + 1) * emj;
            }
            A[0] = Math.sqrt(Math.E);
            B[0] = 0.5D;
            for (int j = 1; j < 12; j++) {
                A[j] = Math.sqrt(A[j - 1]);
                B[j] = B[j - 1] * 0.5D;
            }
        }

        private static double lambert(final double z) {
            if (z >= 0) {
                // Argument z out of range
                return z == 0 ? Double.NEGATIVE_INFINITY : Double.NaN;
            }
            if (z < -0.35) {
                final double p2 = 2 * (Math.E * z + 1);
                if (p2 > 0) {
                    return lambertWSeries(-Math.sqrt(p2));
                }
                if (p2 == 0) {
                    return -1;
                }
                // Argument z out of range
                return Double.NaN;
            }
            int n = 2;
            if (G[n - 1] <= z) {
                for (n = 4; n <= 64; n *= 2) {
                    if (G[n - 1] > z) {
                        break;
                    }
                }
                if (n > 64) {
                    // Argument z too small
                    double lz = Math.log(-z);
                    double llz = Math.log(-lz);
                    double guess = lz - llz + (llz / lz); // Chapeau-Blondeau equation 20, page 2162.
                    return halley_iterate(guess, z);
                }

                int nh = n / 2;
                for (int j = 1; j <= 5; ++j) {
                    nh /= 2;
                    if (nh <= 0) {
                        break;
                    }
                    if (G[n - nh - 1] > z) {
                        n -= nh;
                    }
                }
            }
            --n;
            int jmax = 11;
            if (n >= 8) {
                jmax = 8;
            } else if (n >= 3) {
                jmax = 9;
            } else if (n == 2) {
                jmax = 10;
            }
            double w = -n;
            double y = z * E[n - 1];
            for (int j = 0; j < jmax; ++j) {
                final double wj = w - B[j];
                final double yj = y * A[j];
                if (wj < yj) {
                    w = wj;
                    y = yj;
                }
            }
            return finalResult(w, y);
        }
    }
}
