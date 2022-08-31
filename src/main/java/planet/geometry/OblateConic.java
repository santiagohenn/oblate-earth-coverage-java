package planet.geometry;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import planet.coordinates.Converter;
import planet.ellipsoids.Ellipsoid;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeoutException;

import static java.lang.Math.*;

public class OblateConic {

    public static final double WGS84_F = 1 / 298.257223563;
    public static final double WGS84_E = 0.081819190842613;
    public static final double a_avg = 6371.009D;
    public static final double a = 6378.137D;
    public static final double a2 = Math.pow(a, 2);
    public static final double a4 = Math.pow(a, 4);
    public static final double b = a * (1 - WGS84_F);
    public static final double b2 = Math.pow(b, 2);
    public static final double b4 = Math.pow(b, 4);
    public static final double E2 = Math.pow(WGS84_E, 2);

    public static final Ellipsoid e = new Ellipsoid(6378.137,1 / 298.257223563,3.986004418E14,7.2921150E-5);
    public static final Converter conv = new Converter(e);

    public static List<double[]> drawLLAConic(double x, double y, double z, double eta, int segments) throws Exception {
        List<double[]> coordinates = drawConic(x, y, z, eta, 0, 0, segments, 0);
        return conv.ecef2lla(coordinates);
    }

    public static List<double[]> drawLLAConic(double x, double y, double z, double epsilon, double tol, int segments) throws Exception {
        List<double[]> coordinates = drawConic(x, y, z, 0, epsilon, tol, segments, 1);
        return conv.ecef2lla(coordinates);
    }

    public static List<double[]> drawConic(double x, double y, double z, double eta, int segments) throws Exception {
        return drawConic(x, y, z, eta, 0, 0, segments, 0);
    }

    public static List<double[]> drawConic(double x, double y, double z, double epsilon, double tol, int segments) throws Exception {
        return drawConic(x, y, z, 0, epsilon, tol, segments, 1);
    }

    /**
     * This method returns a list of coordinates depicting the access area of a satellite located at the Earth Centered
     * Earth Fixed coordinates (x,y,z). The two possible methods are: (0) the cone of access is computed using a conic
     * on the satellite with eta = half-angle of visibility or (1) the cone of access is computed using the elevation
     * threshold epsilon.
     *
     * @param x the x component of the satellite's position (Km)
     * @param y the y component of the satellite's position (Km)
     * @param z the z component of the satellite's position (Km)
     * @param eta the satellite's cone half-angle of visibility (degrees)
     * @param epsilon elevation threshold (degrees)
     * @param tol tolerance when obtaining the half-angle from a given elevation threshold (degrees)
     * @param segments the amount of segments for a half-cone of visibility
     * @return a List of double[] containing the polygons coordinates in degrees
     **/
    public static List<double[]> drawConic(double x, double y, double z, double eta, double epsilon, double tol,
                                           int segments, int type) throws Exception {

        Vector3D pos = new Vector3D(x, y, z);

        // n
        Vector3D gcPos = pos.scalarMultiply(-1.0D / pos.getNorm());

        // Line-of-sight projection onto the Earth surface (nadir)
        Vector3D nadir = getNadir(pos, gcPos);

        // Definition of a new reference frame
        // The new reference frame is shifted w.r.t. the Geocentric frame:
        // - origin in the intersection of the line of sight with the Equatorial Plane
        // - axes parallel to the Geocentric Inertial frame

        // Opposite direction to the line of sight
        Vector3D o = gcPos.scalarMultiply(-1D);

        // Determination of the geographic coordinates
        double lon_nadir = Math.atan2(o.getY(), o.getX());
        double lat_nadir = Math.asin(o.getZ());

        // Definition of the angles of the rotation matrix 321
        double phi = lon_nadir;
        double theta = -lat_nadir;

        List<double[]> coordinates1 = new ArrayList<>(segments);
        List<double[]> coordinates2 = new ArrayList<>(2 * segments);

        double step = Math.PI / (segments - 1);
        for (double psi = 0; psi <= Math.PI * 1.00001; psi += step) {

            if (psi > Math.PI) {
                psi = Math.PI;
            }

            // Rotation matrix from Geocentric to horizon frame
            double[][] A_321 = {{cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta)},
                    {-cos(psi) * sin(phi) + sin(psi) * sin(theta) * cos(phi), cos(psi) * cos(phi) + sin(psi) * sin(theta) * sin(phi),
                            sin(psi) * cos(theta)},
                    {sin(psi) * sin(phi) + sin(theta) * cos(phi) * cos(psi), -sin(psi) * cos(phi) + sin(theta) * sin(phi) * cos(psi),
                            cos(theta) * cos(psi)}};

            Vector3D n_prime = getNPrime(A_321);

            double d = n_prime.dotProduct(pos);

            // Set the components of the normal to the ellipse plane
            double n1_prime = n_prime.getX();
            double n2_prime = n_prime.getY();
            double n3_prime = n_prime.getZ();

            // Ellipse's centre in the Geocentric inertial frame
            // (Transposed)
            Vector3D s = new Vector3D((n1_prime * d) / (1 - E2 * Math.pow(n3_prime, 2)),
                    (n2_prime * d) / (1 - E2 * Math.pow(n3_prime, 2)),
                    ((b2 / a2) * n3_prime * d) / (1 - E2 * Math.pow(n3_prime, 2)));

            // Ellipse's semi-major axis
            double a_tilde = a * sqrt(1 - Math.pow(d, 2.0) / (a2 * (1 - E2 * Math.pow(n3_prime, 2))));

            // Ellipse's semi-minor axis
            double b_tilde = b * (sqrt(1 - (Math.pow(d, 2.0) / a2) - E2 * Math.pow(n3_prime, 2))) / (1 - E2 * Math.pow(n3_prime, 2.0));

            Vector3D e;
            Vector3D u;

            if (n1_prime == 0 && n2_prime == 0 && n3_prime == 1) {
                // Ellipse 's semi-major axis direction in the Geocentric frame
                e = new Vector3D(1, 0, 0);
                // Ellipse 's semi-major axis direction in the Geocentric frame
                u = new Vector3D(0, 1, 0);
            } else if (n1_prime == 0 && n2_prime == 0 && n3_prime == -1) {
                // Ellipse 's semi-major axis direction in the Geocentric frame
                e = new Vector3D(1, 0, 0);
                // Ellipse 's semi-major axis direction in the Geocentric frame
                u = new Vector3D(0, -1, 0);
            } else {
                var arg = Math.pow(n1_prime, 2.0) + Math.pow(n2_prime, 2.0);

                // Ellipse 's semi-major axis direction in the Geocentric frame
                e = new Vector3D(n2_prime, -n1_prime, 0);
                e = e.scalarMultiply(1 / sqrt(arg));

                // Ellipse 's semi-minor axis direction in the Geocentric frame
                u = new Vector3D(-n1_prime * n3_prime, -n2_prime * n3_prime, arg);
                u = u.scalarMultiply(-1.00 / sqrt(arg));

            }

            //  Position vector of the S/C w.r.to the ellipse's centre
            Vector3D r_line = pos.subtract(s);

            // Rotation matrix from Geocentric to the Local frame of the ellipse
            double[][] aRot = {{e.getX(), e.getY(), e.getZ()},
                    {u.getX(), u.getY(), u.getZ()},
                    {n_prime.getX(), n_prime.getY(), n_prime.getZ()}};

            // Position vector of the S/C w.r.to the ellipse centre aRot*r_line;
            Vector3D r_line_local =
                    new Vector3D(aRot[0][0] * r_line.getX() + aRot[0][1] * r_line.getY() + aRot[0][2] * r_line.getZ(),
                            aRot[1][0] * r_line.getX() + aRot[1][1] * r_line.getY() + aRot[1][2] * r_line.getZ(),
                            aRot[2][0] * r_line.getX() + aRot[2][1] * r_line.getY() + aRot[2][2] * r_line.getZ());

            // Definition of the S/C coordinates in the local frame
            double e_sc = r_line_local.getX();
            double u_sc = r_line_local.getY();

            // Coefficients derived from the normalisation of the second-degree
            // equation to reduce te numerical error
            double coeff_1 = (2 * e_sc * u_sc) / (Math.pow(a_tilde, 2) - Math.pow(e_sc, 2));
            double coeff_2 = (Math.pow(b_tilde, 2) - Math.pow(u_sc, 2)) / (Math.pow(a_tilde, 2) - Math.pow(e_sc, 2));

            // Discrimant of the second degree-equation
            double Delta_hor = Math.pow(coeff_1, 2) - 4 * coeff_2;

            // Slopes of the tangents: condition to assign the proper P1 and P2
            double m_T1, m_T2;

            if (e_sc * u_sc >= 0) {
                m_T1 = (-coeff_1 + sqrt(Delta_hor)) / 2;
                m_T2 = (-coeff_1 - sqrt(Delta_hor)) / 2;
            } else {
                m_T1 = (-coeff_1 - sqrt(Delta_hor)) / 2;
                m_T2 = (-coeff_1 + sqrt(Delta_hor)) / 2;
            }

            // Vertical intercept of the tangents
            double q_T1 = -m_T1 * e_sc + u_sc;
            double q_T2 = -m_T2 * e_sc + u_sc;

            // Horizon Points
            double e_T1 = (m_T1 * (-q_T1) * Math.pow(a_tilde, 2.0)) / (Math.pow(b_tilde, 2.0)
                    + Math.pow(a_tilde, 2.0) * Math.pow(m_T1, 2.0));
            double e_T2 = (m_T2 * (-q_T2) * Math.pow(a_tilde, 2.0)) / (Math.pow(b_tilde, 2.0)
                    + Math.pow(a_tilde, 2.0) * Math.pow(m_T2, 2.0));

            double u_T1 = m_T1 * e_T1 + q_T1;
            double u_T2 = m_T2 * e_T2 + q_T2;
            Vector3D r_T1 = new Vector3D(e_T1, u_T1, 0);
            Vector3D r_T2 = new Vector3D(e_T2, u_T2, 0);

            // Assign the position vector P1 and P2 to the forward and backward points
            Vector3D cross_prod = r_T2.crossProduct(r_T1);

            if (cross_prod.getZ() < 0) {
                r_T1 = new Vector3D(e_T2, u_T2, 0);
                r_T2 = new Vector3D(e_T1, u_T1, 0);
            }

            // Boresight-angle and ground range-angle
            double arg1 = (r_line_local.subtract(r_T1)).dotProduct(r_line_local);
            double eta_hor_1 = Math.acos((arg1) / (r_line_local.subtract(r_T1).getNorm() * r_line_local.getNorm()));
            double arg2 = (r_line_local.subtract(r_T2)).dotProduct(r_line_local);
            double eta_hor_2 = Math.acos((arg2) / (r_line_local.subtract(r_T2).getNorm() * r_line_local.getNorm()));

            double lambda_hor = Math.acos(r_T1.dotProduct(r_T2) / (r_T1.getNorm() * r_T2.getNorm()));

            // Angle between the S/C line-of-sight and the the semi-major axis direction e
            double alpha_SC = Math.atan2(u_sc, e_sc);

            if (Math.abs(alpha_SC - PI / 2) < 10E-15) {
                alpha_SC = PI / 2;
            } else if (Math.abs(alpha_SC + PI / 2) < 10E-15) {
                alpha_SC = -PI / 2;
            } else if (Math.abs(alpha_SC - PI) < 10E-15) {
                alpha_SC = PI;
            } else if (Math.abs(alpha_SC + PI) < 10E-15) {
                alpha_SC = -PI;
            }

            // Total horizon boresight angle
            double eta_hor = eta_hor_1 + eta_hor_2;

            // Computation of the effective horizon
            List<double[]> vectors = new ArrayList<>();
            // Compute the points in the local frame according to the method

            if (type == 0) {    // Using eta
                vectors = computeHalfAperture(a_tilde, b_tilde, alpha_SC, eta, eta_hor_1, eta_hor_2, r_line_local);
            } else if (type == 1) {     // Using th
//                vectors = computeWithElevationGD(a_tilde, b_tilde, alpha_SC, eta_hor_1, eta_hor_2, epsilon, tol, r_line_local);
                vectors = computeWithElevationBisect(a_tilde, b_tilde, alpha_SC, eta_hor_1, eta_hor_2, epsilon, tol, r_line_local);

            }

            // Inverse transformation from local to Geocentric inertial frame
            double[] P1 = vectors.get(0);
            double[] P2 = vectors.get(1);

            Vector3D P1_3D = new Vector3D(aRot[0][0] * P1[0] + aRot[1][0] * P1[1] + aRot[2][0] * P1[2],
                    aRot[0][1] * P1[0] + aRot[1][1] * P1[1] + aRot[2][1] * P1[2],
                    aRot[0][2] * P1[0] + aRot[1][2] * P1[1] + aRot[2][2] * P1[2]);

            Vector3D P2_3D = new Vector3D(aRot[0][0] * P2[0] + aRot[1][0] * P2[1] + aRot[2][0] * P2[2],
                    aRot[0][1] * P2[0] + aRot[1][1] * P2[1] + aRot[2][1] * P2[2],
                    aRot[0][2] * P2[0] + aRot[1][2] * P2[1] + aRot[2][2] * P2[2]);

            // Relative reference systems formulation to align the origin
            P1_3D.add(s);
            P2_3D.add(s);

            coordinates1.add(new double[]{P1_3D.getX(), P1_3D.getY(), P1_3D.getZ()});
            coordinates2.add(new double[]{P2_3D.getX(), P2_3D.getY(), P2_3D.getZ()});

        }

        coordinates1.addAll(coordinates2);
        return coordinates1;

    }


    private static List<double[]> computeHalfAperture(double a_tilde, double b_tilde, double alpha_SC, double etaDeg,
                                                      double eta_hor_1, double eta_hor_2, Vector3D r_line_local) throws Exception {

        List<double[]> coordinates = new ArrayList<>(2);

        // Definition of the S/C coordinates in the local frame
        double e_sc = r_line_local.getX();
        double u_sc = r_line_local.getY();

        // Aperture angles initialisation
        double eta = Math.toRadians(etaDeg);

        // Check if your aperture angle is bigger than the horizon one
        // if eta > eta_hor_1 || eta > eta_hor_2
        // Log.error("The half-aperture angle is greater than the horizon-boresight angle. Reduce the value of the half-aperture angle.");

        double eta_1, eta_2;

        if (eta > eta_hor_1) {
            eta_1 = eta_hor_1 - 0.0001 * Math.PI / 180;
        } else {
            eta_1 = eta;
        }

        if (eta > eta_hor_2) {
            eta_2 = eta_hor_2 - 0.0001 * Math.PI / 180;
        } else {
            eta_2 = eta;
        }

        // Angle of the secants w.r.to the semi-major axis direction
        double alpha_P1 = (alpha_SC - eta_1);
        double alpha_P2 = (alpha_SC + eta_2);

        //Slopes of the two secants
        double m_P1 = Math.tan(alpha_P1);
        double m_P2 = Math.tan(alpha_P2);

        double[] eComponents = getEComponents(alpha_SC, alpha_P1, alpha_P2, a_tilde, b_tilde, e_sc, u_sc);

        double e_P1 = eComponents[0];
        double e_P2 = eComponents[1];
        double u_P1 = m_P1 * e_P1 - m_P1 * e_sc + u_sc;
        double u_P2 = m_P2 * e_P2 - m_P2 * e_sc + u_sc;

        //// Determination of epsilons
        double[] epsilons = getEpsilons(e_P1, u_P1, e_P2, u_P2, m_P1, m_P2, a_tilde, b_tilde);

        coordinates.add(new double[]{e_P1, u_P1, 0});
        coordinates.add(new double[]{e_P2, u_P2, 0});

        return coordinates;

    }

    public static List<double[]> computeWithElevationGD(double a_tilde, double b_tilde, double alpha_SC,
                                                        double eta_hor_1, double eta_hor_2, double epsilon,
                                                        double tol, Vector3D r_line_local) throws Exception {

        List<double[]> coordinates = new ArrayList<>();

        // Definition of the S/C coordinates in the local frame
        double e_sc = r_line_local.getX();
        double u_sc = r_line_local.getY();

        // Aperture angles initialisation
        double eta_1 = eta_hor_1 - Math.toRadians(0.0001);
        double eta_2 = eta_hor_2 - Math.toRadians(0.0001);

        // Errors and tolerance initialisation
        double err_1 = 1;
        double err_2 = 1;

        double e_P1 = 0, e_P2 = 0, u_P1 = 0, u_P2 = 0;
        double alpha_P1 = 0, alpha_P2 = 0;
        int wdt = -1;

        while ((err_1 > tol || err_2 > tol) && wdt < 1000) {
            wdt++;

            // Angle of the secants w.r.to the semi-major axis direction
            alpha_P1 = (alpha_SC - eta_1);
            alpha_P2 = (alpha_SC + eta_2);

            //Slopes of the two secants
            double m_P1 = Math.tan(alpha_P1);
            double m_P2 = Math.tan(alpha_P2);

            double[] eComponents = getEComponents(alpha_SC, alpha_P1, alpha_P2, a_tilde, b_tilde, e_sc, u_sc);

            e_P1 = eComponents[0];
            e_P2 = eComponents[1];
            u_P1 = m_P1 * e_P1 - m_P1 * e_sc + u_sc;
            u_P2 = m_P2 * e_P2 - m_P2 * e_sc + u_sc;

            double[] epsilons = getEpsilons(e_P1, u_P1, e_P2, u_P2, m_P1, m_P2, a_tilde, b_tilde);

            // Error update
            double epsilon_1 = epsilons[0];
            double epsilon_2 = epsilons[1];

            // Error update
            err_1 = epsilon - epsilon_1;
            err_2 = epsilon - epsilon_2;

            // Decrease the initial guess for the half-aperture angle
            if (epsilon_1 < epsilon) {
                eta_1 = eta_1 - 0.0001;
            }

            if (epsilon_2 < epsilon) {
                eta_2 = eta_2 - 0.0001;
            }

        }

        if (wdt >= 1000) {
            throw new TimeoutException("The algorithm was unable to converge to an eta for the given epsilon");
        }

        coordinates.add(new double[]{e_P1, u_P1, 0});
        coordinates.add(new double[]{e_P2, u_P2, 0});

        return coordinates;

    }

    public static List<double[]> computeWithElevationBisect(double a_tilde, double b_tilde, double alpha_SC,
                                                            double eta_hor_1, double eta_hor_2, double epsilon,
                                                            double tol, Vector3D r_line_local) throws Exception {

        List<double[]> coordinates = new ArrayList<>(2);

        // Definition of the S/C coordinates in the local frame
        double e_sc = r_line_local.getX();
        double u_sc = r_line_local.getY();

        // Aperture angles initialization
        double etaSpherical = getSphericalEta(r_line_local, epsilon);
        // Double eta_1_sup = eta_hor_1 - 0.0001 * PI / 180;
        double eta_1_sup = eta_hor_1 - 0.0001 * PI / 180;
        double eta_2_sup = eta_hor_2 - 0.0001 * PI / 180;

        // 1 deg = 0.017453292519943 rad maximum error from spherical (according to paper)
        double eta_1_inf = etaSpherical - Math.toRadians(1);
        double eta_2_inf = etaSpherical - Math.toRadians(1);

        double eta_1 = (eta_1_sup + eta_1_inf) / 2; // Initial guess
        double eta_2 = (eta_2_sup + eta_2_inf) / 2; // Initial guess

        // Vectors init
        double e_P1 = 0, e_P2 = 0, u_P1 = 0, u_P2 = 0;
        double eps1 = 0, eps2 = 0;
        int it = 0;

        while ((Math.abs(epsilon - eps1) > tol || Math.abs(epsilon - eps2) > tol) && it < 1000) {

            // Angle of the secants w.r.to the semi-major axis direction
            double alpha_P1 = (alpha_SC - eta_1);
            double alpha_P2 = (alpha_SC + eta_2);

            //Slopes of the two secants
            double m_P1 = Math.tan(alpha_P1);
            double m_P2 = Math.tan(alpha_P2);

            double[] eComponents = getEComponents(alpha_SC, alpha_P1, alpha_P2, a_tilde, b_tilde, e_sc, u_sc);

            e_P1 = eComponents[0];
            e_P2 = eComponents[1];
            u_P1 = m_P1 * e_P1 - m_P1 * e_sc + u_sc;
            u_P2 = m_P2 * e_P2 - m_P2 * e_sc + u_sc;

            double[] epsilons = getEpsilons(e_P1, u_P1, e_P2, u_P2, m_P1, m_P2, a_tilde, b_tilde);

            // Error update
            eps1 = epsilons[0];
            eps2 = epsilons[1];

            if (Math.abs(epsilon - eps1) < tol && Math.abs(epsilon - eps2) < tol) {
                break;
            }

            if (Math.abs(epsilon - eps1) > tol) {
                if (eps1 < epsilon) {
                    eta_1_sup = eta_1;
                } else {
                    eta_1_inf = eta_1;
                }
                // Update eta
                eta_1 = (eta_1_sup + eta_1_inf) / 2;
            }

            if (Math.abs(epsilon - eps2) > tol) {
                if (eps2 > epsilon) {
                    eta_2_inf = eta_2;
                } else {
                    eta_2_sup = eta_2;
                }
                // Update eta
                eta_2 = (eta_2_sup + eta_2_inf) / 2;
            }

            it++;

        }

        coordinates.add(new double[]{e_P1, u_P1, 0});
        coordinates.add(new double[]{e_P2, u_P2, 0});

        return coordinates;

    }


    private static double[] getEpsilons(double e_P1, double u_P1, double e_P2, double u_P2, double m_P1, double m_P2, double a_tilde, double b_tilde) {

        // Determination of the tangents in P1 and P2
        double m_t_P1, m_t_P2;

        // Tangent slope
        var v0 = pow(a_tilde, 2.0) * sqrt(1 - pow((e_P1 / a_tilde), 2.0));
        var v1 = pow(a_tilde, 2.0) * sqrt(1 - pow((e_P2 / a_tilde), 2.0));

        if (u_P1 >= 0) {     // Selection of the semi-ellipse
            m_t_P1 = (-b_tilde * e_P1) / v0;
        } else {
            m_t_P1 = (b_tilde * e_P1) / v0;
        }

        if (u_P2 >= 0) {
            m_t_P2 = (-b_tilde * e_P2) / v1;
        } else {
            m_t_P2 = (b_tilde * e_P2) / v1;
        }

        // Elevation angle computation: angle between two lines
        double epsilon_1 = Math.atan((-m_t_P1 + m_P1) / (1 + m_t_P1 * m_P1));
        double epsilon_2 = Math.atan((m_t_P2 - m_P2) / (1 + m_t_P2 * m_P2));

        // Conversion into degrees
        epsilon_1 = Math.toDegrees(epsilon_1);
        epsilon_2 = Math.toDegrees(epsilon_2);

        return new double[]{epsilon_1, epsilon_2};

    }

    private static double[] getEComponents(double alpha_SC, double alpha_P1, double alpha_P2, double a_tilde, double b_tilde, double e_sc, double u_sc) throws Exception {

        //Slopes of the two secants
        double m_P1 = Math.tan(alpha_P1);
        double m_P2 = Math.tan(alpha_P2);

        //Coefficients to normalise the equation and reduce the error
        double coeff2_1 = (2.0 * m_P1 * pow(a_tilde, 2.0) * (u_sc - m_P1 * e_sc)) / (pow(a_tilde, 2.0) * pow(m_P1, 2.0) + pow(b_tilde, 2.0));
        double coeff3_1 = (-pow(a_tilde, 2.0) * pow(b_tilde, 2.0) + pow(a_tilde, 2.0) * pow((u_sc - m_P1 * e_sc), 2)) / (pow(a_tilde, 2.0) * pow(m_P1, 2.0) + pow(b_tilde, 2.0));

        double coeff2_2 = (2.0 * m_P2 * pow(a_tilde, 2.0) * (u_sc - m_P2 * e_sc)) / (pow(a_tilde, 2.0) * pow(m_P2, 2.0) + pow(b_tilde, 2.0));
        double coeff3_2 = (-pow(a_tilde, 2.0) * pow(b_tilde, 2.0) + pow(a_tilde, 2.0) * pow((u_sc - m_P2 * e_sc), 2)) / (pow(a_tilde, 2.0) * pow(m_P2, 2.0) + pow(b_tilde, 2.0));

        //Discriminant of the second-degree equation
        double Delta_P1 = pow(coeff2_1, 2.0) - 4.0 * coeff3_1;
        double Delta_P2 = pow(coeff2_2, 2.0) - 4.0 * coeff3_2;

        double e_P1 = 0, e_P2 = 0;
        // First case
        if (alpha_SC >= 0 && alpha_SC <= Math.PI / 2) { // First Quadrant
            e_P1 = (-coeff2_1 + sqrt(Delta_P1)) / 2;

            //Selection of the two right solutions among the possible four
            if (alpha_P2 <= Math.PI / 2) {
                e_P2 = (-coeff2_2 + sqrt(Delta_P2)) / 2;
            } else {
                e_P2 = (-coeff2_2 - sqrt(Delta_P2)) / 2;
            }

        } else if (alpha_SC > Math.PI / 2 && alpha_SC <= Math.PI) {   // Second Quadrant
            //Selection of the two right solutions among the possible four
            if (alpha_P1 <= Math.PI / 2) {
                e_P1 = (-coeff2_1 + sqrt(Delta_P1)) / 2;
            } else {
                e_P1 = (-coeff2_1 - sqrt(Delta_P1)) / 2;
            }

            e_P2 = (-coeff2_2 - sqrt(Delta_P2)) / 2;

        } else if (alpha_SC >= -Math.PI && alpha_SC <= -Math.PI / 2) {   // Third Quadrant
            e_P1 = (-coeff2_1 - sqrt(Delta_P1)) / 2;

            //Selection of the two right solutions among the possible four
            if (alpha_P2 <= -Math.PI / 2) {
                e_P2 = (-coeff2_2 - sqrt(Delta_P2)) / 2;
            } else {
                e_P2 = (-coeff2_2 + sqrt(Delta_P2)) / 2;
            }

        } else if (alpha_SC > -Math.PI / 2.0 && alpha_SC < 0.0) {   // Fourth Quadrant
            //Selection of the two right solutions among the possible four
            if (alpha_P1 <= Math.PI / 2.0) {
                e_P1 = (-coeff2_1 - sqrt(Delta_P1)) / 2;
            } else {
                e_P1 = (-coeff2_1 + sqrt(Delta_P1)) / 2;
            }

            e_P2 = (-coeff2_2 + sqrt(Delta_P2)) / 2;

        } else {
            throw new Exception("Not considered alpha_SC case");
        }

        return new double[]{e_P1, e_P2};

    }

    private static Vector3D getNPrime(double[][] A_321) {
        return new Vector3D(A_321[2][0], A_321[2][1], A_321[2][2]);
    }

    private static Vector3D getNadir(Vector3D pos, Vector3D gcPos) {

        // Coordinates of the S/C in the Geocentric Inertial Frame
        double XX = pos.getX();
        double YY = pos.getY();
        double ZZ = pos.getZ();

        double nx = gcPos.getX();
        double ny = gcPos.getY();
        double nz = gcPos.getZ();

        double deltaTan = 2 * b4 * nx * ny * XX * YY + 2 * a2 * b2 * nx * nz * XX * ZZ + 2 * a2 * b2 * ny * nz * YY * ZZ
                - b4 * Math.pow(nx, 2) * Math.pow(YY, 2) - a2 * b2 * Math.pow(nx, 2) * Math.pow(ZZ, 2) + a2 * b4 * Math.pow(nx, 2) - b4 * Math.pow(ny, 2) * Math.pow(XX, 2)
                - a2 * b2 * Math.pow(ny, 2) * Math.pow(ZZ, 2) + a2 * b4 * Math.pow(ny, 2) - a2 * b2 * Math.pow(nz, 2) * Math.pow(XX, 2)
                - a2 * b2 * Math.pow(YY, 2) * Math.pow(nz, 2) + a4 * b2 * Math.pow(nz, 2);

        // Solution of the parametric equation
        var v = -(b2 * nx * XX + b2 * ny * YY + a2 * nz * ZZ);
        var v1 = b2 * Math.pow(nx, 2) + b2 * Math.pow(ny, 2) + a2 * Math.pow(nz, 2);
        double t_1 = (v + Math.sqrt(deltaTan)) / v1;
        double t_2 = (v - Math.sqrt(deltaTan)) / v1;

        // Selection of the minimum value corresponding to the real projection
        double t = Math.min(t_1, t_2);

        // Computation of the Intersection point coordinates
        double nadirX = XX + nx * t;
        double nadirY = YY + ny * t;
        double nadirZ = ZZ + nz * t;

        return new Vector3D(nadirX, nadirY, nadirZ);

    }

    private static double getSphericalEta(Vector3D ecefPos, double epsilon) {
        double height = ecefPos.getNorm() - a_avg;
        return Math.asin((a_avg * Math.cos(Math.toRadians(epsilon))) / (a_avg + height));
    }





}
