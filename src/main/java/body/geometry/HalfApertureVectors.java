package body.geometry;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.*;
import static java.lang.Math.sqrt;

public class HalfApertureVectors {


    public static List<double[]> computeHalfAperture(double a_tilde, double b_tilde, double alpha_SC, double etaDeg,
                                                     double eta_hor_1, double eta_hor_2, Vector3D r_line_local) {

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

        // Determination of epsilons
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
            throw new Exception("The algorithm was unable to converge to an eta for the given epsilon");
        }

        coordinates.add(new double[]{e_P1, u_P1, 0});
        coordinates.add(new double[]{e_P2, u_P2, 0});

        return coordinates;

    }

    public static List<double[]> computeWithElevationBisect(double a_tilde, double b_tilde, double alpha_SC,
                                                            double eta_hor_1, double eta_hor_2, double epsilon,
                                                            double tol, Vector3D r_line_local) {

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

    private static double[] getEComponents(double alpha_SC, double alpha_P1, double alpha_P2, double a_tilde, double b_tilde, double e_sc, double u_sc) {

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
            throw new ArithmeticException("Not considered alpha_SC case");
        }

        return new double[]{e_P1, e_P2};

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

    private static double getSphericalEta(Vector3D ecefPos, double epsilon) {
        double a_avg = 6371.009D;
        double height = ecefPos.getNorm() - a_avg;
        return Math.asin((a_avg * Math.cos(Math.toRadians(epsilon))) / (a_avg + height));
    }

}
