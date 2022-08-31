package planet.geometry;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import planet.coordinates.Converter;
import planet.ellipsoids.Ellipsoid;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.*;

public class OblateConic {

    public Ellipsoid e = new Ellipsoid(6378.137,1 / 298.257223563,3.986004418E14,7.2921150E-5);
    public Converter conv = new Converter(e);

    private int drawingMethod;
    private double tolerance;
    private double halfApertureOfSensor;
    private double epsilon;
    private int segments;

    public OblateConic(Ellipsoid e) {
        this.e = e;
        this.drawingMethod = 0;
    }

    public void setDrawingMethod(int drawingMethod) {
        this.drawingMethod = drawingMethod;
    }

    /**
     * This method ...
     *
     * @param epsilon elevation threshold (degrees)
     * @param tolerance tolerance when obtaining the half-angle from a given elevation threshold (degrees)
     * @param segments the amount of segments for a half-cone of visibility
     **/
    public void setElevationParams(double epsilon, double tolerance, int segments) {
        this.epsilon = epsilon;
        this.tolerance = tolerance;
        this.segments = segments;
    }

    /**
     * This method ...
     *
     * @param halfApertureOfSensor the satellite's cone half-angle of visibility (degrees)
     * @param segments the amount of segments for a half-cone of visibility
     **/
    public void setSensorParams(double halfApertureOfSensor, int segments) {
        this.halfApertureOfSensor = halfApertureOfSensor;
        this.segments = segments;
    }

    public List<double[]> drawLLAConic(double x, double y, double z) throws Exception {
        List<double[]> coordinates = drawConic(x, y, z);
        return conv.ecef2lla(coordinates);
    }

    /**
     * This method returns a list of coordinates depicting the access area of a satellite located at the Earth Centered
     * Earth Fixed coordinates (x,y,z). The two possible methods are: (0) the cone of access is computed using a conic
     * sensor on the satellite with the half-angle of visibility or (1) the cone of access is computed using the elevation
     * threshold epsilon.
     *
     * @param x the x component of the satellite's position (Km)
     * @param y the y component of the satellite's position (Km)
     * @param z the z component of the satellite's position (Km)
     * @return a List of double[] containing the polygons coordinates in degrees
     **/
    public List<double[]> drawConic(double x, double y, double z) throws Exception {

        Vector3D pos = new Vector3D(x, y, z);

        // n
        Vector3D gcPos = pos.scalarMultiply(-1.0D / pos.getNorm());

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

            Vector3D n_prime = AuxiliaryVectors.getNPrime(A_321);

            double d = n_prime.dotProduct(pos);

            // Set the components of the normal to the ellipse plane
            double n1_prime = n_prime.getX();
            double n2_prime = n_prime.getY();
            double n3_prime = n_prime.getZ();

            // Ellipse's centre in the Geocentric inertial frame
            // (Transposed)
            Vector3D s = new Vector3D((n1_prime * d) / (1 - e.getE2() * Math.pow(n3_prime, 2)),
                    (n2_prime * d) / (1 - e.getE2() * Math.pow(n3_prime, 2)),
                    ((e.getB2() / e.getA2()) * n3_prime * d) / (1 - e.getE2() * Math.pow(n3_prime, 2)));

            // Ellipse's semi-major axis
            double a_tilde = e.getA() * sqrt(1 - Math.pow(d, 2.0) / (e.getA2() * (1 - e.getE2() * Math.pow(n3_prime, 2))));

            // Ellipse's semi-minor axis
            double b_tilde = e.getB() * (sqrt(1 - (Math.pow(d, 2.0) / e.getA2()) - e.getE2() * Math.pow(n3_prime, 2))) / (1 - e.getE2() * Math.pow(n3_prime, 2.0));

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
            double v = pow(a_tilde, 2) - pow(e_sc, 2);
            double coeff_1 = (2 * e_sc * u_sc) / v;
            double coeff_2 = (Math.pow(b_tilde, 2) - Math.pow(u_sc, 2)) / v;

            // Discriminant of the second degree-equation
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

            // Compute the points in the local frame according to the chosen method
            if (drawingMethod == 0) {    // Using eta
                vectors = HalfApertureVectors.computeHalfAperture(a_tilde, b_tilde, alpha_SC, halfApertureOfSensor, eta_hor_1, eta_hor_2, r_line_local);
            } else if (drawingMethod == 1) {     // Using th
                vectors = HalfApertureVectors.computeWithElevationBisect(a_tilde, b_tilde, alpha_SC, eta_hor_1, eta_hor_2, epsilon, tolerance, r_line_local);
            } else if (drawingMethod == 2) {
                vectors = HalfApertureVectors.computeWithElevationGD(a_tilde, b_tilde, alpha_SC, eta_hor_1, eta_hor_2, epsilon, tolerance, r_line_local);
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

}
