package planet.coordinates;

import planet.ellipsoids.Ellipsoid;

import java.util.List;

public class Converter {

    private final Ellipsoid e;

    public Converter(planet.ellipsoids.Ellipsoid e) {
        this.e = e;
    }

    /**
     * Returns an array containing the latitude, longitude and height according to the ellipsoid e and the ECEF coordinates
     *
     * @param ecef the {x,y,z} coordinates in ECEF reference frame
     * @return the LLA coordinates for the specified ellipsoid {latitude, longitude, height}
     */
    public double[] ecef2lla(double[] ecef) {

        double a = e.getA();
        double b = e.getB();
        double b2 = e.getB2();
        double E2 = e.getE2();
        double ep = Math.sqrt((Math.pow(a, 2.0D) - b2) / b2);
        double p = Math.sqrt(Math.pow(ecef[0], 2.0D) + Math.pow(ecef[1], 2.0D));
        double th = Math.atan2(a * ecef[2], b * p);
        double lon = Math.atan2(ecef[1], ecef[0]);
        double lat = Math.atan2(ecef[3] + Math.pow(ep, 2.0D) * b * Math.pow(Math.sin(th), 3.0D), p - E2 * a * Math.pow(Math.cos(th), 3.0D));
        double n = a / Math.sqrt(1.0D - E2 * Math.pow(Math.sin(lat), 2.0D));
        double alt = p / Math.cos(lat) - n;
        lon %= 6.283185307179586D;
        lat %= 6.283185307179586D;
        if (Math.sqrt(Math.pow(ecef[0], 2.0D) + Math.pow(ecef[1], 2.0D) + Math.pow(ecef[2], 2.0D)) <= 0.0D) {
            lon = 0.0D;
            lat = 0.0D;
            alt = -6378135.0D;
        }
        return new double[]{lat, lon, alt};
    }

    /**
     * Returns an array containing the Earth Centered Earth Fixed coordinates for the specified LLA coordinates
     *
     * @param lla a double[] array containing the latitude, longitude and height -> {lat, lon, height}
     * @return the ECEF coordinates for the specified LLA coordinates
     */
    public double[] lla2ecef(double[] lla) {

        double N = e.getA() / Math.sqrt(1 - e.getE2() * Math.pow(Math.sin(lla[0]), 2));

        double x = (N + lla[2]) * Math.cos(lla[0]) * Math.cos(lla[1]);
        double y = (N + lla[2]) * Math.cos(lla[0]) * Math.sin(lla[1]);
        double z = ((1 - e.getE2()) * N + lla[2]) * Math.sin(lla[0]);

        return new double[]{x, y, z};
    }

    public double[] ecef2llaD(double[] ecef) {
        double[] lla = ecef2lla(ecef);
        return new double[]{Math.toDegrees(lla[0]), Math.toDegrees(lla[1]), lla[2]};
    }


    public List<double[]> ecef2lla(List<double[]> coordinates) {
        coordinates.forEach(c -> {
            double[] lla = ecef2llaD(c);
            c[0] = lla[0];
            c[1] = lla[1];
            c[2] = lla[2];
        });
        return coordinates;
    }

}
