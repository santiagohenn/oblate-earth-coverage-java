package planet.geometry;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import planet.ellipsoids.Ellipsoid;

public class AuxiliaryVectors {

    public static Vector3D getNPrime(double[][] A_321) {
        return new Vector3D(A_321[2][0], A_321[2][1], A_321[2][2]);
    }

    public static Vector3D getNadir(double x, double y, double z, Ellipsoid e) {

        Vector3D pos = new Vector3D(x, y, z);   // n
        Vector3D gcPos = pos.scalarMultiply(-1.0D / pos.getNorm());

        // Line-of-sight projection onto the Earth surface (nadir)
        return AuxiliaryVectors.getNadir(pos, gcPos, e);

    }

    public static Vector3D getNadir(Vector3D pos, Vector3D gcPos, Ellipsoid e) {

        // Coordinates of the S/C in the Geocentric Inertial Frame
        double XX = pos.getX();
        double YY = pos.getY();
        double ZZ = pos.getZ();

        double nx = gcPos.getX();
        double ny = gcPos.getY();
        double nz = gcPos.getZ();

        double deltaTan = 2 * e.getB4() * nx * ny * XX * YY + 2 * e.getA2() * e.getB2() * nx * nz * XX * ZZ + 2 * e.getA2() * e.getB2() * ny * nz * YY * ZZ
                - e.getB4() * Math.pow(nx, 2) * Math.pow(YY, 2) - e.getA2() * e.getB2() * Math.pow(nx, 2) * Math.pow(ZZ, 2) + e.getA2() * e.getB4() * Math.pow(nx, 2) - e.getB4() * Math.pow(ny, 2) * Math.pow(XX, 2)
                - e.getA2() * e.getB2() * Math.pow(ny, 2) * Math.pow(ZZ, 2) + e.getA2() * e.getB4() * Math.pow(ny, 2) - e.getA2() * e.getB2() * Math.pow(nz, 2) * Math.pow(XX, 2)
                - e.getA2() * e.getB2() * Math.pow(YY, 2) * Math.pow(nz, 2) + Math.pow(e.getA2(), 2) * e.getB2() * Math.pow(nz, 2);

        // Solution of the parametric equation
        var v = -(e.getB2() * nx * XX + e.getB2() * ny * YY + e.getA2() * nz * ZZ);
        var v1 = e.getB2() * Math.pow(nx, 2) + e.getB2() * Math.pow(ny, 2) + e.getA2() * Math.pow(nz, 2);
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

}
