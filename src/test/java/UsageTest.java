import org.junit.Test;
import body.ellipsoids.Ellipsoid;
import body.geometry.OblateConic;

import java.util.List;

public class UsageTest {

    @Test
    public void UsageExampleTest() {

        // Define an ellipsoid (WGS84 in this example)
        Ellipsoid e = new Ellipsoid(6378.137, 1 / 298.257223563, 3.986004418E14, 7.2921150E-5);

        // Define an OblateConic object associated with the ellipsoid e
        OblateConic conic = new OblateConic(e);

        // Set the conic drawing method:
        // 0: Sensor attached to the satellite
        // 1: Horizon elevation threshold - bi-section method
        // 2: Horizon elevation threshold - gradient descent method (original)
        conic.setDrawingMethod(1);

        // Set the parameters.
        // Epsilon: elevation threshold over the horizon in degrees
        // Tolerance: the tolerance to find epsilon
        // Segments: the amount of segments of the polygon drawn on the surface of the ellipsoid
        conic.setElevationParams(5, 1E-3, 50);

        // Obtain the list of lat,lon coordinates for a polygon on the ellipsoid, for a given x,y,z
        // position in space, in Kilometers
        List<double[]> coordinates = conic.drawConic(-990.945443, -5817.571039, 3334.217811);

    }

    @Test
    public void ConicDrawingMethodComparisonTest() {

        // Define an ellipsoid (WGS84 in this example)
        Ellipsoid e = new Ellipsoid(6378.137, 1 / 298.257223563, 3.986004418E14, 7.2921150E-5);
        OblateConic conic = new OblateConic(e);

        // Draw with bisection
        System.out.println("Bisection");
        conic.setDrawingMethod(1);
        conic.setElevationParams(5, 1E-3, 4);
        List<double[]> coordinates = conic.drawConic(-990.945443, -5817.571039, 3334.217811);
        coordinates.forEach(c -> System.out.println(c[0] + "," + c[1] + "," + c[2] + "," + c[4] + "," + c[5]));

        // Draw with original method
        System.out.println("Original");
        conic.setDrawingMethod(2);
        coordinates = conic.drawConic(-990.945443, -5817.571039, 3334.217811);
        coordinates.forEach(c -> System.out.println(c[0] + "," + c[1] + "," + c[2] + "," + c[4] + "," + c[5]));

    }

}
