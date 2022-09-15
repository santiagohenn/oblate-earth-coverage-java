import org.junit.Test;
import planet.ellipsoids.Ellipsoid;
import planet.geometry.OblateConic;

import java.util.ArrayList;
import java.util.List;

public class UsageTest {

    @Test
    public void UsageExample() {

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

        try {
            // Obtain the list of lat,lon coordinates for a polygon on the ellipsoid, for a given x,y,z
            // position in space, in KM
            List<double[]> coordinates = conic.drawConic(-990.945443, -5817.571039, 3334.217811);
        } catch (Exception exc) {
            exc.printStackTrace();
        }

    }

}
