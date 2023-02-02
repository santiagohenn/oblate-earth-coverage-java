package body.geometry;

import body.ellipsoids.Ellipsoid;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class OblateConicTest {

    @Test
    public void ConicTestCartesian() {

        Ellipsoid e = new Ellipsoid(6378.137, 1 / 298.257223563, 3.986004418E14, 7.2921150E-5);
        OblateConic conic = new OblateConic(e);
        conic.setDrawingMethod(1);
        conic.setElevationParams(5, 1E-3, 50);
        List<double[]> coordinates = new ArrayList<>();
        try {
            coordinates = conic.drawConic(-990.945443, -5817.571039, 3334.217811);
        } catch (Exception exc) {
            exc.printStackTrace();
        }
        double delta = 1;
        assertEquals(coordinates.get(0)[0], 781.944574031817, delta);
        assertEquals(coordinates.get(0)[1], -5557.11486252831, delta);
        assertEquals(coordinates.get(0)[2], 3020.9539773423, delta);
        assertEquals(coordinates.get(49)[0], -2577.62844776305, delta);
        assertEquals(coordinates.get(49)[1], -4984.85651724034, delta);
        assertEquals(coordinates.get(49)[2], 3020.9539773423, delta);
        assertEquals(coordinates.get(99)[0], 781.944574031817, delta);
        assertEquals(coordinates.get(99)[1], -5557.11486252831, delta);
        assertEquals(coordinates.get(99)[2], 3020.9539773423, delta);

    }

    @Test
    public void ConicTestCartesian2() {

        Ellipsoid e = new Ellipsoid(6378.137, 1 / 298.257223563, 3.986004418E14, 7.2921150E-5);
        OblateConic conic = new OblateConic(e);
        conic.setDrawingMethod(1);
        conic.setElevationParams(5, 1E-3, 50);
        List<double[]> coordinates = new ArrayList<>();
        try {
            coordinates = conic.drawConic(-2338.7552, -3781.5960, -5125.0175);
        } catch (Exception exc) {
            exc.printStackTrace();
        }
        double delta = 1;
        assertEquals(coordinates.get(0)[0], -632.500818849155, delta);
        assertEquals(coordinates.get(0)[1], -4328.16963282358, delta);
        assertEquals(coordinates.get(0)[2], -4626.36732596462, delta);
        assertEquals(coordinates.get(49)[0], -3589.9006229839, delta);
        assertEquals(coordinates.get(49)[1], -2499.14448830811, delta);
        assertEquals(coordinates.get(49)[2], -4626.36732596462, delta);
        assertEquals(coordinates.get(99)[0], -632.500818849155, delta);
        assertEquals(coordinates.get(99)[1], -4328.16963282358, delta);
        assertEquals(coordinates.get(99)[2], -4626.36732596462, delta);

    }

    @Test
    public void ConicTestLLA() {

        Ellipsoid e = new Ellipsoid(6378.137, 1 / 298.257223563, 3.986004418E14, 7.2921150E-5);
        OblateConic conic = new OblateConic(e);
        conic.setDrawingMethod(1);
        conic.setElevationParams(5, 1E-3, 50);
        List<double[]> coordinates = new ArrayList<>();

        try {
            coordinates = conic.drawLLAConic(-2338.7552, -3781.5960, -5125.0175);
        } catch (Exception exc) {
            exc.printStackTrace();
        }
        double delta = 1;
        assertEquals(coordinates.get(0)[0], -46.7972968355143, delta);
        assertEquals(coordinates.get(0)[1], -98.3141171304491, delta);
        assertEquals(coordinates.get(0)[2], 0, delta);
        assertEquals(coordinates.get(24)[0], -33.4204855504455, delta);
        assertEquals(coordinates.get(24)[1], -121.136471925589, delta);
        assertEquals(coordinates.get(24)[2], 0, delta);
        assertEquals(coordinates.get(49)[0], -46.7972968355143, delta);
        assertEquals(coordinates.get(49)[1], -145.155915787073, delta);
        assertEquals(coordinates.get(49)[2], 0, delta);
        assertEquals(coordinates.get(99)[0], -46.7972968355143, delta);
        assertEquals(coordinates.get(99)[1], -98.3141171304491, delta);
        assertEquals(coordinates.get(99)[2], 0, delta);

    }

}
