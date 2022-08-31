import org.junit.Test;
import planet.ellipsoids.Ellipsoid;
import planet.geometry.OblateConic;

import java.util.ArrayList;
import java.util.List;

public class OblateConicTest {

    @Test
    public void ConicTest(){

        Ellipsoid e = new Ellipsoid(6378.137,1 / 298.257223563,3.986004418E14,7.2921150E-5);
        OblateConic conic = new OblateConic(e);
        conic.setDrawingMethod(1);
        conic.setElevationParams(5, 1E-3, 25);
        List<double[]> coordinates = new ArrayList<>();
        try {
            coordinates = conic.drawLLAConic(2000.0, 2000.0, 2000.0);
        } catch (Exception exc) {
            exc.printStackTrace();
        }
        coordinates.forEach(c -> System.out.println(c[0] + "," + c[1] + "," + c[2]));

    }

}
