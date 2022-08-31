import org.junit.Test;
import planet.ellipsoids.Ellipsoid;
import static org.junit.Assert.assertEquals;

public class EllipsoidTest {

    @Test
    public void WGS84Test() {

        Ellipsoid wgs84 = new Ellipsoid(6378.137,1 / 298.257223563,3.986004418E14,7.2921150E-5);
        assertEquals(wgs84.getE(), 0.081819190842613, 1E-16);
        assertEquals(wgs84.getB(), 6356.7523142, 1E-16);

    }


}
