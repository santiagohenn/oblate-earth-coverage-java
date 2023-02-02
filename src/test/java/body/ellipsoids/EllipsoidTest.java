package body.ellipsoids;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class EllipsoidTest {

    @Test
    public void WGS84Test() {

        Ellipsoid wgs84 = new Ellipsoid(6378.137,1 / 298.257223563,3.986004418E14,7.2921150E-5);
        assertEquals(wgs84.getE(), 0.081819190842613, 1E-12);
        assertEquals(wgs84.getB(), 6356.752314245179, 1E-12);

    }

    @Test
    public void Test(){

    }


}
