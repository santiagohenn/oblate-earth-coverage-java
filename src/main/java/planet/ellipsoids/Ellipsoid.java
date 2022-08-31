package planet.ellipsoids;

/**
 * This class holds the parameters of an ellipsoid
 */
public class Ellipsoid {

    private final double a;
    private final double a2;
    private final double b;
    private final double b2;
    private final double b4;
    private final double f;
    private final double e;
    private final double e2;
    private final double gravitationalParameter;
    private final double spinRate;

    /**
     * Class constructor
     */
    public Ellipsoid(double semiMajorAxis, double flattening, double gravitationalParameter, double spinRate) {
        this.a = semiMajorAxis;
        this.a2 = a * a;
        this.b = a * (1 - flattening);
        this.b2 = b * b;
        this.b4 = b2 * b2;
        this.f = flattening;
        this.e = Math.sqrt(2 * f - Math.pow(f, 2));
        this.e2 = e * e;
        this.gravitationalParameter = gravitationalParameter;
        this.spinRate = spinRate;
    }

    /**
     * <p> Returns the semi-major axis of the ellipsoid
     * </p>
     * @return the ellipsoid semi-major axis
     */
    public double getA() {
        return a;
    }

    /**
     * <p> Returns the square of the semi-major axis of the ellipsoid
     * </p>
     * @return the ellipsoid's semi-major axis squared
     */
    public double getA2() {
        return a2;
    }

    /**
     * <p> Returns the semi-minor axis of the ellipsoid
     * </p>
     * @return the ellipsoid semi-minor axis
     */
    public double getB() {
        return b;
    }

    /**
     * <p> Returns the square of the semi-minor axis of the ellipsoid
     * </p>
     * @return the ellipsoid's semi-minor axis squared
     */
    public double getB2() {
        return b2;
    }

    /**
     * <p> Returns the 4th power of the semi-minor axis of the ellipsoid
     * </p>
     * @return the fourth power of the semi-minor axis of the ellipsoid
     */
    public double getB4() {
        return b4;
    }

    /**
     * <p> Returns the ellipsoid's flattening
     * </p>
     * @return the flattening of the ellipsoid
     */
    public double getF() {
        return f;
    }

    /**
     * <p> Returns the ellipsoid's eccentricity
     * </p>
     * @return the eccentricity of the ellipsoid
     */
    public double getE() {
        return e;
    }

    /**
     * <p> Returns the square of the ellipsoid's eccentricity
     * </p>
     * @return the eccentricity of the ellipsoid squared
     */
    public double getE2() {
        return e2;
    }

    /**
     * <p> Returns the ellipsoid's gravitational parameter
     * </p>
     * @return the gravitational parameter of the ellipsoid
     */
    public double getGravitationalParameter() {
        return gravitationalParameter;
    }

    /**
     * <p> Returns the ellipsoid's spin rate
     * </p>
     * @return the spin rate of the ellipsoid
     */
    public double getSpinRate() {
        return spinRate;
    }

}
