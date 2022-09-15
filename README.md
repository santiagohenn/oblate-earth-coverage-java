# oblate-body-coverage
This is a JAVA port for the Nunges/Colombo Algorithm to find intersecting cones for a satellite's Field of View. The original MATLAB algorithm can be found [here](https://de.mathworks.com/matlabcentral/fileexchange/81848-coverage-area-determination-considering-an-oblate-earth) and the paper [here](https://doi.org/10.48550/arXiv.1906.12318). 

# Usage

Define an ellipsoid (WGS84 in this example)
```java
Ellipsoid e = new Ellipsoid(6378.137, 1 / 298.257223563, 3.986004418E14, 7.2921150E-5);
```

Define an OblateConic object associated with the ellipsoid e
```java
OblateConic conic = new OblateConic(e);
```

Set the conic drawing method:

- 0: Sensor attached to the satellite
- 1: Horizon elevation threshold - bi-section method
- 2: Horizon elevation threshold - gradient descent method (original from the paper)

Define an OblateConic object associated with the ellipsoid e
```java
conic.setDrawingMethod(1);
```

Set the parameters for drawing method 0.
- HalfAperture: the sensor's conic half aperture in degrees
- Segments: the amount of segments of the polygon drawn on the surface of the ellipsoid
     
```java
conic.setSensorParams(65, 50);
```

Set the parameters for drawing methods 1 and 2.
- Epsilon: elevation threshold over the horizon in degrees
- Tolerance: the tolerance to find epsilon
- Segments: the amount of segments of the polygon drawn on the surface of the ellipsoid
     
```java
conic.setElevationParams(5, 1E-3, 50);
```

Obtain the list of x, y, z ECEF coordinates for a polygon on the ellipsoid, for a given x, y, z ECEF position in space, in Kilometers
```java
List<double[]> coordinates = conic.drawConic(-990.945443, -5817.571039, 3334.217811);
```

Obtain the list of lat,lon coordinates in degrees for a polygon on the ellipsoid, for a given x, y, z ECEF position in space, in Kilometers
```java
List<double[]> coordinates = conic.drawLLAConic(-990.945443, -5817.571039, 3334.217811);
``` 