# Apollonius

*2D Apollonius Graphs*

___

An Apollonius graph is also known as an *additively weighted Voronoï diagram*. 
Roughly speaking, this is the Voronoï diagram of a series of weighted points, 
that can be seen as circles (whose radii are the weights). When all the 
weights are equal, this is the ordinary Voronoï diagram.


![](https://raw.githubusercontent.com/stla/Apollonius/main/inst/screenshots/agraph03.gif)

![](https://raw.githubusercontent.com/stla/Apollonius/main/inst/screenshots/enclosedSquare.gif)

___

## Installation

```r
remotes::install_github("stla/gyro")
remotes::install_github("stla/Apollonius")
```

___

## License

This package is provided under the GPL-3 license but it uses the C++ library 
[CGAL](https://www.cgal.org/) which requires a license from the 
[GeometryFactory](https://geometryfactory.com) if you wish to use it for 
commercial purposes.
