(omega-user-data-types)=

## Data Types and Precision

Omega supports all standard data types and uses some specific defined
types to guarantee a specific level of precision. There is only one
user-configurable option for precision. When a specific floating point
precision is not required, we use a Real data type that is, by default,
double precision (8 bytes/64-bit) but if the code is built with a
`-DSINGLE_PRECISION` (see insert link to build system) preprocessor flag,
the default Real becomes single precision (4-byte/32-bit). Users are
encouraged to use the default double precision unless exploring the
performance or accuracy characteristics of single precision.
