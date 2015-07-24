# special-fun

[![Status][status-img]][status-url]

Special functions for Rust by binding to the [Cephes library][1].

The following families of functions currently have Rust bindings for `f32` and
`f64`:

* Bessel functions
* Beta functions
* Error functions
* Gamma functions
* Hypergeometric functions
* Zeta functions
* Normal probability distribution

Cephes implements a lot more functions that are not yet exposed in the Rust
interface.


## Installing

Cargo is used to build the included Cephes library (which is written in C) and
to create a Rust library that statically links to Cephes.


## License

The bindings and the Cephes library use the BSD license.
The author of Cephes (Stephen Moshier) has agreed to this.


## Related Projects

* [special][2]: Special functions implemented in pure Rust. Has less functions
  implemented and only supports `f64`.


[1]: http://www.moshier.net/#Cephes
[2]: https://github.com/stainless-steel/special

[status-img]: https://travis-ci.org/vks/special-fun.svg?branch=master
[status-url]: https://travis-ci.org/vks/special-fun
