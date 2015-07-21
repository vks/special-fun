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

While it is not explicitly stated [on the website of the author][1], Cephes has
been BSD licensed in the past (see [here][2] and [here][3]).

## Related Projects

* [special][4]: Special functions implemented in pure Rust. Has less functions
  implemented and only supports `f64`.



[1]: http://www.moshier.net/#Cephes
[2]: https://lists.debian.org/debian-legal/2004/12/msg00295.html
[3]: https://github.com/jucor/torch-cephes/blob/master/LICENSE.txt
[4]: https://github.com/stainless-steel/special

[status-img]: https://travis-ci.org/vks/special-fun.svg?branch=master
[status-url]: https://travis-ci.org/vks/special-fun
