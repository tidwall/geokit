// https://github.com/tidwall/geokit.c
// 
// Copyright 2023 Joshua J Baker. All rights reserved.
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file.

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include "geokit.h"

#define EARTH 6378137                // WGS84 equatorial radius in meters
#define FLATTENING (1/298.257223563) // WGS84 flattening

#define RADIANS (M_PI / 180)
#define DEGREES (180 / M_PI)

#define MAX_ITERS 100  // limit number of iterations. sanity.
#define ACCURACY 1e-3  // desired accuracy in meters.

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wold-style-declaration"
#endif

#pragma region geodesic.h
#define GEOD_DLL static
/**
 * \file geodesic.h
 * \brief API for the geodesic routines in C
 *
 * These routines are a simple transcription of the corresponding C++ classes
 * in <a href="https://geographiclib.sourceforge.io"> GeographicLib</a>.  The
 * "class data" is represented by the structs geod_geodesic, geod_geodesicline,
 * geod_polygon and pointers to these objects are passed as initial arguments
 * to the member functions.  Most of the internal comments have been retained.
 * However, in the process of transcription some documentation has been lost
 * and the documentation for the C++ classes, GeographicLib::Geodesic,
 * GeographicLib::GeodesicLine, and GeographicLib::PolygonAreaT, should be
 * consulted.  The C++ code remains the "reference implementation".  Think
 * twice about restructuring the internals of the C code since this may make
 * porting fixes from the C++ code more difficult.
 *
 * Copyright (c) Charles Karney (2012-2022) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 **********************************************************************/

#if !defined(GEODESIC_H)
#define GEODESIC_H 1

/**
 * The major version of the geodesic library.  (This tracks the version of
 * GeographicLib.)
 **********************************************************************/
#define GEODESIC_VERSION_MAJOR 2
/**
 * The minor version of the geodesic library.  (This tracks the version of
 * GeographicLib.)
 **********************************************************************/
#define GEODESIC_VERSION_MINOR 1
/**
 * The patch level of the geodesic library.  (This tracks the version of
 * GeographicLib.)
 **********************************************************************/
#define GEODESIC_VERSION_PATCH 0

/**
 * Pack the version components into a single integer.  Users should not rely on
 * this particular packing of the components of the version number; see the
 * documentation for ::GEODESIC_VERSION, below.
 **********************************************************************/
#define GEODESIC_VERSION_NUM(a,b,c) ((((a) * 10000 + (b)) * 100) + (c))

/**
 * The version of the geodesic library as a single integer, packed as MMmmmmpp
 * where MM is the major version, mmmm is the minor version, and pp is the
 * patch level.  Users should not rely on this particular packing of the
 * components of the version number.  Instead they should use a test such as
 * @code{.c}
   #if GEODESIC_VERSION >= GEODESIC_VERSION_NUM(1,40,0)
   ...
   #endif
 * @endcode
 **********************************************************************/
#define GEODESIC_VERSION \
 GEODESIC_VERSION_NUM(GEODESIC_VERSION_MAJOR, \
                      GEODESIC_VERSION_MINOR, \
                      GEODESIC_VERSION_PATCH)

#if !defined(GEOD_DLL)
#if defined(_MSC_VER) && defined(PROJ_MSVC_DLL_EXPORT)
#define GEOD_DLL __declspec(dllexport)
#elif defined(__GNUC__)
#define GEOD_DLL __attribute__ ((visibility("default")))
#else
#define GEOD_DLL
#endif
#endif

#if defined(PROJ_RENAME_SYMBOLS)
#include "proj_symbol_rename.h"
#endif

#if defined(__cplusplus)
extern "C" {
#endif

  /**
   * The struct containing information about the ellipsoid.  This must be
   * initialized by geod_init() before use.
   **********************************************************************/
  struct geod_geodesic {
    double a;                   /**< the equatorial radius */
    double f;                   /**< the flattening */
    /**< @cond SKIP */
    double f1, e2, ep2, n, b, c2, etol2;
    double A3x[6], C3x[15], C4x[21];
    /**< @endcond */
  };

  /**
   * The struct containing information about a single geodesic.  This must be
   * initialized by geod_lineinit(), geod_directline(), geodx_gendirectline(),
   * or geod_inverseline() before use.
   **********************************************************************/
  struct geod_geodesicline {
    double lat1;                /**< the starting latitude */
    double lon1;                /**< the starting longitude */
    double azi1;                /**< the starting azimuth */
    double a;                   /**< the equatorial radius */
    double f;                   /**< the flattening */
    double salp1;               /**< sine of \e azi1 */
    double calp1;               /**< cosine of \e azi1 */
    double a13;                 /**< arc length to reference point */
    double s13;                 /**< distance to reference point */
    /**< @cond SKIP */
    double b, c2, f1, salp0, calp0, k2,
      ssig1, csig1, dn1, stau1, ctau1, somg1, comg1,
      A1m1, A2m1, A3c, B11, B21, B31, A4, B41;
    double C1a[6+1], C1pa[6+1], C2a[6+1], C3a[6], C4a[6];
    /**< @endcond */
    unsigned caps;              /**< the capabilities */
  };

  /**
   * The struct for accumulating information about a geodesic polygon.  This is
   * used for computing the perimeter and area of a polygon.  This must be
   * initialized by geod_polygon_init() before use.
   **********************************************************************/
  struct geod_polygon {
    double lat;                 /**< the current latitude */
    double lon;                 /**< the current longitude */
    /**< @cond SKIP */
    double lat0;
    double lon0;
    double A[2];
    double P[2];
    int polyline;
    int crossings;
    /**< @endcond */
    unsigned num;               /**< the number of points so far */
  };

  /**
   * Initialize a geod_geodesic object.
   *
   * @param[out] g a pointer to the object to be initialized.
   * @param[in] a the equatorial radius (meters).
   * @param[in] f the flattening.
   **********************************************************************/
  void GEOD_DLL geod_init(struct geod_geodesic* g, double a, double f);

  /**
   * Solve the direct geodesic problem.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] lat1 latitude of point 1 (degrees).
   * @param[in] lon1 longitude of point 1 (degrees).
   * @param[in] azi1 azimuth at point 1 (degrees).
   * @param[in] s12 distance from point 1 to point 2 (meters); it can be
   *   negative.
   * @param[out] plat2 pointer to the latitude of point 2 (degrees).
   * @param[out] plon2 pointer to the longitude of point 2 (degrees).
   * @param[out] pazi2 pointer to the (forward) azimuth at point 2 (degrees).
   *
   * \e g must have been initialized with a call to geod_init().  \e lat1
   * should be in the range [&minus;90&deg;, 90&deg;].  The values of \e lon2
   * and \e azi2 returned are in the range [&minus;180&deg;, 180&deg;].  Any of
   * the "return" arguments \e plat2, etc., may be replaced by 0, if you do not
   * need some quantities computed.
   *
   * If either point is at a pole, the azimuth is defined by keeping the
   * longitude fixed, writing \e lat = &plusmn;(90&deg; &minus; &epsilon;), and
   * taking the limit &epsilon; &rarr; 0+.  An arc length greater that 180&deg;
   * signifies a geodesic which is not a shortest path.  (For a prolate
   * ellipsoid, an additional condition is necessary for a shortest path: the
   * longitudinal extent must not exceed of 180&deg;.)
   *
   * Example, determine the point 10000 km NE of JFK:
   @code{.c}
   struct geod_geodesic g;
   double lat, lon;
   geod_init(&g, 6378137, 1/298.257223563);
   geod_direct(&g, 40.64, -73.78, 45.0, 10e6, &lat, &lon, 0);
   printf("%.5f %.5f\n", lat, lon);
   @endcode
   **********************************************************************/
  void GEOD_DLL geod_direct(const struct geod_geodesic* g,
                            double lat1, double lon1, double azi1, double s12,
                            double* plat2, double* plon2, double* pazi2);

  /**
   * The general direct geodesic problem.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] lat1 latitude of point 1 (degrees).
   * @param[in] lon1 longitude of point 1 (degrees).
   * @param[in] azi1 azimuth at point 1 (degrees).
   * @param[in] flags bitor'ed combination of ::geod_flags; \e flags &
   *   ::GEOD_ARCMODE determines the meaning of \e s12_a12 and \e flags &
   *   ::GEOD_LONG_UNROLL "unrolls" \e lon2.
   * @param[in] s12_a12 if \e flags & ::GEOD_ARCMODE is 0, this is the distance
   *   from point 1 to point 2 (meters); otherwise it is the arc length
   *   from point 1 to point 2 (degrees); it can be negative.
   * @param[out] plat2 pointer to the latitude of point 2 (degrees).
   * @param[out] plon2 pointer to the longitude of point 2 (degrees).
   * @param[out] pazi2 pointer to the (forward) azimuth at point 2 (degrees).
   * @param[out] ps12 pointer to the distance from point 1 to point 2
   *   (meters).
   * @param[out] pm12 pointer to the reduced length of geodesic (meters).
   * @param[out] pM12 pointer to the geodesic scale of point 2 relative to
   *   point 1 (dimensionless).
   * @param[out] pM21 pointer to the geodesic scale of point 1 relative to
   *   point 2 (dimensionless).
   * @param[out] pS12 pointer to the area under the geodesic
   *   (meters<sup>2</sup>).
   * @return \e a12 arc length from point 1 to point 2 (degrees).
   *
   * \e g must have been initialized with a call to geod_init().  \e lat1
   * should be in the range [&minus;90&deg;, 90&deg;].  The function value \e
   * a12 equals \e s12_a12 if \e flags & ::GEOD_ARCMODE.  Any of the "return"
   * arguments, \e plat2, etc., may be replaced by 0, if you do not need some
   * quantities computed.
   *
   * With \e flags & ::GEOD_LONG_UNROLL bit set, the longitude is "unrolled" so
   * that the quantity \e lon2 &minus; \e lon1 indicates how many times and in
   * what sense the geodesic encircles the ellipsoid.
   **********************************************************************/
  double GEOD_DLL geod_gendirect(const struct geod_geodesic* g,
                                 double lat1, double lon1, double azi1,
                                 unsigned flags, double s12_a12,
                                 double* plat2, double* plon2, double* pazi2,
                                 double* ps12, double* pm12,
                                 double* pM12, double* pM21,
                                 double* pS12);

  /**
   * Solve the inverse geodesic problem.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] lat1 latitude of point 1 (degrees).
   * @param[in] lon1 longitude of point 1 (degrees).
   * @param[in] lat2 latitude of point 2 (degrees).
   * @param[in] lon2 longitude of point 2 (degrees).
   * @param[out] ps12 pointer to the distance from point 1 to point 2
   *   (meters).
   * @param[out] pazi1 pointer to the azimuth at point 1 (degrees).
   * @param[out] pazi2 pointer to the (forward) azimuth at point 2 (degrees).
   *
   * \e g must have been initialized with a call to geod_init().  \e lat1 and
   * \e lat2 should be in the range [&minus;90&deg;, 90&deg;].  The values of
   * \e azi1 and \e azi2 returned are in the range [&minus;180&deg;, 180&deg;].
   * Any of the "return" arguments, \e ps12, etc., may be replaced by 0, if you
   * do not need some quantities computed.
   *
   * If either point is at a pole, the azimuth is defined by keeping the
   * longitude fixed, writing \e lat = &plusmn;(90&deg; &minus; &epsilon;), and
   * taking the limit &epsilon; &rarr; 0+.
   *
   * The solution to the inverse problem is found using Newton's method.  If
   * this fails to converge (this is very unlikely in geodetic applications
   * but does occur for very eccentric ellipsoids), then the bisection method
   * is used to refine the solution.
   *
   * Example, determine the distance between JFK and Singapore Changi Airport:
   @code{.c}
   struct geod_geodesic g;
   double s12;
   geod_init(&g, 6378137, 1/298.257223563);
   geod_inverse(&g, 40.64, -73.78, 1.36, 103.99, &s12, 0, 0);
   printf("%.3f\n", s12);
   @endcode
   **********************************************************************/
  void GEOD_DLL geod_inverse(const struct geod_geodesic* g,
                             double lat1, double lon1,
                             double lat2, double lon2,
                             double* ps12, double* pazi1, double* pazi2);

  /**
   * The general inverse geodesic calculation.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] lat1 latitude of point 1 (degrees).
   * @param[in] lon1 longitude of point 1 (degrees).
   * @param[in] lat2 latitude of point 2 (degrees).
   * @param[in] lon2 longitude of point 2 (degrees).
   * @param[out] ps12 pointer to the distance from point 1 to point 2
   *  (meters).
   * @param[out] pazi1 pointer to the azimuth at point 1 (degrees).
   * @param[out] pazi2 pointer to the (forward) azimuth at point 2 (degrees).
   * @param[out] pm12 pointer to the reduced length of geodesic (meters).
   * @param[out] pM12 pointer to the geodesic scale of point 2 relative to
   *   point 1 (dimensionless).
   * @param[out] pM21 pointer to the geodesic scale of point 1 relative to
   *   point 2 (dimensionless).
   * @param[out] pS12 pointer to the area under the geodesic
   *   (meters<sup>2</sup>).
   * @return \e a12 arc length from point 1 to point 2 (degrees).
   *
   * \e g must have been initialized with a call to geod_init().  \e lat1 and
   * \e lat2 should be in the range [&minus;90&deg;, 90&deg;].  Any of the
   * "return" arguments \e ps12, etc., may be replaced by 0, if you do not need
   * some quantities computed.
   **********************************************************************/
  double GEOD_DLL geod_geninverse(const struct geod_geodesic* g,
                                  double lat1, double lon1,
                                  double lat2, double lon2,
                                  double* ps12, double* pazi1, double* pazi2,
                                  double* pm12, double* pM12, double* pM21,
                                  double* pS12);

  /**
   * Initialize a geod_geodesicline object.
   *
   * @param[out] l a pointer to the object to be initialized.
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] lat1 latitude of point 1 (degrees).
   * @param[in] lon1 longitude of point 1 (degrees).
   * @param[in] azi1 azimuth at point 1 (degrees).
   * @param[in] caps bitor'ed combination of ::geod_mask values specifying the
   *   capabilities the geod_geodesicline object should possess, i.e., which
   *   quantities can be returned in calls to geod_position() and
   *   geod_genposition().
   *
   * \e g must have been initialized with a call to geod_init().  \e lat1
   * should be in the range [&minus;90&deg;, 90&deg;].
   *
   * The ::geod_mask values are:
   * - \e caps |= ::GEOD_LATITUDE for the latitude \e lat2; this is
   *   added automatically,
   * - \e caps |= ::GEOD_LONGITUDE for the latitude \e lon2,
   * - \e caps |= ::GEOD_AZIMUTH for the latitude \e azi2; this is
   *   added automatically,
   * - \e caps |= ::GEOD_DISTANCE for the distance \e s12,
   * - \e caps |= ::GEOD_REDUCEDLENGTH for the reduced length \e m12,
   * - \e caps |= ::GEOD_GEODESICSCALE for the geodesic scales \e M12
   *   and \e M21,
   * - \e caps |= ::GEOD_AREA for the area \e S12,
   * - \e caps |= ::GEOD_DISTANCE_IN permits the length of the
   *   geodesic to be given in terms of \e s12; without this capability the
   *   length can only be specified in terms of arc length.
   * .
   * A value of \e caps = 0 is treated as ::GEOD_LATITUDE | ::GEOD_LONGITUDE |
   * ::GEOD_AZIMUTH | ::GEOD_DISTANCE_IN (to support the solution of the
   * "standard" direct problem).
   *
   * When initialized by this function, point 3 is undefined (l->s13 = l->a13 =
   * NaN).
   **********************************************************************/
  void GEOD_DLL geod_lineinit(struct geod_geodesicline* l,
                              const struct geod_geodesic* g,
                              double lat1, double lon1, double azi1,
                              unsigned caps);

  /**
   * Initialize a geod_geodesicline object in terms of the direct geodesic
   * problem.
   *
   * @param[out] l a pointer to the object to be initialized.
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] lat1 latitude of point 1 (degrees).
   * @param[in] lon1 longitude of point 1 (degrees).
   * @param[in] azi1 azimuth at point 1 (degrees).
   * @param[in] s12 distance from point 1 to point 2 (meters); it can be
   *   negative.
   * @param[in] caps bitor'ed combination of ::geod_mask values specifying the
   *   capabilities the geod_geodesicline object should possess, i.e., which
   *   quantities can be returned in calls to geod_position() and
   *   geod_genposition().
   *
   * This function sets point 3 of the geod_geodesicline to correspond to point
   * 2 of the direct geodesic problem.  See geod_lineinit() for more
   * information.
   **********************************************************************/
  void GEOD_DLL geod_directline(struct geod_geodesicline* l,
                                const struct geod_geodesic* g,
                                double lat1, double lon1,
                                double azi1, double s12,
                                unsigned caps);

  /**
   * Initialize a geod_geodesicline object in terms of the direct geodesic
   * problem specified in terms of either distance or arc length.
   *
   * @param[out] l a pointer to the object to be initialized.
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] lat1 latitude of point 1 (degrees).
   * @param[in] lon1 longitude of point 1 (degrees).
   * @param[in] azi1 azimuth at point 1 (degrees).
   * @param[in] flags either ::GEOD_NOFLAGS or ::GEOD_ARCMODE to determining
   *   the meaning of the \e s12_a12.
   * @param[in] s12_a12 if \e flags = ::GEOD_NOFLAGS, this is the distance
   *   from point 1 to point 2 (meters); if \e flags = ::GEOD_ARCMODE, it is
   *   the arc length from point 1 to point 2 (degrees); it can be
   *   negative.
   * @param[in] caps bitor'ed combination of ::geod_mask values specifying the
   *   capabilities the geod_geodesicline object should possess, i.e., which
   *   quantities can be returned in calls to geod_position() and
   *   geod_genposition().
   *
   * This function sets point 3 of the geod_geodesicline to correspond to point
   * 2 of the direct geodesic problem.  See geod_lineinit() for more
   * information.
   **********************************************************************/
  void GEOD_DLL geod_gendirectline(struct geod_geodesicline* l,
                                   const struct geod_geodesic* g,
                                   double lat1, double lon1, double azi1,
                                   unsigned flags, double s12_a12,
                                   unsigned caps);

  /**
   * Initialize a geod_geodesicline object in terms of the inverse geodesic
   * problem.
   *
   * @param[out] l a pointer to the object to be initialized.
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] lat1 latitude of point 1 (degrees).
   * @param[in] lon1 longitude of point 1 (degrees).
   * @param[in] lat2 latitude of point 2 (degrees).
   * @param[in] lon2 longitude of point 2 (degrees).
   * @param[in] caps bitor'ed combination of ::geod_mask values specifying the
   *   capabilities the geod_geodesicline object should possess, i.e., which
   *   quantities can be returned in calls to geod_position() and
   *   geod_genposition().
   *
   * This function sets point 3 of the geod_geodesicline to correspond to point
   * 2 of the inverse geodesic problem.  See geod_lineinit() for more
   * information.
   **********************************************************************/
  void GEOD_DLL geod_inverseline(struct geod_geodesicline* l,
                                 const struct geod_geodesic* g,
                                 double lat1, double lon1,
                                 double lat2, double lon2,
                                 unsigned caps);

  /**
   * Compute the position along a geod_geodesicline.
   *
   * @param[in] l a pointer to the geod_geodesicline object specifying the
   *   geodesic line.
   * @param[in] s12 distance from point 1 to point 2 (meters); it can be
   *   negative.
   * @param[out] plat2 pointer to the latitude of point 2 (degrees).
   * @param[out] plon2 pointer to the longitude of point 2 (degrees); requires
   *   that \e l was initialized with \e caps |= ::GEOD_LONGITUDE.
   * @param[out] pazi2 pointer to the (forward) azimuth at point 2 (degrees).
   *
   * \e l must have been initialized with a call, e.g., to geod_lineinit(),
   * with \e caps |= ::GEOD_DISTANCE_IN (or \e caps = 0).  The values of \e
   * lon2 and \e azi2 returned are in the range [&minus;180&deg;, 180&deg;].
   * Any of the "return" arguments \e plat2, etc., may be replaced by 0, if you
   * do not need some quantities computed.
   *
   * Example, compute way points between JFK and Singapore Changi Airport
   * the "obvious" way using geod_direct():
   @code{.c}
   struct geod_geodesic g;
   double s12, azi1, lat[101], lon[101];
   int i;
   geod_init(&g, 6378137, 1/298.257223563);
   geod_inverse(&g, 40.64, -73.78, 1.36, 103.99, &s12, &azi1, 0);
   for (i = 0; i < 101; ++i) {
     geod_direct(&g, 40.64, -73.78, azi1, i * s12 * 0.01, lat + i, lon + i, 0);
     printf("%.5f %.5f\n", lat[i], lon[i]);
   }
   @endcode
   * A faster way using geod_position():
   @code{.c}
   struct geod_geodesic g;
   struct geod_geodesicline l;
   double lat[101], lon[101];
   int i;
   geod_init(&g, 6378137, 1/298.257223563);
   geod_inverseline(&l, &g, 40.64, -73.78, 1.36, 103.99, 0);
   for (i = 0; i <= 100; ++i) {
     geod_position(&l, i * l.s13 * 0.01, lat + i, lon + i, 0);
     printf("%.5f %.5f\n", lat[i], lon[i]);
   }
   @endcode
   **********************************************************************/
  void GEOD_DLL geod_position(const struct geod_geodesicline* l, double s12,
                              double* plat2, double* plon2, double* pazi2);

  /**
   * The general position function.
   *
   * @param[in] l a pointer to the geod_geodesicline object specifying the
   *   geodesic line.
   * @param[in] flags bitor'ed combination of ::geod_flags; \e flags &
   *   ::GEOD_ARCMODE determines the meaning of \e s12_a12 and \e flags &
   *   ::GEOD_LONG_UNROLL "unrolls" \e lon2; if \e flags & ::GEOD_ARCMODE is 0,
   *   then \e l must have been initialized with \e caps |= ::GEOD_DISTANCE_IN.
   * @param[in] s12_a12 if \e flags & ::GEOD_ARCMODE is 0, this is the
   *   distance from point 1 to point 2 (meters); otherwise it is the
   *   arc length from point 1 to point 2 (degrees); it can be
   *   negative.
   * @param[out] plat2 pointer to the latitude of point 2 (degrees).
   * @param[out] plon2 pointer to the longitude of point 2 (degrees); requires
   *   that \e l was initialized with \e caps |= ::GEOD_LONGITUDE.
   * @param[out] pazi2 pointer to the (forward) azimuth at point 2 (degrees).
   * @param[out] ps12 pointer to the distance from point 1 to point 2
   *   (meters); requires that \e l was initialized with \e caps |=
   *   ::GEOD_DISTANCE.
   * @param[out] pm12 pointer to the reduced length of geodesic (meters);
   *   requires that \e l was initialized with \e caps |= ::GEOD_REDUCEDLENGTH.
   * @param[out] pM12 pointer to the geodesic scale of point 2 relative to
   *   point 1 (dimensionless); requires that \e l was initialized with \e caps
   *   |= ::GEOD_GEODESICSCALE.
   * @param[out] pM21 pointer to the geodesic scale of point 1 relative to
   *   point 2 (dimensionless); requires that \e l was initialized with \e caps
   *   |= ::GEOD_GEODESICSCALE.
   * @param[out] pS12 pointer to the area under the geodesic
   *   (meters<sup>2</sup>); requires that \e l was initialized with \e caps |=
   *   ::GEOD_AREA.
   * @return \e a12 arc length from point 1 to point 2 (degrees).
   *
   * \e l must have been initialized with a call to geod_lineinit() with \e
   * caps |= ::GEOD_DISTANCE_IN.  The value \e azi2 returned is in the range
   * [&minus;180&deg;, 180&deg;].  Any of the "return" arguments \e plat2,
   * etc., may be replaced by 0, if you do not need some quantities
   * computed.  Requesting a value which \e l is not capable of computing
   * is not an error; the corresponding argument will not be altered.
   *
   * With \e flags & ::GEOD_LONG_UNROLL bit set, the longitude is "unrolled" so
   * that the quantity \e lon2 &minus; \e lon1 indicates how many times and in
   * what sense the geodesic encircles the ellipsoid.
   *
   * Example, compute way points between JFK and Singapore Changi Airport using
   * geod_genposition().  In this example, the points are evenly spaced in arc
   * length (and so only approximately equally spaced in distance).  This is
   * faster than using geod_position() and would be appropriate if drawing the
   * path on a map.
   @code{.c}
   struct geod_geodesic g;
   struct geod_geodesicline l;
   double lat[101], lon[101];
   int i;
   geod_init(&g, 6378137, 1/298.257223563);
   geod_inverseline(&l, &g, 40.64, -73.78, 1.36, 103.99,
                    GEOD_LATITUDE | GEOD_LONGITUDE);
   for (i = 0; i <= 100; ++i) {
     geod_genposition(&l, GEOD_ARCMODE, i * l.a13 * 0.01,
                      lat + i, lon + i, 0, 0, 0, 0, 0, 0);
     printf("%.5f %.5f\n", lat[i], lon[i]);
   }
   @endcode
   **********************************************************************/
  double GEOD_DLL geod_genposition(const struct geod_geodesicline* l,
                                   unsigned flags, double s12_a12,
                                   double* plat2, double* plon2, double* pazi2,
                                   double* ps12, double* pm12,
                                   double* pM12, double* pM21,
                                   double* pS12);

  /**
   * Specify position of point 3 in terms of distance.
   *
   * @param[in,out] l a pointer to the geod_geodesicline object.
   * @param[in] s13 the distance from point 1 to point 3 (meters); it
   *   can be negative.
   *
   * This is only useful if the geod_geodesicline object has been constructed
   * with \e caps |= ::GEOD_DISTANCE_IN.
   **********************************************************************/
  void GEOD_DLL geod_setdistance(struct geod_geodesicline* l, double s13);

  /**
   * Specify position of point 3 in terms of either distance or arc length.
   *
   * @param[in,out] l a pointer to the geod_geodesicline object.
   * @param[in] flags either ::GEOD_NOFLAGS or ::GEOD_ARCMODE to determining
   *   the meaning of the \e s13_a13.
   * @param[in] s13_a13 if \e flags = ::GEOD_NOFLAGS, this is the distance
   *   from point 1 to point 3 (meters); if \e flags = ::GEOD_ARCMODE, it is
   *   the arc length from point 1 to point 3 (degrees); it can be
   *   negative.
   *
   * If flags = ::GEOD_NOFLAGS, this calls geod_setdistance().  If flags =
   * ::GEOD_ARCMODE, the \e s13 is only set if the geod_geodesicline object has
   * been constructed with \e caps |= ::GEOD_DISTANCE.
   **********************************************************************/
  void GEOD_DLL geod_gensetdistance(struct geod_geodesicline* l,
                                    unsigned flags, double s13_a13);

  /**
   * Initialize a geod_polygon object.
   *
   * @param[out] p a pointer to the object to be initialized.
   * @param[in] polylinep non-zero if a polyline instead of a polygon.
   *
   * If \e polylinep is zero, then the sequence of vertices and edges added by
   * geod_polygon_addpoint() and geod_polygon_addedge() define a polygon and
   * the perimeter and area are returned by geod_polygon_compute().  If \e
   * polylinep is non-zero, then the vertices and edges define a polyline and
   * only the perimeter is returned by geod_polygon_compute().
   *
   * The area and perimeter are accumulated at two times the standard floating
   * point precision to guard against the loss of accuracy with many-sided
   * polygons.  At any point you can ask for the perimeter and area so far.
   *
   * An example of the use of this function is given in the documentation for
   * geod_polygon_compute().
   **********************************************************************/
  void GEOD_DLL geod_polygon_init(struct geod_polygon* p, int polylinep);

  /**
   * Clear the polygon, allowing a new polygon to be started.
   *
   * @param[in,out] p a pointer to the object to be cleared.
   **********************************************************************/
  void GEOD_DLL geod_polygon_clear(struct geod_polygon* p);

  /**
   * Add a point to the polygon or polyline.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in,out] p a pointer to the geod_polygon object specifying the
   *   polygon.
   * @param[in] lat the latitude of the point (degrees).
   * @param[in] lon the longitude of the point (degrees).
   *
   * \e g and \e p must have been initialized with calls to geod_init() and
   * geod_polygon_init(), respectively.  The same \e g must be used for all the
   * points and edges in a polygon.  \e lat should be in the range
   * [&minus;90&deg;, 90&deg;].
   *
   * An example of the use of this function is given in the documentation for
   * geod_polygon_compute().
   **********************************************************************/
  void GEOD_DLL geod_polygon_addpoint(const struct geod_geodesic* g,
                                      struct geod_polygon* p,
                                      double lat, double lon);

  /**
   * Add an edge to the polygon or polyline.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in,out] p a pointer to the geod_polygon object specifying the
   *   polygon.
   * @param[in] azi azimuth at current point (degrees).
   * @param[in] s distance from current point to next point (meters).
   *
   * \e g and \e p must have been initialized with calls to geod_init() and
   * geod_polygon_init(), respectively.  The same \e g must be used for all the
   * points and edges in a polygon.  This does nothing if no points have been
   * added yet.  The \e lat and \e lon fields of \e p give the location of the
   * new vertex.
   **********************************************************************/
  void GEOD_DLL geod_polygon_addedge(const struct geod_geodesic* g,
                                     struct geod_polygon* p,
                                     double azi, double s);

  /**
   * Return the results for a polygon.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] p a pointer to the geod_polygon object specifying the polygon.
   * @param[in] reverse if non-zero then clockwise (instead of
   *   counter-clockwise) traversal counts as a positive area.
   * @param[in] sign if non-zero then return a signed result for the area if
   *   the polygon is traversed in the "wrong" direction instead of returning
   *   the area for the rest of the earth.
   * @param[out] pA pointer to the area of the polygon (meters<sup>2</sup>);
   *   only set if \e polyline is non-zero in the call to geod_polygon_init().
   * @param[out] pP pointer to the perimeter of the polygon or length of the
   *   polyline (meters).
   * @return the number of points.
   *
   * The area and perimeter are accumulated at two times the standard floating
   * point precision to guard against the loss of accuracy with many-sided
   * polygons.  Arbitrarily complex polygons are allowed.  In the case of
   * self-intersecting polygons the area is accumulated "algebraically", e.g.,
   * the areas of the 2 loops in a figure-8 polygon will partially cancel.
   * There's no need to "close" the polygon by repeating the first vertex.  Set
   * \e pA or \e pP to zero, if you do not want the corresponding quantity
   * returned.
   *
   * More points can be added to the polygon after this call.
   *
   * Example, compute the perimeter and area of the geodesic triangle with
   * vertices (0&deg;N,0&deg;E), (0&deg;N,90&deg;E), (90&deg;N,0&deg;E).
   @code{.c}
   double A, P;
   int n;
   struct geod_geodesic g;
   struct geod_polygon p;
   geod_init(&g, 6378137, 1/298.257223563);
   geod_polygon_init(&p, 0);

   geod_polygon_addpoint(&g, &p,  0,  0);
   geod_polygon_addpoint(&g, &p,  0, 90);
   geod_polygon_addpoint(&g, &p, 90,  0);
   n = geod_polygon_compute(&g, &p, 0, 1, &A, &P);
   printf("%d %.8f %.3f\n", n, P, A);
   @endcode
   **********************************************************************/
  unsigned GEOD_DLL geod_polygon_compute(const struct geod_geodesic* g,
                                         const struct geod_polygon* p,
                                         int reverse, int sign,
                                         double* pA, double* pP);

  /**
   * Return the results assuming a tentative final test point is added;
   * however, the data for the test point is not saved.  This lets you report a
   * running result for the perimeter and area as the user moves the mouse
   * cursor.  Ordinary floating point arithmetic is used to accumulate the data
   * for the test point; thus the area and perimeter returned are less accurate
   * than if geod_polygon_addpoint() and geod_polygon_compute() are used.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] p a pointer to the geod_polygon object specifying the polygon.
   * @param[in] lat the latitude of the test point (degrees).
   * @param[in] lon the longitude of the test point (degrees).
   * @param[in] reverse if non-zero then clockwise (instead of
   *   counter-clockwise) traversal counts as a positive area.
   * @param[in] sign if non-zero then return a signed result for the area if
   *   the polygon is traversed in the "wrong" direction instead of returning
   *   the area for the rest of the earth.
   * @param[out] pA pointer to the area of the polygon (meters<sup>2</sup>);
   *   only set if \e polyline is non-zero in the call to geod_polygon_init().
   * @param[out] pP pointer to the perimeter of the polygon or length of the
   *   polyline (meters).
   * @return the number of points.
   *
   * \e lat should be in the range [&minus;90&deg;, 90&deg;].
   **********************************************************************/
  unsigned GEOD_DLL geod_polygon_testpoint(const struct geod_geodesic* g,
                                           const struct geod_polygon* p,
                                           double lat, double lon,
                                           int reverse, int sign,
                                           double* pA, double* pP);

  /**
   * Return the results assuming a tentative final test point is added via an
   * azimuth and distance; however, the data for the test point is not saved.
   * This lets you report a running result for the perimeter and area as the
   * user moves the mouse cursor.  Ordinary floating point arithmetic is used
   * to accumulate the data for the test point; thus the area and perimeter
   * returned are less accurate than if geod_polygon_addedge() and
   * geod_polygon_compute() are used.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] p a pointer to the geod_polygon object specifying the polygon.
   * @param[in] azi azimuth at current point (degrees).
   * @param[in] s distance from current point to final test point (meters).
   * @param[in] reverse if non-zero then clockwise (instead of
   *   counter-clockwise) traversal counts as a positive area.
   * @param[in] sign if non-zero then return a signed result for the area if
   *   the polygon is traversed in the "wrong" direction instead of returning
   *   the area for the rest of the earth.
   * @param[out] pA pointer to the area of the polygon (meters<sup>2</sup>);
   *   only set if \e polyline is non-zero in the call to geod_polygon_init().
   * @param[out] pP pointer to the perimeter of the polygon or length of the
   *   polyline (meters).
   * @return the number of points.
   **********************************************************************/
  unsigned GEOD_DLL geod_polygon_testedge(const struct geod_geodesic* g,
                                          const struct geod_polygon* p,
                                          double azi, double s,
                                          int reverse, int sign,
                                          double* pA, double* pP);

  /**
   * A simple interface for computing the area of a geodesic polygon.
   *
   * @param[in] g a pointer to the geod_geodesic object specifying the
   *   ellipsoid.
   * @param[in] lats an array of latitudes of the polygon vertices (degrees).
   * @param[in] lons an array of longitudes of the polygon vertices (degrees).
   * @param[in] n the number of vertices.
   * @param[out] pA pointer to the area of the polygon (meters<sup>2</sup>).
   * @param[out] pP pointer to the perimeter of the polygon (meters).
   *
   * \e lats should be in the range [&minus;90&deg;, 90&deg;].
   *
   * Arbitrarily complex polygons are allowed.  In the case self-intersecting
   * of polygons the area is accumulated "algebraically", e.g., the areas of
   * the 2 loops in a figure-8 polygon will partially cancel.  There's no need
   * to "close" the polygon by repeating the first vertex.  The area returned
   * is signed with counter-clockwise traversal being treated as positive.
   *
   * Example, compute the area of Antarctica:
   @code{.c}
   double
     lats[] = {-72.9, -71.9, -74.9, -74.3, -77.5, -77.4, -71.7, -65.9, -65.7,
               -66.6, -66.9, -69.8, -70.0, -71.0, -77.3, -77.9, -74.7},
     lons[] = {-74, -102, -102, -131, -163, 163, 172, 140, 113,
                88, 59, 25, -4, -14, -33, -46, -61};
   struct geod_geodesic g;
   double A, P;
   geod_init(&g, 6378137, 1/298.257223563);
   geod_polygonarea(&g, lats, lons, (sizeof lats) / (sizeof lats[0]), &A, &P);
   printf("%.0f %.2f\n", A, P);
   @endcode
   **********************************************************************/
  void GEOD_DLL geod_polygonarea(const struct geod_geodesic* g,
                                 double lats[], double lons[], int n,
                                 double* pA, double* pP);

  /**
   * mask values for the \e caps argument to geod_lineinit().
   **********************************************************************/
  enum geod_mask {
    GEOD_NONE         = 0U,                    /**< Calculate nothing */
    GEOD_LATITUDE     = 1U<<7  | 0U,           /**< Calculate latitude */
    GEOD_LONGITUDE    = 1U<<8  | 1U<<3,        /**< Calculate longitude */
    GEOD_AZIMUTH      = 1U<<9  | 0U,           /**< Calculate azimuth */
    GEOD_DISTANCE     = 1U<<10 | 1U<<0,        /**< Calculate distance */
    GEOD_DISTANCE_IN  = 1U<<11 | 1U<<0 | 1U<<1,/**< Allow distance as input  */
    GEOD_REDUCEDLENGTH= 1U<<12 | 1U<<0 | 1U<<2,/**< Calculate reduced length */
    GEOD_GEODESICSCALE= 1U<<13 | 1U<<0 | 1U<<2,/**< Calculate geodesic scale */
    GEOD_AREA         = 1U<<14 | 1U<<4,        /**< Calculate reduced length */
    GEOD_ALL          = 0x7F80U| 0x1FU         /**< Calculate everything */
  };

  /**
   * flag values for the \e flags argument to geod_gendirect() and
   * geod_genposition()
   **********************************************************************/
  enum geod_flags {
    GEOD_NOFLAGS      = 0U,     /**< No flags */
    GEOD_ARCMODE      = 1U<<0,  /**< Position given in terms of arc distance */
    GEOD_LONG_UNROLL  = 1U<<15  /**< Unroll the longitude */
  };

#if defined(__cplusplus)
}
#endif

#endif
#pragma endregion geodesic.h
#pragma region geodesic.c
// Using v2.1 (2023-01-04).

/**
 * \file geodesic.c
 * \brief Implementation of the geodesic routines in C
 *
 * For the full documentation see geodesic.h.
 **********************************************************************/

/** @cond SKIP */

/*
 * This is a C implementation of the geodesic algorithms described in
 *
 *   C. F. F. Karney,
 *   Algorithms for geodesics,
 *   J. Geodesy <b>87</b>, 43--55 (2013);
 *   https://doi.org/10.1007/s00190-012-0578-z
 *   Addenda: https://geographiclib.sourceforge.io/geod-addenda.html
 *
 * See the comments in geodesic.h for documentation.
 *
 * Copyright (c) Charles Karney (2012-2022) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * https://geographiclib.sourceforge.io/
 */

#include <math.h>
#include <float.h>

#if !defined(__cplusplus)
#define nullptr 0
#endif

#define GEOGRAPHICLIB_GEODESIC_ORDER 6
#define nA1   GEOGRAPHICLIB_GEODESIC_ORDER
#define nC1   GEOGRAPHICLIB_GEODESIC_ORDER
#define nC1p  GEOGRAPHICLIB_GEODESIC_ORDER
#define nA2   GEOGRAPHICLIB_GEODESIC_ORDER
#define nC2   GEOGRAPHICLIB_GEODESIC_ORDER
#define nA3   GEOGRAPHICLIB_GEODESIC_ORDER
#define nA3x  nA3
#define nC3   GEOGRAPHICLIB_GEODESIC_ORDER
#define nC3x  ((nC3 * (nC3 - 1)) / 2)
#define nC4   GEOGRAPHICLIB_GEODESIC_ORDER
#define nC4x  ((nC4 * (nC4 + 1)) / 2)
#define nC    (GEOGRAPHICLIB_GEODESIC_ORDER + 1)

typedef int boolx;
enum booly { FALSE = 0, TRUE = 1 };
/* qd = quarter turn / degree
 * hd = half turn / degree
 * td = full turn / degree */
enum dms { qd = 90, hd = 2 * qd, td = 2 * hd };

static unsigned init = 0;
static unsigned digits, maxit1, maxit2;
static double epsilon, realmin, pi, degree, NaN,
  tiny, tol0, tol1, tol2, tolb, xthresh;

static void Init(void) {
  if (!init) {
    digits = DBL_MANT_DIG;
    epsilon = DBL_EPSILON;
    realmin = DBL_MIN;
#if defined(M_PI)
    pi = M_PI;
#else
    pi = atan2(0.0, -1.0);
#endif
    maxit1 = 20;
    maxit2 = maxit1 + digits + 10;
    tiny = sqrt(realmin);
    tol0 = epsilon;
    /* Increase multiplier in defn of tol1 from 100 to 200 to fix inverse case
     * 52.784459512564 0 -52.784459512563990912 179.634407464943777557
     * which otherwise failed for Visual Studio 10 (Release and Debug) */
    tol1 = 200 * tol0;
    tol2 = sqrt(tol0);
    /* Check on bisection interval */
    tolb = tol0;
    xthresh = 1000 * tol2;
    degree = pi/hd;
    NaN = nan("0");
    init = 1;
  }
}

enum captype {
  CAP_NONE = 0U,
  CAP_C1   = 1U<<0,
  CAP_C1p  = 1U<<1,
  CAP_C2   = 1U<<2,
  CAP_C3   = 1U<<3,
  CAP_C4   = 1U<<4,
  CAP_ALL  = 0x1FU,
  OUT_ALL  = 0x7F80U
};

static double sq(double x) { return x * x; }

static double sumx(double u, double v, double* t) {
  volatile double s = u + v;
  volatile double up = s - v;
  volatile double vpp = s - up;
  up -= u;
  vpp -= v;
  if (t) *t = s != 0 ? 0 - (up + vpp) : s;
  /* error-free sum:
   * u + v =       s      + t
   *       = round(u + v) + t */
  return s;
}

static double polyvalx(int N, const double p[], double x) {
  double y = N < 0 ? 0 : *p++;
  while (--N >= 0) y = y * x + *p++;
  return y;
}

static void swapx(double* x, double* y)
{ double t = *x; *x = *y; *y = t; }

static void norm2(double* sinx, double* cosx) {
#if defined(_MSC_VER) && defined(_M_IX86)
  /* hypot for Visual Studio (A=win32) fails monotonicity, e.g., with
   *   x  = 0.6102683302836215
   *   y1 = 0.7906090004346522
   *   y2 = y1 + 1e-16
   * the test
   *   hypot(x, y2) >= hypot(x, y1)
   * fails.  See also
   *   https://bugs.python.org/issue43088 */
  double r = sqrt(*sinx * *sinx + *cosx * *cosx);
#else
  double r = hypot(*sinx, *cosx);
#endif
  *sinx /= r;
  *cosx /= r;
}

static double AngNormalize(double x) {
  double y = remainder(x, (double)td);
  return fabs(y) == hd ? copysign((double)hd, x) : y;
}

static double LatFix(double x)
{ return fabs(x) > qd ? NaN : x; }

static double AngDiff(double x, double y, double* e) {
  /* Use remainder instead of AngNormalize, since we treat boundary cases
   * later taking account of the error */
  double t, d = sumx(remainder(-x, (double)td), remainder( y, (double)td), &t);
  /* This second sum can only change d if abs(d) < 128, so don't need to
   * apply remainder yet again. */
  d = sumx(remainder(d, (double)td), t, &t);
  /* Fix the sign if d = -180, 0, 180. */
  if (d == 0 || fabs(d) == hd)
    /* If t == 0, take sign from y - x
     * else (t != 0, implies d = +/-180), d and t must have opposite signs */
    d = copysign(d, t == 0 ? y - x : -t);
  if (e) *e = t;
  return d;
}

static double AngRound(double x) {
  /* False positive in cppcheck requires "1.0" instead of "1" */
  const double z = 1.0/16.0;
  volatile double y = fabs(x);
  volatile double w = z - y;
  /* The compiler mustn't "simplify" z - (z - y) to y */
  y = w > 0 ? z - w : y;
  return copysign(y, x);
}

static void sincosdx(double x, double* sinx, double* cosx) {
  /* In order to minimize round-off errors, this function exactly reduces
   * the argument to the range [-45, 45] before converting it to radians. */
  double r, s, c; int q = 0;
  r = remquo(x, (double)qd, &q);
  /* now abs(r) <= 45 */
  r *= degree;
  /* Possibly could call the gnu extension sincos */
  s = sin(r); c = cos(r);
  switch ((unsigned)q & 3U) {
  case 0U: *sinx =  s; *cosx =  c; break;
  case 1U: *sinx =  c; *cosx = -s; break;
  case 2U: *sinx = -s; *cosx = -c; break;
  default: *sinx = -c; *cosx =  s; break; /* case 3U */
  }
  /* http://www.open-std.org/jtc1/sc22/wg14/www/docs/n1950.pdf */
  *cosx += 0;                   /* special values from F.10.1.12 */
  /* special values from F.10.1.13 */
  if (*sinx == 0) *sinx = copysign(*sinx, x);
}

static void sincosde(double x, double t, double* sinx, double* cosx) {
  /* In order to minimize round-off errors, this function exactly reduces
   * the argument to the range [-45, 45] before converting it to radians. */
  double r, s, c; int q = 0;
  r = AngRound(remquo(x, (double)qd, &q) + t);
  /* now abs(r) <= 45 */
  r *= degree;
  /* Possibly could call the gnu extension sincos */
  s = sin(r); c = cos(r);
  switch ((unsigned)q & 3U) {
  case 0U: *sinx =  s; *cosx =  c; break;
  case 1U: *sinx =  c; *cosx = -s; break;
  case 2U: *sinx = -s; *cosx = -c; break;
  default: *sinx = -c; *cosx =  s; break; /* case 3U */
  }
  /* http://www.open-std.org/jtc1/sc22/wg14/www/docs/n1950.pdf */
  *cosx += 0;                   /* special values from F.10.1.12 */
  /* special values from F.10.1.13 */
  if (*sinx == 0) *sinx = copysign(*sinx, x);
}

static double atan2dx(double y, double x) {
  /* In order to minimize round-off errors, this function rearranges the
   * arguments so that result of atan2 is in the range [-pi/4, pi/4] before
   * converting it to degrees and mapping the result to the correct
   * quadrant. */
  int q = 0; double ang;
  if (fabs(y) > fabs(x)) { swapx(&x, &y); q = 2; }
  if (signbit(x)) { x = -x; ++q; }
  /* here x >= 0 and x >= abs(y), so angle is in [-pi/4, pi/4] */
  ang = atan2(y, x) / degree;
  switch (q) {
  case 1: ang = copysign((double)hd, y) - ang; break;
  case 2: ang =                  qd       - ang; break;
  case 3: ang =                 -qd       + ang; break;
  default: break;
  }
  return ang;
}

static void A3coeff(struct geod_geodesic* g);
static void C3coeff(struct geod_geodesic* g);
static void C4coeff(struct geod_geodesic* g);
static double SinCosSeries(boolx sinp,
                           double sinx, double cosx,
                           const double c[], int n);
static void Lengths(const struct geod_geodesic* g,
                    double eps, double sig12,
                    double ssig1, double csig1, double dn1,
                    double ssig2, double csig2, double dn2,
                    double cbet1, double cbet2,
                    double* ps12b, double* pm12b, double* pm0,
                    double* pM12, double* pM21,
                    /* Scratch area of the right size */
                    double Ca[]);
static double Astroid(double x, double y);
static double InverseStart(const struct geod_geodesic* g,
                           double sbet1, double cbet1, double dn1,
                           double sbet2, double cbet2, double dn2,
                           double lam12, double slam12, double clam12,
                           double* psalp1, double* pcalp1,
                           /* Only updated if return val >= 0 */
                           double* psalp2, double* pcalp2,
                           /* Only updated for short lines */
                           double* pdnm,
                           /* Scratch area of the right size */
                           double Ca[]);
static double Lambda12(const struct geod_geodesic* g,
                       double sbet1, double cbet1, double dn1,
                       double sbet2, double cbet2, double dn2,
                       double salp1, double calp1,
                       double slam120, double clam120,
                       double* psalp2, double* pcalp2,
                       double* psig12,
                       double* pssig1, double* pcsig1,
                       double* pssig2, double* pcsig2,
                       double* peps,
                       double* pdomg12,
                       boolx diffp, double* pdlam12,
                       /* Scratch area of the right size */
                       double Ca[]);
static double A3f(const struct geod_geodesic* g, double eps);
static void C3f(const struct geod_geodesic* g, double eps, double c[]);
static void C4f(const struct geod_geodesic* g, double eps, double c[]);
static double A1m1f(double eps);
static void C1f(double eps, double c[]);
static void C1pf(double eps, double c[]);
static double A2m1f(double eps);
static void C2f(double eps, double c[]);
static int transit(double lon1, double lon2);
static int transitdirect(double lon1, double lon2);
static void accini(double s[]);
static void acccopy(const double s[], double t[]);
static void accadd(double s[], double y);
static double accsum(const double s[], double y);
static void accneg(double s[]);
static void accrem(double s[], double y);
static double areareduceA(double area[], double area0,
                          int crossings, boolx reverse, boolx sign);
static double areareduceB(double area, double area0,
                          int crossings, boolx reverse, boolx sign);

void geod_init(struct geod_geodesic* g, double a, double f) {
  if (!init) Init();
  g->a = a;
  g->f = f;
  g->f1 = 1 - g->f;
  g->e2 = g->f * (2 - g->f);
  g->ep2 = g->e2 / sq(g->f1);   /* e2 / (1 - e2) */
  g->n = g->f / ( 2 - g->f);
  g->b = g->a * g->f1;
  g->c2 = (sq(g->a) + sq(g->b) *
           (g->e2 == 0 ? 1 :
            (g->e2 > 0 ? atanh(sqrt(g->e2)) : atan(sqrt(-g->e2))) /
            sqrt(fabs(g->e2))))/2; /* authalic radius squared */
  /* The sig12 threshold for "really short".  Using the auxiliary sphere
   * solution with dnm computed at (bet1 + bet2) / 2, the relative error in the
   * azimuth consistency check is sig12^2 * abs(f) * min(1, 1-f/2) / 2.  (Error
   * measured for 1/100 < b/a < 100 and abs(f) >= 1/1000.  For a given f and
   * sig12, the max error occurs for lines near the pole.  If the old rule for
   * computing dnm = (dn1 + dn2)/2 is used, then the error increases by a
   * factor of 2.)  Setting this equal to epsilon gives sig12 = etol2.  Here
   * 0.1 is a safety factor (error decreased by 100) and max(0.001, abs(f))
   * stops etol2 getting too large in the nearly spherical case. */
  g->etol2 = 0.1 * tol2 /
    sqrt( fmax(0.001, fabs(g->f)) * fmin(1.0, 1 - g->f/2) / 2 );

  A3coeff(g);
  C3coeff(g);
  C4coeff(g);
}

static void geod_lineinit_int(struct geod_geodesicline* l,
                              const struct geod_geodesic* g,
                              double lat1, double lon1,
                              double azi1, double salp1, double calp1,
                              unsigned caps) {
  double cbet1, sbet1, eps;
  l->a = g->a;
  l->f = g->f;
  l->b = g->b;
  l->c2 = g->c2;
  l->f1 = g->f1;
  /* If caps is 0 assume the standard direct calculation */
  l->caps = (caps ? caps : GEOD_DISTANCE_IN | GEOD_LONGITUDE) |
    /* always allow latitude and azimuth and unrolling of longitude */
    GEOD_LATITUDE | GEOD_AZIMUTH | GEOD_LONG_UNROLL;

  l->lat1 = LatFix(lat1);
  l->lon1 = lon1;
  l->azi1 = azi1;
  l->salp1 = salp1;
  l->calp1 = calp1;

  sincosdx(AngRound(l->lat1), &sbet1, &cbet1); sbet1 *= l->f1;
  /* Ensure cbet1 = +epsilon at poles */
  norm2(&sbet1, &cbet1); cbet1 = fmax(tiny, cbet1);
  l->dn1 = sqrt(1 + g->ep2 * sq(sbet1));

  /* Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0), */
  l->salp0 = l->salp1 * cbet1; /* alp0 in [0, pi/2 - |bet1|] */
  /* Alt: calp0 = hypot(sbet1, calp1 * cbet1).  The following
   * is slightly better (consider the case salp1 = 0). */
  l->calp0 = hypot(l->calp1, l->salp1 * sbet1);
  /* Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
   * sig = 0 is nearest northward crossing of equator.
   * With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
   * With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
   * With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
   * Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
   * With alp0 in (0, pi/2], quadrants for sig and omg coincide.
   * No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
   * With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi. */
  l->ssig1 = sbet1; l->somg1 = l->salp0 * sbet1;
  l->csig1 = l->comg1 = sbet1 != 0 || l->calp1 != 0 ? cbet1 * l->calp1 : 1;
  norm2(&l->ssig1, &l->csig1); /* sig1 in (-pi, pi] */
  /* norm2(somg1, comg1); -- don't need to normalize! */

  l->k2 = sq(l->calp0) * g->ep2;
  eps = l->k2 / (2 * (1 + sqrt(1 + l->k2)) + l->k2);

  if (l->caps & CAP_C1) {
    double s, c;
    l->A1m1 = A1m1f(eps);
    C1f(eps, l->C1a);
    l->B11 = SinCosSeries(TRUE, l->ssig1, l->csig1, l->C1a, nC1);
    s = sin(l->B11); c = cos(l->B11);
    /* tau1 = sig1 + B11 */
    l->stau1 = l->ssig1 * c + l->csig1 * s;
    l->ctau1 = l->csig1 * c - l->ssig1 * s;
    /* Not necessary because C1pa reverts C1a
     *    B11 = -SinCosSeries(TRUE, stau1, ctau1, C1pa, nC1p); */
  }

  if (l->caps & CAP_C1p)
    C1pf(eps, l->C1pa);

  if (l->caps & CAP_C2) {
    l->A2m1 = A2m1f(eps);
    C2f(eps, l->C2a);
    l->B21 = SinCosSeries(TRUE, l->ssig1, l->csig1, l->C2a, nC2);
  }

  if (l->caps & CAP_C3) {
    C3f(g, eps, l->C3a);
    l->A3c = -l->f * l->salp0 * A3f(g, eps);
    l->B31 = SinCosSeries(TRUE, l->ssig1, l->csig1, l->C3a, nC3-1);
  }

  if (l->caps & CAP_C4) {
    C4f(g, eps, l->C4a);
    /* Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0) */
    l->A4 = sq(l->a) * l->calp0 * l->salp0 * g->e2;
    l->B41 = SinCosSeries(FALSE, l->ssig1, l->csig1, l->C4a, nC4);
  }

  l->a13 = l->s13 = NaN;
}

void geod_lineinit(struct geod_geodesicline* l,
                   const struct geod_geodesic* g,
                   double lat1, double lon1, double azi1, unsigned caps) {
  double salp1, calp1;
  azi1 = AngNormalize(azi1);
  /* Guard against underflow in salp0 */
  sincosdx(AngRound(azi1), &salp1, &calp1);
  geod_lineinit_int(l, g, lat1, lon1, azi1, salp1, calp1, caps);
}

void geod_gendirectline(struct geod_geodesicline* l,
                        const struct geod_geodesic* g,
                        double lat1, double lon1, double azi1,
                        unsigned flags, double s12_a12,
                        unsigned caps) {
  geod_lineinit(l, g, lat1, lon1, azi1, caps);
  geod_gensetdistance(l, flags, s12_a12);
}

void geod_directline(struct geod_geodesicline* l,
                        const struct geod_geodesic* g,
                        double lat1, double lon1, double azi1,
                        double s12, unsigned caps) {
  geod_gendirectline(l, g, lat1, lon1, azi1, GEOD_NOFLAGS, s12, caps);
}

double geod_genposition(const struct geod_geodesicline* l,
                        unsigned flags, double s12_a12,
                        double* plat2, double* plon2, double* pazi2,
                        double* ps12, double* pm12,
                        double* pM12, double* pM21,
                        double* pS12) {
  double lat2 = 0, lon2 = 0, azi2 = 0, s12 = 0,
    m12 = 0, M12 = 0, M21 = 0, S12 = 0;
  /* Avoid warning about uninitialized B12. */
  double sig12, ssig12, csig12, B12 = 0, AB1 = 0;
  double omg12, lam12, lon12;
  double ssig2, csig2, sbet2, cbet2, somg2, comg2, salp2, calp2, dn2;
  unsigned outmask =
    (plat2 ? GEOD_LATITUDE : GEOD_NONE) |
    (plon2 ? GEOD_LONGITUDE : GEOD_NONE) |
    (pazi2 ? GEOD_AZIMUTH : GEOD_NONE) |
    (ps12 ? GEOD_DISTANCE : GEOD_NONE) |
    (pm12 ? GEOD_REDUCEDLENGTH : GEOD_NONE) |
    (pM12 || pM21 ? GEOD_GEODESICSCALE : GEOD_NONE) |
    (pS12 ? GEOD_AREA : GEOD_NONE);

  outmask &= l->caps & OUT_ALL;
  if (!( (flags & GEOD_ARCMODE || (l->caps & (GEOD_DISTANCE_IN & OUT_ALL))) ))
    /* Impossible distance calculation requested */
    return NaN;

  if (flags & GEOD_ARCMODE) {
    /* Interpret s12_a12 as spherical arc length */
    sig12 = s12_a12 * degree;
    sincosdx(s12_a12, &ssig12, &csig12);
  } else {
    /* Interpret s12_a12 as distance */
    double
      tau12 = s12_a12 / (l->b * (1 + l->A1m1)),
      s = sin(tau12),
      c = cos(tau12);
    /* tau2 = tau1 + tau12 */
    B12 = - SinCosSeries(TRUE,
                         l->stau1 * c + l->ctau1 * s,
                         l->ctau1 * c - l->stau1 * s,
                         l->C1pa, nC1p);
    sig12 = tau12 - (B12 - l->B11);
    ssig12 = sin(sig12); csig12 = cos(sig12);
    if (fabs(l->f) > 0.01) {
      /* Reverted distance series is inaccurate for |f| > 1/100, so correct
       * sig12 with 1 Newton iteration.  The following table shows the
       * approximate maximum error for a = WGS_a() and various f relative to
       * GeodesicExact.
       *     erri = the error in the inverse solution (nm)
       *     errd = the error in the direct solution (series only) (nm)
       *     errda = the error in the direct solution (series + 1 Newton) (nm)
       *
       *       f     erri  errd errda
       *     -1/5    12e6 1.2e9  69e6
       *     -1/10  123e3  12e6 765e3
       *     -1/20   1110 108e3  7155
       *     -1/50  18.63 200.9 27.12
       *     -1/100 18.63 23.78 23.37
       *     -1/150 18.63 21.05 20.26
       *      1/150 22.35 24.73 25.83
       *      1/100 22.35 25.03 25.31
       *      1/50  29.80 231.9 30.44
       *      1/20   5376 146e3  10e3
       *      1/10  829e3  22e6 1.5e6
       *      1/5   157e6 3.8e9 280e6 */
      double serr;
      ssig2 = l->ssig1 * csig12 + l->csig1 * ssig12;
      csig2 = l->csig1 * csig12 - l->ssig1 * ssig12;
      B12 = SinCosSeries(TRUE, ssig2, csig2, l->C1a, nC1);
      serr = (1 + l->A1m1) * (sig12 + (B12 - l->B11)) - s12_a12 / l->b;
      sig12 = sig12 - serr / sqrt(1 + l->k2 * sq(ssig2));
      ssig12 = sin(sig12); csig12 = cos(sig12);
      /* Update B12 below */
    }
  }

  /* sig2 = sig1 + sig12 */
  ssig2 = l->ssig1 * csig12 + l->csig1 * ssig12;
  csig2 = l->csig1 * csig12 - l->ssig1 * ssig12;
  dn2 = sqrt(1 + l->k2 * sq(ssig2));
  if (outmask & (GEOD_DISTANCE | GEOD_REDUCEDLENGTH | GEOD_GEODESICSCALE)) {
    if (flags & GEOD_ARCMODE || fabs(l->f) > 0.01)
      B12 = SinCosSeries(TRUE, ssig2, csig2, l->C1a, nC1);
    AB1 = (1 + l->A1m1) * (B12 - l->B11);
  }
  /* sin(bet2) = cos(alp0) * sin(sig2) */
  sbet2 = l->calp0 * ssig2;
  /* Alt: cbet2 = hypot(csig2, salp0 * ssig2); */
  cbet2 = hypot(l->salp0, l->calp0 * csig2);
  if (cbet2 == 0)
    /* I.e., salp0 = 0, csig2 = 0.  Break the degeneracy in this case */
    cbet2 = csig2 = tiny;
  /* tan(alp0) = cos(sig2)*tan(alp2) */
  salp2 = l->salp0; calp2 = l->calp0 * csig2; /* No need to normalize */

  if (outmask & GEOD_DISTANCE)
    s12 = (flags & GEOD_ARCMODE) ?
      l->b * ((1 + l->A1m1) * sig12 + AB1) :
      s12_a12;

  if (outmask & GEOD_LONGITUDE) {
    double E = copysign(1, l->salp0); /* east or west going? */
    /* tan(omg2) = sin(alp0) * tan(sig2) */
    somg2 = l->salp0 * ssig2; comg2 = csig2;  /* No need to normalize */
    /* omg12 = omg2 - omg1 */
    omg12 = (flags & GEOD_LONG_UNROLL)
      ? E * (sig12
             - (atan2(    ssig2, csig2) - atan2(    l->ssig1, l->csig1))
             + (atan2(E * somg2, comg2) - atan2(E * l->somg1, l->comg1)))
      : atan2(somg2 * l->comg1 - comg2 * l->somg1,
              comg2 * l->comg1 + somg2 * l->somg1);
    lam12 = omg12 + l->A3c *
      ( sig12 + (SinCosSeries(TRUE, ssig2, csig2, l->C3a, nC3-1)
                 - l->B31));
    lon12 = lam12 / degree;
    lon2 = (flags & GEOD_LONG_UNROLL) ? l->lon1 + lon12 :
      AngNormalize(AngNormalize(l->lon1) + AngNormalize(lon12));
  }

  if (outmask & GEOD_LATITUDE)
    lat2 = atan2dx(sbet2, l->f1 * cbet2);

  if (outmask & GEOD_AZIMUTH)
    azi2 = atan2dx(salp2, calp2);

  if (outmask & (GEOD_REDUCEDLENGTH | GEOD_GEODESICSCALE)) {
    double
      B22 = SinCosSeries(TRUE, ssig2, csig2, l->C2a, nC2),
      AB2 = (1 + l->A2m1) * (B22 - l->B21),
      J12 = (l->A1m1 - l->A2m1) * sig12 + (AB1 - AB2);
    if (outmask & GEOD_REDUCEDLENGTH)
      /* Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
       * accurate cancellation in the case of coincident points. */
      m12 = l->b * ((dn2 * (l->csig1 * ssig2) - l->dn1 * (l->ssig1 * csig2))
                    - l->csig1 * csig2 * J12);
    if (outmask & GEOD_GEODESICSCALE) {
      double t = l->k2 * (ssig2 - l->ssig1) * (ssig2 + l->ssig1) /
        (l->dn1 + dn2);
      M12 = csig12 + (t *  ssig2 -  csig2 * J12) * l->ssig1 / l->dn1;
      M21 = csig12 - (t * l->ssig1 - l->csig1 * J12) *  ssig2 /  dn2;
    }
  }

  if (outmask & GEOD_AREA) {
    double
      B42 = SinCosSeries(FALSE, ssig2, csig2, l->C4a, nC4);
    double salp12, calp12;
    if (l->calp0 == 0 || l->salp0 == 0) {
      /* alp12 = alp2 - alp1, used in atan2 so no need to normalize */
      salp12 = salp2 * l->calp1 - calp2 * l->salp1;
      calp12 = calp2 * l->calp1 + salp2 * l->salp1;
    } else {
      /* tan(alp) = tan(alp0) * sec(sig)
       * tan(alp2-alp1) = (tan(alp2) -tan(alp1)) / (tan(alp2)*tan(alp1)+1)
       * = calp0 * salp0 * (csig1-csig2) / (salp0^2 + calp0^2 * csig1*csig2)
       * If csig12 > 0, write
       *   csig1 - csig2 = ssig12 * (csig1 * ssig12 / (1 + csig12) + ssig1)
       * else
       *   csig1 - csig2 = csig1 * (1 - csig12) + ssig12 * ssig1
       * No need to normalize */
      salp12 = l->calp0 * l->salp0 *
        (csig12 <= 0 ? l->csig1 * (1 - csig12) + ssig12 * l->ssig1 :
         ssig12 * (l->csig1 * ssig12 / (1 + csig12) + l->ssig1));
      calp12 = sq(l->salp0) + sq(l->calp0) * l->csig1 * csig2;
    }
    S12 = l->c2 * atan2(salp12, calp12) + l->A4 * (B42 - l->B41);
  }

  /* In the pattern
   *
   *   if ((outmask & GEOD_XX) && pYY)
   *     *pYY = YY;
   *
   * the second check "&& pYY" is redundant.  It's there to make the CLang
   * static analyzer happy.
   */
  if ((outmask & GEOD_LATITUDE) && plat2)
    *plat2 = lat2;
  if ((outmask & GEOD_LONGITUDE) && plon2)
    *plon2 = lon2;
  if ((outmask & GEOD_AZIMUTH) && pazi2)
    *pazi2 = azi2;
  if ((outmask & GEOD_DISTANCE) && ps12)
    *ps12 = s12;
  if ((outmask & GEOD_REDUCEDLENGTH) && pm12)
    *pm12 = m12;
  if (outmask & GEOD_GEODESICSCALE) {
    if (pM12) *pM12 = M12;
    if (pM21) *pM21 = M21;
  }
  if ((outmask & GEOD_AREA) && pS12)
    *pS12 = S12;

  return (flags & GEOD_ARCMODE) ? s12_a12 : sig12 / degree;
}

void geod_setdistance(struct geod_geodesicline* l, double s13) {
  l->s13 = s13;
  l->a13 = geod_genposition(l, GEOD_NOFLAGS, l->s13, nullptr, nullptr, nullptr,
                            nullptr, nullptr, nullptr, nullptr, nullptr);
}

static void geod_setarc(struct geod_geodesicline* l, double a13) {
  l->a13 = a13; l->s13 = NaN;
  geod_genposition(l, GEOD_ARCMODE, l->a13, nullptr, nullptr, nullptr, &l->s13,
                   nullptr, nullptr, nullptr, nullptr);
}

void geod_gensetdistance(struct geod_geodesicline* l,
 unsigned flags, double s13_a13) {
  (flags & GEOD_ARCMODE) ?
    geod_setarc(l, s13_a13) :
    geod_setdistance(l, s13_a13);
}

void geod_position(const struct geod_geodesicline* l, double s12,
                   double* plat2, double* plon2, double* pazi2) {
  geod_genposition(l, FALSE, s12, plat2, plon2, pazi2,
                   nullptr, nullptr, nullptr, nullptr, nullptr);
}

double geod_gendirect(const struct geod_geodesic* g,
                      double lat1, double lon1, double azi1,
                      unsigned flags, double s12_a12,
                      double* plat2, double* plon2, double* pazi2,
                      double* ps12, double* pm12, double* pM12, double* pM21,
                      double* pS12) {
  struct geod_geodesicline l;
  unsigned outmask =
    (plat2 ? GEOD_LATITUDE : GEOD_NONE) |
    (plon2 ? GEOD_LONGITUDE : GEOD_NONE) |
    (pazi2 ? GEOD_AZIMUTH : GEOD_NONE) |
    (ps12 ? GEOD_DISTANCE : GEOD_NONE) |
    (pm12 ? GEOD_REDUCEDLENGTH : GEOD_NONE) |
    (pM12 || pM21 ? GEOD_GEODESICSCALE : GEOD_NONE) |
    (pS12 ? GEOD_AREA : GEOD_NONE);

  geod_lineinit(&l, g, lat1, lon1, azi1,
                /* Automatically supply GEOD_DISTANCE_IN if necessary */
                outmask |
                ((flags & GEOD_ARCMODE) ? GEOD_NONE : GEOD_DISTANCE_IN));
  return geod_genposition(&l, flags, s12_a12,
                          plat2, plon2, pazi2, ps12, pm12, pM12, pM21, pS12);
}

void geod_direct(const struct geod_geodesic* g,
                 double lat1, double lon1, double azi1,
                 double s12,
                 double* plat2, double* plon2, double* pazi2) {
  geod_gendirect(g, lat1, lon1, azi1, GEOD_NOFLAGS, s12, plat2, plon2, pazi2,
                 nullptr, nullptr, nullptr, nullptr, nullptr);
}

static double geod_geninverse_int(const struct geod_geodesic* g,
                                  double lat1, double lon1,
                                  double lat2, double lon2,
                                  double* ps12,
                                  double* psalp1, double* pcalp1,
                                  double* psalp2, double* pcalp2,
                                  double* pm12, double* pM12, double* pM21,
                                  double* pS12) {
  double s12 = 0, m12 = 0, M12 = 0, M21 = 0, S12 = 0;
  double lon12, lon12s;
  int latsign, lonsign, swapp;
  double sbet1, cbet1, sbet2, cbet2, s12x = 0, m12x = 0;
  double dn1, dn2, lam12, slam12, clam12;
  double a12 = 0, sig12, calp1 = 0, salp1 = 0, calp2 = 0, salp2 = 0;
  double Ca[nC];
  boolx meridian;
  /* somg12 == 2 marks that it needs to be calculated */
  double omg12 = 0, somg12 = 2, comg12 = 0;

  unsigned outmask =
    (ps12 ? GEOD_DISTANCE : GEOD_NONE) |
    (pm12 ? GEOD_REDUCEDLENGTH : GEOD_NONE) |
    (pM12 || pM21 ? GEOD_GEODESICSCALE : GEOD_NONE) |
    (pS12 ? GEOD_AREA : GEOD_NONE);

  outmask &= OUT_ALL;
  /* Compute longitude difference (AngDiff does this carefully).  Result is
   * in [-180, 180] but -180 is only for west-going geodesics.  180 is for
   * east-going and meridional geodesics. */
  lon12 = AngDiff(lon1, lon2, &lon12s);
  /* Make longitude difference positive. */
  lonsign = signbit(lon12) ? -1 : 1;
  lon12 *= lonsign; lon12s *= lonsign;
  lam12 = lon12 * degree;
  /* Calculate sincos of lon12 + error (this applies AngRound internally). */
  sincosde(lon12, lon12s, &slam12, &clam12);
  lon12s = (hd - lon12) - lon12s; /* the supplementary longitude difference */

  /* If really close to the equator, treat as on equator. */
  lat1 = AngRound(LatFix(lat1));
  lat2 = AngRound(LatFix(lat2));
  /* Swap points so that point with higher (abs) latitude is point 1
   * If one latitude is a nan, then it becomes lat1. */
  swapp = fabs(lat1) < fabs(lat2) || lat2 != lat2 ? -1 : 1;
  if (swapp < 0) {
    lonsign *= -1;
    swapx(&lat1, &lat2);
  }
  /* Make lat1 <= -0 */
  latsign = signbit(lat1) ? 1 : -1;
  lat1 *= latsign;
  lat2 *= latsign;
  /* Now we have
   *
   *     0 <= lon12 <= 180
   *     -90 <= lat1 <= -0
   *     lat1 <= lat2 <= -lat1
   *
   * longsign, swapp, latsign register the transformation to bring the
   * coordinates to this canonical form.  In all cases, 1 means no change was
   * made.  We make these transformations so that there are few cases to
   * check, e.g., on verifying quadrants in atan2.  In addition, this
   * enforces some symmetries in the results returned. */

  sincosdx(lat1, &sbet1, &cbet1); sbet1 *= g->f1;
  /* Ensure cbet1 = +epsilon at poles */
  norm2(&sbet1, &cbet1); cbet1 = fmax(tiny, cbet1);

  sincosdx(lat2, &sbet2, &cbet2); sbet2 *= g->f1;
  /* Ensure cbet2 = +epsilon at poles */
  norm2(&sbet2, &cbet2); cbet2 = fmax(tiny, cbet2);

  /* If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the
   * |bet1| - |bet2|.  Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1 is
   * a better measure.  This logic is used in assigning calp2 in Lambda12.
   * Sometimes these quantities vanish and in that case we force bet2 = +/-
   * bet1 exactly.  An example where is is necessary is the inverse problem
   * 48.522876735459 0 -48.52287673545898293 179.599720456223079643
   * which failed with Visual Studio 10 (Release and Debug) */

  if (cbet1 < -sbet1) {
    if (cbet2 == cbet1)
      sbet2 = copysign(sbet1, sbet2);
  } else {
    if (fabs(sbet2) == -sbet1)
      cbet2 = cbet1;
  }

  dn1 = sqrt(1 + g->ep2 * sq(sbet1));
  dn2 = sqrt(1 + g->ep2 * sq(sbet2));

  meridian = lat1 == -qd || slam12 == 0;

  if (meridian) {

    /* Endpoints are on a single full meridian, so the geodesic might lie on
     * a meridian. */

    double ssig1, csig1, ssig2, csig2;
    calp1 = clam12; salp1 = slam12; /* Head to the target longitude */
    calp2 = 1; salp2 = 0;           /* At the target we're heading north */

    /* tan(bet) = tan(sig) * cos(alp) */
    ssig1 = sbet1; csig1 = calp1 * cbet1;
    ssig2 = sbet2; csig2 = calp2 * cbet2;

    /* sig12 = sig2 - sig1 */
    sig12 = atan2(fmax(0.0, csig1 * ssig2 - ssig1 * csig2) + 0,
                            csig1 * csig2 + ssig1 * ssig2);
    Lengths(g, g->n, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
            cbet1, cbet2, &s12x, &m12x, nullptr,
            (outmask & GEOD_GEODESICSCALE) ? &M12 : nullptr,
            (outmask & GEOD_GEODESICSCALE) ? &M21 : nullptr,
            Ca);
    /* Add the check for sig12 since zero length geodesics might yield m12 <
     * 0.  Test case was
     *
     *    echo 20.001 0 20.001 0 | GeodSolve -i
     *
     * In fact, we will have sig12 > pi/2 for meridional geodesic which is
     * not a shortest path. */
    if (sig12 < 1 || m12x >= 0) {
      /* Need at least 2, to handle 90 0 90 180 */
      if (sig12 < 3 * tiny ||
          /* Prevent negative s12 or m12 for short lines */
          (sig12 < tol0 && (s12x < 0 || m12x < 0)))
        sig12 = m12x = s12x = 0;
      m12x *= g->b;
      s12x *= g->b;
      a12 = sig12 / degree;
    } else
      /* m12 < 0, i.e., prolate and too close to anti-podal */
      meridian = FALSE;
  }

  if (!meridian &&
      sbet1 == 0 &&           /* and sbet2 == 0 */
      /* Mimic the way Lambda12 works with calp1 = 0 */
      (g->f <= 0 || lon12s >= g->f * hd)) {

    /* Geodesic runs along equator */
    calp1 = calp2 = 0; salp1 = salp2 = 1;
    s12x = g->a * lam12;
    sig12 = omg12 = lam12 / g->f1;
    m12x = g->b * sin(sig12);
    if (outmask & GEOD_GEODESICSCALE)
      M12 = M21 = cos(sig12);
    a12 = lon12 / g->f1;

  } else if (!meridian) {

    /* Now point1 and point2 belong within a hemisphere bounded by a
     * meridian and geodesic is neither meridional or equatorial. */

    /* Figure a starting point for Newton's method */
    double dnm = 0;
    sig12 = InverseStart(g, sbet1, cbet1, dn1, sbet2, cbet2, dn2,
                         lam12, slam12, clam12,
                         &salp1, &calp1, &salp2, &calp2, &dnm,
                         Ca);

    if (sig12 >= 0) {
      /* Short lines (InverseStart sets salp2, calp2, dnm) */
      s12x = sig12 * g->b * dnm;
      m12x = sq(dnm) * g->b * sin(sig12 / dnm);
      if (outmask & GEOD_GEODESICSCALE)
        M12 = M21 = cos(sig12 / dnm);
      a12 = sig12 / degree;
      omg12 = lam12 / (g->f1 * dnm);
    } else {

      /* Newton's method.  This is a straightforward solution of f(alp1) =
       * lambda12(alp1) - lam12 = 0 with one wrinkle.  f(alp) has exactly one
       * root in the interval (0, pi) and its derivative is positive at the
       * root.  Thus f(alp) is positive for alp > alp1 and negative for alp <
       * alp1.  During the course of the iteration, a range (alp1a, alp1b) is
       * maintained which brackets the root and with each evaluation of
       * f(alp) the range is shrunk, if possible.  Newton's method is
       * restarted whenever the derivative of f is negative (because the new
       * value of alp1 is then further from the solution) or if the new
       * estimate of alp1 lies outside (0,pi); in this case, the new starting
       * guess is taken to be (alp1a + alp1b) / 2. */
      double ssig1 = 0, csig1 = 0, ssig2 = 0, csig2 = 0, eps = 0, domg12 = 0;
      unsigned numit = 0;
      /* Bracketing range */
      double salp1a = tiny, calp1a = 1, salp1b = tiny, calp1b = -1;
      boolx tripn = FALSE;
      boolx tripb = FALSE;
      for (;; ++numit) {
        /* the WGS84 test set: mean = 1.47, sd = 1.25, max = 16
         * WGS84 and random input: mean = 2.85, sd = 0.60 */
        double dv = 0,
          v = Lambda12(g, sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1,
                        slam12, clam12,
                        &salp2, &calp2, &sig12, &ssig1, &csig1, &ssig2, &csig2,
                        &eps, &domg12, numit < maxit1, &dv, Ca);
        if (tripb ||
            /* Reversed test to allow escape with NaNs */
            !(fabs(v) >= (tripn ? 8 : 1) * tol0) ||
            /* Enough bisections to get accurate result */
            numit == maxit2)
          break;
        /* Update bracketing values */
        if (v > 0 && (numit > maxit1 || calp1/salp1 > calp1b/salp1b))
          { salp1b = salp1; calp1b = calp1; }
        else if (v < 0 && (numit > maxit1 || calp1/salp1 < calp1a/salp1a))
          { salp1a = salp1; calp1a = calp1; }
        if (numit < maxit1 && dv > 0) {
          double
            dalp1 = -v/dv;
          if (fabs(dalp1) < pi) {
            double
              sdalp1 = sin(dalp1), cdalp1 = cos(dalp1),
              nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
            if (nsalp1 > 0) {
              calp1 = calp1 * cdalp1 - salp1 * sdalp1;
              salp1 = nsalp1;
              norm2(&salp1, &calp1);
              /* In some regimes we don't get quadratic convergence because
               * slope -> 0.  So use convergence conditions based on epsilon
               * instead of sqrt(epsilon). */
              tripn = fabs(v) <= 16 * tol0;
              continue;
            }
          }
        }
        /* Either dv was not positive or updated value was outside legal
         * range.  Use the midpoint of the bracket as the next estimate.
         * This mechanism is not needed for the WGS84 ellipsoid, but it does
         * catch problems with more eccentric ellipsoids.  Its efficacy is
         * such for the WGS84 test set with the starting guess set to alp1 =
         * 90deg:
         * the WGS84 test set: mean = 5.21, sd = 3.93, max = 24
         * WGS84 and random input: mean = 4.74, sd = 0.99 */
        salp1 = (salp1a + salp1b)/2;
        calp1 = (calp1a + calp1b)/2;
        norm2(&salp1, &calp1);
        tripn = FALSE;
        tripb = (fabs(salp1a - salp1) + (calp1a - calp1) < tolb ||
                 fabs(salp1 - salp1b) + (calp1 - calp1b) < tolb);
      }
      Lengths(g, eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
              cbet1, cbet2, &s12x, &m12x, nullptr,
              (outmask & GEOD_GEODESICSCALE) ? &M12 : nullptr,
              (outmask & GEOD_GEODESICSCALE) ? &M21 : nullptr, Ca);
      m12x *= g->b;
      s12x *= g->b;
      a12 = sig12 / degree;
      if (outmask & GEOD_AREA) {
        /* omg12 = lam12 - domg12 */
        double sdomg12 = sin(domg12), cdomg12 = cos(domg12);
        somg12 = slam12 * cdomg12 - clam12 * sdomg12;
        comg12 = clam12 * cdomg12 + slam12 * sdomg12;
      }
    }
  }

  if (outmask & GEOD_DISTANCE)
    s12 = 0 + s12x;             /* Convert -0 to 0 */

  if (outmask & GEOD_REDUCEDLENGTH)
    m12 = 0 + m12x;             /* Convert -0 to 0 */

  if (outmask & GEOD_AREA) {
    double
      /* From Lambda12: sin(alp1) * cos(bet1) = sin(alp0) */
      salp0 = salp1 * cbet1,
      calp0 = hypot(calp1, salp1 * sbet1); /* calp0 > 0 */
    double alp12;
    if (calp0 != 0 && salp0 != 0) {
      double
        /* From Lambda12: tan(bet) = tan(sig) * cos(alp) */
        ssig1 = sbet1, csig1 = calp1 * cbet1,
        ssig2 = sbet2, csig2 = calp2 * cbet2,
        k2 = sq(calp0) * g->ep2,
        eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2),
        /* Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0). */
        A4 = sq(g->a) * calp0 * salp0 * g->e2;
      double B41, B42;
      norm2(&ssig1, &csig1);
      norm2(&ssig2, &csig2);
      C4f(g, eps, Ca);
      B41 = SinCosSeries(FALSE, ssig1, csig1, Ca, nC4);
      B42 = SinCosSeries(FALSE, ssig2, csig2, Ca, nC4);
      S12 = A4 * (B42 - B41);
    } else
      /* Avoid problems with indeterminate sig1, sig2 on equator */
      S12 = 0;

    if (!meridian && somg12 == 2) {
      somg12 = sin(omg12); comg12 = cos(omg12);
    }

    if (!meridian &&
        /* omg12 < 3/4 * pi */
        comg12 > -0.7071 &&     /* Long difference not too big */
        sbet2 - sbet1 < 1.75) { /* Lat difference not too big */
      /* Use tan(Gamma/2) = tan(omg12/2)
       * * (tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
       * with tan(x/2) = sin(x)/(1+cos(x)) */
      double
        domg12 = 1 + comg12, dbet1 = 1 + cbet1, dbet2 = 1 + cbet2;
      alp12 = 2 * atan2( somg12 * ( sbet1 * dbet2 + sbet2 * dbet1 ),
                         domg12 * ( sbet1 * sbet2 + dbet1 * dbet2 ) );
    } else {
      /* alp12 = alp2 - alp1, used in atan2 so no need to normalize */
      double
        salp12 = salp2 * calp1 - calp2 * salp1,
        calp12 = calp2 * calp1 + salp2 * salp1;
      /* The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
       * salp12 = -0 and alp12 = -180.  However this depends on the sign
       * being attached to 0 correctly.  The following ensures the correct
       * behavior. */
      if (salp12 == 0 && calp12 < 0) {
        salp12 = tiny * calp1;
        calp12 = -1;
      }
      alp12 = atan2(salp12, calp12);
    }
    S12 += g->c2 * alp12;
    S12 *= swapp * lonsign * latsign;
    /* Convert -0 to 0 */
    S12 += 0;
  }

  /* Convert calp, salp to azimuth accounting for lonsign, swapp, latsign. */
  if (swapp < 0) {
    swapx(&salp1, &salp2);
    swapx(&calp1, &calp2);
    if (outmask & GEOD_GEODESICSCALE)
      swapx(&M12, &M21);
  }

  salp1 *= swapp * lonsign; calp1 *= swapp * latsign;
  salp2 *= swapp * lonsign; calp2 *= swapp * latsign;

  if (psalp1) *psalp1 = salp1;
  if (pcalp1) *pcalp1 = calp1;
  if (psalp2) *psalp2 = salp2;
  if (pcalp2) *pcalp2 = calp2;

  if (outmask & GEOD_DISTANCE)
    *ps12 = s12;
  if (outmask & GEOD_REDUCEDLENGTH)
    *pm12 = m12;
  if (outmask & GEOD_GEODESICSCALE) {
    if (pM12) *pM12 = M12;
    if (pM21) *pM21 = M21;
  }
  if (outmask & GEOD_AREA)
    *pS12 = S12;

  /* Returned value in [0, 180] */
  return a12;
}

double geod_geninverse(const struct geod_geodesic* g,
                       double lat1, double lon1, double lat2, double lon2,
                       double* ps12, double* pazi1, double* pazi2,
                       double* pm12, double* pM12, double* pM21,
                       double* pS12) {
  double salp1, calp1, salp2, calp2,
    a12 = geod_geninverse_int(g, lat1, lon1, lat2, lon2, ps12,
                              &salp1, &calp1, &salp2, &calp2,
                              pm12, pM12, pM21, pS12);
  if (pazi1) *pazi1 = atan2dx(salp1, calp1);
  if (pazi2) *pazi2 = atan2dx(salp2, calp2);
  return a12;
}

void geod_inverseline(struct geod_geodesicline* l,
                      const struct geod_geodesic* g,
                      double lat1, double lon1, double lat2, double lon2,
                      unsigned caps) {
  double salp1, calp1,
    a12 = geod_geninverse_int(g, lat1, lon1, lat2, lon2, nullptr,
                              &salp1, &calp1, nullptr, nullptr,
                              nullptr, nullptr, nullptr, nullptr),
    azi1 = atan2dx(salp1, calp1);
  caps = caps ? caps : GEOD_DISTANCE_IN | GEOD_LONGITUDE;
  /* Ensure that a12 can be converted to a distance */
  if (caps & (OUT_ALL & GEOD_DISTANCE_IN)) caps |= GEOD_DISTANCE;
  geod_lineinit_int(l, g, lat1, lon1, azi1, salp1, calp1, caps);
  geod_setarc(l, a12);
}

void geod_inverse(const struct geod_geodesic* g,
                  double lat1, double lon1, double lat2, double lon2,
                  double* ps12, double* pazi1, double* pazi2) {
  geod_geninverse(g, lat1, lon1, lat2, lon2, ps12, pazi1, pazi2,
                  nullptr, nullptr, nullptr, nullptr);
}

double SinCosSeries(boolx sinp, double sinx, double cosx,
                    const double c[], int n) {
  /* Evaluate
   * y = sinp ? sum(c[i] * sin( 2*i    * x), i, 1, n) :
   *            sum(c[i] * cos((2*i+1) * x), i, 0, n-1)
   * using Clenshaw summation.  N.B. c[0] is unused for sin series
   * Approx operation count = (n + 5) mult and (2 * n + 2) add */
  double ar, y0, y1;
  c += (n + sinp);              /* Point to one beyond last element */
  ar = 2 * (cosx - sinx) * (cosx + sinx); /* 2 * cos(2 * x) */
  y0 = (n & 1) ? *--c : 0; y1 = 0;        /* accumulators for sum */
  /* Now n is even */
  n /= 2;
  while (n--) {
    /* Unroll loop x 2, so accumulators return to their original role */
    y1 = ar * y0 - y1 + *--c;
    y0 = ar * y1 - y0 + *--c;
  }
  return sinp
    ? 2 * sinx * cosx * y0      /* sin(2 * x) * y0 */
    : cosx * (y0 - y1);         /* cos(x) * (y0 - y1) */
}

void Lengths(const struct geod_geodesic* g,
             double eps, double sig12,
             double ssig1, double csig1, double dn1,
             double ssig2, double csig2, double dn2,
             double cbet1, double cbet2,
             double* ps12b, double* pm12b, double* pm0,
             double* pM12, double* pM21,
             /* Scratch area of the right size */
             double Ca[]) {
  double m0 = 0, J12 = 0, A1 = 0, A2 = 0;
  double Cb[nC];

  /* Return m12b = (reduced length)/b; also calculate s12b = distance/b,
   * and m0 = coefficient of secular term in expression for reduced length. */
  boolx redlp = pm12b || pm0 || pM12 || pM21;
  if (ps12b || redlp) {
    A1 = A1m1f(eps);
    C1f(eps, Ca);
    if (redlp) {
      A2 = A2m1f(eps);
      C2f(eps, Cb);
      m0 = A1 - A2;
      A2 = 1 + A2;
    }
    A1 = 1 + A1;
  }
  if (ps12b) {
    double B1 = SinCosSeries(TRUE, ssig2, csig2, Ca, nC1) -
      SinCosSeries(TRUE, ssig1, csig1, Ca, nC1);
    /* Missing a factor of b */
    *ps12b = A1 * (sig12 + B1);
    if (redlp) {
      double B2 = SinCosSeries(TRUE, ssig2, csig2, Cb, nC2) -
        SinCosSeries(TRUE, ssig1, csig1, Cb, nC2);
      J12 = m0 * sig12 + (A1 * B1 - A2 * B2);
    }
  } else if (redlp) {
    /* Assume here that nC1 >= nC2 */
    int l;
    for (l = 1; l <= nC2; ++l)
      Cb[l] = A1 * Ca[l] - A2 * Cb[l];
    J12 = m0 * sig12 + (SinCosSeries(TRUE, ssig2, csig2, Cb, nC2) -
                        SinCosSeries(TRUE, ssig1, csig1, Cb, nC2));
  }
  if (pm0) *pm0 = m0;
  if (pm12b)
    /* Missing a factor of b.
     * Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
     * accurate cancellation in the case of coincident points. */
    *pm12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) -
      csig1 * csig2 * J12;
  if (pM12 || pM21) {
    double csig12 = csig1 * csig2 + ssig1 * ssig2;
    double t = g->ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2);
    if (pM12)
      *pM12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1;
    if (pM21)
      *pM21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2;
  }
}

double Astroid(double x, double y) {
  /* Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive root k.
   * This solution is adapted from Geocentric::Reverse. */
  double k;
  double
    p = sq(x),
    q = sq(y),
    r = (p + q - 1) / 6;
  if ( !(q == 0 && r <= 0) ) {
    double
      /* Avoid possible division by zero when r = 0 by multiplying equations
       * for s and t by r^3 and r, resp. */
      S = p * q / 4,            /* S = r^3 * s */
      r2 = sq(r),
      r3 = r * r2,
      /* The discriminant of the quadratic equation for T3.  This is zero on
       * the evolute curve p^(1/3)+q^(1/3) = 1 */
      disc = S * (S + 2 * r3);
    double u = r;
    double v, uv, w;
    if (disc >= 0) {
      double T3 = S + r3, T;
      /* Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss
       * of precision due to cancellation.  The result is unchanged because
       * of the way the T is used in definition of u. */
      T3 += T3 < 0 ? -sqrt(disc) : sqrt(disc); /* T3 = (r * t)^3 */
      /* N.B. cbrt always returns the double root.  cbrt(-8) = -2. */
      T = cbrt(T3);            /* T = r * t */
      /* T can be zero; but then r2 / T -> 0. */
      u += T + (T != 0 ? r2 / T : 0);
    } else {
      /* T is complex, but the way u is defined the result is double. */
      double ang = atan2(sqrt(-disc), -(S + r3));
      /* There are three possible cube roots.  We choose the root which
       * avoids cancellation.  Note that disc < 0 implies that r < 0. */
      u += 2 * r * cos(ang / 3);
    }
    v = sqrt(sq(u) + q);              /* guaranteed positive */
    /* Avoid loss of accuracy when u < 0. */
    uv = u < 0 ? q / (v - u) : u + v; /* u+v, guaranteed positive */
    w = (uv - q) / (2 * v);           /* positive? */
    /* Rearrange expression for k to avoid loss of accuracy due to
     * subtraction.  Division by 0 not possible because uv > 0, w >= 0. */
    k = uv / (sqrt(uv + sq(w)) + w);   /* guaranteed positive */
  } else {               /* q == 0 && r <= 0 */
    /* y = 0 with |x| <= 1.  Handle this case directly.
     * for y small, positive root is k = abs(y)/sqrt(1-x^2) */
    k = 0;
  }
  return k;
}

double InverseStart(const struct geod_geodesic* g,
                    double sbet1, double cbet1, double dn1,
                    double sbet2, double cbet2, double dn2,
                    double lam12, double slam12, double clam12,
                    double* psalp1, double* pcalp1,
                    /* Only updated if return val >= 0 */
                    double* psalp2, double* pcalp2,
                    /* Only updated for short lines */
                    double* pdnm,
                    /* Scratch area of the right size */
                    double Ca[]) {
  double salp1 = 0, calp1 = 0, salp2 = 0, calp2 = 0, dnm = 0;

  /* Return a starting point for Newton's method in salp1 and calp1 (function
   * value is -1).  If Newton's method doesn't need to be used, return also
   * salp2 and calp2 and function value is sig12. */
  double
    sig12 = -1,               /* Return value */
    /* bet12 = bet2 - bet1 in [0, pi); bet12a = bet2 + bet1 in (-pi, 0] */
    sbet12 = sbet2 * cbet1 - cbet2 * sbet1,
    cbet12 = cbet2 * cbet1 + sbet2 * sbet1;
  double sbet12a;
  boolx shortline = cbet12 >= 0 && sbet12 < 0.5 && cbet2 * lam12 < 0.5;
  double somg12, comg12, ssig12, csig12;
  sbet12a = sbet2 * cbet1 + cbet2 * sbet1;
  if (shortline) {
    double sbetm2 = sq(sbet1 + sbet2), omg12;
    /* sin((bet1+bet2)/2)^2
     * =  (sbet1 + sbet2)^2 / ((sbet1 + sbet2)^2 + (cbet1 + cbet2)^2) */
    sbetm2 /= sbetm2 + sq(cbet1 + cbet2);
    dnm = sqrt(1 + g->ep2 * sbetm2);
    omg12 = lam12 / (g->f1 * dnm);
    somg12 = sin(omg12); comg12 = cos(omg12);
  } else {
    somg12 = slam12; comg12 = clam12;
  }

  salp1 = cbet2 * somg12;
  calp1 = comg12 >= 0 ?
    sbet12 + cbet2 * sbet1 * sq(somg12) / (1 + comg12) :
    sbet12a - cbet2 * sbet1 * sq(somg12) / (1 - comg12);

  ssig12 = hypot(salp1, calp1);
  csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;

  if (shortline && ssig12 < g->etol2) {
    /* really short lines */
    salp2 = cbet1 * somg12;
    calp2 = sbet12 - cbet1 * sbet2 *
      (comg12 >= 0 ? sq(somg12) / (1 + comg12) : 1 - comg12);
    norm2(&salp2, &calp2);
    /* Set return value */
    sig12 = atan2(ssig12, csig12);
  } else if (fabs(g->n) > 0.1 || /* No astroid calc if too eccentric */
             csig12 >= 0 ||
             ssig12 >= 6 * fabs(g->n) * pi * sq(cbet1)) {
    /* Nothing to do, zeroth order spherical approximation is OK */
  } else {
    /* Scale lam12 and bet2 to x, y coordinate system where antipodal point
     * is at origin and singular point is at y = 0, x = -1. */
    double x, y, lamscale, betscale;
    double lam12x = atan2(-slam12, -clam12); /* lam12 - pi */
    if (g->f >= 0) {            /* In fact f == 0 does not get here */
      /* x = dlong, y = dlat */
      {
        double
          k2 = sq(sbet1) * g->ep2,
          eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2);
        lamscale = g->f * cbet1 * A3f(g, eps) * pi;
      }
      betscale = lamscale * cbet1;

      x = lam12x / lamscale;
      y = sbet12a / betscale;
    } else {                    /* f < 0 */
      /* x = dlat, y = dlong */
      double
        cbet12a = cbet2 * cbet1 - sbet2 * sbet1,
        bet12a = atan2(sbet12a, cbet12a);
      double m12b, m0;
      /* In the case of lon12 = 180, this repeats a calculation made in
       * Inverse. */
      Lengths(g, g->n, pi + bet12a,
              sbet1, -cbet1, dn1, sbet2, cbet2, dn2,
              cbet1, cbet2, nullptr, &m12b, &m0, nullptr, nullptr, Ca);
      x = -1 + m12b / (cbet1 * cbet2 * m0 * pi);
      betscale = x < -0.01 ? sbet12a / x :
        -g->f * sq(cbet1) * pi;
      lamscale = betscale / cbet1;
      y = lam12x / lamscale;
    }

    if (y > -tol1 && x > -1 - xthresh) {
      /* strip near cut */
      if (g->f >= 0) {
        salp1 = fmin(1.0, -x); calp1 = - sqrt(1 - sq(salp1));
      } else {
        calp1 = fmax(x > -tol1 ? 0.0 : -1.0, x);
        salp1 = sqrt(1 - sq(calp1));
      }
    } else {
      /* Estimate alp1, by solving the astroid problem.
       *
       * Could estimate alpha1 = theta + pi/2, directly, i.e.,
       *   calp1 = y/k; salp1 = -x/(1+k);  for f >= 0
       *   calp1 = x/(1+k); salp1 = -y/k;  for f < 0 (need to check)
       *
       * However, it's better to estimate omg12 from astroid and use
       * spherical formula to compute alp1.  This reduces the mean number of
       * Newton iterations for astroid cases from 2.24 (min 0, max 6) to 2.12
       * (min 0 max 5).  The changes in the number of iterations are as
       * follows:
       *
       * change percent
       *    1       5
       *    0      78
       *   -1      16
       *   -2       0.6
       *   -3       0.04
       *   -4       0.002
       *
       * The histogram of iterations is (m = number of iterations estimating
       * alp1 directly, n = number of iterations estimating via omg12, total
       * number of trials = 148605):
       *
       *  iter    m      n
       *    0   148    186
       *    1 13046  13845
       *    2 93315 102225
       *    3 36189  32341
       *    4  5396      7
       *    5   455      1
       *    6    56      0
       *
       * Because omg12 is near pi, estimate work with omg12a = pi - omg12 */
      double k = Astroid(x, y);
      double
        omg12a = lamscale * ( g->f >= 0 ? -x * k/(1 + k) : -y * (1 + k)/k );
      somg12 = sin(omg12a); comg12 = -cos(omg12a);
      /* Update spherical estimate of alp1 using omg12 instead of lam12 */
      salp1 = cbet2 * somg12;
      calp1 = sbet12a - cbet2 * sbet1 * sq(somg12) / (1 - comg12);
    }
  }
  /* Sanity check on starting guess.  Backwards check allows NaN through. */
  if (!(salp1 <= 0))
    norm2(&salp1, &calp1);
  else {
    salp1 = 1; calp1 = 0;
  }

  *psalp1 = salp1;
  *pcalp1 = calp1;
  if (shortline)
    *pdnm = dnm;
  if (sig12 >= 0) {
    *psalp2 = salp2;
    *pcalp2 = calp2;
  }
  return sig12;
}

double Lambda12(const struct geod_geodesic* g,
                double sbet1, double cbet1, double dn1,
                double sbet2, double cbet2, double dn2,
                double salp1, double calp1,
                double slam120, double clam120,
                double* psalp2, double* pcalp2,
                double* psig12,
                double* pssig1, double* pcsig1,
                double* pssig2, double* pcsig2,
                double* peps,
                double* pdomg12,
                boolx diffp, double* pdlam12,
                /* Scratch area of the right size */
                double Ca[]) {
  double salp2 = 0, calp2 = 0, sig12 = 0,
    ssig1 = 0, csig1 = 0, ssig2 = 0, csig2 = 0, eps = 0,
    domg12 = 0, dlam12 = 0;
  double salp0, calp0;
  double somg1, comg1, somg2, comg2, somg12, comg12, lam12;
  double B312, eta, k2;

  if (sbet1 == 0 && calp1 == 0)
    /* Break degeneracy of equatorial line.  This case has already been
     * handled. */
    calp1 = -tiny;

  /* sin(alp1) * cos(bet1) = sin(alp0) */
  salp0 = salp1 * cbet1;
  calp0 = hypot(calp1, salp1 * sbet1); /* calp0 > 0 */

  /* tan(bet1) = tan(sig1) * cos(alp1)
   * tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1) */
  ssig1 = sbet1; somg1 = salp0 * sbet1;
  csig1 = comg1 = calp1 * cbet1;
  norm2(&ssig1, &csig1);
  /* norm2(&somg1, &comg1); -- don't need to normalize! */

  /* Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
   * about this case, since this can yield singularities in the Newton
   * iteration.
   * sin(alp2) * cos(bet2) = sin(alp0) */
  salp2 = cbet2 != cbet1 ? salp0 / cbet2 : salp1;
  /* calp2 = sqrt(1 - sq(salp2))
   *       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
   * and subst for calp0 and rearrange to give (choose positive sqrt
   * to give alp2 in [0, pi/2]). */
  calp2 = cbet2 != cbet1 || fabs(sbet2) != -sbet1 ?
    sqrt(sq(calp1 * cbet1) +
         (cbet1 < -sbet1 ?
          (cbet2 - cbet1) * (cbet1 + cbet2) :
          (sbet1 - sbet2) * (sbet1 + sbet2))) / cbet2 :
    fabs(calp1);
  /* tan(bet2) = tan(sig2) * cos(alp2)
   * tan(omg2) = sin(alp0) * tan(sig2). */
  ssig2 = sbet2; somg2 = salp0 * sbet2;
  csig2 = comg2 = calp2 * cbet2;
  norm2(&ssig2, &csig2);
  /* norm2(&somg2, &comg2); -- don't need to normalize! */

  /* sig12 = sig2 - sig1, limit to [0, pi] */
  sig12 = atan2(fmax(0.0, csig1 * ssig2 - ssig1 * csig2) + 0,
                          csig1 * csig2 + ssig1 * ssig2);

  /* omg12 = omg2 - omg1, limit to [0, pi] */
  somg12 = fmax(0.0, comg1 * somg2 - somg1 * comg2) + 0;
  comg12 =           comg1 * comg2 + somg1 * somg2;
  /* eta = omg12 - lam120 */
  eta = atan2(somg12 * clam120 - comg12 * slam120,
              comg12 * clam120 + somg12 * slam120);
  k2 = sq(calp0) * g->ep2;
  eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2);
  C3f(g, eps, Ca);
  B312 = (SinCosSeries(TRUE, ssig2, csig2, Ca, nC3-1) -
          SinCosSeries(TRUE, ssig1, csig1, Ca, nC3-1));
  domg12 = -g->f * A3f(g, eps) * salp0 * (sig12 + B312);
  lam12 = eta + domg12;

  if (diffp) {
    if (calp2 == 0)
      dlam12 = - 2 * g->f1 * dn1 / sbet1;
    else {
      Lengths(g, eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
              cbet1, cbet2, nullptr, &dlam12, nullptr, nullptr, nullptr, Ca);
      dlam12 *= g->f1 / (calp2 * cbet2);
    }
  }

  *psalp2 = salp2;
  *pcalp2 = calp2;
  *psig12 = sig12;
  *pssig1 = ssig1;
  *pcsig1 = csig1;
  *pssig2 = ssig2;
  *pcsig2 = csig2;
  *peps = eps;
  *pdomg12 = domg12;
  if (diffp)
    *pdlam12 = dlam12;

  return lam12;
}

double A3f(const struct geod_geodesic* g, double eps) {
  /* Evaluate A3 */
  return polyvalx(nA3 - 1, g->A3x, eps);
}

void C3f(const struct geod_geodesic* g, double eps, double c[]) {
  /* Evaluate C3 coeffs
   * Elements c[1] through c[nC3 - 1] are set */
  double mult = 1;
  int o = 0, l;
  for (l = 1; l < nC3; ++l) {   /* l is index of C3[l] */
    int m = nC3 - l - 1;        /* order of polynomial in eps */
    mult *= eps;
    c[l] = mult * polyvalx(m, g->C3x + o, eps);
    o += m + 1;
  }
}

void C4f(const struct geod_geodesic* g, double eps, double c[]) {
  /* Evaluate C4 coeffs
   * Elements c[0] through c[nC4 - 1] are set */
  double mult = 1;
  int o = 0, l;
  for (l = 0; l < nC4; ++l) {   /* l is index of C4[l] */
    int m = nC4 - l - 1;        /* order of polynomial in eps */
    c[l] = mult * polyvalx(m, g->C4x + o, eps);
    o += m + 1;
    mult *= eps;
  }
}

/* The scale factor A1-1 = mean value of (d/dsigma)I1 - 1 */
double A1m1f(double eps)  {
  static const double coeff[] = {
    /* (1-eps)*A1-1, polynomial in eps2 of order 3 */
    1, 4, 64, 0, 256,
  };
  int m = nA1/2;
  double t = polyvalx(m, coeff, sq(eps)) / coeff[m + 1];
  return (t + eps) / (1 - eps);
}

/* The coefficients C1[l] in the Fourier expansion of B1 */
void C1f(double eps, double c[])  {
  static const double coeff[] = {
    /* C1[1]/eps^1, polynomial in eps2 of order 2 */
    -1, 6, -16, 32,
    /* C1[2]/eps^2, polynomial in eps2 of order 2 */
    -9, 64, -128, 2048,
    /* C1[3]/eps^3, polynomial in eps2 of order 1 */
    9, -16, 768,
    /* C1[4]/eps^4, polynomial in eps2 of order 1 */
    3, -5, 512,
    /* C1[5]/eps^5, polynomial in eps2 of order 0 */
    -7, 1280,
    /* C1[6]/eps^6, polynomial in eps2 of order 0 */
    -7, 2048,
  };
  double
    eps2 = sq(eps),
    d = eps;
  int o = 0, l;
  for (l = 1; l <= nC1; ++l) {  /* l is index of C1p[l] */
    int m = (nC1 - l) / 2;      /* order of polynomial in eps^2 */
    c[l] = d * polyvalx(m, coeff + o, eps2) / coeff[o + m + 1];
    o += m + 2;
    d *= eps;
  }
}

/* The coefficients C1p[l] in the Fourier expansion of B1p */
void C1pf(double eps, double c[])  {
  static const double coeff[] = {
    /* C1p[1]/eps^1, polynomial in eps2 of order 2 */
    205, -432, 768, 1536,
    /* C1p[2]/eps^2, polynomial in eps2 of order 2 */
    4005, -4736, 3840, 12288,
    /* C1p[3]/eps^3, polynomial in eps2 of order 1 */
    -225, 116, 384,
    /* C1p[4]/eps^4, polynomial in eps2 of order 1 */
    -7173, 2695, 7680,
    /* C1p[5]/eps^5, polynomial in eps2 of order 0 */
    3467, 7680,
    /* C1p[6]/eps^6, polynomial in eps2 of order 0 */
    38081, 61440,
  };
  double
    eps2 = sq(eps),
    d = eps;
  int o = 0, l;
  for (l = 1; l <= nC1p; ++l) { /* l is index of C1p[l] */
    int m = (nC1p - l) / 2;     /* order of polynomial in eps^2 */
    c[l] = d * polyvalx(m, coeff + o, eps2) / coeff[o + m + 1];
    o += m + 2;
    d *= eps;
  }
}

/* The scale factor A2-1 = mean value of (d/dsigma)I2 - 1 */
double A2m1f(double eps)  {
  static const double coeff[] = {
    /* (eps+1)*A2-1, polynomial in eps2 of order 3 */
    -11, -28, -192, 0, 256,
  };
  int m = nA2/2;
  double t = polyvalx(m, coeff, sq(eps)) / coeff[m + 1];
  return (t - eps) / (1 + eps);
}

/* The coefficients C2[l] in the Fourier expansion of B2 */
void C2f(double eps, double c[])  {
  static const double coeff[] = {
    /* C2[1]/eps^1, polynomial in eps2 of order 2 */
    1, 2, 16, 32,
    /* C2[2]/eps^2, polynomial in eps2 of order 2 */
    35, 64, 384, 2048,
    /* C2[3]/eps^3, polynomial in eps2 of order 1 */
    15, 80, 768,
    /* C2[4]/eps^4, polynomial in eps2 of order 1 */
    7, 35, 512,
    /* C2[5]/eps^5, polynomial in eps2 of order 0 */
    63, 1280,
    /* C2[6]/eps^6, polynomial in eps2 of order 0 */
    77, 2048,
  };
  double
    eps2 = sq(eps),
    d = eps;
  int o = 0, l;
  for (l = 1; l <= nC2; ++l) { /* l is index of C2[l] */
    int m = (nC2 - l) / 2;     /* order of polynomial in eps^2 */
    c[l] = d * polyvalx(m, coeff + o, eps2) / coeff[o + m + 1];
    o += m + 2;
    d *= eps;
  }
}

/* The scale factor A3 = mean value of (d/dsigma)I3 */
void A3coeff(struct geod_geodesic* g) {
  static const double coeff[] = {
    /* A3, coeff of eps^5, polynomial in n of order 0 */
    -3, 128,
    /* A3, coeff of eps^4, polynomial in n of order 1 */
    -2, -3, 64,
    /* A3, coeff of eps^3, polynomial in n of order 2 */
    -1, -3, -1, 16,
    /* A3, coeff of eps^2, polynomial in n of order 2 */
    3, -1, -2, 8,
    /* A3, coeff of eps^1, polynomial in n of order 1 */
    1, -1, 2,
    /* A3, coeff of eps^0, polynomial in n of order 0 */
    1, 1,
  };
  int o = 0, k = 0, j;
  for (j = nA3 - 1; j >= 0; --j) {             /* coeff of eps^j */
    int m = nA3 - j - 1 < j ? nA3 - j - 1 : j; /* order of polynomial in n */
    g->A3x[k++] = polyvalx(m, coeff + o, g->n) / coeff[o + m + 1];
    o += m + 2;
  }
}

/* The coefficients C3[l] in the Fourier expansion of B3 */
void C3coeff(struct geod_geodesic* g) {
  static const double coeff[] = {
    /* C3[1], coeff of eps^5, polynomial in n of order 0 */
    3, 128,
    /* C3[1], coeff of eps^4, polynomial in n of order 1 */
    2, 5, 128,
    /* C3[1], coeff of eps^3, polynomial in n of order 2 */
    -1, 3, 3, 64,
    /* C3[1], coeff of eps^2, polynomial in n of order 2 */
    -1, 0, 1, 8,
    /* C3[1], coeff of eps^1, polynomial in n of order 1 */
    -1, 1, 4,
    /* C3[2], coeff of eps^5, polynomial in n of order 0 */
    5, 256,
    /* C3[2], coeff of eps^4, polynomial in n of order 1 */
    1, 3, 128,
    /* C3[2], coeff of eps^3, polynomial in n of order 2 */
    -3, -2, 3, 64,
    /* C3[2], coeff of eps^2, polynomial in n of order 2 */
    1, -3, 2, 32,
    /* C3[3], coeff of eps^5, polynomial in n of order 0 */
    7, 512,
    /* C3[3], coeff of eps^4, polynomial in n of order 1 */
    -10, 9, 384,
    /* C3[3], coeff of eps^3, polynomial in n of order 2 */
    5, -9, 5, 192,
    /* C3[4], coeff of eps^5, polynomial in n of order 0 */
    7, 512,
    /* C3[4], coeff of eps^4, polynomial in n of order 1 */
    -14, 7, 512,
    /* C3[5], coeff of eps^5, polynomial in n of order 0 */
    21, 2560,
  };
  int o = 0, k = 0, l, j;
  for (l = 1; l < nC3; ++l) {                    /* l is index of C3[l] */
    for (j = nC3 - 1; j >= l; --j) {             /* coeff of eps^j */
      int m = nC3 - j - 1 < j ? nC3 - j - 1 : j; /* order of polynomial in n */
      g->C3x[k++] = polyvalx(m, coeff + o, g->n) / coeff[o + m + 1];
      o += m + 2;
    }
  }
}

/* The coefficients C4[l] in the Fourier expansion of I4 */
void C4coeff(struct geod_geodesic* g) {
  static const double coeff[] = {
    /* C4[0], coeff of eps^5, polynomial in n of order 0 */
    97, 15015,
    /* C4[0], coeff of eps^4, polynomial in n of order 1 */
    1088, 156, 45045,
    /* C4[0], coeff of eps^3, polynomial in n of order 2 */
    -224, -4784, 1573, 45045,
    /* C4[0], coeff of eps^2, polynomial in n of order 3 */
    -10656, 14144, -4576, -858, 45045,
    /* C4[0], coeff of eps^1, polynomial in n of order 4 */
    64, 624, -4576, 6864, -3003, 15015,
    /* C4[0], coeff of eps^0, polynomial in n of order 5 */
    100, 208, 572, 3432, -12012, 30030, 45045,
    /* C4[1], coeff of eps^5, polynomial in n of order 0 */
    1, 9009,
    /* C4[1], coeff of eps^4, polynomial in n of order 1 */
    -2944, 468, 135135,
    /* C4[1], coeff of eps^3, polynomial in n of order 2 */
    5792, 1040, -1287, 135135,
    /* C4[1], coeff of eps^2, polynomial in n of order 3 */
    5952, -11648, 9152, -2574, 135135,
    /* C4[1], coeff of eps^1, polynomial in n of order 4 */
    -64, -624, 4576, -6864, 3003, 135135,
    /* C4[2], coeff of eps^5, polynomial in n of order 0 */
    8, 10725,
    /* C4[2], coeff of eps^4, polynomial in n of order 1 */
    1856, -936, 225225,
    /* C4[2], coeff of eps^3, polynomial in n of order 2 */
    -8448, 4992, -1144, 225225,
    /* C4[2], coeff of eps^2, polynomial in n of order 3 */
    -1440, 4160, -4576, 1716, 225225,
    /* C4[3], coeff of eps^5, polynomial in n of order 0 */
    -136, 63063,
    /* C4[3], coeff of eps^4, polynomial in n of order 1 */
    1024, -208, 105105,
    /* C4[3], coeff of eps^3, polynomial in n of order 2 */
    3584, -3328, 1144, 315315,
    /* C4[4], coeff of eps^5, polynomial in n of order 0 */
    -128, 135135,
    /* C4[4], coeff of eps^4, polynomial in n of order 1 */
    -2560, 832, 405405,
    /* C4[5], coeff of eps^5, polynomial in n of order 0 */
    128, 99099,
  };
  int o = 0, k = 0, l, j;
  for (l = 0; l < nC4; ++l) {        /* l is index of C4[l] */
    for (j = nC4 - 1; j >= l; --j) { /* coeff of eps^j */
      int m = nC4 - j - 1;           /* order of polynomial in n */
      g->C4x[k++] = polyvalx(m, coeff + o, g->n) / coeff[o + m + 1];
      o += m + 2;
    }
  }
}

int transit(double lon1, double lon2) {
  double lon12;
  /* Return 1 or -1 if crossing prime meridian in east or west direction.
   * Otherwise return zero. */
  /* Compute lon12 the same way as Geodesic::Inverse. */
  lon12 = AngDiff(lon1, lon2, nullptr);
  lon1 = AngNormalize(lon1);
  lon2 = AngNormalize(lon2);
  return
    lon12 > 0 && ((lon1 < 0 && lon2 >= 0) ||
                  (lon1 > 0 && lon2 == 0)) ? 1 :
    (lon12 < 0 && lon1 >= 0 && lon2 < 0 ? -1 : 0);
}

int transitdirect(double lon1, double lon2) {
  /* Compute exactly the parity of
   *   int(floor(lon2 / 360)) - int(floor(lon1 / 360)) */
  lon1 = remainder(lon1, 2.0 * td); lon2 = remainder(lon2, 2.0 * td);
  return ( (lon2 >= 0 && lon2 < td ? 0 : 1) -
           (lon1 >= 0 && lon1 < td ? 0 : 1) );
}

void accini(double s[]) {
  /* Initialize an accumulator; this is an array with two elements. */
  s[0] = s[1] = 0;
}

void acccopy(const double s[], double t[]) {
  /* Copy an accumulator; t = s. */
  t[0] = s[0]; t[1] = s[1];
}

void accadd(double s[], double y) {
  /* Add y to an accumulator. */
  double u, z = sumx(y, s[1], &u);
  s[0] = sumx(z, s[0], &s[1]);
  if (s[0] == 0)
    s[0] = u;
  else
    s[1] = s[1] + u;
}

double accsum(const double s[], double y) {
  /* Return accumulator + y (but don't add to accumulator). */
  double t[2];
  acccopy(s, t);
  accadd(t, y);
  return t[0];
}

void accneg(double s[]) {
  /* Negate an accumulator. */
  s[0] = -s[0]; s[1] = -s[1];
}

void accrem(double s[], double y) {
  /* Reduce to [-y/2, y/2]. */
  s[0] = remainder(s[0], y);
  accadd(s, 0.0);
}

void geod_polygon_init(struct geod_polygon* p, boolx polylinep) {
  p->polyline = (polylinep != 0);
  geod_polygon_clear(p);
}

void geod_polygon_clear(struct geod_polygon* p) {
  p->lat0 = p->lon0 = p->lat = p->lon = NaN;
  accini(p->P);
  accini(p->A);
  p->num = p->crossings = 0;
}

void geod_polygon_addpoint(const struct geod_geodesic* g,
                           struct geod_polygon* p,
                           double lat, double lon) {
  if (p->num == 0) {
    p->lat0 = p->lat = lat;
    p->lon0 = p->lon = lon;
  } else {
    double s12, S12 = 0;     /* Initialize S12 to stop Visual Studio warning */
    geod_geninverse(g, p->lat, p->lon, lat, lon,
                    &s12, nullptr, nullptr, nullptr, nullptr, nullptr,
                    p->polyline ? nullptr : &S12);
    accadd(p->P, s12);
    if (!p->polyline) {
      accadd(p->A, S12);
      p->crossings += transit(p->lon, lon);
    }
    p->lat = lat; p->lon = lon;
  }
  ++p->num;
}

void geod_polygon_addedge(const struct geod_geodesic* g,
                          struct geod_polygon* p,
                          double azi, double s) {
  if (p->num) {              /* Do nothing is num is zero */
    /* Initialize S12 to stop Visual Studio warning.  Initialization of lat and
     * lon is to make CLang static analyzer happy. */
    double lat = 0, lon = 0, S12 = 0;
    geod_gendirect(g, p->lat, p->lon, azi, GEOD_LONG_UNROLL, s,
                   &lat, &lon, nullptr,
                   nullptr, nullptr, nullptr, nullptr,
                   p->polyline ? nullptr : &S12);
    accadd(p->P, s);
    if (!p->polyline) {
      accadd(p->A, S12);
      p->crossings += transitdirect(p->lon, lon);
    }
    p->lat = lat; p->lon = lon;
    ++p->num;
  }
}

unsigned geod_polygon_compute(const struct geod_geodesic* g,
                              const struct geod_polygon* p,
                              boolx reverse, boolx sign,
                              double* pA, double* pP) {
  double s12, S12, t[2];
  if (p->num < 2) {
    if (pP) *pP = 0;
    if (!p->polyline && pA) *pA = 0;
    return p->num;
  }
  if (p->polyline) {
    if (pP) *pP = p->P[0];
    return p->num;
  }
  geod_geninverse(g, p->lat, p->lon, p->lat0, p->lon0,
                  &s12, nullptr, nullptr, nullptr, nullptr, nullptr, &S12);
  if (pP) *pP = accsum(p->P, s12);
  acccopy(p->A, t);
  accadd(t, S12);
  if (pA) *pA = areareduceA(t, 4 * pi * g->c2,
                            p->crossings + transit(p->lon, p->lon0),
                            reverse, sign);
  return p->num;
}

unsigned geod_polygon_testpoint(const struct geod_geodesic* g,
                                const struct geod_polygon* p,
                                double lat, double lon,
                                boolx reverse, boolx sign,
                                double* pA, double* pP) {
  double perimeter, tempsum;
  int crossings, i;
  unsigned num = p->num + 1;
  if (num == 1) {
    if (pP) *pP = 0;
    if (!p->polyline && pA) *pA = 0;
    return num;
  }
  perimeter = p->P[0];
  tempsum = p->polyline ? 0 : p->A[0];
  crossings = p->crossings;
  for (i = 0; i < (p->polyline ? 1 : 2); ++i) {
    double s12, S12 = 0;     /* Initialize S12 to stop Visual Studio warning */
    geod_geninverse(g,
                    i == 0 ? p->lat  : lat, i == 0 ? p->lon  : lon,
                    i != 0 ? p->lat0 : lat, i != 0 ? p->lon0 : lon,
                    &s12, nullptr, nullptr, nullptr, nullptr, nullptr,
                    p->polyline ? nullptr : &S12);
    perimeter += s12;
    if (!p->polyline) {
      tempsum += S12;
      crossings += transit(i == 0 ? p->lon  : lon,
                           i != 0 ? p->lon0 : lon);
    }
  }

  if (pP) *pP = perimeter;
  if (p->polyline)
    return num;

  if (pA) *pA = areareduceB(tempsum, 4 * pi * g->c2, crossings, reverse, sign);
  return num;
}

unsigned geod_polygon_testedge(const struct geod_geodesic* g,
                               const struct geod_polygon* p,
                               double azi, double s,
                               boolx reverse, boolx sign,
                               double* pA, double* pP) {
  double perimeter, tempsum;
  int crossings;
  unsigned num = p->num + 1;
  if (num == 1) {               /* we don't have a starting point! */
    if (pP) *pP = NaN;
    if (!p->polyline && pA) *pA = NaN;
    return 0;
  }
  perimeter = p->P[0] + s;
  if (p->polyline) {
    if (pP) *pP = perimeter;
    return num;
  }

  tempsum = p->A[0];
  crossings = p->crossings;
  {
    /* Initialization of lat, lon, and S12 is to make CLang static analyzer
     * happy. */
    double lat = 0, lon = 0, s12, S12 = 0;
    geod_gendirect(g, p->lat, p->lon, azi, GEOD_LONG_UNROLL, s,
                   &lat, &lon, nullptr,
                   nullptr, nullptr, nullptr, nullptr, &S12);
    tempsum += S12;
    crossings += transitdirect(p->lon, lon);
    geod_geninverse(g, lat,  lon, p->lat0,  p->lon0,
                    &s12, nullptr, nullptr, nullptr, nullptr, nullptr, &S12);
    perimeter += s12;
    tempsum += S12;
    crossings += transit(lon, p->lon0);
  }

  if (pP) *pP = perimeter;
  if (pA) *pA = areareduceB(tempsum, 4 * pi * g->c2, crossings, reverse, sign);
  return num;
}

void geod_polygonarea(const struct geod_geodesic* g,
                      double lats[], double lons[], int n,
                      double* pA, double* pP) {
  int i;
  struct geod_polygon p;
  geod_polygon_init(&p, FALSE);
  for (i = 0; i < n; ++i)
    geod_polygon_addpoint(g, &p, lats[i], lons[i]);
  geod_polygon_compute(g, &p, FALSE, TRUE, pA, pP);
}

double areareduceA(double area[], double area0,
                   int crossings, boolx reverse, boolx sign) {
  accrem(area, area0);
  if (crossings & 1)
    accadd(area, (area[0] < 0 ? 1 : -1) * area0/2);
  /* area is with the clockwise sense.  If !reverse convert to
   * counter-clockwise convention. */
  if (!reverse)
    accneg(area);
  /* If sign put area in (-area0/2, area0/2], else put area in [0, area0) */
  if (sign) {
    if (area[0] > area0/2)
      accadd(area, -area0);
    else if (area[0] <= -area0/2)
      accadd(area, +area0);
  } else {
    if (area[0] >= area0)
      accadd(area, -area0);
    else if (area[0] < 0)
      accadd(area, +area0);
  }
  return 0 + area[0];
}

double areareduceB(double area, double area0,
                   int crossings, boolx reverse, boolx sign) {
  area = remainder(area, area0);
    if (crossings & 1)
    area += (area < 0 ? 1 : -1) * area0/2;
  /* area is with the clockwise sense.  If !reverse convert to
   * counter-clockwise convention. */
  if (!reverse)
    area *= -1;
  /* If sign put area in (-area0/2, area0/2], else put area in [0, area0) */
  if (sign) {
    if (area > area0/2)
      area -= area0;
    else if (area <= -area0/2)
      area += area0;
  } else {
    if (area >= area0)
      area -= area0;
    else if (area < 0)
      area += area0;
  }
  return 0 + area;
}

/** @endcond */
#pragma endregion geodesic.c

static double sign(double x) {
    return x < 0 ? -1 : x > 0 ? 1 : 0;
}

static double wrap_360(double degs) {
    if (degs < 0 || degs >= 360) {
        degs = fmod(degs, 360);
        if (degs < 0) {
            degs += 360;
        } else if (degs >= 360) {
            degs -= 360;
        }
    }
    return degs;
}

static double sphere_dist(double lat1, double lon1, double lat2, double lon2) {
    // haversine
    double rlat1 = lat1 * RADIANS;
    double rlat2 = lat2 * RADIANS;
    double rlat = (lat2-lat1) * RADIANS;
    double rlon = (lon2-lon1) * RADIANS;
    double a = sin(rlat/2) * sin(rlat/2) + 
        cos(rlat1) * cos(rlat2) * sin(rlon/2) * sin(rlon/2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    return EARTH * c;
}

static double sphere_initial_azi(double lat1, double lon1, double lat2,
    double lon2)
{
    double rlat1 = lat1 * RADIANS;
    double rlat2 = lat2 * RADIANS;
    double rlon = (lon2 - lon1) * RADIANS;
    double y = sin(rlon) * cos(rlat2);
    double x = cos(rlat1) * sin(rlat2) - sin(rlat1) * cos(rlat2) * cos(rlon);
    double a = atan2(y, x);
    double bearing = a * DEGREES;
    return wrap_360(bearing);
}

static double sphere_final_azi(double lat1, double lon1, double lat2, 
    double lon2)
{
    return wrap_360(sphere_initial_azi(lat2, lon2, lat1, lon1) + 180.0);
}

static void sphere_dest(double lat1, double lon1, double azi1, double s12,
    double *lat2, double *lon2)
{
    double dr = s12 / EARTH;
    double ar = azi1 * RADIANS;
    double rlat1 = lat1 * RADIANS;
    double rlon1 = lon1 * RADIANS;
    double phi = asin(sin(rlat1) * cos(dr) + cos(rlat1) * sin(dr) * cos(ar));
    double lambda = rlon1 + atan2(sin(ar) * sin(dr) * cos(rlat1), cos(dr) - 
        sin(rlat1) * sin(phi));
    lambda = fmod(lambda+3*M_PI, 2*M_PI) - M_PI; // normalise to -180..+180°
    *lat2 = phi * DEGREES;
    *lon2 = lambda * DEGREES;
}

static void geokit_unused(void) {
    (void)geokit_unused;
    (void)geod_polygonarea;
    (void)geod_polygon_testpoint;
    (void)geod_polygon_testedge;
    (void)geod_polygon_addedge;
    (void)geod_position;
    (void)geod_inverseline;
    (void)geod_directline;
    (void)geod_inverse;
}

///////////////////////////////////////////////////////////////////////////////

static _Thread_local bool geodesic_initd = false;
static _Thread_local struct geod_geodesic geodesic;

static struct geod_geodesic *local_geod(void) {
    if (!geodesic_initd) {
        geod_init(&geodesic, EARTH, FLATTENING);
        geodesic_initd = true;
    }
    return &geodesic;
}

static double inverse_ellipsoid(double lat1, double lon1, double lat2, 
    double lon2, double* ps12, double* pazi1, double* pazi2)
{
    double arc = geod_geninverse(local_geod(), lat1, lon1, lat2, lon2, ps12,
        pazi1, pazi2, 0, 0, 0, 0);
    if (pazi1) *pazi1 = wrap_360(*pazi1);
    if (pazi2) *pazi2 = wrap_360(*pazi2);
    return arc;
}

static void direct_ellipsoid(double lat1, double lon1, double azi1, double s12,
    double* plat2, double* plon2, double* pazi2)
{
    geod_direct(local_geod(), lat1, lon1, azi1, s12, plat2, plon2, pazi2);
}

static void inverse_sphere(double lat1, double lon1, double lat2, double lon2,
    double* ps12, double* pazi1, double* pazi2)
{
    if (ps12) *ps12 = sphere_dist(lat1, lon1, lat2, lon2);
    if (pazi1) *pazi1 = sphere_initial_azi(lat1, lon1, lat2, lon2);
    if (pazi2) *pazi2 = sphere_final_azi(lat1, lon1, lat2, lon2);
}

static void direct_sphere(double lat1, double lon1, double azi1, double s12,
    double* plat2, double* plon2, double* pazi2)
{
    double lat2, lon2;
    sphere_dest(lat1, lon1, azi1, s12, &lat2, &lon2);
    if (plat2) *plat2 = lat2;
    if (plon2) *plon2 = lon2;
    if (pazi2) *pazi2 = sphere_final_azi(lat1, lon1, lat2, lon2);
}

void geokit_inverse(enum geokit_shape shape, double lat1, double lon1, double lat2, 
    double lon2, double* ps12, double* pazi1, double* pazi2)
{
    (void)geokit_unused;
    if (shape == GEOKIT_ELLIPSOID) {
        inverse_ellipsoid(lat1, lon1, lat2, lon2, ps12, pazi1, pazi2);
    } else {
        inverse_sphere(lat1, lon1, lat2, lon2, ps12, pazi1, pazi2);
    }
}

void geokit_direct(enum geokit_shape shape, double lat1, double lon1, double azi1, 
    double s12, double* plat2, double* plon2, double* pazi2)
{
    if (shape == GEOKIT_ELLIPSOID) {
        direct_ellipsoid(lat1, lon1, azi1, s12, plat2, plon2, pazi2);
    } else {
        direct_sphere(lat1, lon1, azi1, s12, plat2, plon2, pazi2);
    }
}

static void track_ellipsoid_correction(double lat1, double lon1, double azi1, 
    double lat2, double lon2, double *psat, double *psxt, double *plat3, 
    double *plon3)
{
    // This corrects the geokit_track sat and sxt values using the bisection
    // method, terminating when the sat window is within the desired accuracy.
    // The initial window is approximated through probing.
    const double INIT_WINDOW = 2;      // initial sat window in meters.
    const double PROBE_MULTIPLIER = 2; // window multiplier when probing.
    double sat = *psat;
    double sxt = *psxt;
    double lat3, lon3, initial_dist;
    double initial_sat = sat;
    direct_ellipsoid(lat1, lon1, azi1, initial_sat, &lat3, &lon3, 0);
    inverse_ellipsoid(lat2, lon2, lat3, lon3, &initial_dist, 0, 0);
    double min_sat, max_sat;
    double min_dist, max_dist;
    double probe;
    probe = INIT_WINDOW/2;
    min_sat = sat-probe;
    direct_ellipsoid(lat1, lon1, azi1, min_sat, &lat3, &lon3, 0);
    inverse_ellipsoid(lat2, lon2, lat3, lon3, &min_dist, 0, 0);
    while (min_dist < initial_dist) {
        probe *= PROBE_MULTIPLIER;
        min_sat = sat-probe;
        direct_ellipsoid(lat1, lon1, azi1, min_sat, &lat3, &lon3, 0);
        inverse_ellipsoid(lat2, lon2, lat3, lon3, &min_dist, 0, 0);
    }
    probe = INIT_WINDOW/2;
    max_sat = sat+probe;
    direct_ellipsoid(lat1, lon1, azi1, max_sat, &lat3, &lon3, 0);
    inverse_ellipsoid(lat2, lon2, lat3, lon3, &max_dist, 0, 0);
    while (max_dist < initial_dist) {
        probe *= PROBE_MULTIPLIER;
        max_sat = sat+probe;
        direct_ellipsoid(lat1, lon1, azi1, max_sat, &lat3, &lon3, 0);
        inverse_ellipsoid(lat2, lon2, lat3, lon3, &max_dist, 0, 0);
    }
    double mid_sat = sat;
    double mid_dist = initial_dist;
    int iters = 0;
    while (iters < MAX_ITERS && max_sat - min_sat > ACCURACY) {
        mid_sat = (max_sat+min_sat) / 2;
        direct_ellipsoid(lat1, lon1, azi1, mid_sat, &lat3, &lon3, 0);
        inverse_ellipsoid(lat2, lon2, lat3, lon3, &mid_dist, 0, 0);
        if (min_dist < max_dist) {
            max_sat = mid_sat;
            max_dist = mid_dist;
        } else {
            min_sat = mid_sat;
            min_dist = mid_dist;
        }
        iters++;
    }
    sat = mid_sat;
    sxt = mid_dist;
    if (iters == 0) {
        direct_ellipsoid(lat1, lon1, azi1, sat, &lat3, &lon3, 0);
    }
    *psat = sat;
    *psxt = sxt;
    *plat3 = lat3;
    *plon3 = lon3;
}

static void geokit_track0(enum geokit_shape shape, double lat1, double lon1,
    double azi1, double lat2, double lon2, double *psat, double *psxt,
    double *plat3, double *plon3)
{
    if (!psat && !psxt) {
        return;
    }
    if (lat1 == lat2 && lon1 == lon2) {
        if (psat) *psat = 0;
        if (psxt) *psxt = 0;
        if (plat3) *plat3 = lat1;
        if (plon3) *plon3 = lon1;
        return;
    }

    // Perform the cross track formula and along track formula
    // See: http://www.movable-type.co.uk/scripts/latlong.html#cross-track
    double s13, a13, δ13, θ13, θ12, δxt, δat;
    geokit_inverse(shape, lat1, lon1, lat2, lon2, &s13, &a13, 0);
    δ13 = s13 / EARTH;
    θ13 = a13 * RADIANS;
    θ12 = azi1 * RADIANS;
    δxt = asin(sin(δ13) * sin(θ13-θ12));
    δat = acos(cos(δ13) / fabs(cos(δxt)));
    double sxt = δxt * EARTH;
    double sat = δat * sign(cos(θ12-θ13)) * EARTH;

    // The sxt and sat now hold the cross track and along track distances.
    double lat3, lon3;
    if (shape == GEOKIT_ELLIPSOID) {
        // A spherical formula was used to get δxt and δat, resulting in the
        // xat and sat values being slightly off for WGS84 ellipsoid.
        // The values need correcting.
        track_ellipsoid_correction(lat1, lon1, azi1, lat2, lon2, &sat, &sxt, 
            &lat3, &lon3);
    } else if (plat3 || plon3) {
        geokit_direct(shape, lat1, lon1, azi1, sat, &lat3, &lon3, 0);
    }
    sxt = fabs(sxt);
    if (psat) *psat = sat;
    if (psxt) *psxt = sxt;
    if (plat3) *plat3 = lat3;
    if (plon3) *plon3 = lon3;
}


void geokit_track(enum geokit_shape shape, double lat1, double lon1, double azi1,
    double lat2, double lon2, double *psat, double *psxt)
{
    geokit_track0(shape, lat1, lon1, azi1, lat2, lon2, psat, psxt, 0, 0);
}

// azi_to_perc returns the position, as a percentage, of azi with the range
// of aazi and bazi. Thus is dazi == aazi then the result is 0.0, and if 
// dazi == bazi then the result is 1.0.
static double azi_to_perc(double azi, double aazi, double bazi) {
    if (fabs((aazi+360)-(bazi+360)) > 180) {
        if (azi < 180) {
            azi += 360;
        }
        if (aazi > bazi) {
            return (azi-aazi) / (bazi+360-aazi);
        } else {
            return 1 - (azi-bazi) / (aazi+360-bazi);
        }
    } else {
        return ((azi-aazi) / (bazi-aazi));
    }
}

// geokit_intersect returns true if line AB intersects line CD.
// Also returns the point of interection, if true.
static bool geokit_intersect0(enum geokit_shape shape,
    double alat, double alon, double blat, double blon,
    double clat, double clon, double dlat, double dlon,
    double *plat, double *plon,
    /* extra fields for composite formulas */
    double *pab_s12, double *pab_azi1, double *pab_azi2,
    double *pcd_s12, double *pcd_azi1, double *pcd_azi2)
{
    // Check if the AB line is directly connected to CD. If so, then the lines
    // are considered interecting and the result lat/lon position is the
    // connected point.
    if ((alat == clat && alon == clon) || (alat == dlat && alon == dlon)) {
        if (plat) *plat = alat;
        if (plon) *plon = alon;
        if (pab_s12 || pab_azi1 || pab_azi2) {
            geokit_inverse(shape, alat, alon, blat, blon, pab_s12, pab_azi1, 
                pab_azi2);
        }
        if (pcd_s12 || pcd_azi1 || pcd_azi2) {
            geokit_inverse(shape, clat, clon, dlat, dlon, pcd_s12, pcd_azi1, 
                pcd_azi2);
        }
        return true;
    }
    if ((blat == clat && blon == clon) || (blat == dlat && blon == dlon)) {
        if (plat) *plat = blat;
        if (plon) *plon = blon;
        if (pab_s12 || pab_azi1 || pab_azi2) {
            geokit_inverse(shape, alat, alon, blat, blon, pab_s12, pab_azi1, 
                pab_azi2);
        }
        if (pcd_s12 || pcd_azi1 || pcd_azi2) {
            geokit_inverse(shape, clat, clon, dlat, dlon, pcd_s12, pcd_azi1, 
                pcd_azi2);
        }
        return true;
    }


    double ab_s12, ab_azi1;
    geokit_inverse(shape, alat, alon, blat, blon, &ab_s12, &ab_azi1, pab_azi2);
    if (pab_s12) *pab_s12 = ab_s12;
    if (pab_azi1) *pab_azi1 = ab_azi1;
    if (clat == dlat && clon == dlon) {
        // CD is a point, not a line, and has no real direction. Treat
        // this condition like an along track operation.
        if (pcd_s12 || pcd_azi1 || pcd_azi2) {
            geokit_inverse(shape, clat, clon, dlat, dlon, pcd_s12, pcd_azi1, 
                pcd_azi2);
        }
        double sat, sxt;
        geokit_track(shape, alat, alon, ab_azi1, clat, clon, &sat, &sxt);
        if (sxt > ACCURACY || sat < 0 || sat > ab_s12) {
            if (plat) *plat = 0;
            if (plon) *plon = 0;
            return false;
        }
        geokit_direct(shape, alat, alon, ab_azi1, sat, plat, plon, 0);
        return true;
    }

    // The intersecting algorithm is below.
    double cd_s12, cd_azi1;
    geokit_inverse(shape, clat, clon, dlat, dlon, &cd_s12, &cd_azi1, pcd_azi2);
    if (pcd_s12) *pcd_s12 = cd_s12;
    if (pcd_azi1) *pcd_azi1 = cd_azi1;
    double ca_s12, ca_azi1;
    geokit_inverse(shape, clat, clon, alat, alon, &ca_s12, &ca_azi1, 0);
    double cb_s12, cb_azi1;
    geokit_inverse(shape, clat, clon, blat, blon, &cb_s12, &cb_azi1, 0);
    double perc = azi_to_perc(cd_azi1, ca_azi1, cb_azi1);
    bool inside = !(perc < 0 || perc > 1);
    if (!inside) {
        if (plat) *plat = 0;
        if (plon) *plon = 0;
        return false;
    }
    double elat = alat;
    double elon = alon;
    double ce_s12 = ca_s12;
    double ce_azi1 = ca_azi1;
    double min = 0;
    double max = ab_s12;
    int iters = 0;
    while (iters < MAX_ITERS && max-min > ACCURACY) {
        double mid = (max+min)/2;
        geokit_direct(shape, alat, alon, ab_azi1, mid, &elat, &elon, 0);
        geokit_inverse(shape, clat, clon, elat, elon, &ce_s12, &ce_azi1, 0);
        double perc2 = azi_to_perc(ce_azi1, ca_azi1, cb_azi1);
        if (perc2 < perc) {
            min = mid;
        } else {
            max = mid;
        }
        iters++;
    }
    if (cd_s12 < ce_s12) {
        if (plat) *plat = 0;
        if (plon) *plon = 0;
        return false;
    }
    if (plat) *plat = elat;
    if (plon) *plon = elon;
    return true;
}

bool geokit_intersect(enum geokit_shape shape,
    double alat, double alon, double blat, double blon, 
    double clat, double clon, double dlat, double dlon,
    double *plat, double *plon)
{
    return geokit_intersect0(shape, alat, alon, blat, blon, clat, clon, dlat, dlon, 
        plat, plon, 0, 0, 0, 0, 0, 0);
}


static void track_nearest(enum geokit_shape shape, double alat, double alon,
    double blat, double blon, double ab_s12,double ab_azi1, double clat, 
    double clon, double *plat, double *plon, double *ps12, double *pazi1,
    double *pazi2)
{
    double lat, lon, sat, sxt;
    geokit_track0(shape, alat, alon, ab_azi1, clat, clon, &sat, &sxt, &lat, &lon);
    if (sat < 0) {
        lat = alat;
        lon = alon;
    } else if (sat > ab_s12) {
        lat = blat;
        lon = blon;
    } else {
        geokit_direct(shape, alat, alon, ab_azi1, sat, &lat, &lon, 0);
    }
    double s12, azi1, azi2;
    geokit_inverse(shape, clat, clon, lat, lon, &s12, &azi1, &azi2);
    if (plat) *plat = lat;
    if (plon) *plon = lon;
    if (ps12) *ps12 = s12;
    if (pazi1) *pazi1 = azi1;
    if (pazi2) *pazi2 = azi2;
}

void geokit_nearest(enum geokit_shape shape,
    double alat, double alon, double blat, double blon, 
    double clat, double clon, double dlat, double dlon,
    double *plat1, double *plon1, double *plat2, double *plon2, 
    double *ps12, double *pazi1, double *pazi2)
{
    double lat, lon;
    double ab_s12, ab_azi1, ab_azi2;
    double cd_s12, cd_azi1, cd_azi2;
    bool isect = geokit_intersect0(shape, alat, alon, blat, blon, clat, clon, 
        dlat, dlon, &lat, &lon, &ab_s12, &ab_azi1, &ab_azi2,
        &cd_s12, &cd_azi1, &cd_azi2);
    if (isect) {
        // intersecting
        if (plat1) *plat1 = lat;
        if (plon1) *plon1 = lon;
        if (plat2) *plat2 = lat;
        if (plon2) *plon2 = lon;
        if (ps12) *ps12 = 0;
        if (pazi1) *pazi1 = 0;
        if (pazi2) *pazi2 = 0;
        return;
    }
    double lat1, lon1, lat2, lon2, s12, azi1, azi2;
    double ts12, tazi1, tazi2;
    track_nearest(shape, clat, clon, dlat, dlon, cd_s12, cd_azi1, alat, alon,
        &lat, &lon, &ts12, &tazi1, &tazi2);
    if (1) { 
        lat1 = alat;
        lon1 = alon;
        lat2 = lat;
        lon2 = lon;
        s12 = ts12;
        azi1 = tazi1;
        azi2 = tazi2;
    }
    track_nearest(shape, clat, clon, dlat, dlon, cd_s12, cd_azi1, blat, blon,
        &lat, &lon, &ts12, &tazi1, &tazi2);
    if (ts12 < s12) { 
        lat1 = blat;
        lon1 = blon;
        lat2 = lat;
        lon2 = lon;
        s12 = ts12;
        azi1 = tazi1;
        azi2 = tazi2;
    }
    track_nearest(shape, alat, alon, blat, blon, ab_s12, ab_azi1, clat, clon,
        &lat, &lon, &ts12, &tazi1, &tazi2);
    if (ts12 < s12) { 
        lat1 = lat;
        lon1 = lon;
        lat2 = clat;
        lon2 = clon;
        s12 = ts12;
        azi1 = wrap_360(tazi1+180);
        azi2 = wrap_360(tazi2+180);
    }
    track_nearest(shape, alat, alon, blat, blon, ab_s12, ab_azi1, dlat, dlon,
        &lat, &lon, &ts12, &tazi1, &tazi2);
    if (ts12 < s12) { 
        lat1 = lat;
        lon1 = lon;
        lat2 = dlat;
        lon2 = dlon;
        s12 = ts12;
        azi1 = wrap_360(tazi1+180);
        azi2 = wrap_360(tazi2+180);
    }
    if (plat1) *plat1 = lat1;
    if (plon1) *plon1 = lon1;
    if (plat2) *plat2 = lat2;
    if (plon2) *plon2 = lon2;
    if (ps12) *ps12 = s12;
    if (pazi1) *pazi1 = azi1;
    if (pazi2) *pazi2 = azi2;
}

struct seg {
    bool filled;
    double sat0;
    double lat0;
    double lon0;
    double sat1;
    double lat1;
    double lon1;  
};

struct res {
    double sat;
    double sxt;
    double lat;
    double lon;
};

static struct seg *fill_seg(enum geokit_shape shape, struct seg segs[], int j, 
    double lat0, double lon0, double azi0)
{
    if (!segs[j].filled) {
        const double EQUATOR = 2 * M_PI * EARTH;
        const double SEGMENT_DIST = EQUATOR / 3;
        double lat1 = lat0;
        double lon1 = lon0;
        double sat1 = 0;
        for (int i = 0; i <= j; i++) {
            double sat0 = sat1;
            sat1 += SEGMENT_DIST;
            if (!segs[i].filled) {
                if (i == 2) {
                    sat1 += SEGMENT_DIST*0.10;
                }
                double lat2, lon2;
                geokit_direct(shape, lat0, lon0, azi0, sat1, &lat2, &lon2, 0);
                
                segs[i].sat0 = sat0;
                segs[i].lat0 = lat1;
                segs[i].lon0 = lon1;
                segs[i].sat1 = sat1;
                segs[i].lat1 = lat2;
                segs[i].lon1 = lon2;
                segs[i].filled = true;
            }
            lat1 = segs[i].lat1;
            lon1 = segs[i].lon1;
        }
    }
    return &segs[j];
}

static void geokit_intercept0(enum geokit_shape shape, double lat1, double lon1, 
    double azi1, double lat2, double lon2, double azi2, double *psat,
    double *psxt, double *plat3, double *plon3)
{
    // intersect test order: 00 01 10 11 21 12 22 20 02
    struct seg segs[2][3] = {{{0},{0},{0}},{{0},{0},{0}}};
    int order[9][2] = {{0,0},{0,1},{1,0},{1,1},{2,1},{1,2},{2,2},{2,0},{0,2}};
    int nresults = 0;
    struct res results[2] = { 0 };
    for (int i = 0; i < 9; i++) {
        int j = order[i][0];
        int k = order[i][1];
        struct seg *seg0 = fill_seg(shape, segs[0], j, lat1, lon1, azi1);
        struct seg *seg1 = fill_seg(shape, segs[1], k, lat2, lon2, azi2);
        double lat, lon;
        bool isect = geokit_intersect(shape, seg0->lat0, seg0->lon0, seg0->lat1, 
            seg0->lon1, seg1->lat0, seg1->lon0, seg1->lat1, seg1->lon1,
            &lat, &lon);
        if (isect) {
            double s12, s13;
            geokit_inverse(shape, seg0->lat0, seg0->lon0, lat, lon, &s12, 0, 0);
            geokit_inverse(shape, seg1->lat0, seg1->lon0, lat, lon, &s13, 0, 0);
            results[nresults++] = (struct res) {
                .sat = seg0->sat0 + s12,
                .sxt = seg1->sat0 + s13,
                .lat = lat,
                .lon = lon,
            };
            if (nresults == 2) {
                break;
            }
        }
    }
    struct res *a = &results[0];
    if (nresults > 0) {
        struct res *b = &results[1];
        if (MAX(a->sat, a->sxt) - MAX(b->sat, b->sxt) > 0) {
            a = b;
        }
    }
    if (psat) *psat = a->sat;
    if (psxt) *psxt = a->sxt;
    if (plat3) *plat3 = a->lat;
    if (plon3) *plon3 = a->lon;
}

void geokit_intercept(enum geokit_shape shape,
    double lat1, double lon1, double azi1,
    double lat2, double lon2, double azi2,
    double *psat, double *psxt)
{
    geokit_intercept0(shape, lat1, lon1, azi1, lat2, lon2, azi2, psat, psxt, 0, 0);
}

void geokit_equibox(enum geokit_shape shape, 
    double alat, double alon, double blat, double blon, 
    double *pminlat, double *pminlon, double *pmaxlat, double *pmaxlon)
{
    double minlat = fmin(alat, blat);
    double minlon = fmin(alon, blon);
    double maxlat = fmax(alat, blat);
    double maxlon = fmax(alon, blon);
    double azi1, s12;
    if (maxlat > 0) {
        geokit_inverse(shape, maxlat, minlon, maxlat, maxlon, &s12, &azi1, 0);
        geokit_direct(shape, maxlat, minlon, azi1, s12/2, &maxlat, 0, 0);
    }
    if (minlat < 0) {
        geokit_inverse(shape, minlat, minlon, minlat, maxlon, &s12, &azi1, 0);
        geokit_direct(shape, minlat, minlon, azi1, s12/2, &minlat, 0, 0);
    }
    if (maxlon-minlon > 180.0) {
        minlon = -180.0;
        maxlon = 180.0;
    }
    if (pminlat) *pminlat = minlat;
    if (pminlon) *pminlon = minlon;
    if (pmaxlat) *pmaxlat = maxlat;
    if (pmaxlon) *pmaxlon = maxlon;
}

// void geokit_


struct sphere_polygon {
    bool polyline;
    double lat_set;
    double lon_set;
    double lat0;
    double lon0;
    double lat1;
    double lon1;
    size_t npoints;
};

struct polygon {
    enum geokit_shape shape;
    bool computed;
    double area;
    double perimeter;
    union {
       struct geod_polygon epoly;
       struct sphere_polygon spoly;
    };
};

static_assert(
    sizeof(struct geokit_polygon) >= sizeof(struct polygon), 
    "geokit_polygon is too small");

static void polygon_init_ellipsoid(struct polygon *poly, bool polyline) {
    geod_polygon_init(&poly->epoly, polyline);
}

static void polygon_init_sphere(struct polygon *poly, bool polyline) {
    poly->spoly.polyline = polyline;
    poly->spoly.lat_set = 0;
    poly->spoly.lon_set = 0;
    poly->spoly.lat0 = 0;
    poly->spoly.lon0 = 0;
    poly->spoly.lat1 = 0;
    poly->spoly.lon1 = 0;
    poly->spoly.npoints = 0;
}
static void polygon_add_ellipsoid(struct polygon *poly, double lat,
    double lon)
{
    geod_polygon_addpoint(local_geod(), &poly->epoly, lat, lon);
}
static void polygon_add_sphere(struct polygon *poly, double lat,
    double lon)
{
    poly->spoly.lat1 = poly->spoly.lat0;
    poly->spoly.lon1 = poly->spoly.lon0;
    poly->spoly.lat0 = lat;
    poly->spoly.lon0 = lon;
    poly->spoly.npoints++;
    if (poly->spoly.npoints == 1) {
        poly->spoly.lat_set = lat;
        poly->spoly.lon_set = lon;
    } else {
        if (!poly->spoly.polyline) {
            // Karney's method: Area of a Spherical Polygon
            double φ1 = poly->spoly.lat0 * RADIANS;
            double φ2 = poly->spoly.lat1 * RADIANS;
            double Δλ = (poly->spoly.lon1 - poly->spoly.lon0) * RADIANS;
            double φ1T2 = tan(φ1/2);
            double φ2T2 = tan(φ2/2);
            double E = 2 * atan2(tan(Δλ/2) * (φ1T2+φ2T2), 1 + φ1T2*φ2T2);
            poly->area += E;
        }
        poly->perimeter += sphere_dist(poly->spoly.lat0, poly->spoly.lon0, 
            poly->spoly.lat1, poly->spoly.lon1);
    }
}

static void polygon_compute_ellipsoid(struct polygon *poly) {
    if (poly->computed) {
        return;
    }
    geod_polygon_compute(local_geod(), &poly->epoly, 0, 1, &poly->area, 
        &poly->perimeter);
    poly->area = fabs(poly->area);
    poly->computed = true;
}
static void polygon_compute_sphere(struct polygon *poly) {
    if (poly->computed) {
        return;
    }
    if (poly->spoly.npoints > 0) {
        if (!poly->spoly.polyline) {
            polygon_add_sphere(poly, poly->spoly.lat_set, poly->spoly.lon_set);
            poly->area = fabs(poly->area * EARTH * EARTH);
        }
    }
    poly->computed = true;
}

static double polygon_area_ellipsoid(struct polygon *poly) {
    polygon_compute_ellipsoid(poly);
    return poly->area;
}

static double polygon_area_sphere(struct polygon *poly) { 
    polygon_compute_sphere(poly);
    return poly->area;
}

static double polygon_perimeter_ellipsoid(struct polygon *poly) {
    polygon_compute_ellipsoid(poly);
    return poly->perimeter;
}

static double polygon_perimeter_sphere(struct polygon *poly) {
    polygon_compute_sphere(poly);
    return poly->perimeter;
}

void geokit_polygon_init(struct geokit_polygon *polygon, enum geokit_shape shape, 
    bool polyline)
{
    struct polygon *poly = (struct polygon*)polygon;
    poly->shape = shape;
    poly->computed = false;
    poly->area = 0;
    poly->perimeter = 0;
    if (poly->shape == GEOKIT_ELLIPSOID) {
        polygon_init_ellipsoid(poly, polyline);
    } else {
        polygon_init_sphere(poly, polyline);
    }
}

void geokit_polygon_add(struct geokit_polygon *polygon, double lat, double lon) {
    struct polygon *poly = (struct polygon*)polygon;
    if (poly->shape == GEOKIT_ELLIPSOID) {
        polygon_add_ellipsoid(poly, lat, lon);
    } else {
        polygon_add_sphere(poly, lat, lon);
    }
}

double geokit_polygon_area(struct geokit_polygon *polygon) {
    struct polygon *poly = (struct polygon*)polygon;
    if (poly->shape == GEOKIT_ELLIPSOID) {
        return polygon_area_ellipsoid(poly);
    } else {
        return polygon_area_sphere(poly);
    }
}

double geokit_polygon_perimeter(struct geokit_polygon *polygon) {
    struct polygon *poly = (struct polygon*)polygon;
    if (poly->shape == GEOKIT_ELLIPSOID) {
        return polygon_perimeter_ellipsoid(poly);
    } else {
        return polygon_perimeter_sphere(poly);
    }
}

static void geokit_flatlat(enum geokit_shape shape,
    double alat, double alon1, double alon2, 
    double blat, double blon1, double blon2,
    double *plat1, double *plon1, double *plat2, double *plon2,
    double *ps12, double *pazi1, double *pazi2)
{
    if (alon2 < alon1) {
        double tmp = alon1;
        alon1 = alon2;
        alon2 = tmp;
    }
    if (blon2 < blon1) {
        double tmp = blon1;
        blon1 = blon2;
        blon2 = tmp;
    }
    double lat1, lon1, lat2, lon2, s12, azi1, azi2;
    if (!(blon1 > alon2 || blon2 < alon1)) {
        // intersecting
        lat1 = alat;
        lat2 = blat;
        lon1 = alon1 < blon1 ? alon1 : blon1;
        lon2 = alon2 > blon2 ? alon2 : blon2;
        lon1 = (lon1 + lon2) / 2;
        lon1 = lon1 < alon1 ? alon1 : lon1;
        lon1 = lon1 > alon2 ? alon2 : lon1;
        lon1 = lon1 < blon1 ? blon1 : lon1;
        lon1 = lon1 > blon2 ? blon2 : lon1;
        lon2 = lon1;
        if (ps12 || pazi1 || pazi2) {
            geokit_inverse(shape, lat1, lon1, lat2, lon2, &s12, &azi1, &azi2);
        }
    } else {
        double s122, azi12, azi22;
        lat1 = alat, lon1 = alon1;
        lat2 = blat, lon2 = blon1;
        geokit_inverse(shape, alat, alon1, blat, blon1, &s12, &azi1, &azi2);
        geokit_inverse(shape, alat, alon1, blat, blon2, &s122, &azi12, &azi22); 
        if (s122 < s12) {
            s12 = s122, azi1 = azi12, azi2 = azi22;
            lat1 = alat, lon1 = alon1;
            lat2 = blat, lon2 = blon2;
        }
        geokit_inverse(shape, alat, alon2, blat, blon1, &s122, &azi12, &azi22); 
        if (s122 < s12) {
            s12 = s122, azi1 = azi12, azi2 = azi22;
            lat1 = alat, lon1 = alon2;
            lat2 = blat, lon2 = blon1;
        }
        geokit_inverse(shape, alat, alon2, blat, blon2, &s122, &azi12, &azi22); 
        if (s122 < s12) {
            s12 = s122, azi1 = azi12, azi2 = azi22;
            lat1 = alat, lon1 = alon2;
            lat2 = blat, lon2 = blon2;
        }
    }
    if (plat1) *plat1 = lat1;
    if (plon1) *plon1 = lon1;
    if (plat2) *plat2 = lat2;
    if (plon2) *plon2 = lon2;
    if (ps12) *ps12 = s12;
    if (pazi1) *pazi1 = azi1;
    if (pazi2) *pazi2 = azi2;
}


static void flatmidpoint(double lat1, double lon1, double lat2, double lon2,
    double *plat3, double *plon3)
{
    *plat3 = (lat1 + lat2) / 2.0;
    *plon3 = (lon1 + lon2) / 2.0;
}


static void point_to_flat_lon(enum geokit_shape shape, double lat, double lon,
    double lat1, double lon1, double lat2, 
    double *plat, double *plon, double *ps12)
{


    double alat = lat1;
    double alon = lon1;
    double blat;
    double blon;
    double clat = lat2;
    double clon = lon1;
    double as12, bs12, cs12;
    double dist;
    geokit_inverse(shape, alat, alon, clat, clon, &dist, 0, 0);
    int iters = 0;
    while (iters < MAX_ITERS && dist > ACCURACY) {
        flatmidpoint(alat, alon, clat, clon, &blat, &blon);
        geokit_inverse(shape, lat, lon, alat, alon, &as12, 0, 0);
        geokit_inverse(shape, lat, lon, blat, blon, &bs12, 0, 0);
        geokit_inverse(shape, lat, lon, clat, clon, &cs12, 0, 0);
        if (as12 < cs12) {
            clat = blat;
            clon = blon;
            cs12 = bs12;
        } else {
            alat = blat;
            alon = blon;
            as12 = bs12;
        }
        geokit_inverse(shape, alat, alon, clat, clon, &dist, 0, 0);
        iters++;
    }
    if (iters == 0) {
        geokit_inverse(shape, lat, lon, alat, alon, &as12, 0, 0);
    }
    *plat = alat;
    *plon = alon;
    *ps12 = as12;
}

static void geokit_flatlon(enum geokit_shape shape,
    double alat1, double alat2, double alon, 
    double blat1, double blat2, double blon,
    double *plat1, double *plon1, double *plat2, double *plon2,
    double *ps12, double *pazi1, double *pazi2)
{

    double lat1, lon1, lat2, lon2, lat3, lon3, lat4, lon4;
    double s121, s122, s123, s124;

    point_to_flat_lon(shape, alat1, alon, blat1, blon, blat2, 
        &lat1, &lon1, &s121);
    point_to_flat_lon(shape, alat2, alon, blat1, blon, blat2, 
        &lat2, &lon2, &s122);
    point_to_flat_lon(shape, blat1, blon, alat1, alon, alat2, 
        &lat3, &lon3, &s123);
    point_to_flat_lon(shape, blat2, blon, alat1, alon, alat2, 
        &lat4, &lon4, &s124);
    double rlat1 = alat1;
    double rlon1 = alon;
    double rlat2 = lat1;
    double rlon2 = lon1;
    double rs12 = s121;
    if (s122 < rs12) {
        rlat1 = alat2;
        rlon1 = alon;
        rlat2 = lat2;
        rlon2 = lon2;
        rs12 = s122;
    }
    if (s123 < rs12) {
        rlat1 = blat1;
        rlon1 = blon;
        rlat2 = lat3;
        rlon2 = lon3;
        rs12 = s123;
    }
    if (s124 < rs12) {
        rlat1 = blat2;
        rlon1 = blon;
        rlat2 = lat4;
        rlon2 = lon4;
        rs12 = s124;
    }

    if (plat1) *plat1 = rlat1;
    if (plon1) *plon1 = rlon1;
    if (plat2) *plat2 = rlat2;
    if (plon2) *plon2 = rlon2;
    if (ps12) *ps12 = rs12;
    if (pazi1 || pazi2) {
        geokit_inverse(shape, rlat1, rlon1, rlat2, rlon2, 0, pazi1, pazi2);
    }
}

static bool equi_intersects(double aminlat, double aminlon, double amaxlat, 
    double amaxlon, double bminlat, double bminlon, double bmaxlat,
    double bmaxlon)
{
    return !((bminlat > amaxlat) | (bmaxlat < aminlat) |
        (bminlon > amaxlon) | (bmaxlon < aminlon));
}

void geokit_equipath(enum geokit_shape shape, 
    double aminlat, double aminlon, double amaxlat, double amaxlon, 
    double bminlat, double bminlon, double bmaxlat, double bmaxlon, 
    double *plat1, double *plon1, double *plat2, double *plon2, 
    double *ps12, double *pazi1, double *pazi2)
{
    if (equi_intersects(aminlat, aminlon, amaxlat, amaxlon, bminlat, bminlon,
        bmaxlat, bmaxlon))
    {
        // intersecting
        if (plat1) *plat1 = 0;
        if (plon1) *plon1 = 0;
        if (plat2) *plat2 = 0;
        if (plon2) *plon2 = 0;
        if (ps12) *ps12 = 0;
        if (pazi1) *pazi1 = 0;
        if (pazi2) *pazi2 = 0;
        return;
    }

    double lat11, lon11, lat21, lon21, s121;
    double lat12, lon12, lat22, lon22, s122;
    double lat13, lon13, lat23, lon23, s123;
    double lat14, lon14, lat24, lon24, s124;
    double lat15, lon15, lat25, lon25, s125;
    double lat16, lon16, lat26, lon26, s126;
    double lat17, lon17, lat27, lon27, s127;
    double lat18, lon18, lat28, lon28, s128;

    geokit_flatlat(shape, 
        aminlat, aminlon, amaxlon,
        bminlat, bminlon, bmaxlon,
        &lat11, &lon11, &lat21, &lon21, &s121, 0, 0);

    geokit_flatlat(shape, 
        aminlat, aminlon, amaxlon,
        bmaxlat, bminlon, bmaxlon,
        &lat12, &lon12, &lat22, &lon22, &s122, 0, 0);

    geokit_flatlat(shape, 
        amaxlat, aminlon, amaxlon,
        bminlat, bminlon, bmaxlon,
        &lat13, &lon13, &lat23, &lon23, &s123, 0, 0);

    geokit_flatlat(shape, 
        amaxlat, aminlon, amaxlon,
        bmaxlat, bminlon, bmaxlon,
        &lat14, &lon14, &lat24, &lon24, &s124, 0, 0);

    geokit_flatlon(shape, 
        aminlat, amaxlat, aminlon,
        bminlat, bmaxlat, bminlon,
        &lat15, &lon15, &lat25, &lon25, &s125, 0, 0);

    geokit_flatlon(shape, 
        aminlat, amaxlat, aminlon,
        bminlat, bmaxlat, bmaxlon,
        &lat16, &lon16, &lat26, &lon26, &s126, 0, 0);

    geokit_flatlon(shape, 
        aminlat, amaxlat, amaxlon,
        bminlat, bmaxlat, bminlon,
        &lat17, &lon17, &lat27, &lon27, &s127, 0, 0);

    geokit_flatlon(shape, 
        aminlat, amaxlat, amaxlon,
        bminlat, bmaxlat, bmaxlon,
        &lat18, &lon18, &lat28, &lon28, &s128, 0, 0);

    double lat1 = lat11;
    double lon1 = lon11;
    double lat2 = lat21;
    double lon2 = lon21;
    double s12 = s121;
    if (s122 < s12) {
        lat1 = lat12;
        lon1 = lon12;
        lat2 = lat22;
        lon2 = lon22;
        s12 = s122;
    }
    if (s123 < s12) {
        lat1 = lat13;
        lon1 = lon13;
        lat2 = lat23;
        lon2 = lon23;
        s12 = s123;
    }
    if (s124 < s12) {
        lat1 = lat14;
        lon1 = lon14;
        lat2 = lat24;
        lon2 = lon24;
        s12 = s124;
    }
    if (s125 < s12) {
        lat1 = lat15;
        lon1 = lon15;
        lat2 = lat25;
        lon2 = lon25;
        s12 = s125;
    }
    if (s126 < s12) {
        lat1 = lat16;
        lon1 = lon16;
        lat2 = lat26;
        lon2 = lon26;
        s12 = s126;
    }
    if (s127 < s12) {
        lat1 = lat17;
        lon1 = lon17;
        lat2 = lat27;
        lon2 = lon27;
        s12 = s127;
    }
    if (s128 < s12) {
        lat1 = lat18;
        lon1 = lon18;
        lat2 = lat28;
        lon2 = lon28;
        s12 = s128;
    }
    if (plat1) *plat1 = lat1;
    if (plon1) *plon1 = lon1;
    if (plat2) *plat2 = lat2;
    if (plon2) *plon2 = lon2;
    if (ps12) *ps12 = s12;
    if (pazi1 || pazi2) {
        geokit_inverse(shape, lat1, lon1, lat2, lon2, ps12, pazi1, pazi2);
    }
}
