// https://github.com/tidwall/geokit.c
// 
// Copyright 2023 Joshua J Baker. All rights reserved.
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file.

#ifndef GEOKIT_H
#define GEOKIT_H

#include <stdbool.h>

enum geokit_shape { GEOKIT_SPHERE, GEOKIT_ELLIPSOID };

void geokit_inverse(enum geokit_shape shape, 
    double lat1, double lon1, double lat2, double lon2, 
    double* ps12, double* pazi1, double* pazi2);

void geokit_direct(enum geokit_shape shape, 
    double lat1, double lon1, double azi1, double s12, 
    double* plat2, double* plon2, double* pazi2);

void geokit_track(enum geokit_shape shape, 
    double lat1, double lon1, double azi1,
    double lat2, double lon2,
    double *psat, double *psxt);

void geokit_intercept(enum geokit_shape shape,
    double lat1, double lon1, double azi1,
    double lat2, double lon2, double azi2,
    double *psat, double *psxt);

bool geokit_intersect(enum geokit_shape shape,
    double alat, double alon, double blat, double blon, 
    double clat, double clon, double dlat, double dlon,
    double *plat, double *plon);

void geokit_nearest(enum geokit_shape shape,
    double alat, double alon, double blat, double blon, 
    double clat, double clon, double dlat, double dlon,
    double *plat1, double *plon1, double *plat2, double *plon2, 
    double *ps12, double *pazi1, double *pazi2);

void geokit_equibox(enum geokit_shape shape, 
    double alat, double alon, double blat, double blon, 
    double *minlat, double *minlon, double *maxlat, double *maxlon);

void geokit_equipath(enum geokit_shape shape, 
    double aminlat, double aminlon, double amaxlat, double amaxlon, 
    double bminlat, double bminlon, double bmaxlat, double bmaxlon, 
    double *plat1, double *plon1, double *plat2, double *plon2, 
    double *ps12, double *pazi1, double *pazi2);


struct geokit_polygon { char priv[128]; };

void geokit_polygon_init(struct geokit_polygon *poly, enum geokit_shape shape,
    bool polyline);
void geokit_polygon_add(struct geokit_polygon *poly, double lat, double lon);
double geokit_polygon_area(struct geokit_polygon *poly);
double geokit_polygon_perimeter(struct geokit_polygon *poly);

#endif
