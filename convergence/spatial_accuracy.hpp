#ifndef UGPLUGIN_XBRAIDFORUG4_CORE_SPATIAL_ACCURACY_HPP
#define UGPLUGIN_XBRAIDFORUG4_CORE_SPATIAL_ACCURACY_HPP

#import <cmath>

int GetSpatialAccuracy(double loose_tol,
                             double tight_tol,
                             double previous_tolerance,
                             double stopping_tolerance,
                             int num_rnorms,
                             int iteration,
                             double * rnorms) {
    double result;
    double interpolation;

    double log_rel_current_tol, log_loose_tol, log_tight_tol, log_rel_stopping_tol;




    double rnorm0 = rnorms[0];
    double rnorm = rnorms[iteration - 1];

    if ((iteration == 0)) {
        return loose_tol;
    }

    log_rel_current_tol = -log10(rnorm / rnorm0);
    log_rel_stopping_tol = -log10(stopping_tolerance / rnorm0);
    log_loose_tol = -log10(loose_tol);
    log_tight_tol = -log10(tight_tol);

    if (log_rel_current_tol >= (9.0 / 10.0) * log_rel_stopping_tol) {
        result = tight_tol;
    } else {
        interpolation = (log_rel_current_tol / log_rel_stopping_tol) * (log_tight_tol - log_loose_tol) + log_loose_tol;
        result = pow(10, -interpolation);

        if ((result > previous_tolerance) && (previous_tolerance > 0)) {
            result = previous_tolerance;
        }
    }

    return result;
}


#endif