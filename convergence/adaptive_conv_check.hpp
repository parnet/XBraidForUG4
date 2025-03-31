#ifndef ADAPTIVE_CONV_CHECK_HPP
#define ADAPTIVE_CONV_CHECK_HPP

namespace ug { namespace xbraid {
class AdaptiveConvCheck {
public:

    void set_tol(double loose = 1e-3, double tight = 1e-14) {
        this->loose_tol_ = loose;
        this->tight_tol_ = tight;
    }

    void set_max_iterations(size_t iter = 100) {
        this->max_iter_ = iter;
    }

private:
    size_t max_iter_;

    double loose_tol_;

    double tight_tol_;
};
}}
#endif
