#include "data_structure.hpp"

/*double xp_dot{0.0};  // longitudinal speed
    double yp_dot{0.0};  // lateral speed
    double psi_dot{0.0}; //
    double epsi{0.0};    //
    double ey{0.0};      // lateral position
    double s{0.0};       // longitudinal position
    double steer{0.0};
    double acc{0.0};*/
FB_state::FB_state(double a, double b, double c, double d, double e, double f, double g, double h):
    xp_dot(a), yp_dot(b), psi_dot(c), epsi(d),
    ey(e), s(f), steer(g), acc(h)
{

}
