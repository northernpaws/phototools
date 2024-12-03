

#include "color.h"

float core::hue2rgb(float p, float q, float t) {
    if (t < 0)
        t += 1;
    if (t > 1)
        t -= 1;
    if (t < 1.0/6.0f)
        return p + (q - p) * 6 * t;
    if (t < 1.0/2.0f)
        return q;
    if (t < 2.0/3.0f)
        return p + (q - p) * (2.0f/3.0f - t) * 6.0f;
    return p;
}