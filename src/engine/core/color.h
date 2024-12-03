

#ifndef PHOTOTOOLS_COLOR_H
#define PHOTOTOOLS_COLOR_H

#include <cstddef>
#include <cstdint>
#include <array>
#include <type_traits>
#include <vector>

namespace core {

    /**
     *
     * @author Maciej Ciemborowicz
     * @see https://gist.github.com/ciembor/1494530#file-gistfile1-c-L61
     * @note Modified to fix some implicit casting between doubles and floats.
     *
     * @param p
     * @param q
     * @param t
     * @return
     */
    float hue2rgb(float p, float q, float t);

    /**
     *
     * @note I was debating between having an all-encompassing Color class vs RGBA, HSL, etc. classes, but ultimately
     *  decided that a singular RGBA-based color class had more sense due to camera sensors and other digital images
     *  only working in RGB data to begin with. Intermediate conversion to other formats can still be done, but having
     *  Color be RGBA helps implicitly enforce working in the RGBA format where ever possible for consistency.
     */
    struct Color {
        // Store the color components in a union of a couple naming conventions plus an array
        // so that they can be accessed in a few different standard methods for convenience.
        union {
            struct {
                float red;
                float green;
                float blue;
                float alpha;
            };

            struct {
                float r = 0.0f;
                float g = 0.0f;
                float b = 0.0f;
                float a = 1.0f;
            };

            std::array<float, 4> components;
        };

        Color() = default;

        /**
         * Construct a Color from the provided RGB floating values.
         *
         * It is assumed that the values are already normalized for a 0.0 to 1.0 floating range.
         *
         * @tparam T
         * @param rgb
         */
        template<typename T, std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
        Color(T p_r, T p_g, T p_b) {
            r = p_r;
            g = p_g;
            b = p_b;
            a = 1.0f;
        }

        /**
         * Construct a Color from the provided RGBA floating values.
         *
         * It is assumed that the values are already normalized for a 0.0 to 1.0 floating range.
         *
         * @tparam T
         * @param rgb
         */
        template<class T, std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
        Color(T p_r, T p_g, T p_b, T p_a) {
            r = p_r;
            g = p_g;
            b = p_b;
            a = p_a;
        }

        /**
         * Construct a Color from an array of RGB floating values.
         *
         * It is assumed that the values are already normalized for a 0.0 to 1.0 floating range.
         *
         * @tparam T
         * @param rgb
         */
        template<typename T, std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
        Color(std::array<T, 3> rgb) {
            r = rgb[0];
            g = rgb[1];
            b = rgb[2];
            a = 1.0f;
        }

        /**
         * Construct a Color from an array of RGBA floating values.
         *
         * It is assumed that the values are already normalized for a 0.0 to 1.0 floating range.
         *
         * @tparam T
         * @param rgb
         */
        template<typename T, std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
        Color(std::array<T, 4> rgba) {
            r = rgba[0];
            g = rgba[1];
            b = rgba[2];
            a = rgba[3];
        }

        /**
         * Construct a Color from an array of RGB integer values.
         *
         * This constructor assumes that the provided color components are using the correct integer type for their
         * resolution. The components will be converted to a 0.0 to 1.0 floating range using the numeric limit for
         * the provided integer type, so using the correct integer size is important.
         *
         * @tparam T
         * @param rgb
         */
        template<typename T, std::enable_if_t<std::is_integral<T>::value, bool> = true>
        Color(T p_r, T p_g, T p_b) {
            // TODO: We don't handle negatives here right now, but should add a type that does
            static_assert(std::is_unsigned<T>::value, "T should be an unsigned integer type");

            r = static_cast<float>(p_r) / std::numeric_limits<T>::max();
            g = static_cast<float>(p_g) / std::numeric_limits<T>::max();
            b = static_cast<float>(p_b) / std::numeric_limits<T>::max();
            a = 1.0f;
        }

        /**
         * Construct a Color from an array of RGBA integer values.
         *
         * This constructor assumes that the provided color components are using the correct integer type for their
         * resolution. The components will be converted to a 0.0 to 1.0 floating range using the numeric limit for
         * the provided integer type, so using the correct integer size is important.
         *
         * @tparam T
         * @param rgb
         */
        template<typename T, std::enable_if_t<std::is_integral<T>::value, bool> = true>
        Color(T p_r, T p_g, T p_b, T p_a) {
            // TODO: We don't handle negatives here right now, but should add a type that does
            static_assert(std::is_unsigned<T>::value, "T should be an unsigned integer type");

            r = static_cast<float>(p_r) / std::numeric_limits<T>::max();
            g = static_cast<float>(p_g) / std::numeric_limits<T>::max();
            b = static_cast<float>(p_b) / std::numeric_limits<T>::max();
            a = static_cast<float>(p_a) / std::numeric_limits<T>::max();
        }

        /**
         * Construct a Color from an array of RGB integer values.
         *
         * This constructor assumes that the provided color components are using the correct integer type for their
         * resolution. The components will be converted to a 0.0 to 1.0 floating range using the numeric limit for
         * the provided integer type, so using the correct integer size is important.
         *
         * @tparam T
         * @param rgb
         */
        template<typename T, std::enable_if_t<std::is_integral<T>::value, bool> = true>
        Color(std::array<T, 3> rgb) {
            // TODO: We don't handle negatives here right now, but should add a type that does
            static_assert(std::is_unsigned<T>::value, "T should be an unsigned integer type");

            r = static_cast<float>(rgb[0]) / std::numeric_limits<T>::max();
            g = static_cast<float>(rgb[1]) / std::numeric_limits<T>::max();
            b = static_cast<float>(rgb[2]) / std::numeric_limits<T>::max();
            a = 1.0f;
        }

        /**
         * Construct a Color from an array of RGBA integer values.
         *
         * This constructor assumes that the provided color components are using the correct integer type for their
         * resolution. The components will be converted to a 0.0 to 1.0 floating range using the numeric limit for
         * the provided integer type, so using the correct integer size is important.
         *
         * @tparam T
         * @param rgb
         */
        template<typename T, std::enable_if_t<std::is_integral<T>::value, bool> = true>
        Color(std::array<T, 4> rgb) {
            // TODO: We don't handle negatives here right now, but should add a type that does
            static_assert(std::is_unsigned<T>::value, "T should be an unsigned integer type");

            r = static_cast<float>(rgb[0]) / std::numeric_limits<T>::max();
            g = static_cast<float>(rgb[1]) / std::numeric_limits<T>::max();
            b = static_cast<float>(rgb[2]) / std::numeric_limits<T>::max();
            a = static_cast<float>(rgb[3]) / std::numeric_limits<T>::max();
        }

        /**
         * Construct a color from a RRGGBBAA formatted uint32.
         * @param color
         */
        Color(uint32_t color) {
            uint8_t red = (color & 0xFF000000) >> 24;
            uint8_t green = (color & 0x00FF0000) >> 16;
            uint8_t blue = (color & 0x0000FF00) >> 8;
            uint8_t alpha = (color & 0x000000FF);

            r = static_cast<float>(red) / 255.0f;
            g = static_cast<float>(green) / 255.0f;
            b = static_cast<float>(blue) / 255.0f;
            a = static_cast<float>(alpha) / 255.0f;
        }

        /**
         * Retrieve the RGB components of the color.
         * @return
         */
        inline std::array<float, 3> rgb() const {
            return {r, g, b};
        }

        /**
         * Retrieve the RGBA components of the color.
         * @return
         */
        inline std::array<float, 4> rgba() const {
            return {r, g, b, a};
        }

        /**
         * Converts an RGB color value to HSL. Conversion formula
         * adapted from http://en.wikipedia.org/wiki/HSL_color_space.
         *
         * @author Maciej Ciemborowicz
         * @see https://gist.github.com/ciembor/1494530#file-gistfile1-c-L18
         * @image https://camo.githubusercontent.com/2f9d65bfb0118b67f0e86810098b1635104bb5a05f31e78e707b5098c0b56003/68747470733a2f2f692e6b796d2d63646e2e636f6d2f656e74726965732f69636f6e732f6f726967696e616c2f3030302f3032362f3438392f637279696e672e6a7067
         *
         * @return HSL in the set [0, 1]
         */
        inline std::array<float, 3> hsl () const {
            std::array<float, 3> hsl{};

            float max = std::max(std::max(r,g),b);
            float min = std::min(std::min(r,g),b);

            hsl[0] = hsl[1] = hsl[2] = (max + min) / 2;

            if (max == min) {
                hsl[0] = hsl[1] = 0; // achromatic
            } else {
                float d = max - min;
                hsl[1] = (hsl[2] > 0.5) ? d / (2 - max - min) : d / (max + min);

                if (max == r) {
                    hsl[0] = (g - b) / d + (g < b ? 6.0f : 0.0f);
                } else if (max == g) {
                    hsl[0] = (b - r) / d + 2;
                } else if (max == b) {
                    hsl[0] = (r - g) / d + 4;
                }

                hsl[0] /= 6;
            }

            return hsl;
        }

        /**
         * Converts an HSL color value to RGB. Conversion formula
         * adapted from http://en.wikipedia.org/wiki/HSL_color_space.
         *
         * Assumes h, s, and l are contained in the set [0, 1].
         *
         * @author Maciej Ciemborowicz
         * @see https://gist.github.com/ciembor/1494530#file-gistfile1-c-L18
         * @image https://camo.githubusercontent.com/2f9d65bfb0118b67f0e86810098b1635104bb5a05f31e78e707b5098c0b56003/68747470733a2f2f692e6b796d2d63646e2e636f6d2f656e74726965732f69636f6e732f6f726967696e616c2f3030302f3032362f3438392f637279696e672e6a7067
         * @note Modified to fix some implicit casting between doubles and floats.
         *
         * @return
         */
        inline void set_hsl(float h, float s, float l) {
            if (0 == s) {
                // modification by alpheratz0
                // see: https://gist.github.com/ciembor/1494530?permalink_comment_id=4494915#gistcomment-4494915
                r = g = b = l * 255; // achromatic
            } else {
                float q = l < 0.5 ? l * (1 + s) : l + s - l * s;
                float p = 2 * l - q;

                r = hue2rgb(p, q, h + 1.0f/3) * 255;
                g = hue2rgb(p, q, h) * 255;
                b = hue2rgb(p, q, h - 1.0f/3) * 255;
            }
        }

        inline Color& operator*(const float v) {
            for(float & component : components) {
                component *= v;
            }

            return *this;
        }

        inline Color operator*(const float v) const {
            return {r * v, g * v, b * v, a * v};
        }

        /**
         * Returns a pointer to the R component in the RGBA list.
         *
         * @return
         */
        inline explicit operator const float *() const {
            return components.data();
        };
    };

    // Some extra helpers for working with specific color formats

    struct RG {
        union {
            struct {
                float r;
                float g;
            };

            std::array<float, 2> components;
        };
    };

    struct RGB {
        union {
            struct {
                float r;
                float g;
                float b;
            };

            std::array<float, 3> components;
        };
    };

    struct RGBA {
        union {
            struct {
                float r;
                float g;
                float b;
                float a;
            };

            std::array<float, 4> components;
        };
    };

    //using Pixel = std::variant<float, RG, RGB, RGBA, std::vector<float>>;

    // TODO: probably not going to use this anywhere, but it was a neat idea
    struct Pixel : std::variant<float, RG, RGB, RGBA, std::vector<float>> {
        using variant::variant;

        /**
         * Size of the pixel in color components.
         *
         * @return
         */
        inline std::size_t size() {
             if (index() < 4) {
                 // first 4 types in the variant have their index correspond to their size.
                 return index();
             }

             return std::get<std::vector<float>>(*this).size();
        }

        inline bool is_r() const {
            return std::holds_alternative<float>(*this);
        }

        inline bool is_rg() const {
            return std::holds_alternative<RG>(*this);
        }

        inline bool is_rgb() const {
            return std::holds_alternative<RGB>(*this);
        }

        inline bool is_rgba() const {
            return std::holds_alternative<RGBA>(*this);
        }

        inline bool is_vector() const {
            return std::holds_alternative<std::vector<float>>(*this);
        }

        /**
         * @note For vector types, up to the first 4 elements are converted as RGBA and the rest discarded.
         *
         * @return
         */
        operator Color() const {
            Color c;

            switch(index()) {
                case 0: { // r
                    c.r = std::get<float>(*this);
                } break;
                case 1: { // rg
                    RG rg = std::get<RG>(*this);;
                    c.r = rg.r;
                    c.g = rg.g;
                } break;
                case 2: { // rgb
                    RGB rg = std::get<RGB>(*this);;
                    c.r = rg.r;
                    c.g = rg.g;
                    c.b = rg.b;
                } break;
                case 3: { // rgba
                    RGBA rg = std::get<RGBA>(*this);;
                    c.r = rg.r;
                    c.g = rg.g;
                    c.b = rg.b;
                    c.a = rg.a;
                } break;
                case 4: { // vector<float>
                    std::vector<float> vec = std::get<std::vector<float>>(*this);;

                    for (int e = 0; e < std::min(4, (int)vec.size()); e++) {
                        c.components[e] = vec[e];
                    }
                } break;
            }

            return c;
        }
    };

}
/**
 * Represents a single RGBA color.
 *
 * Note that the base color type is locked to floats as the majority of the engine
 * operates on floating-point arithmetic for the image processing operations.
 */
/*struct Color {
    union {
        struct {
            float red;
            float green;
            float blue;
            float alpha;
        };

        struct {
            float r;
            float g;
            float b;
            float a;
        };

        std::array<float, 4> components;
    };

    Color() = default;

    Color(uint8_t p_r, uint8_t p_g, uint8_t p_b) {
        r = static_cast<float>(p_r) / 255.0f;
        g = static_cast<float>(p_g) / 255.0f;
        g = static_cast<float>(p_b) / 255.0f;
    }

    Color(uint8_t p_r, uint8_t p_g, uint8_t p_b, uint8_t p_a) {
        r = static_cast<float>(p_r) / 255.0f;
        g = static_cast<float>(p_g) / 255.0f;
        g = static_cast<float>(p_b) / 255.0f;
        a = static_cast<float>(p_a) / 255.0f;
    }
//Construct a color from a RRGGBBAA formatted uint32.
    Color(uint32_t color) {
        uint8_t red = (color & 0xFF000000) >> 24;
        uint8_t green = (color & 0x00FF0000) >> 16;
        uint8_t blue = (color & 0x0000FF00) >> 8;
        uint8_t alpha = (color & 0x000000FF);

        r = static_cast<float>(red) / 255.0f;
        g = static_cast<float>(green) / 255.0f;
        g = static_cast<float>(blue) / 255.0f;
        a = static_cast<float>(alpha) / 255.0f;
    }
};*/

// template specialization
/*template<size_t C>
struct Color {};

template<> struct Color<3> {

};

template<> struct Color<4> {

};*/

/*struct RGB {
    int red;
    int green;
    int blue;
};

struct RGBA {
    int red;
    int green;
    int blue;
    int alpha;
};

using Color = std::variant<RGB, RGBA>;
 // Example usage
    Color myColor = RGB{255, 0, 127};

    // Accessing the color components
    if (auto rgb = std::get_if<RGB>(&myColor)) {
        std::cout << "RGB color: (" << rgb->red << ", " << rgb->green << ", " << rgb->blue << ")\n";
    } else if (auto rgba = std::get_if<RGBA>(&myColor)) {
        std::cout << "RGBA color: (" << rgba->red << ", " << rgba->green << ", " << rgba->blue << ", " << rgba->alpha << ")\n";
    } else {
        std::cout << "Unknown color type\n";
    }*/

//struct RGB {
//    union {
//        struct {
//            float red;
//            float green;
//            float blue;
//        };
//
//        std::array<float, 3> components;
//    };
//};
//
//struct RGBA {
//    union {
//        struct {
//            float red;
//            float green;
//            float blue;
//            float alpha;
//        };
//
//        std::array<float, 4> components;
//    };
//};
//
//using Color = std::variant<RGB, RGBA>;

/*using ColorBase = std::variant<RGB, RGBA>;

// ref: https://stackoverflow.com/a/70566152
template<template<typename...> class V, typename... Vargs>
struct map_inner {
    map_inner(const V<Vargs...>&){}

    template<template<typename> class Wrap>
    using type = V<Wrap<Vargs>...>;
};

template<typename T>
struct my_wrapper {};

// ref: https://stackoverflow.com/a/70566152
using Color = typename decltype(map_inner{std::declval<ColorBase>()})::type<my_wrapper>; // c++17 version

static_assert(std::is_same_v<std::variant<my_wrapper<RGB>, my_wrapper<RGBA>>, Color>);*/

//struct Color: public std::variant<RGB, RGBA> {
//    RGB rgb() const {
//        return std::get_if<RGB>;
//    }
//};

#endif //PHOTOTOOLS_COLOR_H
