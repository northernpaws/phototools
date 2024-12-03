
#ifndef PHOTOTOOLS_IMAGE_H
#define PHOTOTOOLS_IMAGE_H

#include <type_traits>
#include <memory>
#include <cassert>
#include <vector>

#include "src/engine/core/color.h"
#include "src/engine/core/math/matrix.h"

/**
 * A catch-all class for representing a single image and it's data.
 */
struct Image {
protected:
    // TODO: replace with 2D dimension/size type
    uint16_t width = 0;
    uint16_t height = 0;

    /**
     * Count of color channels in the image.
     *
     * @see https://en.wikipedia.org/wiki/Channel_(digital_image)
     */
    uint8_t channels = 0;

//    std::unique_ptr<float[]> data = nullptr;
    // TODO: use some kind of global allocator for image mipmap data so the memory is centrally managed.
    std::vector<float> data;

    // TODO: add colorspace indicator
public:
    Image(uint16_t p_width, uint16_t p_height, uint8_t p_channels) {
        width = p_width;
        height = p_height;
        channels = p_channels;

        // Ensure the vector starts with sufficient memory to contain all the image
        // data so that the address doesn't change as we access the data.
        data.reserve(width * height * channels);
//        data = std::make_unique<float[]>(width * height * channels);
    }

    inline uint16_t get_width() const noexcept {
        return width;
    }

    inline uint16_t get_height() const noexcept {
        return height;
    }

    /**
     * Returns the count of channels in the image.
     * @return
     */
    inline uint8_t get_channels() const noexcept {
        return channels;
    }

    /**
     * @return A pointer to the image data, structured as row-major data
     *   with as many color components as specified by get_channels().
     */
    inline float *const get_data() noexcept {
        return data.data();
    };

    inline const float *const get_data() const noexcept {
        return data.data();
    };

    /**
     * Apply the provided 3x3 matrix as a color transformer to each pixel in the image.
     *
     * @tparam T
     * @param mat
     */
    template<typename T>
    inline void color_transform(const Matrix3x3<T>& mat) {
        if (channels < 3) {
            throw std::runtime_error("cannot 3x3 transform images with less then 3 color components (RGB)");
        }

        // TODO: credit below
        // method adapted from: https://github.com/LibRaw/LibRaw/blob/83bf3ad5e0744e78b82fa3c34e7542b2d62fcfa4/src/postprocessing/postprocessing_utils.cpp#L101C5-L118C6
        for (int row = 0; row < height; row++) {
            for (int col = 0; col < width; col++) {
                auto idx = pixel_index(col, row); // (row * stride) + (col * target.get_channels());

                auto r = mat[0][0] * data[idx+0] + mat[0][1] * data[idx+1] +
                        mat[0][2] * data[idx+2];
                auto g = mat[1][0] * data[idx+0] + mat[1][1] * data[idx+1] +
                        mat[1][2] * data[idx+2];
                auto b = mat[2][0] * data[idx+0] + mat[2][1] * data[idx+1] +
                        mat[2][2] * data[idx+2];

                // TODO: do we want to use clipping to 0,1 here, or allow unbounded values (i.e. HDR)?
                data[idx+0] = r; // CLIPF(out[0]);
                data[idx+1] = g; // CLIPF(out[1]);
                data[idx+2] = b; // CLIPF(out[2]);
            }
        }
    }


    /**
     * Retrieves the pixel at the specified position.
     *
     * @param index
     * @return
     */
     // TODO: doubt i'll use this so probably remove it later
    Pixel operator()(unsigned int x, unsigned int y) {
        assert(x < width);
        assert(y < height);

        auto index = pixel_index(x, y);

        Pixel pix;

        switch (channels) {
            case 1: {
                pix = data[index];
            } break;
            case 2: {
                pix = RG{data[index], data[index+1]};
            } break;
            case 3: {
                pix = RGB{data[index], data[index+1], data[index+2]};
            } break;
            case 4: {
                pix = RGBA{data[index], data[index+1], data[index+2], data[index+3]};
            } break;
            default: {
                std::vector<float> components;
                for (int c = 0; c < channels; c++) {
                    components[c] = data[index+c];
                }
                pix = components;
            } break;
        }

        return pix;
    }

    /**
     * Returns a reference to the index's place in the image data.
     *
     * @param index
     * @return
     */
    float& operator[] (unsigned int index) noexcept {
        assert(index < width * height * channels);

        return data[index];
    }

    /**
     * Returns a reference to the index's place in the image data.
     *
     * @param index
     * @return
     */
    const float& operator[] (unsigned int index) const noexcept {
        assert(index < width * height * channels);

        return data[index];
    }

    /**
     * Constant indexing operator for image data.
     *
     * @param index
     * @return
     */
//    float operator[] (int index) const noexcept {
//        assert(index < width * height * channels);
//
//        return data[index];
//    }

    inline std::size_t size() const {
        return width * height * channels * sizeof(float);
    }

    inline std::size_t pixels() const {
        return width * height * channels;
    }

    /**
     * Returns the index in the data array of the pixel at the provided position.
     *
     * @param x
     * @param y
     * @return
     */
    inline unsigned int pixel_index(unsigned int x, unsigned int y) const {
        assert(x < width);
        assert(y < height);

        auto idx = y * (width * channels) + (x * channels);

        assert(idx < (width * height * channels));

        return idx;
    }
};


//template <typename T>  struct Image {
//    static_assert(std::is_arithmetic<T>::value, "T must be numeric");
//
//protected:
//    // TODO: either use unique pointers or some kind of global allocator for image mipmaps.
//    T* data;
//
//    /**
//     * Count of color channels in the image.
//     *
//     * @see https://en.wikipedia.org/wiki/Channel_(digital_image)
//     */
//    uint8_t channels;
//
//public:
//    /**
//     *
//     * @return A pointer to the image data, structured as row-major data with the type defined by T, with as many color components as specified by get_channels();
//     */
//    virtual T* get_data() const = 0;
//};


#endif //PHOTOTOOLS_IMAGE_H
