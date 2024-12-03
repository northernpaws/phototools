
#ifndef PHOTOTOOLS_RAW_IMAGE_H
#define PHOTOTOOLS_RAW_IMAGE_H

#include <cstdint>
#include <memory>
#include <functional>
#include <cassert>

#include "src/engine/core/math/matrix.h"

struct ColorFilterArray {
    uint8_t width = 0;
    uint8_t height = 0;

    std::vector<unsigned char> layout;

    /**
     * @return count of the filters in the CFA.
     */
    inline int count() const {
        return width*height;
    }

    inline bool is_xtrans() const {
        return width == 6 && height == 6;
    }

    inline bool is_bayer() const {
        // TODO: double check this is correct.
        return width == 4 && height == 4;
    }

    inline unsigned char operator[](int index) const {
        return layout[index];
    }
};

// TODO: Sensor class (with types, bayer, xtrans, etc) that encapsulates stuff like the CFA, to account for "sensors" without a CFA

enum RawImageFormat: int {
    UINT16,
    FLOAT32
};

/**
 * Holds data about a camera RAW image, including the un-mosaiced RAW data, CFA (Color Filter Array) layout, camera
 *  colorspace conversion matrix(es), and the camera metadata.
 *
 * Note that the RawImage class works with RAW image data in it's original format - typically uint16_t or float. This is
 *  because some RAW handling libraries and algorithms expect the raw data to be in those formats, so maintaining those
 *  formats until we actually need floats for the develop pipeline is useful.
 *
 * TODO: At the moment this is only designed for CFA images (i.e. X-trans or Bayer) but could probably easily be
 *   expanded to other raw types supported by Libraw and Rawspeed - i.e. grayscale scanner images, IR cameras, etc.
 */
class RawImage {
public:
protected:
    uint16_t width = 0;
    uint16_t height = 0;

    // describes the datatype in the data field.
    RawImageFormat format;

    ColorFilterArray cfa;

    // TODO: store data as a uint8_t for a catch-all, and provide intermediate casting methods.

    // note: most camera raw data is either uint16 or float32, but we'll use the smallest of the two for compatability and provide utility functions for reinterpreting as the correct types.
    // uses a unique ptr to help ensure the class always maintains ownership over its data.
    // TODO: ensure the deleter is specified correctly and actually deleting all the data.
    std::unique_ptr<uint16_t[]> data = nullptr;
public:
    // TODO: this can be probably handled better, maybe as a generic method of attaching camera-embedded color matrixes.
    // TODO: need to support 4-color sensors (i.e. RGGB to RGB) in the future via 4x3 matrixes
    //   https://github.com/LibRaw/LibRaw/blob/83bf3ad5e0744e78b82fa3c34e7542b2d62fcfa4/src/postprocessing/postprocessing_utils.cpp#L120
    /**
     * A row-major 3x3 matrix that converts from CIE XYZ colorspace values to the camera's native RGB colorspace.
     *
     * This value is sourced either directly from the camera RAW files (for example, TIFF ColorMatrix2)or from a
     * table of user-calibrated conversion matrices for specific cameras, such as from dcraw, rawspeed, or libraw.
     *
     * The inverse of this matrix is used by processing steps to convert from the camera reference colorspace into
     * other colorspace (i.e. sRGB, CIELAB, Linear Rec2020, etc.) that are required for processing the image data.
     *
     * @see https://www.awaresystems.be/imaging/tiff/tifftags/colormatrix2.html
     *
     * @note (rawspeed) under calibration illuminant 21 (D65)
     */
    Matrix3x3d xyz_to_cam;

    // whitebalance coffieciants
    // TODO: move protected and define better
    std::array<float, 4> whitebalance;

    // note: assumes uint16 data
    int white_point;
    // note: assumes uint16 as well
    int black_level;

    RawImage() = default;

    RawImage(RawImageFormat p_format, uint16_t p_width, uint16_t p_height, ColorFilterArray p_cfa, std::unique_ptr<uint16_t[]> p_data) {
        format = p_format;
        width = p_width;
        height = p_height;
        cfa = p_cfa;
        data = std::move(p_data);
    }

    RawImage(RawImageFormat p_format, uint16_t p_width, uint16_t p_height, ColorFilterArray p_cfa, uint16_t* p_data) : RawImage(p_format, p_width, p_height, p_cfa, std::unique_ptr<uint16_t[]>(p_data)) {
        format = p_format;
        width = p_width;
        height = p_height;
        cfa = p_cfa;
    }

    inline uint16_t get_width() const noexcept {
        return width;
    }

    inline uint16_t get_height() const noexcept {
        return height;
    }

    /**
     * Count of bytes per color component.
     * @return
     */
    inline std::size_t bytes_per_component() const {
        switch (format) {
            case UINT16:
                return sizeof(uint16_t);
            case FLOAT32:
                return sizeof(float);
        }
    }

    /**
     * Returns the size of the raw image in bytes.
     * @return
     */
    inline std::size_t size() const {
        return width * height * bytes_per_component();
    }

    inline const ColorFilterArray& get_cfa() const noexcept {
        return cfa;
    }

    inline const bool is_valid() const noexcept {
        return width != 0 && height != 0 && data != nullptr;
    }

    /**
     * @return RAW image data as a uint16_t array.
     */
    inline uint16_t* uint16() const {
        assert(format == RawImageFormat::UINT16 &&
               "Attempting to access floating-point buffer as uint16_t.");
        assert(data != nullptr && "Data not yet allocated.");

        return data.get();
    }

    /**
     * @return RAW image data as a float32 array.
     */
    inline float* float32() const {
        assert(format == RawImageFormat::FLOAT32 &&
               "Attempting to access uint16_t buffer as floating-point.");
        assert(data != nullptr && "Data not yet allocated.");

        return reinterpret_cast<float*>(data.get());
    }

    uint16_t& operator[](unsigned int index) {
        return data.get()[index];
    }

};


#endif //PHOTOTOOLS_RAW_IMAGE_H
