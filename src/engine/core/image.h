

#ifndef PHOTOTOOLS_IMAGE1_H
#define PHOTOTOOLS_IMAGE1_H

#include <cstdint>

namespace core {
    /**
     * Image represents the physical image data in memory and it's associated properties.
     */
    class Image {
        // private as they shouldn't be changed independently of the image data.
        uint16_t m_width = 0;
        uint16_t m_height = 0;
    protected:

    public:

        /**
         * @return The width of the image.
         */
        [[nodiscard]] inline uint16_t get_width() const noexcept {
            return m_width;
        }

        /**
         * @return The height of the image.
         */
        [[nodiscard]] inline uint16_t get_height() const noexcept {
            return m_height;
        }
    };
}


#endif //PHOTOTOOLS_IMAGE_H
